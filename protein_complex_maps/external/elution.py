from __future__ import division
import numpy as np
import operator
import os
import itertools as it
import random
import orth
import pairdict as pd
import score as sc
from Struct import Struct
import utils as ut

def load_elution(fname, getname=True):
    # expected file structure:
    # first col: gene id
    # second col: treat differently if 2nd col header is 'Total' or
    # 'Description'
    # remaining cols: elution profile data
    lines = [l for l in ut.load_tab_file(fname)]
    # final row: total count in msblender output; don't skip in cuihong's data
    skip_final_row = (lines[-1][0][0] == '#')
    rows = lines[1:-1] if skip_final_row else lines[1:]
    fractions = [f for f in lines[0][1:]]
    if fractions[0].lower() in ['total', 'totalcount', 'description']:
        start_data_col = 2
        fractions.remove(fractions[0])
    else:
        start_data_col = 1
    mat = np.matrix([row[start_data_col:] for row in rows],dtype='float32')
    prots = [row[0] for row in rows]
    elut = Struct(mat=mat, prots=prots, fractions=fractions, filename=fname,
                  filename_original=fname)
    if start_data_col == 2:
        col2name_vals = [row[1] for row in rows]
        elut.column2vals = col2name_vals
    if getname: elut.name = os.path.basename(fname).split('.')[0]
    return elut

class NormElut(Struct):

    def __init__(self, filename, sp_base='Hs', norm_rows=False, norm_cols=False):
        e = load_elution(filename)
        self.prots = e.prots
        self.filename = e.filename
        self.normarr = ut.normalize_fracs(e.mat, norm_rows=norm_rows,
                norm_cols=norm_cols)
        self.pinv = ut.list_inv_to_dict(e.prots)
        sp_target = ut.shortname(e.filename)[:2]
        self.baseid2inds = sc.orth_indices(sp_base, sp_target, e.prots, False)



def all_filtered_pairs(fnames, score_keys, cutoff=0.5, sp_base=None,
                       verbose=True, allow_singles=True):
    allpairs = pd.PairDict([])
    for skey,f in it.product(score_keys,fnames):
        if verbose: print skey, cutoff, ut.shortname(f)
        elut = load_elution(f)
        newpairs = passing_pairs(elut, skey, cutoff, allow_singles)
        newpairs = translate_pairs(newpairs, sp_base, file_sp(f))
        allpairs = pd.pd_union_novals(allpairs, newpairs)
    return allpairs

def passing_pairs(elut, score_key, cutoff, allow_singles):
    pairs = sc.pairs_exceeding(elut, score_key, thresh=cutoff)
    singles = sc.prots_singles(elut)
    if allow_singles:
        newpairs = pd.PairDict(pairs)
    else:
        newpairs = pd.PairDict(((p1,p2) for p1,p2 in pairs 
            if p1 not in singles and p2 not in singles))
    return newpairs

def supporting_ppis(ppis, fnames, score_keys, sp_base, cutoff=0.5, verbose=True):
    ppis_support = [pd.PairDict([]) for p in ppis]
    eluts = [load_elution(f) for f in fnames]
    for elut,skey in it.product(eluts, score_keys):
        if verbose: print skey, ut.shortname(elut.filename)
        od = orth.odict(sp_base, file_sp(elut.filename))
        new_pairs = passing_pairs(elut, skey, cutoff)
        for p,pdsupport in zip(ppis,ppis_support):
            for opair in orth.orth_pairs(p[:2], od):
                opair = tuple(opair)
                if new_pairs.contains(opair):
                    pdsupport.set(opair,None)
    return [list(p) + [s.d.keys()] for p,s in zip(ppis, ppis_support)]

def translate_pairs(pairs, sp_base, sp_target):
    t2b = targ2base(sp_base, sp_target)
    if t2b:
        pairs = pd.PairDict(((base1,base2)
                          for t1,t2 in pairs.d
                          for base1,base2 in
                          it.product(t2b.get(t1,[]), t2b.get(t2,[]))))
    return pairs
    
def targ2base(sp_base, sp_target):
    t2b = None
    if sp_base:
        if sp_base != sp_target:
            t2b = orth.odict(sp_target, sp_base)
    return t2b

def file_sp(filename):
    return ut.shortname(filename)[:2]

def all_prots(elut_fs, sp_base=None, min_count=2):
    print 'Loading all proteins from files'
    allprots = set([])
    for f in elut_fs:
        t2b = targ2base(sp_base, file_sp(f))
        elut = load_elution(f)
        if min_count:
            newprots = \
            set(np.array(elut.prots)[np.where(np.max(elut.mat,axis=1) >=
                min_count)[0]][0])
        else:
            newprots = set(elut.prots) 
        if t2b:
            newprots = set((b for t in newprots for b in t2b.get(t,[])))
        allprots = set.union(allprots, newprots)
    return allprots


##################################################
# Special-purpose or single-use
##################################################

def _fraction_elutions(fractions):
    """
    Given a list of fraction names, group them after removing 'FractionXX' from
    the end.  Return a dict of { elutionname: listofindices }.
    Example of fraction name: Orbitrap_HeLaCE_IEF_pH3_to_10_Fraction10
    """
    elution_names = {}
    for i,fname in enumerate(fractions):
        ename = fname[:fname.find('_Fraction')]
        elution_names.setdefault(ename,[]).append(i)
    return elution_names
    
def split_muliple_elutions(big_elut):
    """
    Split an elution into multiple based on use of _fraction_elutions.
    """
    elution_columns = _fraction_elutions(big_elut.fractions)
    eluts = {}
    for elution_name in elution_columns:
        new_elut = Struct()
        new_elut.__dict__ = big_elut.__dict__.copy()
        new_elut.mat = big_elut.mat[:,elution_columns[elution_name]]
        new_elut.fractions = list(np.array(big_elut.fractions)[elution_columns[elution_name]])
        new_elut.filename = big_elut.filename + '__' + elution_name
        eluts[elution_name] = new_elut
    return eluts

def write_elution(elut, fname, forR=False):
    """
    Write out an elution in the spcount format
    $ProtID\tTotalCount\tCol1....
    """
    # First eliminate empty protein rows
    nonzeros = np.sum(np.array(elut.mat),axis=1)>0
    arr = np.array(elut.mat[nonzeros,:])
    prots = list(np.array(elut.prots)[nonzeros])
    if not forR:
        header = "#ProtID TotalCount".split() + elut.fractions
        data = [[prots[i], np.sum(arr[i,:])] + arr[i,:].tolist() for i in
                range(len(prots))]
    else: #R: no column header for first column, and transpose
        header = prots
        data = [[elut.fractions[i]] + arr[:,i].tolist() for i in
                range(len(elut.fractions))]
    ut.write_tab_file([header] + data, fname)

def process_raw_wan(f_source, f_dest=None, first_col_element=1,
                    first_data_col=1, end_description_col=True,
                    first_data_row=1):
    # specific to cuihong's files, and tries to handle the differences seen in them
    # always keeps the first column as variable name
    # processes first column, splitting it and keeping first_col_element
    # for the array, keeps columns [first_data_col:end_data_col]. None works.
    # textlines = [textline for textline in open(f_source)]
    # # handle unix-unreadable linebreaks from excel
    # if len(textlines) == 1:
    #     if textlines[0].find('\r\n') > -1:
    #         textlines = textlines[0].split('\r\n')
    lines = [line.strip().split('\t') for line in open(f_source)if line.strip()!='']
    # simple: one step at a time.
    # column manipulation first.
    if end_description_col:
        lines = [[l[0]] + [l[-1]] + l[first_data_col:-1] for l in lines]
    else:
        lines = [[l[0]] + l[first_data_col:] for l in lines]
    # variable name manipulation
    if first_col_element is not None:
        # manipulate gene name in all but header row. skip anything btw header
        # and first_data_row.
        lines = [lines[0]] + [[l[0].split('|')[first_col_element]] +
                    l[1:] for l in lines[first_data_row:]]
    # rename file
    if f_dest is None:
        split = os.path.splitext(f_source)
        f_dest = split[0] + '_proc' + split[1]
    ut.write_tab_file(lines, f_dest)

def correlate_single(elut1, elut2, prot):
    return np.corrcoef(elut1.mat[elut1.prots.index(prot),:],
        elut2.mat[elut2.prots.index(prot),:])[0][1]

def correlate_matches(elut1, elut2):
    overlap = set.intersection(set(elut1.prots),set(elut2.prots))
    return [correlate_single(elut1, elut2, p) for p in overlap]
    
def correlate_matches_dict(elut1, elut2, pdict_1to2):
    overlap = [p for p in elut1.prots if p in pdict_1to2 and
        list(pdict_1to2[p])[0] in set(elut2.prots)]
    return [(p,np.corrcoef(elut1.mat[elut1.prots.index(p),:],
        elut2.mat[elut2.prots.index(list(pdict_1to2[p])[0]),:])[0][1]) for p in
        overlap]
    
def combine_elutions(e1, e2, combine_corr_func=None):
    # Assumes the fractions from each elution are mutually exclusive; puts them
    # in order of e1fracs+e2fracs.
    # Proteins (rows) are merged.
    allprots = list(set.union(set(e1.prots), set(e2.prots)))
    nprots = len(allprots)
    # use n fractions instead of matrix shape to handle 0-row elutions
    nfracs1 = len(e1.fractions) 
    allfracs = nfracs1 + len(e2.fractions)
    mat = np.matrix(np.zeros((nprots,allfracs)))
    mat[0:len(e1.prots),0:e1.mat.shape[1]] = e1.mat[:,:]
    for elut,(start,stop) in [(e1,(0,nfracs1)),(e2,(nfracs1,None))]:
        for row in range(len(elut.prots)):
            mat[allprots.index(elut.prots[row]), start:stop] = elut.mat[row,:]
    elut = Struct(mat=mat, prots=allprots, fractions=e1.fractions+e2.fractions,
                  filename=e1.filename+e2.filename+str(combine_corr_func))
    if combine_corr_func:
        elut.corr = combine_corrs(e1, e2, allprots, combine_corr_func)
    return elut

def combine_corrs(e1, e2, allprots, combine_func, default_val=None):
    # we combine the symmetric correlation matrices using the specified
    # element-wise function. function examples: max, sum
    # we use the specified ordering of elements in allprots
    default_val = default_val if default_val else -1 if \
        combine_func.__name__.find('max') > -1 else 0
    nprots = len(allprots)
    corr = np.matrix(np.zeros((nprots,nprots)))
    dprots1 = ut.list_inv_to_dict(e1.prots)
    dprots2 = ut.list_inv_to_dict(e2.prots)
    for row,p1 in enumerate(allprots):
        for col,p2 in enumerate(allprots):
            val1 = e1.corr[dprots1[p1], dprots1[p2]] if p1 in dprots1 and p2 in \
                dprots1 else default_val
            val2 = e1.corr[dprots2[p1], dprots2[p2]] if p1 in dprots2 and p2 in \
                dprots2 else default_val
            corr[row,col] = combine_func(val1, val2)
    return corr
    
def test_combined_corrs(eluts, ncomparisons=10):
    # compare at ncomparisons randomly-selected places
    prots_common = list(reduce(set.union,[set(e.prots) for e in eluts]))
    p1s = random.sample(prots_common, ncomparisons)
    p2s = random.sample(prots_common, ncomparisons)
    return [[e.corr[e.prots.index(p1),e.prots.index(p2)] for e in eluts] for
(p1,p2) in zip(p1s, p2s)]

def downsample_elution(elution, downsample, seed=0):
    """
    Return a new elution with every downsample-th fraction.
    """
    down_elut = Struct()
    down_elut.__dict__ = elution.__dict__.copy()
    down_elut.mat = elution.mat[:,seed::2]
    down_elut.fractions = elution.fractions[::2]
    down_elut.name = elution.name + '_down%i' % downsample
    return(down_elut)

def subset_elution(elution, prot_set):
    """ 
    Return an elution only containing the proteins contained in the
    provided prot_set.
    """
    newel = Struct()
    newel.__dict__ = elution.__dict__.copy()
    prot_inds, newprots = zip(*[(i,p) for i,p in enumerate(elution.prots) 
        if p in prot_set])
    newel.mat = elution.mat[prot_inds,:]
    newel.prots = newprots
    print len(newel.prots), 'prots from', elution.filename, 'in set'
    return newel

def filter_matching_elution(edata, efilter, remove_data_ending='.map'):
    """
    Use efilter as a mask on edata, setting edata to 0 based on efilter being 0
    """
    newmat = np.matrix(np.zeros(edata.mat.shape))
    # First create the column-matched array from filter to data
    inv_fracs = ut.list_inv_to_dict(efilter.fractions)
    data_fracs = [f.replace(remove_data_ending,"") for f in edata.fractions]
    arrfilt = np.zeros((efilter.mat.shape[0], edata.mat.shape[1]))
    for i,f in enumerate(data_fracs):
        arrfilt[:,i] = (np.asarray(efilter.mat)[:,inv_fracs[f]] if f in
        inv_fracs else np.zeros(efilter.mat.shape[0]))
    # Then go row-by-row
    filter_map = ut.list_inv_to_dict(efilter.prots)
    for i,g in enumerate(edata.prots):
        if g in filter_map:
            newmat[i,:] = (np.asarray(edata.mat)[i,:] *
                    (arrfilt[filter_map[g],:] > 0).astype(int))
        else:
            newmat[i,:] = np.zeros(edata.mat.shape[1])
    return newmat

def sort_elution(elution):
    newel = ut.struct_copy(elution)
    newel.mat = np.array(newel.mat)
    inds = np.argsort(np.sum(newel.mat, axis=1))[::-1]
    newel.prots = list(np.array(newel.prots)[inds])
    newel.mat = newel.mat[inds]
    return newel

def convert_elution(elution, odict, gt=None, sep='|'):
    newel = ut.struct_copy(elution)
    newel.prots = ['|'.join([gt.id2name.get(g,g) for g in odict.get(xg,[xg])]) 
            for xg in elution.prots]
    return newel
