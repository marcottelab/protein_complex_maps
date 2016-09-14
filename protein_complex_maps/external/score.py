import sys
import itertools
import numpy as np
from scipy import spatial
import os
from scipy import sparse
from collections import defaultdict
import operator
import utils as ut
import elution as el
import orth



def score_array_multi(arr, sp_base, elut_fs, scores, cutoff, verbose=False,
        remove_multi_base=False, gidscheme=None, allow_singles=True):
    """
    - remove_multi_base: This is not the method currently used to filter scores
      in cases of orthogroup fan-outs--this is a stricter earlier version. That
      filter is feature.py: filter_multi_orths(), applied after scoring.
    """
    assert gidscheme=='', "Gidscheme not implemented in scoring."
    current_sp = ''
    if remove_multi_base: 
        print ("Filtering orths: only single base gene in orthogroups.")
    for e,f in [(el.load_elution(f),f) for f in elut_fs]:
        sp_target = ut.shortname(f)[:2]
        if sp_target != current_sp: # Just for status output
            print "Starting first %s file: %s" % (sp_target, ut.shortname(f))
            current_sp = sp_target
        baseid2inds = orth_indices(sp_base, sp_target, e.prots,
                remove_multi_base)
        # singles based on original spec counts
        singles = set([]) if allow_singles else prots_singles(e) 
        for score in scores:
            if verbose: print score, f
            score_array(arr, e, f, score, cutoff, baseid2inds, singles, lambda prots:
                    orth_indices(sp_base, sp_target, prots, remove_multi_base))

def orth_indices(sp_base, sp_target, prot_list, remove_multi_base):
    """
    Using appropriate orthology, take a list of target species gene ids
    (corresponding to rows in the target species score matrix), and
    return a dict mapping base species gene ids to (sets of) indices in that
    list and therefore to (sets of) row/column indices in the square
    interaction score matrix. 
    """
    targ2inds = dict([(k,set([v]))
                      for k,v in ut.list_inv_to_dict(prot_list).items()])
    if sp_base == sp_target:
        return targ2inds
    else:
        base2targ = orth.odict(sp_base, sp_target)
        if remove_multi_base:
            base2targ = remove_multi_keys(base2targ)
        base2inds = ut.compose_dict_sets(base2targ, targ2inds)
        base2inds = dict([(k,v) for k,v in base2inds.items() if len(v)>0])
        return base2inds

def remove_multi_keys(d, max_keys=1):
    """
    Given a dict of key: set(vs), eliminate from the dict any keys that map to
    the same set of vs.
    """
    newd = d.copy()
    dinv = ut.dict_inverse_sets(newd)
    for k,vs in newd.items():
        for v in vs:
            if len(dinv[v]) > max_keys:
                del newd[k]
                break
    return newd

def score_array(arr, elut, fname, score, cutoff, id2inds, singles_exclude,
        recalc_id2inds):
    """
    Use the target species score matrix to get interaction pair in the base
    species array.  Don't score and just leave as default (0 now) cases where
    either: 1) One of the pair is not in this score matrix, or 2) The two base
    ids in the pair map to identical targets, since in that case we also can
    get no information from this data (see notes 2012.08.12).
    Also exclude any proteins with just one total count in this elution.
    - Recalc_id2inds: purpose is for remapping to the right indices in the case
      of swiching out to a new elution file with a differently-ordered matrix.
      This is currently only to handle the ms1 elution data.
    """
    score_mat, new_id2inds, new_prots = scorekey_elution(score, elut, recalc_id2inds)
    id2inds = new_id2inds or id2inds
    prots = new_prots or elut.prots
    score_name = name_score(fname,score)
    for i,row in enumerate(arr):
        id1,id2 = row['id1'],row['id2']
        if id1 in id2inds and id2 in id2inds and id2inds[id1]!=id2inds[id2]: 
            inds1, inds2 = [id2inds[gid] for gid in id1,id2]
            if len(singles_exclude) > 0:
                inds1, inds2 = [remove_labeled(inds, prots, singles_exclude)
                    for inds in inds1,inds2]
            if len(inds1)>0 and len(inds2)>0:
                # Could also check for i!=j but would have no effect here since
                # these mappings come from disjoint orthogroups.
                row[score_name] = max([score_mat[i,j] 
                    for i in inds1 for j in inds2])

def remove_labeled(ids, labels, set_remove):
    return [i for i in ids if labels[i] not in set_remove]

def name_score(fname, score):
    return ut.shortname(fname) + '_' + score 

def prots_singles(elut):
    """
    Using where to find proteins with only one count: messy but fast
    """
    singles_inds = np.array(np.where(elut.mat.sum(axis=1) == 1)[0])[0]
    return set(np.array(elut.prots)[singles_inds])

def scorekey_elution(score, elut, recalc_id2inds):
    new_id2inds = None
    new_prots = None
    if score == 'apex':
        score_mat = ApexScores(elut)
    elif score == 'cosine_old':
        score_mat = CosineLazyScores(elut)
    elif score == 'cosine':
        score_mat = CosineLazyNew(elut)
    elif score == 'euclidean':
        score_mat = pdist_score(elut.mat, norm_rows=True, norm_cols=True,
                metric=score)
    elif score in ('pq_euc', 'pq_unfilt_euc', 'mq_euc'):
        # Use pepquant specific elution file.
        extension = ( '_pqmsb_filtmsb.tab' if score=='pq_euc' else
                '_pqmsb.tab' if score=='pq_unfilt_euc' else
                '.mq_Intensity.tab' if score=='mq_euc' else 0)
        elut = el.load_elution(os.path.splitext(elut.filename)[0] + extension)
        if recalc_id2inds is not None:
            new_id2inds = recalc_id2inds(elut.prots) #cv framework (arrfeats)
        new_prots = elut.prots 
        score_mat = pdist_score(elut.mat, norm_rows=True, norm_cols=True,
                metric='euclidean')
    else:
        fscore = elut.filename + (
                  '.corr_poisson' if score=='poisson' else
                  '.T.wcc_width1' if score=='wcc' else
                  '.corr_euclidean' if score=='euc_poisson' else
                  '.standard' if score=='standard' else # eg elution/testms1
                  0 ) # no score: exception since string and int don't add
        score_mat = precalc_scores(fscore)
    return score_mat, new_id2inds, new_prots

def traver_corr(mat, repeat=1000, norm='columns', verbose=True):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each column to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat, iteration_display, norm):
        if verbose: print iteration_display
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        if norm=='columns': 
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
        elif norm=='rows': # seems to make no performance difference 1/25
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
        corr = np.nan_to_num(np.corrcoef(poisson_mat))
        return corr
    avg_result = (reduce(operator.add, (poisson_corr(mat, i, norm=norm) for i in
                                        range(repeat))) / repeat)
    return avg_result

def pdist_score(mat, metric='euclidean', norm_rows=True,
        norm_cols=True):
    norm_mat = ut.normalize_fracs(mat, norm_rows, norm_cols)
    dists = spatial.distance.pdist(norm_mat, metric=metric)
    dist_mat = spatial.distance.squareform(dists)
    score_mat = 1 - np.nan_to_num(dist_mat)
    return score_mat

def poisson_repeat(mat, repeat=200, **kwargs):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each column to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_dist(mat, iteration_display, metric='cosine', norm_rows=True,
            norm_cols=True, verbose=True):
        if verbose: print iteration_display, metric
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        score_mat = pdist_score(mat, metric=metric,
                norm_rows=norm_rows, norm_cols=norm_cols)
        return score_mat
    avg_result = (reduce(operator.add, (poisson_dist(mat, i, **kwargs) for i in
                                        range(repeat))) / repeat)
    return avg_result

class ApexScores(object):

    def __init__(self, elution):
        self.apex_array = np.argmax(np.array(elution.mat), axis=1)
        self.shape = (len(self.apex_array),len(self.apex_array))

    def __getitem__(self, index):
        return int(self.apex_array[index[0]] == self.apex_array[index[1]])

def precalc_scores(scoref, dtype='f2'):
    """
    Also zero out the diagonal to more efficiently remove all self-interactions
    up-front.
    """
    # NOTE to change dtype you must change it in loadtxt below!!
    save_compact = ut.config()['save_compact_corrs'] 
    compactf = '%s.%s.pyd' % (scoref, dtype)
    if os.path.exists(compactf): 
        mat = ut.loadpy(compactf)
        inds = range(mat.shape[0]) # always square score matrix
        mat[inds, inds] = 0
        return mat
    else:
        ascores = np.loadtxt(scoref, dtype='f2')
        if save_compact:
            print 'saving compact', compactf
            ut.savepy(ascores, compactf)
        return ascores


class CosineLazyNew(object):

    def __init__(self,elution):
        self.norm_mat = np.mat(el.normalize_fracs(elution.mat))
        
    def __getitem__(self, index):
        # Dot product of normed rows
        return coscore(self.norm_mat, index[0], index[1])

def coscore(mat, i, j):
    #return 1 - dotrows(mat,i,j)/(dotrows(mat,i,i)**.5 * dotrows(mat,j,j)**.5)
    return dotrows(mat,i,j)/(dotrows(mat,i,i)**.5 * dotrows(mat,j,j)**.5)

def dotrows(mat, i, j):
    return np.asarray(mat[i,:]*mat[j,:].T)[0][0]


class CosineLazyScores(object):

    def __init__(self,elution):
        mat = elution.mat
        norms = np.apply_along_axis(np.linalg.norm, 1, mat)
        self.mat_rownormed = np.nan_to_num(mat / np.matrix(norms).T)
        assert type(self.mat_rownormed) == type(np.matrix(''))
        self.shape = (mat.shape[0],mat.shape[0])
        
    def __getitem__(self, index):
        # Dot product of normed rows
        return float(self.mat_rownormed[index[0],:] *
                    self.mat_rownormed[index[1],:].T)

def matching_pairs(values, ids):
    """
    Return all pairs of ids for indices in the given list whose values match.
    Will not return identity matches since uses combinations.
    """
    d = defaultdict(list)
    for ind,val in enumerate(values):
        d[val].append(ids[ind])
    return [(i,j) for value in d for i,j in itertools.combinations(d[value],2)]
    
def pairs_exceeding(elut, skey, thresh):
    """
    Doesn't return self-self interactions.
    """
    arr_prots = np.array(elut.prots)
    if skey == 'apex':
        apexes = ApexScores(elut).apex_array
        pairs = matching_pairs(apexes, arr_prots)
    else: # loading precomputed indices is so far massively slower than this
        score_mat, _, new_prots = scorekey_elution(skey, elut, None)
        if new_prots is not None:
            arr_prots = np.array(new_prots)
        rows, cols = np.where(score_mat > thresh)
        p1s, p2s = [arr_prots[ids] for ids in rows, cols]
        pairs =  ut.zip_exact(p1s, p2s)
    return pairs

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 3:
        sys.exit("usage: python score.py filename method(poisson|dotproduct|corrcoef|cov) [argument]") 
    fname = sys.argv[1]
    method = sys.argv[2]
    methodarg = None if nargs < 4 else int(sys.argv[3])
    elut = el.load_elution(fname)
    if method == 'poisson':
        corr = traver_corr(elut.mat, repeat=methodarg) if methodarg else \
            traver_corr(elut.mat)
    elif method in ['cosine_poisson','euclidean_poisson']:
        corr = poisson_repeat(elut.mat, metric=method.split('_')[0],
                repeat=methodarg) if methodarg else poisson_repeat(elut.mat,
                        metric=method)
    elif method in ['euclidean']:
        corr = pdist_score(elut.mat, norm_rows=True, norm_cols=True,
                metric=method)
    #elif method == 'dotproduct':
        #corr = elut.mat * elut.mat.T
    #elif method == 'corrcoef':
        #corr = np.corrcoef(elut.mat)
    #elif method == 'cov':
        #corr = np.cov(elut.mat)
    fileout = fname+'.corr_'+method
    np.savetxt(fileout, corr, delimiter='\t')

