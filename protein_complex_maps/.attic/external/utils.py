from __future__ import division
import cPickle
import datetime
import errno
import itertools
import numpy as np
import os
import operator
import re
import random
import scipy
import shutil
import string
import subprocess
import sys
import warnings
from scipy.misc import comb
from Struct import Struct



######################################################################
## PYTHON OBJECT SAVE/LOAD
######################################################################

def savepy(obj,fname, check_exists=False):
    # If safe is true, we save obj to a temporary file first, then mv the file
    # to its final destination
    #if safe:
        #save(obj, fname+'.partial', safe=False)
        #shutil.move(fname+'.partial', fname)
        #return
    ## Because reloading a module thwarts pickling 
    ## (the new class is not the same as the old)
    #def maybe_reload(obj):
        #if hasattr(obj, '_pickle_reload'):
            #obj = obj._pickle_reload()
        #return obj
    #obj = maybe_reload(obj)
    if check_exists and os.path.exists(fname):
        print "File already exists. Not saved.", fname
    else:
        cPickle.dump(obj, file(fname, 'wb'), protocol=2)

def loadpy(fname):
    obj = cPickle.load(file(fname, 'rb'))
    # if isinstance(obj, bayc.World):
    #     obj.update_version()
    if hasattr(obj, '_pickle_update'):
        # This function is called on loading, and should return the new
        # object, or self if the object has updated itself.
        obj = obj._pickle_update()
    try:
        obj.filename = fname
    except AttributeError:
        pass # Some objects don't have a filename
    return obj

def loadlab(fname, loadfunc=loadpy, copy=True):
    remotef = remoted(fname)
    if os.path.isfile(remotef):
        if copy:
            if (not os.path.isfile(fname) 
                or confirm(prompt="Local file exists.  Overwrite? [y/n]: ")):
                if os.path.isfile(remotef):
                    shutil.copy(remotef, fname)
                    print "Copied locally."
                return loadfunc(fname)
            else:
                print "Not copied."
                return loadfunc(remotef)
        else:
            return loadfunc(remotef)
    elif confirm(prompt="Remote mount not found. Use scp? [y/n]: "):
        source = os.path.join("Dropbox",exp_dirs(),fname)
        run_command("scp libra:%s ." % source)
        print "Copied locally."
        return loadfunc(fname)
    else:
        print "Not found:", remotef

def remoted(fname=''):
    return os.path.join('/Volumes/users/Dropbox/', exp_dirs(), fname)

def exp_dirs():
    return '/'.join(os.path.abspath('.').split('/')[-3:])

def bigd(fname=''):
    bigbase = '/data/' if os.uname()[1]=='libra' else os.path.expanduser('~/bigdata')
    bigpath = os.path.join(bigbase,exp_dirs())
    if not os.path.exists(bigpath):
        print "Dir did not previously exist; ran mkdir" , bigpath
        os.mkdir(bigpath)
    return os.path.join(bigpath,fname)

########################################################################
## COLLECTIONS and math functions
########################################################################

#t_array = type(array([1])) # because type(array([1])) != array

def all_same(f, bag):
    v = f(bag[0])
    for x in bag[1:]:
        if f(x) != v: return False
    return True

def arr_add_feats(arr, names, newtype=None):
    olddtype = arr.dtype.descr[-1][1]
    newtype = newtype or olddtype
    newdescr = arr.dtype.descr + [(n,newtype) for n in names]
    newarr = np.empty(arr.shape, dtype=newdescr)
    for n in arr.dtype.names:
        newarr[n] = arr[n]
    return newarr

def arr_assign_row_values(arr, rowind, colinds, values):
    row = arr[rowind]
    for c,v in zip_exact(colinds,values):
        row[c] = v

def arr_copy(arr):
    newarr = np.empty(arr.shape, dtype=arr.dtype.descr)
    for n in arr.dtype.names:
        newarr[n] = arr[n]
    return newarr

def bin(list, binsize):
    nbins = int(np.ceil(len(list)/binsize))
    return [list[i*binsize:(i+1)*binsize] for i in range(nbins)]

def column_totals(mat):
    return np.array(mat.sum(axis=0)).flatten()

def column_uniques(mat):
    return np.array((mat > 0).sum(axis=0)).flatten()

def count_collect(items, length):
    d = {}
    for i in items:
        init = i[:length]
        d[init] = d.get(init,0) + 1
    return ut.fsort(ut.fsort(d.items(), key=lambda x: x[0]),
                    key=lambda x: len(x[0]), reverse=True)

def ctype(items, newtype):
    return [newtype(i) for i in items]

def every(pred, bag):
    # Like the Common Lisp EVERY
    for x in bag:
        if not pred(x): return False
    return True

def fremove(el, l): # functional remove
    l2 = l[:]
    l2.remove(el)
    return l2

def fsort(l, *args, **kwargs):
    assert(type(l) == list)
    l2 = l[:]
    l2.sort(*args, **kwargs)
    return l2

def hyperg(k, N, m, n):
    exact=1
    return comb(m,k,exact) * comb(N-m,n-k,exact)/ comb(N,n,exact)

def hyperg_cum2(c, N, m, n):
    # Taking 1-sum_from_0_to_c is much, much faster than summing from c
    # to min(m,n) as is done in the paper
    # FIXME: We lose all of our precision substracting from 1
    return 1-sum([hyperg(k,N,m,n) for k in range(0, c)])

def hyperg_cum(c, N, m, n):
    if c == 0:
        return 1
    else:
        s = sum([hyperg(k,int(N),int(m),int(n)) for k in range(c, min(m,n)+1)])
        # Sanity check. Sadly, it has failed before. Calling
        # hyperg(2,5022,22,17) yields a negative number when the 22 is a
        # numpy.int64, but only if comb() is exact (see help(comb))
        assert(s>0)
        return s

def i0(seq):
    return inds(seq,0)

def i1(seq):
    return inds(seq,1)

def inds(seq,ind):
    return map(lambda x: x[ind], seq)

def inset(items, set):
    return [i for i in items if i in set]

def idpr(x): #useful sometimes for printing nested in expressions
    print x
    return x

def list_filter_value(l, index, val):
    return [x for x in l if x[index]==val]

def list_frac(l, frac):
    return l[:int(frac * len(l))]

def list_inds(l, inds):
    return list(np.array(l)[inds])

def least(f, seq, comp=operator.lt):
    # Returns the x in seq for which f(x) is smallest
    l = None
    lval = None
    for el in seq:
        v = f(el)
        if (lval and comp(v, lval)) or (not lval):
            l = el
            lval = v
    return l

def most(f, seq): return least(lambda x: -f(x), seq, comp=operator.gt)

def only(l):
    assert(len(l)==1)
    return l[0]

def onlys(s): # for set
    return only(tuple(s))

def printnow(s):
    print s
    sys.stdout.flush()

def r_squared(a,b):
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(a,b)
    return r_value**2

def ravg(l):
    return rsum(l) / len(l)

def rescale(frac, mult):
    """
    frac: fraction positives
    mult: negative multiplier
    End cases 0 and 1 work here.
    Doing the math: wlog, p+n=1, so frac=p=1-n  (for p=%pos, n=%neg)
    """
    #return p/(p+(1-p)*n) # to show similarity to wrong version below
    return frac / (frac + (1-frac) * mult)

def rescale_dep(p,n):
    """
    This is WRONG. But was used in clusterings etc so I'll keep it for now.
    Rescale posterior probability p according to multiple of negatives n.
    Instead should just be p/n, which hardly needs a function. 
    Think about the edge case of a multiplier of 1--should just return p.
    Must have been thinking additively when I made this, but used
    multiplicatively.
    """
    print "ut.rescale: Deprecated. Wrong."
    return p/(1+(1-p)*n)

def regex_filter(seq, pattern):
    matches = [n for n in seq if len(re.findall(pattern, n))>0]
    return matches

def rsum(l): #reduce sum
    # This is exactly like the Python built-in sum, but scipy overloads it
    # with its own sum function.
    return reduce(operator.add, l)

# def rnd(a,b=None):
#     # excludes b. we return an int if a is an int, a float if it's a float
#     if isinstance(a,int):
#         if b: assert(isinstance(b,int))
#         if b is None:
#             b = a
#             a = 0
#         return random.randint(a, b-1)
#     else:
#         if b is None:
#             b = a
#             a = 0
#         return random.uniform(a,b)

def sample_wr(pop, k):
    n = len(pop)
    _random, _int = random.random, int
    return [pop[_int(_random() * n)] for i in itertools.repeat(None,k)]

def sqrt_shape(k):
    """
    Get a proper-sized (m,n) for m rows, n columns for k items
    """
    srt = np.sqrt(k)
    ans1 = (np.floor(srt), np.ceil(srt))
    if ans1[0]*ans1[1] >= k:
        return ans1
    else:
        return (np.ceil(srt),np.ceil(srt))

def struct_copy(s):
    newstruct = Struct()
    newstruct.__dict__ = s.__dict__.copy()
    return newstruct

def zip_exact(*seqs):
    # Like zip, but generates an error if the seqs are not all the same
    assert(all_same(len, seqs))
    return zip(*seqs)



#####################################################################
##  LOADING DATA
#####################################################################

def convert_list(loi, ctype=float):
    return [ctype(i) for i in loi]

def convert_lol(lol, ctype=float):
    return [convert_list(lst, ctype=ctype) for lst in lol]
    
def pre_ext(fname, preext):
    path, name = os.path.split(fname)
    basename,ext = os.path.splitext(name)
    return os.path.join(path, basename+str(preext)+ext)

def shortname(fname):
    return os.path.basename(fname).split('.')[0]

def load_dict_flat(fname):
    # insted use dict(load_tab_file(fname))
    assert 0==1
    pass

def load_tab_file(fname, sep='\t', use_special_clean=False, dtype=None,
        skip1=False):
    """ Returns a generator of a list of list
    (assumes each element in a line is tab-delimited)
    """
    def clean(x):
        return x.strip() # remove the newlines
    if use_special_clean:
        datagen = (tuple(l.split(sep)[:-1]+[l.split(sep)[-1].strip()]) for l in file(fname, 'r'))
    else:
        datagen = (tuple(clean(l).split(sep)) for l in file(fname, 'r'))
    if skip1:
        datagen.next()
    return datagen

def load_dict_sets(fname, lines='pairs', **kwargs):
    # Returns a key->set(values) mapping
    if lines=='pairs':
        return dict_sets_from_tuples(load_list_of_lists(fname, **kwargs))
    elif lines=='lists':
        return dict([(l[0],set(l[1:])) for l in load_list_of_lists(fname,
            **kwargs)])

def load_list(fname):
    return [row[0] for row in load_tab_file(fname)]

def load_array(fname):
    return np.loadtxt(fname)

def load_list_of_lists(fname, **kwargs):
    return load_list_of_type(fname, list, **kwargs)
load_lol = load_list_of_lists #shortcut

def load_los(fname, **kwargs):
    return load_list_of_type(fname, set, **kwargs)

def load_lot(fname, **kwargs):
    return load_list_of_type(fname, tuple, **kwargs)

def load_list_of_type(fname, argtype, **kwargs):
    return [argtype(l) for l in load_tab_file(fname, **kwargs)]

def write_dict_sets_lines(d,fname, sep='\t'):
    if sep=='\t':
        write_tab_file([[k]+list(v) for k,v in d.items()],fname)
    else:
        write_tab_file([[k] + [sep.join(v)] for k,v in d.items()],fname)
    
def write_dict_sets_tab(d,fname):
    mylist = []
    for k in d:
        mylist.extend([(k,v) for v in d[k]])
    write_tab_file(mylist,fname)

def write_tab_file(ll, fname, formatter='{0}', header=None, islist=False):
    # '{0:10g}' for 10 digit precision
    f = file(fname, 'w')
    if header:
        f.write('\t'.join(formatter.format(x) for x in header) + '\n')
    if islist:
        f.write('\n'.join((formatter.format(l) for l in ll)))
    else:
        for l in ll:
            f.write('\t'.join(formatter.format(x) for x in l) + '\n')

def print_lol(lol):
    for l in lol:
        print '\t'.join([str(i) for i in flatten(l)])

def print_list(l):
    print '\t'.join([str(i) for i in l])

def flatten(lst):
   out = []
   for sublist in lst:
       if not isinstance(sublist,list):
           out.append(sublist)
       else:
           out.extend(sublist)
   return out

def fremove(lst,item):
    newl = list(lst)
    newl.remove(item)
    return newl


#######################################
# DICT OPERATIONS
#######################################

def constant_factory(value):
    return itertools.repeat(value).next

def dict_remove_if_value(d,remove_val):
    emptys = []
    for key in d:
        if d[key] == remove_val:
            emptys.append(key)
    for k in emptys:
        del(d[k])
    return d

def compose_dict_sets(d1,d2,compose_op=set.union,default=set()):
    # d1 must be key1->set(keys2)
    # and d2 be key2->V
    # compose_op must accept two V
    # default is for d2, it should return a V
    return dict(((k,reduce(compose_op, [d2.get(k2,default) for k2 in v],
                           default)) for (k,v) in d1.items()))

def dict_combine_sets(d1,d2):
    # Assume that the values are sets, and combine them pairwise (union)
    return dict([(k, d1.get(k,set()).union(d2.get(k,set())))
                 for k in list(set(d1.keys()+d2.keys()))])

def dict_inverse_sets(d):
    # If d is a K->set(V) mapping, return a V->set(K) dict
    dout = {}
    for key, values in d.items():
        for v in values:
            dout.setdefault(v,set([])).add(key)
    return dout

def dict_inverse(d):
    dout = {}
    for k, v in d.items():
        dout[v] = k
    return dout

def list_inv_to_dict(lst):
    d = {}
    for index,item in enumerate(lst):
        d[item]=index
    return d

def dict_from_lol(lol):
    return dict([(x[0],x[1:]) for x in lol])

def list_trans(items, d):
    return [d[i] for i in items]

def lol_trans(lol, d):
    return [[d[i] for i in items] for items in lol]

def dict_set_defaults(d, defaultd):
    for k,v in defaultd.items():
        if k not in d:
            d[k] = v
    return d

def dict_sets_from_tuples(lot):
    d = {}
    for k,v in [tup for tup in lot if len(tup)>1]:
        d.setdefault(k,set([])).add(v)
    return d

def dict_count_values(d):
    return sum([len(d[k]) for k in d])

def dict_dedup(d):
    # in a dict of sets, for every k:v[i] found, remove v[i]:k
    d_dedup = d.copy()
    for k in d:
        for vi in d[k]:
            if vi in d_dedup:
                if k in d_dedup[vi]:
                    d_dedup[vi].remove(k)
    return d_dedup

def dict_quick_merge(d1,d2):
    """
    Simply combine the two dicts, assuming keys are disjoint.
    """
    return dict([(k,v) for k,v in d1.items()]+[(k,v) for k,v in d2.items()])

def dict_merge_values(d1, d2, default=0):
    """
    Return a new dict of {key: (value1, value2); ...}
    """
    d = {}
    for k,v in d1.items():
        d[k] = (v, default)
    for k,v in d2.items():
        d[k] = (d.get(k, (default,default))[0], v)
    return d

def dict_select(d, keys):
    """
    Return a dict only with keys in the list of keys.
    example:
    clust_kw = 'min_density min_size haircut penalty'.split()
    clust_kwargs = ut.dict_select(kwargs, clust_kw)
    """
    return dict([(k,d[k]) for k in keys if k in d])

def dict_bicliques(d, minlen=1):
    bis = []
    tried = set([])
    ks = [k for k,v in d.items() if len(v)>minlen]
    for k in ks:
        if k not in tried:
            right = d[k]
            invd = dict_inverse_sets(d)
            left = invd[list(right)[0]]
            tried = set.union(tried, left)
            if (all([right==d[l] for l in left]) and all([left==invd[r] for r
                in right])):
                bis.append((k, len(right), len(left), right, left))
    return bis

##########################################
# Arrays
##########################################

def arrnd_retype(arr, newtypes, newnames=None):
    newarr = np.zeros(shape=len(arr), dtype = ','.join(newtypes))
    newarr.dtype.names = tuple(newnames) if newnames else arr.dtype.names
    newarr[:] = arr[:]
    return newarr

def arr_retype(arr, newtype='f2'):
    newarr = np.zeros(arr.shape, dtype=newtype)
    newarr[:,:] = arr[:,:]
    return newarr

def arr_norm(arr, axis=0):
    """
    axis=0: normalize each column; 1: normalize each row
    """
    mat = np.asmatrix(arr)
    return np.asarray(np.nan_to_num(mat / np.sum(mat, axis)))

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
    if norm_cols:
        # Normalize columns first--seems correct for overall elution profile if
        # you're calculating correlation-type scores
        arr = arr_norm(arr, 0)
    if norm_rows:
        arr = arr_norm(arr, 1)
    return arr

##########################################
# Configuration Filenames/Paths/Etc
##########################################

dir_project = "~/Dropbox/complex"
base_path = os.path.expanduser(dir_project)

def proj_path(pathkey,basename=''):
    if basename:
        return os.path.join(base_path, config()[pathkey], basename)
    else:
        return os.path.join(base_path, config()[pathkey])

def config():
    """
    Usage: value = ut.config()['strkey']
    """
    def str_to_bool(string):
        return (True if string.lower() == 'true' else False if string.lower() ==
            'false' else string)
    conf_path = os.path.join(base_path, 'config.txt')
    dconf = dict([(l.split()[0],str_to_bool(l.split()[1])) 
        for l in load_list(conf_path) if len(l)>0 and l[0]!='#'])
    return dconf

##########################################
# File and I/O helper functions
##########################################

def confirm(prompt="Continue? y/n"):
    """
    Usage: if confirm(): take action
    """
    while True:
        ans = raw_input(prompt)
        if ans not in ['y','n']:
            print "answers: y, n"
            continue
        if ans == 'y':
            return True
        if ans == 'n':
            return False
        
def temp_placeholder(path):
    """
    Usage: if temp_placeholder(ftemp): take action; remove ftemp; else: skip.
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True

def run_command(cmd):
    print cmd
    subprocess.call(cmd, shell=True)

def date():
    """
    Returns eg "20130420"
    """
    return date_ymd()

def date_ymd():
    return datetime.date.today().strftime("%Y%m%d")

def date_md():
    return datetime.date.today().strftime("%m%d")

