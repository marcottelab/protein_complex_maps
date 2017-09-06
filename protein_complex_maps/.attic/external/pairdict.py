from __future__ import division
import utils as ut
from numpy import ndarray

class PairDict(object):
    """
    A common storage format for examples and predictions that handles all the
    messiness of merging, deduping, etc.
    """

    def __init__(self, loitems):
        """
        Make this from a list or set of tuples of (id1, id2, val1, val2, ...)
        """
        self.d = {}
        pd_set_loi(self, loitems)
        #self.d = dict([((p[0],p[1]),list(p)[2:]) for p in lopair_vals])


    def set(self, key, val):
        k = self.find(key)
        if k==None:
            k = key
        self.d[k] = val

    def append(self, key, val):
        k = self.find(key)
        if k==None:
            self.d[key] = val
        else:
            if isinstance(self.d[k],set):
                self.d[k].add(*val)
            else:
                self.d[k].append(*val)

    def contains(self, pair):
        return (self.find(pair) is not None)

    def find(self, pair):
        if pair in self.d:
            return pair
        elif pd_flip(pair) in self.d:
            return pd_flip(pair)
        else:
            return None

    def get(self, pair):
        usepair = self.find(pair)
        if usepair:
            return self.d[usepair]

def pd_set_loi(pd, loitems):
    for r in loitems:
        pd.set((r[0],r[1]),list(r)[2:])

def pd_dedupe(pd):
    newpd = PairDict([])
    pd_set_loi(newpd, pd_lol(pd))
    return newpd


def pd_copy(pd):
    newpd = PairDict([])
    newpd.d = pd.d.copy()
    return newpd

def pd_flip(pair):
    return (pair[1],pair[0])

def pd_lol(pd):
    return [[k[0],k[1]] + list(pd.d[k]) for k in pd.d]

def pd_intersect_avals(pda, pdb):
    """
    Merge two PairDicts and return a new one.
    Values in the intersection will be taken from the first argument, pda.
    """
    newpd = PairDict([])
    for pair,val in pda.d.items():
        if pdb.contains(pair):
            newpd.set(pair,val)
    return newpd

def pd_union_novals(a,b):
    """
    Merge two PairDicts and return a new one.
    Any existing values in a will be smashed by b's corresponding values.
    """
    newpd = PairDict([])
    newpd.d = a.d.copy()
    for k,v in b.d.items():
        newpd.set(k,v)
    return newpd

def pd_combine_ppis_matched(a,b,comb_func):
    """
    comb_func: a function used to combine the values, eg [lambda
    x:x[0]] to always just take the first value.
    REQUIRES: a,b contain same set of pairs.  Don't want to deal with it now.
    """
    newpd = PairDict([])
    newpd.d = a.d.copy()
    for p in a.d:
        bp = b.find(p)
        newpd.d[p][0] = comb_func(a.d[p][0],b.d[bp][0])
    return newpd

def pd_combine_ppis(a,b,comb_func):
    """
    comb_func: a function used to combine the values, eg [lambda
    x,y:x] to always just take the first value.
    """
    #a,b = [PairDict(pd_lol(pdx)[:3]) for pdx in a,b] #Get rid of extra columns
    newpd = pd_union_disjoint_vals(a,b,adefaults=[0,-1],bdefaults=[0,-1])
    for pair,(ascore, atrue, bscore, btrue) in newpd.d.items():
        newpd.d[pair] = (comb_func(float(ascore),float(bscore)), max(int(atrue), int(btrue)))
    return newpd

def pd_union_disjoint_vals(a,b,adefaults=None,bdefaults=None):
    """
    Merge two PairDicts and return a new one.
    Values for each pair become a list of avalues+bvalues, with defaults for
    each index provided if desired.
    """
    adefaults = adefaults if adefaults else [None]*len(a.d.values()[0])
    bdefaults = bdefaults if bdefaults else [None]*len(b.d.values()[0])
    newpd = PairDict([])
    def merge_help(from_set,newpd,a,b,adefaults,bdefaults, reverse=False):
        bleftovers = set(b.d.keys())
        for apair in from_set:
            bpair = b.find(apair)
            if bpair:
                newvals = (list(a.d[apair]) + list(b.d[bpair]) if not reverse else
                            list(b.d[bpair]) + list(a.d[apair]))
                bleftovers.remove(bpair)
            else:
                newvals = (list(a.d[apair]) + list(bdefaults) if not reverse
                        else list(bdefaults) + list(a.d[apair]))
            newpd.set(apair,newvals)
        return bleftovers
    bleftovers = merge_help(a.d.keys(),newpd,a,b,adefaults,bdefaults)
    merge_help(bleftovers,newpd,b,a,bdefaults,adefaults,reverse=True)
    return newpd
