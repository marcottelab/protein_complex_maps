
import numpy as np
from scipy.stats.stats import nanmean

#kdrew: functions from Peter Cock: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	try:
		diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	except KeyError:
		return np.nan
	return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two, no_hetero_atoms = True ) :
	"""Returns a matrix of C-alpha distances between two chains"""
	if no_hetero_atoms:
		#kdrew: create new chain_one and chain_two arrays with only residues not hetatoms
		chain_one_tmp = []
		chain_two_tmp = []
		for r1 in chain_one:
			#print r1.get_id()
			if r1.get_id()[0]== ' ':
				chain_one_tmp.append(r1)
		for r2 in chain_two:
			#print r2.get_id()
			if r2.get_id()[0] == ' ':
				chain_two_tmp.append(r2)
		chain1 = chain_one_tmp
		chain2 = chain_two_tmp
	else:
		chain1 = chain_one
		chain2 = chain_two


	answer = np.zeros((len(chain1), len(chain2)), np.float)
	for row, residue_one in enumerate(chain1) :
		for col, residue_two in enumerate(chain2) :
			answer[row, col] = calc_residue_dist(residue_one, residue_two)
	return answer

def min_dist(structure, chain1, chain2, structure2=None):
	return metric_dist( structure, chain1, chain2, structure2=structure2, function=np.nanmin)

def max_dist(structure, chain1, chain2, structure2=None):
	return metric_dist( structure, chain1, chain2, structure2=structure2, function=np.nanmax)

def mean_dist(structure, chain1, chain2, structure2=None):
	return metric_dist( structure, chain1, chain2, structure2=structure2, function=np.nanmean)

#kdrew: structure as parsed by PDBParser and chain letter codes
#kdrew: if chains are in different pdbs (3aja, 3j3b of ribosome) use structure2
def metric_dist(structure, chain1, chain2, structure2=None, function=np.nanmin):
	try:
		chain_one = structure[0][chain1]
	except KeyError:
		print "missing chain: %s" % (chain1,)
		return np.nan

	if structure2 == None:
		try:
			chain_two = structure[0][chain2]
		except KeyError:
			print "missing chain: %s" % (chain2,)
			return np.nan
	else:
		try:
			chain_two = structure2[0][chain2]
		except KeyError:
			print "missing chain: %s" % (chain2,)
			return np.nan

	dmat = calc_dist_matrix(chain_one, chain_two)
	return function(dmat)
	

class PISA_Interfaces(object):

	def __init__(self, pisa_filename):
		self.pisa_filename = pisa_filename
		self.interfaces_area = dict()

		#kdrew: get interface area (angstrom^2) from pisa file for every interacting chain
		pisa_file = open(pisa_filename, "rb")
		for line in pisa_file.readlines():
			line_split = line.split()
			chainA = line_split[0]
			chainB = line_split[1]
			#kdrew: keep compound key in order
			if chainA > chainB:
				chainA, chainB = chainB, chainA
			self.interfaces_area[(chainA, chainB)] = float(line_split[2])

	def surface_area(self, chainA, chainB):
		if chainA > chainB:
			chainA, chainB = chainB, chainA
		
		try:
			sa = self.interfaces_area[(chainA, chainB)]
		except KeyError:
			sa = None

		return sa

