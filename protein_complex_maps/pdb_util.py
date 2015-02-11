
import os.path as op
import numpy as np
from scipy.stats.stats import nanmean

import Bio.PDB
import Bio.SeqIO.PdbIO as pdbio
import Bio.SeqUtils as su

#kdrew: add residues that are missing from pdb atom record (discontinous residue numbering), chain is chain identifier
def add_missing_residues(pdb_filename, chain_id, no_hetero_atoms = True):

	complete_chain = []
	
	#kdrew: parse filename to get pdbid
	pdbid = op.basename(pdb_filename).split('.')[0] 
	try:    
		structure = Bio.PDB.PDBParser().get_structure(pdbid, pdb_filename)
	except IOError as e:
		print "I/O error({0}): {1}".format(e.errno, e.strerror)
		return
	try:
		chain_one = structure[0][chain_id]
	except KeyError:
		print "missing chain: %s" % (chain_id,)
		return 

	f = open(pdb_filename, "rb")
	chain_seqrecord = None
	#kdrew: get sequence of chain's atom record with gaps filled
	for i in pdbio.PdbAtomIterator(f):
		if chain_id == i.id.split(':')[1]:
			chain_seqrecord = i
	#print chain_seqrecord.seq

	#kdrew: create list of residues from chain
	tmp_chain = [r1 for r1 in chain_one]
	
	#kdrew: compare chain residues with gap filled sequence and populate complete_chain list with residues
	curr_resseq = 0
	for j, aa in enumerate(chain_seqrecord):
		#print aa
		if aa == 'X':
			#print "in aa==X"
			curr_resseq=+1
			xres = Bio.PDB.Residue.Residue((' ',curr_resseq, ' '), 'X', ' ')
			complete_chain.append(xres)
		else:
			try:
				r1 = tmp_chain.pop(0)
				#print su.seq1(r1.resname)
				#kdrew: make sure sequences match at this point
				if aa == su.seq1(r1.resname):
					complete_chain.append(r1)
					curr_resseq = r1.get_id()[1]
				else:
					print "WARNING: PDB gap filled sequence does not match pdb atom record: %s  %s" % (aa, su.seq1(r1.resname))
			except IndexError:
				print "WARNING: Extra residues in gap filled sequence but not in pdb atom record: %s" % (aa,)

	if len(tmp_chain) > 0:
		print "WARNING: Extra residues in pdb atom record but not in gap filled sequence: %s" % (tmp_chain,)

	return complete_chain


#kdrew: functions from Peter Cock: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	try:
		diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	except KeyError:
		return np.nan
	return np.sqrt(np.sum(diff_vector * diff_vector))



def calc_dist_matrix(chain_one, chain_two, no_hetero_atoms = True) :
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

#kdrew: calculates distance matrix between all residues in both chains
def calc_dist_matrix_complete(chain_one, chain_two, no_hetero_atoms= True):

	#kdrew: not sure if the direction (topleft, botleft, etc) makes sense with how the hstack and vstack  are used
	pdb_dist_mat_topleft = calc_dist_matrix(chain_one, chain_one)
	print "topleft shape: %s" % (pdb_dist_mat_topleft.shape,)

	pdb_dist_mat_botleft = calc_dist_matrix(chain_one, chain_two)
	print "botleft shape: %s" % (pdb_dist_mat_botleft.shape,)

	pdb_dist_mat_topright = calc_dist_matrix(chain_two, chain_one)
	print "topright shape: %s" % (pdb_dist_mat_topright.shape,)

	pdb_dist_mat_botright = calc_dist_matrix(chain_two, chain_two)
	print "botright shape: %s" % (pdb_dist_mat_botright.shape,)


	pdb_dist_mat_left = np.hstack((pdb_dist_mat_topleft, pdb_dist_mat_botleft))
	print "left shape: %s" % (pdb_dist_mat_left.shape,)

	pdb_dist_mat_right = np.hstack((pdb_dist_mat_topright, pdb_dist_mat_botright))
	print "right shape: %s" % (pdb_dist_mat_right.shape,)

	pdb_dist_mat = np.vstack((pdb_dist_mat_left, pdb_dist_mat_right))
	print "pdb_dist_mat shape: %s" % (pdb_dist_mat.shape,)

	return pdb_dist_mat



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

