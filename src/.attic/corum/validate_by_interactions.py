
import numpy as np
import csv
import MySQLdb
import argparse
import itertools as it
import random

def main():

	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--bicluster_file", action="store", dest="bicluster_filename", required=True, help="")
	parser.add_argument("--corum_file", action="store", dest="corum_file", required=True, help="")
	parser.add_argument("--threshold", action="store", dest="threshold", required=False, help="", default=0.0)
	parser.add_argument("--network_key", action="store", dest="network_key", required=False, help="", default=None)

	args = parser.parse_args()

	print args

	#kdrew: using genemania's physical networks downloaded from http://pages.genemania.org/data/
	db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'genemania_20131024')


	biclusters = dict()
	f = open(args.bicluster_filename,'rb')
	#kdrew: read in biclusters generated from corum_bicluster.py 
	for line in f.readlines():
		#print line
		complex_key = int(line.split()[1])
		zscore = float(line.split()[3])
		gene_str = line.split('genes:')[1]
		genes = gene_str.split()
		try:
			bc_zscore = biclusters[tuple(genes)].zscore
			#kdrew: only store top zscore cluster (lower zscore is better)
			if bc_zscore > zscore:
				biclusters[tuple(genes)] = Cluster(complex_key, zscore, genes)
		except KeyError:
			biclusters[tuple(genes)] = Cluster(complex_key, zscore, genes)


	complexes_dict = dict()
	#kdrew: read in corum complexes to get complete set of protein ids
	with open(args.corum_file, 'rb') as csvfile:
		complexes = csv.reader(csvfile, delimiter=';')
		for row in complexes:
			if row[3] == "Human" and row[4] != '':
				ids = row[4].split(',')
				#kdrew: remove ( ) around extra ids
				#print row
				complexes_dict[int(row[0])] = [x.translate(None, '()') for x in ids]


	#kdrew: create object that stores physical interactions from database
	iset = InteractionSet(db, args.network_key)


	#for complex1 in complexes_dict:
	#	print complex1
	#	for protID in complexes_dict[complex1]:
	#		for protID2 in complexes_dict[complex1]:
	#			print protID, protID2
	#			iset.get_interactions(protID, protID2)
	#			#for i in iset.get_interactions(protID, protID2):
	#			#	print i


	#kdrew: for every bicluster test to see if there are more connections in physical networks than random set from full complex
	for gene_ids in biclusters:
		print gene_ids
		bc = biclusters[gene_ids]
		print bc
		complex_gene_list = complexes_dict[biclusters[gene_ids].complex_key]
		print complex_gene_list
		iset.get_interactions_from_db(complex_gene_list)

		#for protID, protID2 in it.combinations(complexes_dict[biclusters[gene_ids].complex_key], 2):
		#	iset.get_interactions(protID, protID2)

		#kdrew: weight is defined from genemania
		bc_total_weight = iset.get_total_weight(gene_ids, float(args.threshold))
		print "total interaction weight: %s" % (bc_total_weight,)


		#kdrew: pick 1000 random protein sets from parent complex and calculate their physical interaction weight
		bc_count = len(gene_ids)
		random_total_weight_list = []
		for i in xrange(1000):
			random_genes = [ complex_gene_list[i] for i in sorted(random.sample(xrange(len(complex_gene_list)), bc_count)) ]
			random_total_weight = iset.get_total_weight(random_genes, float(args.threshold))
			random_total_weight_list.append(random_total_weight)

		
		rmean = np.array(random_total_weight_list).mean()
		rstd = np.array(random_total_weight_list).std()
		interaction_zscore = (bc_total_weight - rmean)/rstd

		print "rmean: %s rstd: %s" % (rmean, rstd,)
		print "interaction zscore: %s" % (interaction_zscore,)
		print "bc: %s, interaction_zscore: %s" % (bc, interaction_zscore,)
		print ""




#kdrew: class to store bicluster results 
class Cluster(object):
	def __init__(self, complex_key, zscore, gene_list):
		self.complex_key = complex_key
		self.zscore = zscore
		self.gene_list = gene_list
	def __str__(self):
		return "%s %s %s" % (self.complex_key, self.zscore, self.gene_list)

#kdrew: class to store interactions from database
class Interaction(object):
	def __init__(self, gene1, gene2, weight, network_key):
		self.gene1 = gene1
		self.gene2 = gene2
		self.weight = weight
		self.network_key = network_key

	def __str__(self):
		return "%s %s : %s, %s" % (self.gene1, self.gene2, self.weight, self.network_key)

	#kdrew: need for use in set()
	def __eq__(self, other):
		if not isinstance(other, type(self)): 
			return NotImplemented

		case1 = (self.gene1 == other.gene1 and self.gene2 == other.gene2 and self.weight == other.weight and self.network_key == other.network_key)
		case2 = (self.gene1 == other.gene2 and self.gene2 == other.gene1 and self.weight == other.weight and self.network_key == other.network_key)
		return case1 or case2

	#kdrew: need for use in set()
	def __hash__(self):
		return hash((self.gene1, self.gene2, self.weight, self.network_key))


#kdrew: stores set of interactions and caches so not repeatidly hammering database
class InteractionSet(object):
	def __init__(self, db, network_key=None):
		self.db = db
		self.interactions = dict()
		self.network_key = network_key

	#kdrew: gets interactions between two genes from cache if available, from db if not
	#kdrew: returns all found interactions
	def get_interactions(self, gene1, gene2):

		try:
			return self.interactions[(gene1, gene2)]	
		except KeyError:
			try:
				return self.interactions[(gene2, gene1)]
			except KeyError:
				print "non cache hit: %s %s" % (gene1, gene2,)
				self.get_interactions_from_db( [gene1, gene2] )
				return self.get_interactions(gene1, gene2)

	#kdrew: queries db for all interactions between given genes
	#kdrew: stores interactions and returns all interactions found for list of genes
	def get_interactions_from_db( self, genes ):                                          

		cursor = self.db.cursor()
		sql_smt = "select im.name, im2.name, i.* from identifier_mappings as im, identifier_mappings as im2, interactions as i where ((i.gene_A = im.preferred_name and i.gene_B = im2.preferred_name) or (i.gene_B = im.preferred_name and i.gene_A = im2.preferred_name)) and im.name in (%s) and im2.name in (%s)" % (','.join(['%s']*len(genes)), ','.join(['%s']*len(genes)))
		if self.network_key != None:
			sql_smt = sql_smt + " and i.network_key = %s" % (self.network_key)
		#print sql_smt
		cursor.execute(sql_smt, (genes+genes))
		db_interactions = cursor.fetchall()
		#print interactions
		interaction_set = set()

		for gene_pair in it.combinations(genes,2):
			if gene_pair in self.interactions.keys():
				continue
			elif (gene_pair[1], gene_pair[0]) in self.interactions.keys():
				continue
			else:
				self.interactions[gene_pair] = set()

		for i in db_interactions:
			#kdrew: store gene names, weight and network key
			iobj = Interaction(i[0], i[1], i[4], i[5])
			interaction_set.add(iobj)

			#kdrew: store in cache
			try:
				self.interactions[i[0], i[1]].add(iobj)
			except KeyError:
				self.interactions[i[0], i[1]] = set([iobj,])

		#print self.interactions
		return list(interaction_set)

	def get_total_weight( self, gene_list, threshold=0.0 ):

		total_weight = 0.0
		for gene_pair in it.combinations(gene_list, 2):
			for interaction in self.get_interactions(gene_pair[0], gene_pair[1]):
				if interaction.weight > threshold:
					total_weight += interaction.weight

		return total_weight
			


if __name__ == "__main__":
	main()


