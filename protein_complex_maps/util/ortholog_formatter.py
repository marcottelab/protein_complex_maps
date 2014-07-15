

import glob
import cPickle
import argparse

import protein_complex_maps.protein_util as pu

class OrthologEntry(object):
	def __init__(self, cluster_id, guid, species, score, protid):
		self.cluster_id = cluster_id
		self.guid = guid
		self.species = species
		self.score = score
		self.protid = protid

def main():

	parser = argparse.ArgumentParser(description="Tool to read in Mass Spec Data Set (MSDS) from elution files and pickle it")
	parser.add_argument("--filename", action="store", dest="filename", required=True,
						help="Filenames of Blake's inparanoid files")
	parser.add_argument("--output_filename", action="store", dest="out_filename", required=True, 
						help="Output filename ")
	parser.add_argument("--species1", action="store", dest="species1", required=True, 
						help="species listed in OrtoA column")
	parser.add_argument("--species2", action="store", dest="species2", required=True, 
						help="species listed in OrtoB column")
	parser.add_argument("--species1_from_id", action="store", dest="species1_from_id", required=False, default="", 
						help="species1 id type")
	parser.add_argument("--species2_from_id", action="store", dest="species2_from_id", required=False, default="", 
						help="species2 id type")
	parser.add_argument("--to_id", action="store", dest="to_id", required=False, default="ACC", 
				help="convert ids this type, use empty string for original ids (note: automatically converts gene_exact to ACC)")

	args = parser.parse_args()

	#kdrew: read in ms files
	#sample_filenames = "/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/CeDmHsMmSp_ms2_elutions/Hs_*"

	#kdrew: error checking
	if not args.filename:
		print "\nError: Specify --filename \n"
		parser.print_help()
		return

	outfile = open(args.out_filename, "wb")

	ortholog_entry_list = []
	prot_s1_list = []
	prot_s2_list = []
		
	file1 = open(args.filename, 'rb')
	for line in file1.readlines():
		col_split = line.split('\t')
		cluster_id = col_split[0]
		guid = col_split[1]
		species1_list = col_split[2].split()
		species2_list = col_split[3].split()

		for prot, score in zip(species1_list[0::2], species1_list[1::2]):
			#outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cluster_id, guid, args.species1, score, prot,))
			ortholog_entry_list.append(OrthologEntry(cluster_id, guid, args.species1, score, prot))
			prot_s1_list.append(prot)

		for prot, score in zip(species2_list[0::2], species2_list[1::2]):
			#outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cluster_id, guid, args.species2, score, prot,))
			ortholog_entry_list.append(OrthologEntry(cluster_id, guid, args.species2, score, prot))
			prot_s2_list.append(prot)


	if args.to_id != "":
		if args.species1_from_id == "gene_exact":
			s1_map = pu.get_from_uniprot_by_genename(prot_s1_list)
		else:
			s1_map = pu.map_protein_ids(prot_s1_list, from_id=args.species1_from_id, to_id=args.to_id)

		if args.species2_from_id == "gene_exact":
			s2_map = pu.get_from_uniprot_by_genename(prot_s2_list)
		else:
			s2_map = pu.map_protein_ids(prot_s2_list, from_id=args.species2_from_id, to_id=args.to_id)

		ortholog_map = dict(s1_map.items() + s2_map.items())

	for entry in ortholog_entry_list:
		if args.to_id == "":
			outfile.write("%s\t%s\t%s\t%s\t%s\n" % (entry.cluster_id, entry.guid, entry.species, entry.score, entry.protid,))
		else:
			try: 
				for o_id in ortholog_map[entry.protid]:
					outfile.write("%s\t%s\t%s\t%s\t%s\n" % (entry.cluster_id, entry.guid, entry.species, entry.score, o_id,))
			except KeyError:
				outfile.write("%s\t%s\t%s\t%s\t%s\n" % (entry.cluster_id, entry.guid, entry.species, entry.score, entry.protid,))

	file1.close()
	outfile.close()

if __name__ == "__main__":
	main()




