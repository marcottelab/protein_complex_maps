
import logging
import tempfile
import numpy as np
import os.path
import itertools as it
import argparse
import pickle

from Bio import SeqIO
import Bio.Align.Applications as baa

import protein_complex_maps.protein_util as pu


logging.basicConfig(level = logging.INFO,format='%(asctime)s %(levelname)s %(message)s')

def main():

	parser = argparse.ArgumentParser(description="Build MSA ")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--output_msa_file", action="store", dest="output_msa_file", required=True, 
						help="Filename of final MSA")
	parser.add_argument("--orthologs", action="store_true", dest="orthologs", required=False, default=False,
						help=""" Build MSA based on orthologs of given proteins, default=False""")


	args = parser.parse_args()

	tmpdir = tempfile.mkdtemp()
	print "tmpdir: %s" % (tmpdir,)

	#kdrew: dictionary of dictionaries: input_prot_id -> [sequences]
	sequence_map = dict()

	if args.orthologs:

		orthologs_map = pu.get_all_orthologs(args.proteins)

		for key in orthologs_map:

			#sequences = pu.get_from_uniprot(orthologs_map[key], 'sequence')

			primary_sequence = pu.get_sequences_uniprot( [key,], seqrecord=True )
			sequences = pu.get_sequences_uniprot( orthologs_map[key], seqrecord=True )

			print primary_sequence.values()

			sequence_map[key] = primary_sequence.values() + sequences.values()


	#print sequence_map
	for key in sequence_map:
		#kdrew: output sequence_map to fasta in temp dir
		sequences_filename = "%s/%s.fasta" % (tmpdir, key)
		SeqIO.write(sequence_map[key], sequences_filename, "fasta")

		msa_filename = "%s/%s.msa.fasta" % (tmpdir, key)  
		#kdrew: build MSA from fasta
		msa_cline = baa.MuscleCommandline(input=sequences_filename, out=msa_filename)
		msa_cline()


	#kdrew: output MSA to specified outputfile


	

if __name__ == "__main__":
	main()


