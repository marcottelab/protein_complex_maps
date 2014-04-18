

import numpy as np


#kdrew: function to combine peptide_dict and peptide_count_dict
#kdrew: essentially maps protein->fraction->peptide->count
def protein_counts_by_peptide(peptide_dict, peptide_count_dict):
	protein_counts = dict()

	for protein in peptide_dict:
		fraction_dict = dict()
		for peptide in peptide_dict[protein]:
			try:
				for fraction in peptide_count_dict[peptide]:
					count = peptide_count_dict[peptide][fraction]
					try:
						fraction_dict[fraction][peptide] = count
					except:
						fraction_dict[fraction] = dict()
					fraction_dict[fraction][peptide] = count
			except KeyError:
				#print "no counts found for peptide"
				continue

		protein_counts[protein] = fraction_dict

	return protein_counts



#kdrew: read spectral count data by peptide by fraction (example file: test/Hs_test.peplist)
#kdrew: file format is protein_id (ignored) fraction peptide count
def read_peptide_list(peptide_list_file):

	peptide_count_dict = dict()

	for line in peptide_list_file.readlines():
		fraction = line.split()[1]
		peptide = line.split()[2]
		count = float(line.split()[3])

		try:
			count_dict = peptide_count_dict[peptide]
			try:
				count_dict[fraction] = count_dict[fraction] + count
				print "duplicate peptide %s %s %s" % (fraction, peptide, count)
			except KeyError:
				count_dict[fraction] = count

			peptide_count_dict[peptide] = count_dict
		except KeyError:
			peptide_count_dict[peptide] = dict()
			peptide_count_dict[peptide][fraction] = count

	return peptide_count_dict



#kdrew: peptide dict is the file that maps proteins to peptides (example file: test/Hs_test.pepDict)
#kdrew: peptide_file is a file handle, id_list is a list of protein ids, ignore_nonunique flag will not store nonunique peptides found in file
def read_peptide_dict(peptide_file, id_list, ignore_nonunique=True):
	peptide_dict = dict()
	#kdrew: count all the peptides for all the ids in id_dict
	for line in peptide_file.readlines():
		#kdrew: ignore non-unique peptides
		if line.count('|') and ignore_nonunique:
			#print "ignore line"
			continue

		#print line
		#kdrew: count all lines that have protein id in them
		for key in id_list:
			if key in line:
				peptide = line.split()[0]
				try:
					peptide_dict[key].append(peptide)
				except KeyError:
					peptide_dict[key] = [peptide]

	return peptide_dict



#kdrew: reads file in the format of peptide_count protein_id 
#kdrew: used for normalizing, can also be used for length of protein for the same purpose
def read_peptide_counts(peptide_file):
	peptide_dict = dict()
	for line in peptide_file.readlines():
		key = line.split()[1] 
		count = int(line.split()[0])
		print key, count
		peptide_dict[key] = count

	return peptide_dict



