
import numpy as np
import pickle 
import argparse
import os.path


def main(): 

	parser = argparse.ArgumentParser(description="Tool to return populated fraction names and counts of specified proteins")
	parser.add_argument("--input_msds_pickle", action="store", nargs='+', dest="msds_filenames", required=True, 
						help="Filename(s) of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")
	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--threshold", action="store", dest="threshold", type=float, required=False, default=0.0, 
						help="Only report fractions with counts above thresold")

	args = parser.parse_args()


	for p in args.proteins:
		print p
		flist_full = []
		fcount_full = []
		for msds_filename in args.msds_filenames:
			msds = pickle.load( open( msds_filename, "rb" ) )
			try:
				flist, fcounts = fractions_with_counts(msds, p, threshold=args.threshold)
			except KeyError:
				continue
			for i, fname in enumerate(flist):
				print "%s %s %s" % (os.path.basename(msds_filename), fname, fcounts[i])
			flist_full += flist.tolist()
			fcount_full += fcounts.tolist()


		try:
			print "max: %s %s %s" % (p, flist_full[np.argmax(fcount_full)], fcount_full[np.argmax(fcount_full)])
		except ValueError:
			continue


		print "\n"

#kdrew: function to return only the fractions with counts above a specified threshold
def fractions_with_counts( msds, prot_id, threshold = 0.0 ):
	fraction_list = msds.get_fraction_list()

	pindex = msds.get_id_dict()[prot_id]
	fraction_bool_array = msds.get_data_matrix()[pindex] > threshold
	fraction_ids = np.array(np.where(fraction_bool_array)[1])[0]


	fraction_names = np.array(fraction_list)[fraction_ids]
	fraction_counts = np.array(msds.get_data_matrix()[pindex,fraction_ids])[0]

	return fraction_names, fraction_counts

if __name__ == "__main__":
	main()
