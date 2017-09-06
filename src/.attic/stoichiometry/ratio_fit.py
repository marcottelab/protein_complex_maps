
import itertools as it
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os.path
from sklearn import mixture
from scipy.stats import ks_2samp
import scipy
import rpy2.robjects as ro

from multiprocessing import Pool
import multiprocessing as mp
import argparse
import cPickle
import pickle
import MySQLdb
import emd

def main():
	parser = argparse.ArgumentParser(description="Fits GMM to ratios of proteomic profiles for given proteins ")
	parser.add_argument("--input_msds_pickle", action="store", dest="msds_filename", required=True, 
						help="Filename of MSDS pickle: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py")

	parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=False, 
						help="Protein ids in which to anaylze")
	parser.add_argument("--set2_proteins", action="store", dest="set2_proteins", nargs='+', required=False, default=None,
						help="Protein ids in which to anaylze, results in ratio of proteins / set2_proteins")

	parser.add_argument("--complexes", action="store", dest="complexes", nargs='+', required=False, 
						help="Complex names from database in which to analyze")
	parser.add_argument("--input_complex_list", action="store", dest="complex_filename", required=False, 
						help="Filename of complex protein identifiers")
	parser.add_argument("--plot_filename", action="store", dest="plot_filename", required=False, default=None,
						help="Filename of output plot")
	parser.add_argument("--total_occupancy", action="store_true", dest="total_occupancy", required=False, default=False,
						help="Flag to only plot columns where every gene has greater than zero counts")
	parser.add_argument("--log_ratio", action="store_true", dest="log_ratio", required=False, default=False,
						help="Flag to log transform ratios")
	#parser.add_argument("--histogram", action="store_true", dest="histogram", required=False, default=False,
	#					help="Plot ratios in a histogram")
	parser.add_argument("--data_points_threshold", action="store", dest="threshold", required=False, default=10,
						help="Pairs of proteins must have atleast specified number of data points, default=10")
	parser.add_argument("-j", action="store", dest="numOfProcs", required=False, default=1,
						help="Number of processors to use, default=1")
	parser.add_argument("--model_selection", action="store", dest="model_selection", required=False, default='BIC',
						help="Method for choosing GMM model: BIC, AIC or k(int); default=BIC")

	args = parser.parse_args()

	pool = Pool(processes=int(args.numOfProcs))
	#print "setup pool"

	#kdrew: pickle comes from running protein_complex_maps.util.read_ms_elutions_pickle_MSDS.py
	msds = pickle.load( open( args.msds_filename, "rb" ) )
	msds_id_dict = msds.get_id_dict()
	msds_id_set = set(msds_id_dict.keys())

	if args.complex_filename != None:

		complex_file = open(args.complex_filename, 'rb')
		complex_list = []
		for i in complex_file.readlines():
			complex_list.append(i.split())

		print complex_list

		for i, prot_ids in enumerate(complex_list):
			if msds_id_set >= set(prot_ids):                                                                              
				#plot_fname = os.path.splitext(args.plot_filename)[0]+'.'+str(i)+os.path.splitext(args.plot_filename)[1]
				for prot_id_pairs in it.combinations(prot_ids, 2):
					#print "sending job"
					if args.numOfProcs > 1:
						pool.apply_async(fit_ratios, (msds, prot_id_pairs), dict(ids_set2=args.set2_proteins, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, data_threshold=int(args.threshold), model_selection=args.model_selection))
					else:
						print prot_id_pairs
						fit_ratios(msds, prot_id_pairs, ids_set2=args.set2_proteins, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, data_threshold=int(args.threshold), model_selection=args.model_selection)

				#print "sent job"

	elif args.set2_proteins != None:
		fit_ratios(msds, args.proteins, ids_set2=args.set2_proteins, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, data_threshold=int(args.threshold), plot_filename=args.plot_filename, model_selection=args.model_selection, complexes_mode=True)
	elif args.proteins != None:
		fit_ratios(msds, args.proteins, ids_set2=None, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, data_threshold=int(args.threshold), plot_filename=args.plot_filename, model_selection=args.model_selection)
	elif args.complexes != None:
		fit_ratios(msds, args.complexes, ids_set2=None, total_occupancy=args.total_occupancy, log_ratio=args.log_ratio, data_threshold=int(args.threshold), complexes_mode=True, plot_filename=args.plot_filename, model_selection=args.model_selection)
		
	else:
		print "use either --proteins or --complex_filename or --complexes to specify protein ids"

	pool.close()
	pool.join()

def fit_ratios(msds, ids, ids_set2=None, total_occupancy=False, log_ratio=False, data_threshold=10, sql_store=False, complexes_mode=False, plot_filename=None, model_selection='BIC'):

	cursor = None
	if sql_store or complexes_mode:
		db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', 'stoichiometry')
		cursor = db.cursor()
		sql_insert = """INSERT INTO ratio (PROTEIN_KEY1, PROTEIN_KEY2, log_mean, log_median, count, converged, class) 
							VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')"""
							#VALUES ('%d', '%d', '%f', '%f', '%d', '%b', '%d')"""
		sql_complex_insert = """INSERT INTO complex_ratio (COMPLEX_KEY1, COMPLEX_KEY2, log_mean, log_median, count, converged, class) 
							VALUES ('%s', '%s', '%s', '%s', '%s', '%s', '%s')"""


	protein_ids = []

	#kdrew: complexes mode groups sets of proteins together
	#kdrew: for example we can define all of the proteins in the proteasome lid and beta ring and compare stoichiometries between them
	#kdrew: currently a little hacky because I have to have multiple if statements to check for complex mode
	#kdrew: it would be better to treat the single protein vs protein case as a special case of the complex case
	if complexes_mode:
		complex_protein_map = dict()
		if ids_set2 == None:
			for name in ids:
				try:
					cursor.execute("select distinct p.proteinid from complex as c, complex_association as ca, protein as p where p.id = ca.protein_key and ca.complex_key = c.id and c.name = '%s'" % (name) )
					complex_protein_ids = cursor.fetchall()
					print complex_protein_ids
					#kdrew: unpack sql output into list, flattens tuple of tuples
					complex_protein_ids = [ele for tupl in complex_protein_ids for ele in tupl]
					print complex_protein_ids
					complex_protein_map[name] = complex_protein_ids
					protein_ids = protein_ids + complex_protein_ids
				except MySQLdb.Error, e:
					print e
		else:
			complex_protein_map['set1'] = ids
			complex_protein_map['set2'] = ids_set2
			protein_ids = protein_ids + ids
			protein_ids = protein_ids + ids_set2


	else:
		protein_ids = ids

	print protein_ids
	data_set, new_id_map = msds.get_subdata_matrix(protein_ids) 


	#kdrew: only look at columns that have at least one gene with values
	#kdrew: is this worth it? just removes zero columns 
	#kdrew: originally it was suppose to be all genes have values but not implemented that way
	if total_occupancy:
		col_sum = np.sum(data_set,0)
		shared_cols =  np.array(np.nonzero(col_sum)[1]).reshape(-1).tolist()
		#print "shared_cols: %s" % (shared_cols,)
		data_set = data_set[np.ix_(range(0,data_set.shape[0]), shared_cols)]

	#print "printing dataset"
	#print "data_set: %s" % (data_set,)

	ratio_dict = dict()

	#kdrew: code duplication from below with minor tweaks to make work for complex vs complex comparisons
	if complexes_mode:
		for complex1 in complex_protein_map.keys():

			for complex2 in complex_protein_map.keys():
				if complex1 == complex2:
					continue

				ratio_dict[complex1, complex2] = []

				for i, data_r in enumerate(data_set):
					data_row = np.array(data_r.reshape(-1))[0]

					#kdrew: for every other protein
					for j, data_r2 in enumerate(data_set):
						if i == j:
							continue
						#kdrew: make sure proteinid of current indices are members of complex1 and complex2 respectively
						if new_id_map[i] not in complex_protein_map[complex1] or new_id_map[j] not in complex_protein_map[complex2]:
							continue

						data_row2 = np.array(data_r2.reshape(-1))[0]
						#kdrew: for every fraction in data_matrix
						for k in range(len(data_row2)):
							#kdrew: make sure both proteins were seen at this fraction, 
							#kdrew: this might be a problem where floats are sometimes not exactly 0.0
							if data_row[k] != 0.0 and data_row2[k] != 0.0:
								#print "%s : %s : %s" % (data_row[k], data_row2[k], data_row[k]/data_row2[k])
								if log_ratio:
									#kdrew: add ratio to list
									ratio_dict[complex1,complex2].append(np.log(data_row[k]/data_row2[k]))
								else:
									ratio_dict[complex1,complex2].append(data_row[k]/data_row2[k])

	else:
		#kdrew: calculate all (protein) v all ratios of peak intensities
		#kdrew: for every protein
		for i, data_r in enumerate(data_set):
			data_row = np.array(data_r.reshape(-1))[0]

			#kdrew: for every other protein
			for j, data_r2 in enumerate(data_set):
				if i == j:
					continue

				ratio_dict[i,j] = []
				data_row2 = np.array(data_r2.reshape(-1))[0]
				#kdrew: for every fraction in data_matrix
				for k in range(len(data_row2)):
					#kdrew: make sure both proteins were seen at this fraction, 
					#kdrew: this might be a problem where floats are sometimes not exactly 0.0
					if data_row[k] != 0.0 and data_row2[k] != 0.0:
						#print "%s : %s : %s" % (data_row[k], data_row2[k], data_row[k]/data_row2[k])
						if log_ratio:
							#kdrew: add ratio to list
							ratio_dict[i,j].append(np.log(data_row[k]/data_row2[k]))
						else:
							ratio_dict[i,j].append(data_row[k]/data_row2[k])


	#f = open("ratio_dict.out","wb")
	#for r in ratio_dict:
	#	f.write(str(r)+"\n")
	#	for data in ratio_dict[r]:
	#		f.write(str(data)+"\n")
	#	f.write("\n")
	#f.close()

	#print ratio_dict

	#f, data_subplots = plt.subplots(len(ratio_dict.keys())+1,1,sharex='col')


	#pair = ratio_dict.keys()[0]

	#kdrew: for every protein pair
	for i, pair in enumerate(ratio_dict):
		f, data_subplots = plt.subplots(3,1)

		if complexes_mode:
			name0 = pair[0]
			name1 = pair[1]
		else:
			name0 = new_id_map[pair[0]]
			name1 = new_id_map[pair[1]]



		if len(ratio_dict[pair]) == 0:
			print "WARNING: ratio_dict array has zero length for %s : %s" % (name0, name1,)
			return


		maximum = np.max(ratio_dict[pair])
		minimum = np.min(ratio_dict[pair])
		xs = np.linspace(minimum, maximum, 100)

		if complexes_mode:
			print pair
		else:
			print new_id_map[pair[0]], new_id_map[pair[1]]

		#kdrew: check for nans and number fractions between pair is larger than threshold
		if len(ratio_dict[pair]) < data_threshold or np.isnan(np.sum(ratio_dict[pair])):
			print "WARNING: did not pass data_threshold or NANs in ratio calculation"
			return

		#print ratio_dict[pair]
		#print new_id_map[pair[0]], new_id_map[pair[1]]

		#kdrew: set up mixture model functionality
		color_iter = it.cycle(['r', 'g', 'b', 'c', 'm'])
		#clf = mixture.DPGMM(n_components=5, cvtype='diag')


		X = np.array([[x,] for x in ratio_dict[pair]])

		N = np.arange(1, 11)
		models = [None for ii in range(len(N))]

		for ii in range(len(N)):
			models[ii] = mixture.GMM(N[ii], n_iter=10000).fit(X)

		# compute the AIC and the BIC
		AIC = [m.aic(X) for m in models]
		BIC = [m.bic(X) for m in models]

		print "ARGMIN %s" % (np.argmin(BIC),)
		if model_selection == 'BIC':
			model_id = np.argmin(BIC)
		elif model_selection == 'AIC':
			model_id = np.argmin(AIC)
		else:
			try:
				model_id = int(model_selection)-1
			except ValueError, e:
				print "ERROR: problem with model selection input, must be int"
				print e
				return

		clf = models[model_id]

		#clf = mixture.GMM(n_components=3, covariance_type='diag', n_iter=10000)
		#clf.fit(X, n_iter=10000)
		#clf.fit(X)


		#print "weights:"
		#print clf.weights
		print "\n"
		#print clf.converged_
		#if clf.converged_:
		#	print clf

		Y = clf.predict(X)
		#print Y

		#kdrew: add proteins to sql
		protein_key1 = None
		protein_key2 = None
		if sql_store and cursor != None and not complexes_mode:
			print "insert proteins"
			try:
				stmt = "insert ignore into protein (proteinid) values ('%s')" % (new_id_map[pair[0]], )
				cursor.execute(stmt)
				db.commit()
				cursor.execute("select id from protein where proteinid = '%s'" % (new_id_map[pair[0]]) )
				protein_key1 = cursor.fetchone()[0]
				print "protein_key1: %s" % (protein_key1,)
			except:
				print "except protein insert: %s" % (stmt, )
				db.rollback()

			try:
				stmt = "insert ignore into protein (proteinid) values ('%s')" % (new_id_map[pair[1]], )
				cursor.execute(stmt)
				db.commit()
				cursor.execute("select id from protein where proteinid = '%s'" % (new_id_map[pair[1]]) )
				protein_key2 = cursor.fetchone()[0]
				print "protein_key2: %s" % (protein_key2,)
			except:
				db.rollback()


		ratio_dict_pair_array = np.array(ratio_dict[pair])

		#kdrew: call R to do dip test, returns pvalue of multimodal distribution
		ro.r('library(diptest)')
		rdp = np.array(ratio_dict[pair])
		rdp_vec = ro.FloatVector(rdp.transpose())
		diptest = ro.r['dip']
		dippval = diptest(rdp_vec)[0]
		#print "dip test pvalue: %s" % (dippval,)


		#kdrew:  print stats between classes found by GMM
		for j, mean in enumerate(clf.means_):
			class_data = np.array(ratio_dict[pair])[Y == j]

			for k, mean2 in enumerate(clf.means_):
				if j == k:
					continue

				class_data2 = np.array(ratio_dict[pair])[Y == k]

				if len(class_data) < data_threshold or len(class_data2) < data_threshold:
					continue

				#print "class %s (mean: %s) class2 %s (mean: %s) ks: %s" % (j, mean[0], k, mean2[0], ks_2samp(class_data, class_data2), )
				if len(class_data) > data_threshold and len(class_data2) > data_threshold:
					print "class %s (mean: %s) class2 %s (mean: %s) kl: %s, dip_pval: %s" % (j, mean[0], k, mean2[0], KL(class_data, class_data2), dippval )

		#kdrew: print stats of individual classes found by GMM, store in sql if flag 
		#kdrew: also plotting histogram (duplicated code from plot_ratios.py but this works for complexes mode
		for j, (mean, color) in enumerate(zip(clf.means_, color_iter)):

			print "clf.mean: %s clf.covar: %s" % (mean, clf._get_covars()[j][0],)
			class_data = np.array(ratio_dict[pair])[Y == j]

			if len(class_data) > data_threshold:

				if complexes_mode:
					print "%s : %s, converged? %s, mean: %s (%s), median: %s (%s), std: %s (%s), counts: %s" % (pair[0], pair[1], clf.converged_, mean, np.exp(mean), np.median(class_data), np.exp(np.median(class_data)), np.std(class_data), np.std(np.exp(class_data)), len(class_data))

				else:
					print "%s : %s, converged? %s, mean: %s (%s), median: %s (%s), counts: %s" % (new_id_map[pair[0]], new_id_map[pair[1]], clf.converged_, mean, np.exp(mean), np.median(class_data), np.exp(np.median(class_data)), len(class_data))

				if sql_store and cursor != None:
					#print "insert ratios"
					try:
						if complexes_mode:
							#print sql_complex_insert
							sql_insert_loaded = sql_complex_insert % (pair[0], pair[1], mean[0], np.median(class_data), len(class_data), 1 if clf.converged_ else 0, j)
						else:
							#print sql_insert
							sql_insert_loaded = sql_insert % (protein_key1, protein_key2, mean[0], np.median(class_data), len(class_data), 1 if clf.converged_ else 0, j)

						print sql_insert_loaded
						#kdrew: comment out for now
						#cursor.execute(sql_insert_loaded)
						#db.commit()
					except MySQLdb.Error, e:
						db.rollback()
						print e

		data_subplots[1].hist(ratio_dict[pair], bins=30)
				#data_subplots[i].plot(xs, mlab.normpdf(xs,mean,np.sqrt(clf._get_covars()[j][0])), color=color)

		print '\n'
														

		xs_2d = np.array([(x,) for x in list(xs)])

		#logprob, responsibilities = clf.eval(xs_2d)
		logprob, responsibilities = clf.score_samples(xs_2d)
		pdf = np.exp(logprob)
		#print pdf
		#for r in responsibilities:
		#	print r
		pdf_individual = responsibilities * pdf[:, np.newaxis]

		density = scipy.stats.gaussian_kde(ratio_dict[pair])
		#ax3 = data_subplots[i].twinx()
		data_subplots[0].plot(xs, density(xs))

		#ax2 = data_subplots[0].twinx()
		#ax2.plot(xs, pdf, '-k')
		data_subplots[0].plot(xs, pdf_individual, '--k')


		data_subplots[0].set_title("%s : %s, converged? %s" % (name0, name1, clf.converged_, ))

		data_subplots[2].plot(N, AIC, '-k', label='AIC')
		data_subplots[2].plot(N, BIC, '--k', label='BIC')
		data_subplots[2].set_xlabel('n. components')
		data_subplots[2].set_ylabel('information criterion')
		data_subplots[2].legend(loc=2)


		if plot_filename != None:

			plot_filename_base = plot_filename.split('.')[:-1]
			plot_filename_ext = plot_filename.split('.')[-1]
			plot_filename2 = "%s.%s_%s.%s" % ('.'.join(plot_filename_base), name0, name1, plot_filename_ext)
			plt.savefig(plot_filename2)
			plt.close()

		if sql_store or complexes_mode:
			db.close()



		print "Full set:  %s : %s, mean: %s (%s), median: %s (%s), std: %s (%s), counts: %s, dip test pvalue: %s, #ofGaussians: %s" % (name0, name1, ratio_dict_pair_array.mean(), np.exp(ratio_dict_pair_array.mean()), np.median(ratio_dict_pair_array), np.exp(np.median(ratio_dict_pair_array)), np.std(ratio_dict_pair_array), np.std(np.exp(ratio_dict_pair_array)), len(ratio_dict_pair_array), dippval, model_id+1)

def KL( a_in, b_in ):
	TINY_NUM  = 0.00001

	minimum = np.min(a_in.tolist() + b_in.tolist())
	maximum = np.max(a_in.tolist() + b_in.tolist())

	bins = scipy.linspace(minimum, maximum, 100)	

	a_hist = scipy.histogram(a_in, bins=bins, density=True)
	b_hist = scipy.histogram(b_in, bins=bins, density=True)

	a = a_hist[0] + TINY_NUM
	b = b_hist[0] + TINY_NUM

	return np.sum(np.where(a != 0, a * np.log(a / b), 0))

if __name__ == "__main__":
	main()



