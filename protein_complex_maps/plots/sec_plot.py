
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import csv
import protein_complex_maps.sec.sec_ids as si

TINY_NUM = 0.000001
sec_filename = '/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.csv'

#vps35:Q96QK1
#vps26a:O75436
#vps26b:Q4G0F5
#vps29:	Q9UBQ0
#ttc4:	O95801
retromer_ids =  ['Q96QK1','O75436','Q4G0F5','Q9UBQ0','O95801']

#sh3glb1 commd10 commd3  commd5  commd1  ccdc22  ccdc93  commd2
#sh3glb1:Q9Y371:ENSG00000097033 
#commd2:Q86X83:ENSG00000114744:199aa	
#commd10:Q9Y6G5:ENSG00000145781:202aa	 
#commd3:Q9UBI1:ENSG00000148444:195aa	
#commd5:Q9GZQ3:ENSG00000170619:224aa	
#commd1:Q8N668:ENSG00000173163:190aa	
#ccdc22:O60826:ENSG00000101997:627aa
#ccdc93:Q567U6:ENSG00000125633:631aa	

commd_ids = ['Q86X83', 'Q9Y371', 'Q9Y6G5', 'Q9UBI1', 'Q9GZQ3', 'Q8N668',  'O60826',  'Q567U6'] 
commd_id_map = {'Q9Y371':'sh3glb1', 'Q86X83':'commd2', 'Q9Y6G5':'commd10', 'Q9UBI1':'commd3', 'Q9GZQ3':'commd5', 'Q8N668':'commd1',  'O60826':'ccdc22',  'Q567U6':'ccdc93'} 
#kdrew: without sh3glb1
#commd_ids = ['Q86X83', 'Q9Y6G5', 'Q9UBI1', 'Q9GZQ3', 'Q8N668',  'O60826',  'Q567U6'] 
#commd_id_map = {'Q86X83':'commd2', 'Q9Y6G5':'commd10', 'Q9UBI1':'commd3', 'Q9GZQ3':'commd5', 'Q8N668':'commd1',  'O60826':'ccdc22',  'Q567U6':'ccdc93'} 


dimer_ids = ['P35251', 'P12956']
dimer_id_map = {'P12956':'P12956','P35251':'P35251'}

#subcomplex_ids = dimer_ids
subcomplex_ids = commd_ids
id_map = commd_id_map
plot_filename = "/home/kdrew/public_html/test/commd_sec.pdf"
#plot_filename = "dimer_sec_all.pdf"
#plot_filename = "commd_sec_combined.pdf"
combine=False


#kdrew: sec data (fractions 7-40) keyed on ensembl ENSG
sec_data = dict()

with open(sec_filename, 'rb') as sec_file:
	reader = csv.reader(sec_file)
	for row in reader:
		if row[0] in subcomplex_ids:
			sec_data[ row[0] ] = map( float, row[1:34] )

#print sec_data

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax3 = ax1.twinx()


#kdrew: make initial array all equal to 1.0
combined_data = np.repeat(1.0,len(sec_data[sec_data.keys()[0]]))
print combined_data


individual_color = '0.90'

#plt.plot( sec_data['Q16401'] )
#for i in sec_data.keys():
for i in subcomplex_ids:
	if i in subcomplex_ids:
		if combine:
			ax1.plot( sec_data[i], color=individual_color, zorder=0)
			ax2.plot( sec_data[i], color=individual_color, zorder=-1)
		else:
			ax1.plot( sec_data[i])
			ax2.plot( sec_data[i])

		sec_data_noised = np.array(sec_data[i]) + TINY_NUM
		combined_data = combined_data * sec_data_noised
		print "########################"
		print np.array(sec_data[i]).max() 
		print combined_data
		print combined_data.max()
		combined_data = combined_data/combined_data.max()
		print combined_data


if combine:
	ax1.plot( combined_data , linewidth=4, zorder=1)
else:
	#ax1.legend(sec_data.keys())
	ax1.legend([id_map[i] for i in subcomplex_ids])

xlabels = range(10,41,5)

toplabels = [670, 440,130,67,15]

ax1.set_xticks(range(4,33,5))
ax1.set_xticklabels(xlabels)
ax1.set_xlabel('Fraction')

ax2.set_xticks([10,13,15,19,27])
ax2.set_xticklabels(toplabels)
ax2.set_xlabel('kDa')
ax2.set_zorder(-1)

#ax2.set_title("Concensus SEC of Commander Complex subunits")
#ax2.set_title("SEC of Commander Complex subunits")

plt.savefig(plot_filename)


