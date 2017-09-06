import csv

complexes_list = []

#with open('allComplexes.csv','rb') as csvfile:
with open('allComplexesCore.csv','rb') as csvfile:
	complexes = csv.reader(csvfile, delimiter=';')
	for row in complexes:
		if row[3] == "Human" and row[4] != '':
			complexes_list.append(row)

for complex1 in complexes_list:
	#print "%s %s: %s" % (complex1[0], complex1[1], complex1[4],)
	ids = complex1[4].split(',')
	#print "%s" % (complex1[4])
	#print "%s" % (ids)

        trim_ids = []
	for id1 in ids:
		id1 = id1.translate(None, '()')
		print id1

        #print ' '.join(trim_ids)
