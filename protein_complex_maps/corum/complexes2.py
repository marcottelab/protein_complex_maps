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

        if len(ids) <= 1:
            continue

        trim_ids = []
        skip_flag=False
	for id1 in ids:
            if skip_flag:
                continue
            if '(' in id1:
                id1 = id1.translate(None, '()')
                skip_flag = True
		#print id1
            elif ')' in id1:
                skip_flag = False
                continue

            trim_ids.append(id1)

        print ' '.join(trim_ids)
