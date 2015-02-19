import csv

complexes_list = []
complex_dict = {}

with open('allComplexes.csv','rb') as csvfile:
	complexes = csv.reader(csvfile, delimiter=';')
	for row in complexes:
		if row[3] == "Human" and row[4] != '':
			complexes_list.append(row)
			complex_dict[row[0]] = row

subcomplex_dict = {}
for complex1 in complexes_list:
	for complex2 in complexes_list:
		#if complex2[4] in complex1[4] and complex2[0] != complex1[0]:
		if complex2[4] in complex1[4] and complex1[4] not in complex2[4]:
			print "%s %s subcomplex of %s %s" % (complex2[0], complex2[1], complex1[0], complex1[1])
			print "%s : %s" % (complex2[4], complex1[4],)
			print "\n"

			try:
				subcomplex_dict[complex1[0]].append(complex2[0])
			except KeyError:
				subcomplex_dict[complex1[0]] = [complex2[0]]


for c in subcomplex_dict:
	for sc in subcomplex_dict[c]:
		print "%s %s has subcomplex %s %s" % (complex_dict[c][0], complex_dict[c][1], complex_dict[sc][0], complex_dict[sc][1], )
		print "%s : %s" % (' '.join(complex_dict[c][4].split(',')), ' '.join(complex_dict[sc][4].split(',')))
		print "\n"
