
#kdrew: this file is for utilities for querying features about proteins

import urllib, urllib2 
import MySQLdb

ACC_QUERY_LENGTH = 250

##kdrew: queries uniprot for protein sequence length
#def get_length_uniprot( protein_id ):
#	length_dict = dict()
#	report_query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:%s&columns=id,length" % (protein_id,)
#	f = urllib2.urlopen(report_query)
#	for line in f.readlines():
#		if protein_id == line.split()[0]:
#			return int(line.split()[1])
#	return None


def get_ortholog( prot_ids, species1, species2, version="_v8", database='inparanoid', score_threshold = 1.0 ):

	ortholog_map = dict()
	#kdrew: if same species map each protein id to itself
	if species1 == species2:
		for prot in prot_ids:
			ortholog_map[prot] = prot
		return ortholog_map

	db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
	cursor = db.cursor()
	#kdrew: swap names if second is lower alphabetically
	if species1 > species2:
		species2, species1 = species1, species2
	#kdrew: build up query statement
	table_name = "%s_%s%s" % (species1, species2, version)
	query = "select a1.protein_id as pid1, a2.protein_id as pid2 from %s as a1, %s as a2 where a1.cluster = a2.cluster and a1.species <> a2.species" % (table_name, table_name,) 
	query = query + " and a2.score >= %s" % (score_threshold,)
	placeholder = '%s'
	placeholders = ','.join([placeholder] * len(prot_ids))
	query_add = " and a1.protein_id in (%s)" % placeholders
	query = query + query_add
	print query

	cursor.execute(query, prot_ids)

	ortholog_results = cursor.fetchall()
	for pair in ortholog_results:
		ortholog_map[pair[0]] = pair[1]

	return ortholog_map


#kdrew: queries uniprot for protein sequence length
def get_length_uniprot( protein_ids ):
	return get_from_uniprot( protein_ids, "length" )

#kdrew: queries uniprot for protein sequence length
def get_genenames_uniprot( protein_ids ):
	return get_from_uniprot( protein_ids, "genes(PREFERRED)" )

def get_from_uniprot( protein_ids, keyword ):
	return_dict = dict()

	#kdrew: sometimes uniprot accs have added '-1' that does not play well with this webservice
	cleaned_protein_ids = []
	for prot in protein_ids:
		cleaned_protein_ids.append(prot.split('-')[0])

	protein_ids = protein_ids + cleaned_protein_ids

	#kdrew: break up query into chunks so as not to get 414 error
	for i in xrange( (len(protein_ids)/ACC_QUERY_LENGTH)+1 ):
		start_splice = i*ACC_QUERY_LENGTH
		stop_splice = (i+1)*ACC_QUERY_LENGTH
		report_query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:(%s)&columns=id,%s" % ("+or+".join(protein_ids[start_splice:stop_splice]), keyword)
		#print report_query
		f = urllib2.urlopen(report_query)
		for line in f.readlines():
			if line.split()[0] in protein_ids:
				if keyword == "length":
					return_dict[line.split()[0]] = int(line.split()[1])
				else:
					try:
						return_dict[line.split()[0]] = line.split()[1]
					except IndexError:
						return_dict[line.split()[0]] = None


	return return_dict

def get_from_uniprot_by_genename( gene_ids, organism="", gene_prefix="gene_exact", reviewed=False):
	return_dict = dict()
	reviewed_str = ""
	if reviewed:
		reviewed_str = "+and+reviewed:yes"

	#kdrew: break up query into chunks so as not to get 414 error
	for i in xrange( (len(gene_ids)/ACC_QUERY_LENGTH)+1 ):
		start_splice = i*ACC_QUERY_LENGTH
		stop_splice = (i+1)*ACC_QUERY_LENGTH                                                                                 
		genes_formatted = ["%s:%s" % (gene_prefix, id1,) for id1 in gene_ids[start_splice:stop_splice]]
		report_query = "http://www.uniprot.org/uniprot/?format=tab&query=organism:(%s)+and+(%s)%s" % (organism, "+or+".join(genes_formatted), reviewed_str )
		print report_query
		f = urllib2.urlopen(report_query)
		for line in f.readlines():
			#print line
			#kdrew: the returned line has a bunch of extra gene names (and other info), finding the intersection of the gene ids and the line returns the original queried gene id
			intersection = list(set(line.split()).intersection(gene_ids))
			if intersection:
				try:
					return_dict[intersection[0]].append(line.split()[0])
				except:
					return_dict[intersection[0]] = [line.split()[0],] 


	return return_dict


#kdrew: uses uniprot webservice to map ids
#kdrew: from_id and to_id are abbreviations of dbid names which can be found: http://www.uniprot.org/faq/28
def map_protein_ids( id_list, from_id, to_id ):
	url = 'http://www.uniprot.org/mapping/'

	params = {
		'from':'%s' % (from_id,),
		'to':'%s' % (to_id,),
		'format':'tab',
		#kdrew: noticed a TypeError exception in join but not sure why
		'query':' '.join(id_list)
	}

	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "kdrew@utexas.edu" 
	request.add_header('User-Agent', 'Python %s' % contact)
	response = urllib2.urlopen(request)

    #kdrew: put resulting map into dictionary of lists (one to many)
	return_dict = {}
	for line in response.readlines():
		print line
		if line.split()[0] != "From":
			try:
				if len(line.split()) > 2:
					return_dict[line.split()[0]] = return_dict[line.split()[0]] + line.split()[1:]
				else:
					return_dict[line.split()[0]].append(line.split()[1])
			except:
				if len(line.split()) > 2:
					return_dict[line.split()[0]] = line.split()[1:]
				else:
					return_dict[line.split()[0]] = [line.split()[1],] 

	for i in id_list:
		if i not in return_dict.keys():
			print "No id match for %s, adding empty list" % i
			return_dict[i] = []

	return return_dict

