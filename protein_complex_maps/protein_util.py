
#kdrew: this file is for utilities for querying features about proteins

import urllib, urllib2 

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

#kdrew: queries uniprot for protein sequence length
def get_length_uniprot( protein_ids ):
	return get_from_uniprot( protein_ids, "length" )

#kdrew: queries uniprot for protein sequence length
def get_genenames_uniprot( protein_ids ):
	return get_from_uniprot( protein_ids, "genes(PREFERRED)" )

def get_from_uniprot( protein_ids, keyword ):
	return_dict = dict()
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
					return_dict[line.split()[0]] = line.split()[1]

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
		if line.split()[0] != "From":
			try:
				return_dict[line.split()[0]].append(line.split()[1])
			except:
				return_dict[line.split()[0]] = [line.split()[1],] 

	return return_dict

