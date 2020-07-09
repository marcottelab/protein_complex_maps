
#kdrew: this file is for utilities for querying features about proteins
from __future__ import print_function

#import urllib, urllib2 
import urllib.request, urllib.parse, urllib.error, urllib.request, urllib.error, urllib.parse

import requests
import difflib
import Bio.Seq as Seq
import Bio.SeqRecord as SeqRecord
import Bio.Alphabet as ba

ACC_QUERY_LENGTH = 100
FUZZY_MATCH_THRESHOLD = 0.25

##kdrew: queries uniprot for protein sequence length
#def get_length_uniprot( protein_id ):
#   length_dict = dict()
#   report_query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:%s&columns=id,length" % (protein_id,)
#   f = urllib2.urlopen(report_query)
#   for line in f.readlines():
#       if protein_id == line.split()[0]:
#           return int(line.split()[1])
#   return None


def get_organism_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "organism" )

#kdrew: queries uniprot for protein sequence length
def get_length_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "length" )

#kdrew: queries uniprot for protein sequence genename
def get_genenames_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "genes(PREFERRED)" )

def get_sequences_uniprot( protein_ids, seqrecord=False ):
    seq_map = get_from_uniprot(protein_ids, keyword="sequence") 
    if seqrecord:
        return_seq_map = dict()
        organism_map = get_organism_uniprot( protein_ids )
        for key in seq_map:
            if seq_map[key] != None:
                s = Seq.Seq(seq_map[key], ba.generic_protein)
                sr = SeqRecord.SeqRecord(s, id=key, description=organism_map[key])
                #sr = SeqRecord.SeqRecord(Seq.Seq(seq_map[key], ba.generic_protein), id=key, description=organism_map[key])
                #print(sr)
                return_seq_map[key] = sr
        return return_seq_map

    return seq_map

def get_from_uniprot( protein_ids, keyword, return_list=False ):
    return_dict = dict()

    #kdrew: sometimes uniprot accs have added '-1' that does not play well with this webservice
    cleaned_protein_ids = []
    for prot in protein_ids:
        cleaned_protein_ids.append(prot.split('-')[0])

    protein_ids = protein_ids + cleaned_protein_ids

    #kdrew: break up query into chunks so as not to get 414 error
    for i in range( (len(protein_ids)/ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH
        report_query = "https://www.uniprot.org/uniprot/?format=tab&query=accession:(%s)&columns=id,%s" % ("+or+".join(protein_ids[start_splice:stop_splice]), keyword)
        #print(report_query)
        try:
            f = urllib.request.urlopen(report_query)
        except urllib.error.HTTPError as msg:
            print(msg)
            continue

        for line in f.readlines():
            if line.split()[0] in protein_ids:
                if keyword == "length":
                    return_dict[line.split()[0]] = int(line.split()[1])
                elif keyword == "organism":
                    return_dict[line.split()[0]] = line.split('\t')[1]
                elif keyword == "protein+names":
                    return_dict[line.split()[0]] = line.split('\t')[1]
                else:
                    if return_list:
                        try:
                            return_dict[line.split()[0]] = line.split()[1:]
                        except IndexError:
                            return_dict[line.split()[0]] = []
                    else:
                        try:
                            return_dict[line.split()[0]] = line.split()[1]
                        except IndexError:
                            return_dict[line.split()[0]] = None


    return return_dict

def get_from_uniprot_by_genename( gene_ids, organism="", gene_prefix="gene_exact", reviewed=False ):
    return_dict = dict()

    gene_ids_case_incensitive = gene_ids + [x.lower() for x in gene_ids] + [y.upper() for y in gene_ids] 

    reviewed_str = ""
    if reviewed:
        reviewed_str = "+and+reviewed:yes"

    #kdrew: break up query into chunks so as not to get 414 error
    for i in range( (len(gene_ids)/ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH                                                                                 
        genes_formatted = ["%s:%s" % (gene_prefix, id1,) for id1 in gene_ids[start_splice:stop_splice]]
        report_query = "https://www.uniprot.org/uniprot/?format=tab&query=organism:(%s)+and+(%s)%s" % (organism, "+or+".join(genes_formatted), reviewed_str )
        print(report_query)
        try:
            #f = urllib2.urlopen(encoded_query)
            f = requests.get(report_query)
            #print(f.text)
            for line in f.text.split('\n'):
                print(line)
                #kdrew: the returned line has a bunch of extra gene names (and other info), finding the intersection of the gene ids and the line returns the original queried gene id
                intersection = list(set(line.split()).intersection(gene_ids_case_incensitive))
                if intersection:
                    try:
                        return_dict[intersection[0]].append(line.split()[0])
                    except:
                        return_dict[intersection[0]] = [line.split()[0],] 
        except urllib.error.HTTPError as e:
            print(e) #, e.read()
            continue


    return return_dict


#kdrew: uses uniprot webservice to map ids
#kdrew: from_id and to_id are abbreviations of dbid names which can be found: http://www.uniprot.org/faq/28
def map_protein_ids( id_list, from_id, to_id, contact_email, reviewed=False ):
    url = 'https://www.uniprot.org/mapping/'

    query_str = ' '.join(id_list)
    params = {
        'query':query_str,
        'from':'%s' % (from_id,),
        'to':'%s' % (to_id,),
        'format':'tab',
        #kdrew: noticed a TypeError exception in join but not sure why
    }

    print(query_str)

    data = urllib.parse.urlencode(params)
    print(data)
    request = urllib.request.Request(url, data)
    contact = contact_email
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)

    #kdrew: put resulting map into dictionary of lists (one to many)
    return_dict = {}
    for line in response.readlines():
        print(line)
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
            print("No id match for %s, adding empty list" % i)
            return_dict[i] = []

    if reviewed:
        #kdrew: flatten ids into set
        ret_id_list = list(set([item for sublist in return_dict.values() for item in sublist]))
        #print ret_id_list
        is_reviewed_list = get_from_uniprot(ret_id_list, 'reviewed')
        #print is_reviewed_list
        unreviewed_list = [item for item in is_reviewed_list if is_reviewed_list[item] == 'unreviewed']
        #print unreviewed_list

        for i in return_dict:
            return_dict[i] = list(set(return_dict[i]) - set(unreviewed_list))

    return return_dict

