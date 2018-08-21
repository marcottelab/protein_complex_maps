
#kdrew: this file is for utilities for querying features about proteins

import urllib, urllib2 
import requests
import MySQLdb
import difflib
import Bio.Seq as Seq
import Bio.SeqRecord as SeqRecord
import Bio.Alphabet as ba

ACC_QUERY_LENGTH = 500
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

def get_all_orthologs( prot_ids, version="_v8", database='inparanoid', score_threshold = 1.0 ): 
    ortholog_map = dict()
    db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
    cursor = db.cursor()
    query = "show tables"
    cursor.execute(query)
    for tableresult in cursor.fetchall():
        print tableresult
        tablename = tableresult[0]
        omap = get_ortholog_by_table( prot_ids, tablename, database, score_threshold, species=None)
        print omap
        for key in omap:
            try:
                ortholog_map[key].append(omap[key])
            except KeyError:
                ortholog_map[key] = [omap[key],]

    return ortholog_map

#kdrew: wrapper class
#kdrew: species2 is the organism to map to, species1 is the organism of the input prot_ids
def get_ortholog( prot_ids, species1, species2=None, version="_v8", database='inparanoid', score_threshold = 1.0, reversible=False ):
    ortholog_map = dict()
    if species2 == None:
        #kdrew: find all tables in db and call get_ortholog_by_table for each table
        db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
        cursor = db.cursor()
        query = "show tables"
        cursor.execute(query)
        for tableresult in cursor.fetchall():
            tablename = tableresult[0]
            if species1 in tablename:
                omap = get_ortholog_by_table( prot_ids, tablename, database, score_threshold, species=species1 ) 
                ortholog_map = dict( ortholog_map.items() + omap.items() )

        missing_in_orthologmap = set()
        for protid in prot_ids:
            if protid not in ortholog_map:
                missing_in_orthologmap.add(protid)

        if 0 < len(missing_in_orthologmap):
            print "Do fuzzy match, missing: %s" % (missing_in_orthologmap,)
            ids2species = get_organism_uniprot(list(missing_in_orthologmap))
            for id1 in ids2species:
                #kdrew: check to see if the passed in species is the same as id's species
                #kdrew: have to do a fuzzy match because what gets returned by uniprot is the full name 
                #kdrew: and the inparanoid name is an abbreviation ( I guess this is why they invented taxids )
                #kdrew: maybe all of this can be avoided if I create a Hsapiens_Hsapiens table in inparanoid
                print ids2species[id1], species1
                if difflib.SequenceMatcher(None, species1,ids2species[id1]).ratio() > FUZZY_MATCH_THRESHOLD:
                    ortholog_map[id1] = id1

        if reversible:
            ortholog_map.update({ortholog_map[x]:x for x in ortholog_map})

        return ortholog_map
                
    else:
        if species1 == species2:
            for prot in prot_ids:
                ortholog_map[prot] = prot
            return ortholog_map
        else:
            #kdrew: swap names if second is lower alphabetically
            if species1 > species2:
                species2_tmp, species1_tmp = species1, species2
            else:
                species1_tmp, species2_tmp = species1, species2
            #kdrew: build up query statement
            tablename = "%s_%s%s" % (species1_tmp, species2_tmp, version)

            return get_ortholog_by_table( prot_ids, tablename, database, score_threshold, species=species2) 


#kdrew: species is the organism in which to map to
def get_ortholog_by_table( prot_ids, tablename, database='inparanoid', score_threshold = 1.0, species=None ): 
    ortholog_map = dict()
    #kdrew: if same species map each protein id to itself

    db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
    cursor = db.cursor()
    query = "select a1.protein_id as pid1, a2.protein_id as pid2 from %s as a1, %s as a2 where a1.cluster = a2.cluster and a1.species <> a2.species" % (tablename, tablename,) 
    query = query + " and a2.score >= %s" % (score_threshold,)
    placeholder = '%s'
    placeholders = ','.join([placeholder] * len(prot_ids))
    query_add = " and a1.protein_id in (%s)" % placeholders
    query = query + query_add
    if species != None:
        query_add = " and a2.species = '%s' " % (species[0]+'.'+species[1:],)
        query = query + query_add
        #print len(prot_ids)
    #print query

    cursor.execute(query, tuple(prot_ids))

    ortholog_results = cursor.fetchall()
    for pair in ortholog_results:
        ortholog_map[pair[0]] = pair[1]

    return ortholog_map

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
                #print sr
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

    print "len(protein_ids) : %s" % len(protein_ids)
    #kdrew: break up query into chunks so as not to get 414 error
    for i in xrange( (len(protein_ids)/ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH
        report_query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:(%s)&columns=id,%s" % ("+or+".join(protein_ids[start_splice:stop_splice]), keyword)
        print report_query
        try:
            f = urllib2.urlopen(report_query)
        except urllib2.HTTPError, msg:
            print msg
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
    for i in xrange( (len(gene_ids)/ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH                                                                                 
        genes_formatted = ["%s:%s" % (gene_prefix, id1,) for id1 in gene_ids[start_splice:stop_splice]]
        report_query = "https://www.uniprot.org/uniprot/?format=tab&query=organism:(%s)+and+(%s)%s" % (organism, "+or+".join(genes_formatted), reviewed_str )
        print report_query
        try:
            #f = urllib2.urlopen(encoded_query)
            f = requests.get(report_query)
            #print f.text
            for line in f.text.split('\n'):
                print line
                #kdrew: the returned line has a bunch of extra gene names (and other info), finding the intersection of the gene ids and the line returns the original queried gene id
                intersection = list(set(line.split()).intersection(gene_ids_case_incensitive))
                if intersection:
                    try:
                        return_dict[intersection[0]].append(line.split()[0])
                    except:
                        return_dict[intersection[0]] = [line.split()[0],] 
        except urllib2.HTTPError, e:
            print e, e.read()
            continue


    return return_dict


#kdrew: map_protein_ids using uniprot changed their to_pdb formatting so it no longer has chain information
#kdrew: this is from Martin's Lab server at UCL: http://bioinformatics.oxfordjournals.org/content/21/23/4297.full
#kdrew: ids should be uniprot ac 
#kdrew: there is a server but downloading the db and querying it should be faster
def map_protein_ids_to_pdb( id_list, database='pdbsws' ):
    #http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=ac&id=P38764

    id2pdb = dict()
    db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
    cursor = db.cursor()
    #kdrew: build up query statement
    placeholder = '%s'
    placeholders = ','.join([placeholder] * len(id_list))
    query = "select pdbid, chain, uniprot from pdb2uniprot as p2u where uniprot in (%s)" % placeholders
    print query
    
    cursor.execute(query, id_list)
    
    results = cursor.fetchall()
    for mapping in results:
        try:
            id2pdb[mapping[2]].append(mapping[0:2])
        except KeyError:
            id2pdb[mapping[2]] = [mapping[0:2],]
    
    return id2pdb

#kdrew: return dictionary of acc -> list(chainids)
def get_pdb_protein_ids_chainid_map( pdbid, database='sifts_pdbsws'):
#http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=ac&id=P38764

    protid2chain = dict()
    db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
    cursor = db.cursor()
    query = "select pdbid, chain, uniprot from pdb2uniprot as p2u where pdbid = '%s'" % (pdbid,)
    print query

    cursor.execute(query)

    results = cursor.fetchall()
    for mapping in results:
        print mapping
        acc = mapping[2]
        if acc == '?':
            acc = None
        try:
            protid2chain[acc].add(mapping[1])
        except KeyError:
            protid2chain[acc] = set([mapping[1]])

    return protid2chain

#kdrew: return dictionary of chainid -> acc
#def get_pdb_protein_ids( pdbid, database='pdbsws', reversible=False ):
def get_pdb_protein_ids( pdbid, database='sifts_pdbsws', reversible=False ):
#http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=ac&id=P38764

    chain2protid = dict()
    db = MySQLdb.connect("localhost", 'kdrew', 'kdrew_utexas', database)
    cursor = db.cursor()
    query = "select pdbid, chain, uniprot from pdb2uniprot as p2u where pdbid = '%s'" % (pdbid,)
    print query

    cursor.execute(query)

    results = cursor.fetchall()
    for mapping in results:
        acc = mapping[2]
        if acc == '?':
            acc = None
        chain2protid[mapping[1]] = acc

    if reversible:
        ##kdrew: return lists instead, this is getting hacky, probably should be a separate object
        ##kdrew: initialize lists
        #for x in chain2protid:
        #    chain2protid[chain2protid[x]] = list()
        ##kdrew: initialize populate
        #for x in chain2protid:
        #    chain2protid[chain2protid[x]].append(x)
        
        #kdrew: overwrites multiple chains
        chain2protid.update({chain2protid[x]:x for x in chain2protid})

    
    return chain2protid

#kdrew: uses uniprot webservice to map ids
#kdrew: from_id and to_id are abbreviations of dbid names which can be found: http://www.uniprot.org/faq/28
def map_protein_ids( id_list, from_id, to_id, reviewed=False ):
    url = 'https://www.uniprot.org/mapping/'

    query_str = ' '.join(id_list)
    params = {
        'query':query_str,
        'from':'%s' % (from_id,),
        'to':'%s' % (to_id,),
        'format':'tab',
        #kdrew: noticed a TypeError exception in join but not sure why
    }

    print query_str

    data = urllib.urlencode(params)
    print data
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

