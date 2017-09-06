#kdrew: this script is for creating a dictionary of pdbs that have valid stoichiometry information and protein ids
#kdrew: it produces a pickle of the dictionary 


import urllib2, StringIO, csv
import string
import cPickle


search_url = 'http://www.rcsb.org/pdb/rest/search'
report_url="http://www.rcsb.org"

#queryText = """
#<?xml version="1.0" encoding="UTF-8"?>
#<orgPdbQuery>
#<version>B0907</version>
#<queryType>org.pdb.query.simple.ExpTypeQuery</queryType>
#<description>Experimental Method Search : Experimental Method=SOLID-STATE NMR</description>
#<mvStructure.expMethod.value>SOLID-STATE NMR</mvStructure.expMethod.value>
#</orgPdbQuery>
#"""

#kdrew: taken from the pdbs stoichiometry page: http://www.rcsb.org/pdb/browse/stoichiometry.do
entity_counts_list = [[1,1],[2,1],[2,2],[3,1],[3,2],[3,3],[4,1],[4,2],[4,3],[4,4],[5,1],[5,5],[6,1],[6,2],[6,3],[6,6],[7,7],[8,1],[8,2],[8,4],[8,6],[8,8],[9,1],[9,3],[9,9],[10,5],[10,10],[12,1],[12,2],[12,6],[12,12],[14,7],[14,14],[15,12],[16,8],[18,12],[20,20],[24,12],[24,24],[60,48],[60,60],[64,32],[108,108],[120,60],[132,33],[144,36],[168,42],[180,60],[180,120],[180,180],[240,240],[360,180],[420,420],[780,60],[780,120],[1,1,1],[2,1,1],[2,2,1],[2,2,2],[3,1,1],[3,3,1],[3,3,2],[3,3,3],[4,1,1],[4,2,2],[4,3,1],[4,4,1],[4,4,4],[5,1,1],[5,2,1],[5,3,2],[5,5,5],[6,3,3],[6,4,2],[6,4,4],[6,6,3],[8,6,6],[10,10,1],[10,10,10],[12,2,2],[14,14,14],[18,18,12],[24,24,24],[36,12,6],[60,60,60],[120,120,60],[180,120,120],[240,180,180],[240,240,240],[480,240,240],[720,120,60],[780,120,60],[900,60,60],[2,1,1,1],[2,2,1,1],[2,2,2,1],[2,2,2,2],[3,1,1,1],[3,2,1,1],[3,3,1,1],[3,3,2,2],[3,3,3,3],[4,2,2,2],[4,4,1,1],[4,4,4,2],[4,4,4,4],[6,3,3,3],[8,8,8,8],[12,6,6,6],[14,4,4,4],[14,5,5,5],[14,6,5,5],[14,6,6,6],[18,18,6,3],[60,60,60,60],[72,72,36,36],[120,60,60,60],[176,176,176,176],[180,180,120,120],[180,180,180,180],[240,60,60,60],[600,240,60,60],[720,60,60,60],[1200,120,120,60],[2,1,1,1,1],[2,2,2,2,1],[2,2,2,2,2],[3,3,1,1,1],[7,2,2,2,1],[10,3,3,1,1],[10,8,4,4,2],[12,6,6,6,6],[16,8,4,4,4],[60,60,60,60,60],[120,60,60,60,60],[180,60,60,60,60],[240,240,240,240,240],[720,240,120,60,60],[2,1,1,1,1,1],[2,2,1,1,1,1],[3,1,1,1,1,1],[3,3,1,1,1,1],[8,3,3,1,1,1],[10,3,3,1,1,1],[60,60,60,60,60,60],[600,600,120,60,60,60],[600,600,180,120,60,60],[3,3,2,2,1,1,1],[4,2,2,2,2,2,2],[2,2,2,2,2,2,2,1],[2,2,2,2,2,2,2,2],[4,2,2,2,2,2,2,2],[2,2,2,2,2,2,2,2,2],[3,3,1,1,1,1,1,1,1],[4,2,2,2,2,2,2,2,2],[4,2,2,2,2,2,2,2,2,2],[10,6,4,2,2,2,2,2,2,2],[20,6,6,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,1,1],[6,3,3,3,3,3,3,3,3,3,3],[1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,1],[1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,4,2,2,2,2,2,2,2,2,2,2,2,2],[6,4,4,4,2,2,2,2,2,2,2,2,2,2,2],[14,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

#entity_counts_list = [[6,3,3,3],[8,8,8,8],[12,6,6,6],[14,4,4,4],[14,5,5,5],[14,6,5,5]]
#entity_counts_list = [[3,1,1],[3,3,1],[3,3,2],[3,3,3],[4,1,1],[4,2,2]]
#entity_counts_list = [[10,10],[12,1],[12,2],[12,6],[12,12],[14,7]]

identityCutoff = 100

def main():
	pdbid_dict = dict()
	for entity_counts in entity_counts_list:
		print entity_counts
		pdbids = get_pdbids_by_stoichiometry(entity_counts)
		for pdbid in pdbids:
			protids_dict, sanity = map_stoichiometry_to_protids(pdbid, entity_counts)
			print "%s : %s -> %s" % (pdbid, sanity, protids_dict)
			if sanity:
				pdbid_dict[pdbid] = protids_dict
	
	#kdrew: remove redundant pdbs that have the same protein ids and stoichiometry
	pdbid_dict_uniq = dict()
	for key, value in pdbid_dict.items():
		if value not in pdbid_dict_uniq.values():
			pdbid_dict_uniq[key] = value

	print pdbid_dict_uniq
	cPickle.dump( pdbid_dict_uniq, open( "pdb_stoichiometry.p", "wb"))

#entity_counts = (780, 120)
#entity_counts = (2, 1)

def get_pdbids_by_stoichiometry(entity_counts):
	stoichiometry_string = ""
	for i, cnt in enumerate(entity_counts):
		stoichiometry_string += string.uppercase[i]
		#kdrew: pdb stoichiometry format is not explicit about '1' so don't add cnt if 1
		if cnt != 1:
			stoichiometry_string += str(cnt)

	##kdrew: pdb stoichiometry format is not explicit about '1' so change those to empty strings
	#stoichiometry = "A%sB%s" % tuple([x if x != 1 else '' for x in entity_counts])

#	queryText = """
#	<orgPdbQuery>
#	<version>head</version>
#	<queryType>org.pdb.query.simple.StoichiometryQuery</queryType>
#	<description>Stoichiometry in biological assembly</description>
#	<stoichiometry>%s</stoichiometry>
#	</orgPdbQuery>
#	""" % (stoichiometry_string,)


	queryText = """
		<orgPdbCompositeQuery version="1.0">
			<queryRefinement>
				<queryRefinementLevel>0</queryRefinementLevel>
				<orgPdbQuery>
					<version>head</version>
					<queryType>org.pdb.query.simple.StoichiometryQuery</queryType>
					<stoichiometry>%s</stoichiometry>
				</orgPdbQuery>
			</queryRefinement>
			<queryRefinement>
				<queryRefinementLevel>1</queryRefinementLevel>
				<conjunctionType>and</conjunctionType>
				<orgPdbQuery>
					<version>head</version>
					<queryType>org.pdb.query.simple.HomologueReductionQuery</queryType>
					<identityCutoff>%s</identityCutoff>
				</orgPdbQuery>
			</queryRefinement>
		</orgPdbCompositeQuery>
	""" % (stoichiometry_string, identityCutoff)

	print "query:\n", queryText

	print "querying PDB...\n"

	req = urllib2.Request(search_url, data=queryText)

	f = urllib2.urlopen(req)

	result = f.read()


	if result:
			print "Found number of PDB entries:", result.count('\n')
	else:
			print "Failed to retrieve results" 

	result_pdbids = result.split()

	pdbids_str = ','.join(result_pdbids)

	return result_pdbids

def map_stoichiometry_to_protids(pdbid, entity_counts):
	protid_dict = dict()
	stoichiometry_dict = dict()

	#kdrew: query pdb for chains/entities, protein ids, protein id db, and molecular type by pdbid
	report_query="%s/pdb/rest/customReport?pdbids=%s&customReportColumns=entityId,db_id,db_name,entityMacromoleculeType&service=wsdisplay&format=csv&ssa=n" % (report_url, pdbid)

	print report_query

	f2 = urllib2.urlopen(report_query)
	response2 = f2.read()
	response2 = response2.replace('<br />', '\n')
	output2 = StringIO.StringIO(response2)
	result2 = csv.reader(output2)

	#kdrew: count the number of pdb chains by protein id to determine entity counts
	for row in result2:
		print row
		#kdrew: ignore header line
		if row[0] != pdbid or row[5] != 'Polypeptide(L)':
			continue
		try:
			stoichiometry_dict[row[3]] += 1
		except KeyError:
			stoichiometry_dict[row[3]] = 1

	sorted_stoich_dict = sorted(stoichiometry_dict.items(), key=lambda x: x[1], reverse=True)
	print sorted_stoich_dict

	#kdrew: this is a flag for whether the observed entity counts match the reported stoichiometry
	overall_sanity = True
	for i, stoich in enumerate(sorted_stoich_dict):
		try:
			print "%s : %s" % (stoich[0], entity_counts[i])
			protid_dict[stoich[0]] = entity_counts[i]
			#kdrew: this check is needed so we do not check ratio for last element and one that does not exist (i.e. i+1)
			if i+1 < len(sorted_stoich_dict):
				#kdrew: make sure the ratio of counts reported by stoichiometry is the same as the ratio of observed protein ids
				ratio_sanity = (1.0*sorted_stoich_dict[i][1])/sorted_stoich_dict[i+1][1] == 1.0*entity_counts[i]/entity_counts[i+1] 
				print "ratio sanity check: %s" % ( ratio_sanity, )
				if not ratio_sanity:
					overall_sanity = False
		except IndexError:
			overall_sanity = False
	
	#kdrew: make sure the number of entities match the number of protein ids found
	if len(sorted_stoich_dict) != len(entity_counts):
		overall_sanity = False

	return protid_dict, overall_sanity


if __name__ == "__main__":
	main()

