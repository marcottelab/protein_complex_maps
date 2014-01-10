
import cPickle

def main():

	pdbs = cPickle.load(open("pdb_stoichiometry.p"))

	f = open("/home/kdrew/data/protein_complex_maps/shared_complexes/source_data/elutions_protein_counts/ids/pdb_protids_grep.txt")
	ms_protids = []
	for i in f.readlines():
		ms_protids.append(i.replace('\n',''))

	f2 = open("/home/kdrew/data/protein_complex_maps/lamond_sec/mcp.M113.032367-2.uniprot.ids")
	sec_protids = []
	for i in f2.readlines():
		sec_protids.append(i.replace('\n',''))

	
	#kdrew: find pdbids whose protein ids are also in our ms data
	ms_complete = set()
	for id in pdbs:
		pdb_complete = True
		for protid in pdbs[id]:
			if protid not in ms_protids:
				pdb_complete = False

		if pdb_complete:
			ms_complete.add(id)

	for i in ms_complete:
		print i, pdbs[i]

	sec_ms_complete = set()

	for i in ms_complete:
		pdb_complete = True
		for protid in pdbs[i]:
			if protid not in sec_protids:
					pdb_complete = False
		if pdb_complete:
			sec_ms_complete.add(i)
	
	print "\n"
	for i in sec_ms_complete:
		print i, pdbs[i]

	print len(ms_complete)
	print len(sec_ms_complete)

if __name__ == "__main__":
	main()

