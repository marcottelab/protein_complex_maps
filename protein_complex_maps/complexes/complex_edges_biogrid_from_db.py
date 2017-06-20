
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.protein_util as pu

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Outputs table of protein complex edges from database")
    parser.add_argument("--outputfile", action="store", dest="outputfile", required=False, default=None,
                                    help="Filename to store output, default=stdout")
    parser.add_argument("--field_delimiter", action="store", dest="field_delimiter", required=False, default=',',
                                    help="delimiter between fields, default = ,")
    parser.add_argument("--inner_field_delimiter", action="store", dest="inner_field_delimiter", required=False, default='\t',
                                    help="delimiter within fields, default = \\t")
    parser.add_argument("--taxon", action="store", dest="taxon", required=False, default=9606,
                                    help="Global taxon of protein interactions, default = 9606")
    parser.add_argument("--pubmedid", action="store", dest="pubmedid", required=False, default='28596423',
                                    help="PubMedID for publication, default = 28596423 (Drew et al. 2017 MolSysBio)")


    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    edges = db.session.query(cdb.Edge).all()

    proteins = set()
    protein_edges = set()
    output_str = "InteractorA%sInteractorB%sInteractorA_genename%sInteractorB_genename%sInteractorA_taxon%sInteractorB_taxon%sEvidence_Code%sPubmedID%sScore%sModification%sPhenotypes%sComment_Qualifier\n" % tuple([args.field_delimiter]*11)

    for edge in edges:
        prots = edge.get_proteins()

        #kdrew: determine evidence type
        evidence_codes = []
        evidence_list = [ev.evidence_type for ev in edge.evidences]
        if len(set(["bioplex","bioplex_prey","hein_prey","hein"]).intersection(evidence_list)) > 0:
            evidence_codes.append("Affinity Capture-MS")
        if "fraction" in evidence_list:
            evidence_codes.append("Co-fractionation")

        edge_sum = 0
        PMID_set = set()
        edge_class = ""
        if 'bioplex'in evidence_list:
            edge_sum+=1
            edge_class = 'BioPlex'
            PMID_set.add("26186194")
        if 'hein' in evidence_list:
            edge_sum+=1
            if edge_sum > 1:
                edge_class = 'Multiple_Evidences'
            else:
                edge_class = 'Hein et al.'
            PMID_set.add("26496610")
        if 'fraction' in evidence_list:
            edge_sum+=1
            if edge_sum > 1:
                edge_class = 'Multiple_Evidences'
            else:
                edge_class = 'Wan et al.'
            PMID_set.add("26344197")
        if 'bioplex_prey' in evidence_list and 'bioplex' not in evidence_list:
            edge_sum+=1
            if edge_sum > 1:
                edge_class = 'Multiple_Evidences'
            else:
                edge_class = 'Matrix_Model'
            PMID_set.add("26186194")
        if 'hein_prey' in evidence_list and 'hein' not in evidence_list:
            edge_sum+=1
            if edge_sum > 1 and edge_class != 'Matrix_Model':
                edge_class = 'Multiple_Evidences'
            else:
                edge_class = 'Matrix_Model'    
            PMID_set.add("26496610")

        pmid_str = ""
        if len(PMID_set) >0:
            pmid_str = "(PMIDs: %s)" % (' '.join(list(PMID_set))) 

        output_str = output_str + args.field_delimiter.join(["%s"]*12) % (
                                                            prots[0].gene_id,
                                                            prots[1].gene_id,
                                                            prots[0].genename(),
                                                            prots[1].genename(),
                                                            args.taxon,
                                                            args.taxon,
                                                            args.inner_field_delimiter.join(evidence_codes),
                                                            args.pubmedid,  #kdrew: PubmedID
                                                            edge.score, 
                                                            "", #kdrew: Modifications
                                                            "", #kdrew: Phenotypes
                                                            "%s %s" % (edge_class, pmid_str) #Comment
                                                            ) + "\n" 

        proteins.add(prots[0].gene_id)
        proteins.add(prots[1].gene_id)
        protein_edges.add(frozenset([prots[0].gene_id,prots[1].gene_id]))
        
    if args.outputfile != None:
        f = open(args.outputfile, "wb")
        f.write(output_str)
        f.close()
    else:
        print output_str

    print "#ofProteins: %s, #ofEdges: %s" % (len(proteins),len(protein_edges))

if __name__ == "__main__":
    main()


