echo "Loading OrthogroupIDs"
python load_ortho.py --orthogroup_file static/data/all_OrthogroupIDs.txt
#Then load rest. Since Orthogroup and Protein IDs already exist and don't have to be created, the rest goes faster
echo "Loading clusters"
#Each orthogroup belongs in multiple clusters.
python load_clusters.py --cluster_table static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold504.walktrap.csv
echo "Load orthogroup annotation"
python load_orthoannot.py --annotation_file static/data/virNOG_collapse_annotations.txt 
echo "Load CF-MS scores"
python load_scores.py --score_file static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.top100k 

python load_prot.py --conversion_file static/data/plant_virNOG_orthology.csv    

#echo "Loading ProteinID"
#python load_prot.py --protein_file static/data/all_ProteinIDs.txt
#echo "Load conversion between protein and orthogroup"
#python load_conversion.py --conversion_file static/data/plant_virNOG_orthology.csv    
