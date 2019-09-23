
echo "Loading OrthogroupIDs"
gunzip static/data/all_OrthogroupIDs.txt.gz
python scripts/load_ortho.py --orthogroup_file static/data/all_OrthogroupIDs.txt
#Then load rest. Since Orthogroup IDs already exist and don't have to be created, the rest goes faster

echo "Loading clusters"
#Each orthogroup belongs in multiple clusters.
gunzip static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold504.walktrap.csv.gz
python scripts/load_clusters.py --cluster_table static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold504.walktrap.csv

echo "Load orthogroup annotation"
gunzip static/data/virNOG_collapse_annotations.txt.gz 
python scripts/load_orthoannot.py --annotation_file static/data/virNOG_collapse_annotations.txt 

echo "Load CF-MS scores"
#gunzip static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.top100k.gz
python scripts/load_scores.py --score_file static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.top100k 

echo "Load Orthogroup Protein conversion"
#gunzip static/data/complete_orthology_w_atnums.csv
python scripts/load_prot.py --conversion_file static/data/complete_orthology_w_atnums.csv
#python scripts/load_prot.py --conversion_file static/data/plant_virNOG_orthology.csv    
