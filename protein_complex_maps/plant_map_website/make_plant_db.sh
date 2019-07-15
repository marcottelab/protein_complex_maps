echo "Loading OrthogroupIDs"
python scripts/load_ortho.py --orthogroup_file static/data/all_OrthogroupIDs.txt
#Then load rest. Since Orthogroup IDs already exist and don't have to be created, the rest goes faster
echo "Loading clusters"
#Each orthogroup belongs in multiple clusters.
python scripts/load_clusters.py --cluster_table static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold504.walktrap.csv
echo "Load orthogroup annotation"
python scripts/load_orthoannot.py --annotation_file static/data/virNOG_collapse_annotations.txt 
echo "Load CF-MS scores"
gunzip static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.top10k.gz
python scripts/load_scores.py --score_file static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.top10k 
echo "Load Orthogroup Protein conversion"
python scripts/load_prot.py --conversion_file static/data/plant_virNOG_orthology.csv    
