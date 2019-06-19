python load_orthoannot.py --annotation_file static/data/virNOG_collapse_annotations.txt 
python load_clusters.py --cluster_table static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold504.walktrap.csv
python load_conversion.py --conversion_file static/data/plant_virNOG_orthology.csv    
