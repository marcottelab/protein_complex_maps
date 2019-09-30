#Each orthogroup belongs in multiple clusters.
python scripts/setup_images.py --cluster_table static/data/allplants_feature_matrix_missing1.unscaled.top100.edges.3steps.fdr10.threshold509.walktrap.csv

Rscript scripts/image_making.R 

