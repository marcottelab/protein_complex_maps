#Each orthogroup belongs in multiple clusters.
python scripts/setup_images.py --cluster_table static/data/panplant_clusters.csv
Rscript scripts/image_making.R
