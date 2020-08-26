
echo $1 #cluster file
echo $2 #pairwise file

#kdrew: convert clusters to pairwise with cluster ids attached
echo "running cluster2pairwise"
python ~/scripts/protein_complex_maps/protein_complex_maps/postprocessing_util/cluster2pairwise.py --filename $1 --output_filename ${1/.txt/.pairsWclustID.txt} --add_cluster_id

#kdrew: generate node attribute table
echo "running cluster2node_table"
python ~/scripts/protein_complex_maps/protein_complex_maps/postprocessing_util/cluster2node_table.py --cluster_filename $1 --from_id P_ENTREZGENEID --reviewed --output_filename ${1/.txt/.nodeTable.txt}

#kdrew: generate edge attribute table
echo "running pairwise2clusterid"
python ~/scripts/protein_complex_maps/protein_complex_maps/postprocessing_util/pairwise2clusterid.py --pairwise_filename $2 --cluster_filename $1 --output_filename ${1/.txt/.edgeAttributeWClusterid.txt}

