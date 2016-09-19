cd /project/cmcwhite/protein_complex_maps/protein_complex_maps/orthology_proteomics/

echo "consolidate MSblender output"
bash scripts/consolidate_MSblender_output.sh

echo "get elution profiles"
bash scripts/get_elutions.sh

echo "get eggnog output"
bash scripts/get_eggnog_output.sh

echo "get peptides from proteomes"
echo "This step could potentially be done on the fly"
bash scripts/get_peptides.sh








