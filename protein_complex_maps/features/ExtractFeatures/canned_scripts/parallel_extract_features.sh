## Example of how to parallelize feature extraction with xargs ##

# Fill in the path to your copy of extract_features.pyO
SCRIPT='~/Programs/protein_complex_maps/protein_complex_maps/features/ExtractFeatures/canned_scripts/extract_features.py'

# Change the feature, resampling, etc. below
ls *.csv | xargs -I {} -n 1 -P 15 \
python $SCRIPT {} \
--feature canberra \
--resampling poisson_noise \
--iterations 100 \
--threshold 5 \
--as_pickle
