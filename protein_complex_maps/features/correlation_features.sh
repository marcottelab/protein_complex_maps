if [ $# -ne 2 ]; then
    echo usage: bash correlation_features.sh elutions.txt correlation_type
    echo ----elutions.txt: file with 1 elution filename per line
    echo ----correlation_type: one of the following
    echo ---------euclidean
    echo ---------pearson
    echo ---------spearman
    echo ---------cov
    exit 1
fi

ELUTION_FILE=$1
CORR=$2


echo starting correlation

cat $ELUTION_FILE | parallel -j14 python ../../external/score.py {} $CORR ','

echo correlation done

cat $ELUTION_FILE | parallel -j4 python ../../features/convert_correlation.py --input_correlation_matrix {}.corr_${CORR} --input_elution_profile {} --output_file {}.corr_${CORR}.pairs 

echo pairs done
cat $ELUTION_FILE | parallel -j4 python ../../features/alphabetize_pairs.py --feature_pairs {}.corr_${CORR}.pairs --outfile {}.corr_${CORR}.pairs.ordered 





