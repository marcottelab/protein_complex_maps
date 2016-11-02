echo "concatenating experiments"
#To remove

#Fix for this will be making the output of eggnog be comma separated


tail -n +2 elutiongs/*elution.csv | cat > for_mySQL/mySQL_elutions.tmp

#HEADER=`head -1 elutiongs/*elution.csv`

cd for_mySQL
rm -f elutions.tab
echo "ExperimentID,FractionID,Peptide,PeptideCount" > elutions.csv

sort -u mySQL_elutions.tmp >> elutions.csv
sed -i 's/^/,/g' elutions.csv #blank column to hold auto_increment ID

rm mySQL_elutions.tmp




