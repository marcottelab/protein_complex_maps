
for f in proteomes/?????/*.fasta
do 

   echo $f   
   python scripts/trypsin.py --input $f --output ${f%.*}_peptides.csv --miss 2

done

#echo "concatenating peptide lists, removing duplicates"

#echo "ProteinID,Peptide" > for_mySQL/peptides.csv
#cat proteomes/?????/*peptides.csv | sort -u >> for_mySQL/peptides.csv

#sed -i 's/^/,/' for_mySQL/peptides.csv #First column empty to hold auto_increment ID





