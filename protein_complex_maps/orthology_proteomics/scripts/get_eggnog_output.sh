echo "concatenating peptide lists"


#Fix for this will be making the output of eggnog be comma separated


tail -n +2 eggnog_output/*tophit.txt | cat > for_mySQL/orthology.tmp
tail -n +2 eggnog_output/*nonhits.txt | cat >> for_mySQL/orthology.tmp



#HEADER=`head -1 eggnog_output/*tophit.txt`
#


cd for_mySQL
rm -f orthology.tab
echo "GroupID	Rank	Level	Species	ProteinID	evalue	QueryRange	ProteomeID	Hitlength	Sequence	Annotation" > orthology.tab

sort -u orthology.tmp >> orthology.tab
sed -i 's/^/	/g' orthology.tab  #First column must be blank

rm orthology.tmp




