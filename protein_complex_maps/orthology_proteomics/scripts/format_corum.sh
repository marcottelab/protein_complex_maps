#This script takes the two corum (all mammals) files from http://mips.helmholtz-muenchen.de/genre/proj/corum/index.html (downloaded 6/6/2016)
#and the replaces the uniprot IDs with eggnog IDs

level=$1

BASEDIR=$( pwd )


cd $BASEDIR/eggnog_output

awk -F'|' -v OFS='\t' '{print $1, $2}' corum.${level}_orthology.tab > corum.${level}_orthology.tmp 
awk -F'\t' -v OFS='\t' '{print $1, $5}' corum.${level}_orthology.tmp > $BASEDIR/corum/corum_${level}_orthology_dict.tab 


cd $BASEDIR/corum

head -2 corum_${level}_orthology_dict.tab


grep "Human" allComplexesCore.csv > allComplexesCore_human.csv
grep "Human" allComplexes.csv > allComplexes_human.csv

awk -F';' -v OFS=' ' '{print $5}' allComplexesCore_human.csv > allComplexesCore_human_accs.csv
sed -i 's/,/ /g' allComplexesCore_human_accs.csv
cp allComplexesCore_human_accs.csv allComplexesCore_human_${level}.csv

awk -F';' -v OFS=' ' '{print $5}' allComplexes_human.csv > allComplexes_human_accs.csv
sed -i 's/,/ /g' allComplexes_human_accs.csv

cp allComplexes_human_accs.csv allComplexes_human_${level}.csv


#Replace Uniprot Accession IDs with ${level}
while read p; do

    acc=${p#*	}
    group=${p%	*}
    echo $acc replaced with $group
    sed -i "s@$acc@$group@g" allComplexesCore_human_${level}.csv
    sed -i "s@$acc@$group@g" allComplexes_human_${level}.csv
 


done < corum_${level}_orthology_dict.tab

sed -i "s/(//g" allComplexes_human_${level}.csv
sed -i "s/)//g" allComplexes_human_${level}.csv

sed -i "s/(//g" allComplexesCore_human_${level}.csv
sed -i "s/)//g" allComplexesCore_human_${level}.csv





