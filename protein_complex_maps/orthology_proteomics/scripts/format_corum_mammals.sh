#This script takes the two corum (all mammals) files from http://mips.helmholtz-muenchen.de/genre/proj/corum/index.html (downloaded 6/6/2016)
#and the replaces the uniprot IDs with eggnog IDs

level=$1


if [[ $# -eq 0 ]] ; then
    echo 'No argument supplied to format_corum_mammals.sh'
    echo 'Please provide eggNOG level as argument'
    echo 'Ex. bash format_corum_mammals.sh euNOG'
    echo 'exiting'
    exit 1
fi


BASEDIR=$( pwd )


cd $BASEDIR/eggnog_output

awk -F'|' -v OFS='\t' '{print $1, $2}' corumcore.${level}_orthology.tab > corumcore.${level}_orthology.tmp 
awk -F'\t' -v OFS='\t' '{print $1, $6}' corumcore.${level}_orthology.tmp > $BASEDIR/corum/corumcore_${level}_orthology_dict.tab 
rm corumcore.${level}_orthology.tmp 


cd $BASEDIR/corum

head -2 corumcore_${level}_orthology_dict.tab


#grep "Human" allComplexesCore.csv > allComplexesCore_mammals.csv
#grep "Human" allComplexes.csv > allComplexes_mammals.csv

awk -F';' -v OFS=' ' '{print $5}' allComplexesCore.csv > allComplexesCore_mammals_accs.csv
sed -i 's/,/ /g' allComplexesCore_mammals_accs.csv
cp allComplexesCore_mammals_accs.csv allComplexesCore_mammals_${level}.csv



#Replace Uniprot Accession IDs with ${level}
while read p; do

    acc=${p#*	}
    group=${p%	*}
    echo $acc replaced with $group
    sed -i "s@$acc@$group@g" allComplexesCore_mammals_${level}.csv
 


done < corumcore_${level}_orthology_dict.tab


sed -i "s/(//g" allComplexesCore_mammals_${level}.csv
sed -i "s/)//g" allComplexesCore_mammals_${level}.csv



python $BASEDIR/scripts/get_unique_within_corumrows.py allComplexesCore_mammals_${level}.csv


