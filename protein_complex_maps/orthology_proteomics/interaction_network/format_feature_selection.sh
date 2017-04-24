feats=$1
infile=$2
condition=$3

echo $infile 

awk '{print $1}' $feats > ${feats}_num


echo '$1' >  ${feats}_num_inc

while read NUM
do
echo $((NUM+1)) >> ${feats}_num_inc
done < ${feats}_num

sort -n ${feats}_num_inc > ${feats}_num_ord


selection=`cat ${feats}_num_ord | tr '\n' ',' |  sed 's/,$//' |   sed 's/,/, $/g'`

#echo $selection


sed  "s@colselect@$selection@" template_featselect.sh > ${condition}_featselect.sh
sed -i "s@infile@$infile@" ${condition}_featselect.sh
sed -i "s@condition@$condition@" ${condition}_featselect.sh


#cat ${condition}_featselect.sh




#awk -F' ' '{print $1,$31,$36,$39,$51,$53,$58,$64,$67,$71,$72,$78,$82,$84,$85,$89,$93,$94,$100,$102,$105,$106,$107,$109,$112}' $INFILE0 > ${INFILE0%.txt}RTE_top24.txt



