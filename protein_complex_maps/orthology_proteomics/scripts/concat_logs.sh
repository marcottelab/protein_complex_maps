echo "concatenating log files with run statistics into a single tidy data table"



cd logs/

#tail -n +2 ${spec}.${level}_tophit.txt | cat > ${spec}.${level}_orthology.tmp
#tail -n +2 ${spec}.${level}_nonhits.txt | cat >> ${spec}.${level}_orthology.tmp


echo "Species	Experiment	Level	Attribute	Count" > concat_orthoprot_logs.txt


cat *.log >> concat_orthoprot_logs.txt


