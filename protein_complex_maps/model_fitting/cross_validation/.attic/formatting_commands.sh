
full_commands_file=atobsc3_full_libsvm1.scale.poissonplus.formatting_COMMANDS.sh
rm $full_commands_file


#Higher C = more support vectors
#for cvalue in 0.00390625 0.0078125 2 32 128


for cvalue in  32 128 250 500 1000
#for cvalue in 500 1000
do   
      #Higher G = more influence from individual datapoint
      #for gvalue in 0.001953125 0.00390625 0.0078125 0.015625
      for gvalue in 0.015625 0.03125 0.0625 0.125 0.5
      do

          for segnum in 1 2 3 4 5
          do
          for analysis in test train
          do

              basename=atobsc3_full_libsvm1.scale.poissonplus.${analysis}.seg${segnum}
              probs=$basename.c${cvalue}_g${gvalue//.}.resultsWprob
              commands_file=${probs}_formatting_COMMANDS.sh
              rm $commands_file


              echo "tail -n +2 $probs > ${probs}_noheader" >> $commands_file
              echo "paste -d' '  ${basename}_labels  ${probs}_noheader > ${probs}_labels" >> $commands_file
              echo "awk -F' ' '{print \$1, \$2, \$4}' ${probs}_labels > ${probs}_pairs" >> $commands_file
              echo "sort -g -k3 -r ${probs}_pairs > ${probs}_pairs_sort" >> $commands_file
 
              done
          done
      done

done


for f in *formatting_COMMANDS.sh
do 

   echo "bash $f" >> $full_commands_file

done



