
commands_file=atobsc3_full_libsvm1.scale.poissonplus.train_COMMANDS.sh
rm $commands_file

#Higher C = more support vectors
#for cvalue in 0.00390625 0.0078125 2 32 128 
#for cvalue in  32 128 250
for cvalue in 32 128 250 500 1000
do   
      #Higher G = more influence from individual datapoint
      for gvalue in 0.015625 0.03125 0.0625 0.125 0.5
      do

          for segnum in 1 2 3 4 5
          do
  
              echo "/home1/03491/cmcwhite/libsvm/svm-train -m 1000 -h 0 -b 1 -c $cvalue -g $gvalue atobsc3_full_libsvm1.scale.poissonplus.train.seg${segnum} atobsc3_full_libsvm1.scale.poissonplus.train.seg${segnum}.model_c${cvalue}_g${gvalue//.} > log_${cvalue}_${gvalue//.}_seg${segnum}.txt" >> $commands_file 


          done
      done

done

