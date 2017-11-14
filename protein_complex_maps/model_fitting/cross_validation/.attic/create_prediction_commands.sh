



commands_file=atobsc3_u30full_libsvm1.scale.predict_COMMANDS.sh
rm $commands_file


#Higher C = more support vectors
#for cvalue in 0.00390625 0.0078125 2 32 128 



for cvalue in 16 32 64 128 256 512
do
      #Higher G = more influence from individual datapoint
      #for gvalue in 0.015625 0.03125 0.0625 0.125 0.5
      for gvalue in 0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625 0.125 0.5
      do

          for segnum in 1 2 3 4 5
          do

              echo "/home1/03491/cmcwhite/libsvm/svm-predict -b 1 atobsc3_u30full_libsvm1.scale.test.seg${segnum} atobsc3_u30full_libsvm1.scale.train.seg${segnum}.model_c${cvalue}_g${gvalue//.}   atobsc3_u30full_libsvm1.scale.test.seg${segnum}.c${cvalue}_g${gvalue//.}.resultsWprob" >> $commands_file


              echo "/home1/03491/cmcwhite/libsvm/svm-predict -b 1 atobsc3_u30full_libsvm1.scale.train.seg${segnum} atobsc3_u30full_libsvm1.scale.train.seg${segnum}.model_c${cvalue}_g${gvalue//.} atobsc3_u30full_libsvm1.scale.train.seg${segnum}.c${cvalue}_g${gvalue//.}.resultsWprob" >> $commands_file


          done
      done

done

