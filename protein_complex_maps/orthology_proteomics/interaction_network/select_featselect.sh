
PROJECT_DIR=/home/claire/protein_complex_maps/protein_complex_maps
BASEDIR=/home/claire/protein_complex_maps/protein_complex_maps/orthology_proteomics/interaction_network/training_dir
LIBSVMDIR=/home/kdrew/programs/libsvm-3.20
INFILE1=$1
INFILE0=$2

cd $BASEDIR
#echo "Start select cols zero threshold"
#awk -F' ' '{print $1,$29,$31,$33,$36,$39,$43,$46,$49,$50,$51,$52,$53,$54,$55,$57,$58,$59,$60,$63,$64,$67,$68,$70,$71,$72,$73,$75,$76,$77,$78,$79,$82,$83,$84,$85,$89,$91,$93,$94,$99,$100,$102,$105,$106,$107,$109,$112,$115}' $INFILE1 > ${INFILE1%.txt}RTE_abovezero.txt
#echo "End select cols zero threshold"


#echo "Start select cols 24 threshold"
#awk -F' ' '{print $1,$31,$36,$39,$51,$53,$58,$64,$67,$71,$72,$78,$82,$84,$85,$89,$93,$94,$100,$102,$105,$106,$107,$109,$112}' $INFILE1 > ${INFILE1%.txt}RTE_top24.txt
#echo "End select cols 24 threshold"

echo "scale training zero"
$LIBSVMDIR/svm-scale -s ${INFILE1%.txt}RTE_abovezero.scale_parameters ${INFILE1%.txt}RTE_abovezero.txt > ${INFILE1%.txt}RTE_abovezero.scale.txt
echo "scale training zero done"


echo "scale training 24"
$LIBSVMDIR/svm-scale -s ${INFILE1%.txt}RTE_top24.scale_parameters ${INFILE1%.txt}RTE_top24.txt > ${INFILE1%.txt}RTE_top24.scale.txt
echo "scale training 24 done"

#python $LIBSVMDIR/tools/grid.py ${INFILE1%.txt}RTE_abovezero.scale.txt
#python $LIBSVMDIR/tools/grid.py ${INFILE1%.txt}RTE_top24.scale.txt

#echo "grid.pys for libsvm1 done"


#template command
#/home/kdrew/programs/libsvm-3.20/svm-scale -r $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm1.scale_parameters $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.txt > $PROJECT_DIR/interaction_network/atobs_euNOG_corum_train_downsamp_labeled.libsvm0.scaleByTrain.txt



awk -F' ' '{print $1,$29,$31,$33,$36,$39,$43,$46,$49,$50,$51,$52,$53,$54,$55,$57,$58,$59,$60,$63,$64,$67,$68,$70,$71,$72,$73,$75,$76,$77,$78,$79,$82,$83,$84,$85,$89,$91,$93,$94,$99,$100,$102,$105,$106,$107,$109,$112,$115}' $INFILE0 > ${INFILE0%.txt}RTE_abovezero.txt


echo "scale training 24"
awk -F' ' '{print $1,$31,$36,$39,$51,$53,$58,$64,$67,$71,$72,$78,$82,$84,$85,$89,$93,$94,$100,$102,$105,$106,$107,$109,$112}' $INFILE0 > ${INFILE0%.txt}RTE_top24.txt
echo "scale training 24 done"

echo "scale training zero"
$LIBSVMDIR/svm-scale -r ${INFILE1%.txt}RTE_abovezero.scale_parameters ${INFILE0%.txt}RTE_abovezero.txt >  ${INFILE0%.txt}RTE_abovezero.scaleByTrain.txt
echo "scale training zero done"


echo "scale training 24"
$LIBSVMDIR/svm-scale -r ${INFILE%.txt}RTE_top24.scale_parameters ${INFILE0%.txt}RTE_top24.txt >  ${INFILE0%.txt}RTE_top24.scale.txt
echo "scale training 24 done"


#Can't preidct without a model
#$LIBSVMDIR/svm-predict -b 1 ${INFILE0%.txt}RTE_abovezero.scale.txt  ${INFILE1%.txt}RTE_top24.scale_parameters  ${INFILE0}_corumtrain_labeled.libsvm0.scaleByTrain_c8_g05_h0.resultsWprob


#$LIBSVMDIR/svm-predict -b 1 ${INFILE0%.txt}RTE_top24.scale.txt ${INFILE1%.txt}RTE_top24.scale_parameters  ${INFILE0}_corumtrain_labeled.libsvm0.scaleByTrain_c8_g05_h0.resultsWprob
#
#python $LIBSVMDIR/tools/grid.py ${INFILE0%.txt}RTE_top24.scale.txt















