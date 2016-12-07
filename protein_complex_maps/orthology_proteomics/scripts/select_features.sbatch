#!/bin/bash
#SBATCH -J runsvm      # job name
#SBATCH -o runsvm.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -c 16 		# number of CPUs/task
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=claire.mcwhite@utexas.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH -A A-cm10          # Specify allocation to charge against

module load python/2.7.12

PROJECT_DIR=$SCRATCH/protein_complex_maps/protein_complex_maps
BASEDIR=$SCRATCH/protein_complex_maps/protein_complex_maps/orthology_proteomics
LIBSVMDIR=$HOME/libsvm
INFILE1=atobsc_euNOG_corumtrain_labeled.libsvm1.txt
INFILE0=atobsc_euNOG_corumtrain_labeled.libsvm0.txt


echo "Start select cols zero threshold"
awk -F' ' '{print $1,$29,$31,$33,$36,$39,$43,$46,$49,$50,$51,$52,$53,$54,$55,$57,$58,$59,$60,$63,$64,$67,$68,$70,$71,$72,$73,$75,$76,$77,$78,$79,$82,$83,$84,$85,$89,$91,$93,$94,$99,$100,$102,$105,$106,$107,$109,$112,$115}' $INFILE1 > ${INFILE1%.txt}RTE_abovezero.txt
echo "End select cols zero threshold"


echo "Start select cols 24 threshold"
awk -F' ' '{print $1,$31,$36,$39,$51,$53,$58,$64,$67,$71,$72,$78,$82,$84,$85,$89,$93,$94,$100,$102,$105,$106,$107,$109,$112}' $INFILE1 > ${INFILE1%.txt}RTE_top24.txt
echo "End select cols 24 threshold"

echo "scale training zero"
$LIBSVMDIR/svm-scale -s ${INFILE1%.txt}RTE_abovezero.scale_parameters ${INFILE1%.txt}RTE_abovezero.txt > ${INFILE1%.txt}RTE_abovezero.scale.txt
echo "scale training zero done"


echo "scale training 24"
$LIBSVMDIR/svm-scale -s ${INFILE1%.txt}RTE_top24.scale_parameters ${INFILE1%.txt}RTE_top24.txt > ${INFILE1%.txt}RTE_top24.scale.txt
echo "scale training 24 done"

python $LIBSVMDIR/tools/grid.py ${INFILE1%.txt}RTE_abovezero.scale.txt
python $LIBSVMDIR/tools/grid.py ${INFILE1%.txt}RTE_top24.scale.txt

echo "grid.pys for libsvm1 done"


awk -F' ' '{print $1,$29,$31,$33,$36,$39,$43,$46,$49,$50,$51,$52,$53,$54,$55,$57,$58,$59,$60,$63,$64,$67,$68,$70,$71,$72,$73,$75,$76,$77,$78,$79,$82,$83,$84,$85,$89,$91,$93,$94,$99,$100,$102,$105,$106,$107,$109,$112,$115}' $INFILE0 > ${INFILE0%.txt}RTE_abovezero.txt


echo "scale training 24"
awk -F' ' '{print $1,$31,$36,$39,$51,$53,$58,$64,$67,$71,$72,$78,$82,$84,$85,$89,$93,$94,$100,$102,$105,$106,$107,$109,$112}' $INFILE0 > ${INFILE0%.txt}RTE_top24.txt
echo "scale training 24 done"

echo "scale training zero"
$LIBSVMDIR/svm-scale -s ${INFILE0%.txt}RTE_abovezero.scale_parameters ${INFILE0%.txt}RTE_abovezero.txt >  ${INFILE0%.txt}RTE_abovezero.scale.txt
echo "scale training zero done"


echo "scale training 24"
$LIBSVMDIR/svm-scale -s ${INFILE0%.txt}RTE_top24.scale_parameters ${INFILE0%.txt}RTE_top24.txt >  ${INFILE0%.txt}RTE_top24.scale.txt
echo "scale training 24 done"


python $LIBSVMDIR/tools/grid.py ${INFILE0%.txt}RTE_abovezero.scale.txt
python $LIBSVMDIR/tools/grid.py ${INFILE0%.txt}RTE_top24.scale.txt












