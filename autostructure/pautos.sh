#!/bin/bash
# script for running parallelalized AUTOSTRUCTURE as a batch job
AUTOSDIR=$HOME"/hydrocal/autostructure"
FN=$1
WORKDIR=$PWD
SCRATCHDIR=$PWD/$FN
echo " "
echo "starting AUTOSTRUCTURE"
echo "        file name: $FN"
echo "working directory: $WORKDIR"
echo "scratch directory: $SCRATCHDIR"
echo "autostructure dir: $AUTOSDIR"
echo " "


#-----------------------------------------------
#use current directory as working directory
#$ -cwd
# combine std out and std err into one file
#$ -j y
#-----------------------------------------------

# ---------------------------
# set up the mpich version to use
# ---------------------------
MPICH_PROCESS_GROUP=no
#
# load the module
. /etc/profile.d/modules.sh
#module add default-ethernet

#----------------------------
# set up the parameters for qsub
# ---------------------------

#  mail at beginning/end/abort/on suspension
#$ -m beas
#$ -M Stefan.E.Schippers@iamp.physik.uni-giessen.de

#  export these environment variables
#$ -v MPI_HOME

# request between 2 and 8 slots
#$ -pe mpich 8
# gridengine allocates the max number of free slots
# and sets the variable $NSLOTS.
echo "got $NSLOTS slots."

#-----------------------------------------------
# connect files
#-----------------------------------------------
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR
ln -s $WORKDIR/$FN.in das >& /dev/null
ln -s $WORKDIR/$FN.radwin radwin >& /dev/null
ln -s $WORKDIR/$FN.o1 oic >& /dev/null

for ((i=0;i<NSLOTS;i++))
do 
	printf -v ofile "ols%02d" $i
	ln -s /dev/null $ofile >& /dev/null
done

# ---------------------------
# run the job
# ---------------------------
echo "Will run command: $MPI_HOME/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines $AUTOSDIR/PAUTOS.exe"
$MPI_HOME/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines $AUTOSDIR/PAUTOS.exe

#----------------------------
# concatenate output files
#----------------------------
cat oic00 > o1
for ((i=1;i<NSLOTS;i++))
do 
	printf -v ofile "oic%02d" $i
	printf -v ifile "o%d" $i
	mv $ofile $ifile >& /dev/null
done

