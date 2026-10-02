#!/bin/bash
AUTOSDIR=$HOME"/hydrocal/autostructure"
echo " "
if [ $# -lt 1 ]
then
  echo " "
  echo "usage: autoshf <filename> [-r]"
  echo " "
  echo "no input files required "
  echo " "
  echo "option -r: reformat only"
  echo " "
  echo "produces output file <filename>.radwin "
  echo " "
  exit
fi
FN=$1
WORKDIR=$PWD
SCRATCHDIR=$PWD/$FN
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR

if [ $# -eq 1 ] # just one command line parameter
then
  echo " "
  echo "starting FROESE-FISCHER HF calculation"
  echo " "
  echo "       file name : $FN"
  echo "working directory: $WORKDIR"
  echo "scratch directory: $SCRATCHDIR"
  echo "autostructure dir: $AUTOSDIR"
  echo " "
  $AUTOSDIR/AUTOSHF.exe
fi
echo " "
echo " now reformating wavefunctions for AUTOSTRUCTURE"
echo " "
ANSWER=y
TEST=y
while [ $ANSWER == $TEST ]
do
	$AUTOSDIR/AUTOSHFR.exe
	echo " "
	echo -n "repeat with other radius? (y/n): "
	read ANSWER
done
cp -f $SCRATCHDIR/orbital.inp $WORKDIR/$FN.radwin
cd $WORKDIR
