#!/bin/bash
AUTOSDIR=$HOME"/hydrocal/autostructure"
if [ $# -lt 1 ]
then
  echo " "
  echo "usage: autospp <filename> [<fnl-filename>]"
  echo " "
  echo "requires user input file <filename>.ppin"
  echo "         and optionally <fnl-filename>.fnl"
  echo " "
  echo "produces output files <filename>.alpha and <filename>.autos"
  echo " "
  exit
fi	
FN=$1
FN2=$2
WORKDIR=$PWD
SCRATCHDIR=$PWD/$FN
echo " "
echo " starting AUTOSTRUCTURE post processor"
echo "        file name: $FN"
echo "working directory: $WORKDIR"
echo "scratch directory: $SCRATCHDIR"
echo "autostructure dir: $AUTOSDIR"
echo " "
mkdir -p $SCRATCHDIR
touch $FN.autos &> /dev/null # create empty file
touch $FN.alpha &> /dev/null # create empty file
cd $SCRATCHDIR
ln -s $WORKDIR/$FN.ppin ppin &> /dev/null
ln -s $WORKDIR/$FN2.fnl fnl &> /dev/null
ln -s $WORKDIR/$FN.autos ocs &> /dev/null
ln -s $WORKDIR/$FN.alpha doutgnu &> /dev/null
$AUTOSDIR/AUTOSPP.exe < ppin
