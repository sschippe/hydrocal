#!/bin/bash
AUTOSDIR=$HOME"/hydrocal/autostructure"
if [ $# -lt 1 ]
then
  echo " "
  echo "usage: autosrr <filename> [<fnl-filename>]"
  echo " "
  echo "requires user input file <filename>.rrin"
  echo "         and optionally <fnl-filename>.fnl"
  echo " "
  echo "produces output files <filename>.arr and <filename>.brr"
  echo " "
  exit
fi	
FN=$1
FNL=$2
WORKDIR=$PWD
SCRATCHDIR=$PWD/$FN
echo " "
echo " starting AUTOSTRUCTURE RR post processor"
echo "        file name: $FN"
echo "working directory: $WORKDIR"
echo "scratch directory: $SCRATCHDIR"
echo "autostructure dir: $AUTOSDIR"
echo " "
touch $WORKDIR/$FN.rrout &> /dev/null
touch $WORKDIR/$FN.arr &> /dev/null
#touch $WORKDIR/$FN.brr &> /dev/null
cd $SCRATCHDIR
rm -f adasout ocs XRRTOT fnl rrin &> /dev/null
ln -sf $WORKDIR/$FN.rrin rrin &> /dev/null
ln -sf $WORKDIR/$FNL.fnl fnl &> /dev/null
ln -sf $WORKDIR/$FN.rrout adasout &> /dev/null
ln -sf $WORKDIR/$FN.arr XRRTOT &> /dev/null
$AUTOSDIR/AUTOSRR.exe < rrin
mv oic o1 &> /dev/null
mv opic op1 &> /dev/null
cd $WORKDIR
echo " "
echo "general output written to $FN.rrout"
echo "RR cross section written to $FN.arr"
echo " "

