#!/bin/bash
AUTOSDIR=$HOME"/hydrocal/autostructure"
if [ $# -lt 1 ]
then
  echo " "
  echo "usage: autos <filename>"
  echo " "
  echo "requires input files <filename>.in and optionally <filename>.radwin"
  echo " "
  echo "produces output file <filename>.o1 (formatted) or <filename>.o1u (unformatted)"
  echo " "
  exit
fi	
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
mkdir -p $SCRATCHDIR
cd $SCRATCHDIR
ln -sf $WORKDIR/$FN.in das &> /dev/null
ln -sf $WORKDIR/$FN.radwin radwin &> /dev/null
$AUTOSDIR/AUTOS.exe < das
if grep "CUP='IC" das
then
   if grep "PRINT='UNFORM'" das
   then
      mv oicu o1u
   else
      mv oic o1
   fi
else
   if grep "PRINT='UNFORM'" das
   then
      mv olsu o1u
   else
      mv ols o1
   fi
fi
