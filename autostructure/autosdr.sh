#!/bin/bash
AUTOSDIR=$HOME"/hydrocal/autostructure"
if [ $# -lt 1 ]
then
  echo " "
  echo "usage: autosdr <filename> [<fnl-filename>]"
  echo " "
  echo "requires user input file <filename>.drin"
  echo "         and optionally <fnl-filename>.fnl"
  echo " "
  echo "produces output files <filename>.adr and <filename>.bdr"
  echo " "
  exit
fi	
FN=$1
FNL=$2
WORKDIR=$PWD
SCRATCHDIR=$PWD/$FN
echo " "
echo " starting AUTOSTRUCTURE DR post processor"
echo "        file name: $FN"
echo "working directory: $WORKDIR"
echo "scratch directory: $SCRATCHDIR"
echo "autostructure dir: $AUTOSDIR"
echo " "
touch $WORKDIR/$FN.drout &> /dev/null
touch $WORKDIR/$FN.adr &> /dev/null
touch $WORKDIR/$FN.bdr &> /dev/null
cd $SCRATCHDIR
rm -f adasout ocs XDRTOT drin fnl &> /dev/null
ln -sf $WORKDIR/$FN.drin drin &> /dev/null
ln -sf $WORKDIR/$FNL.fnl fnl &> /dev/null
ln -sf $WORKDIR/$FN.drout adasout &> /dev/null
ln -sf $WORKDIR/$FN.bdr ocs &> /dev/null
ln -sf $WORKDIR/$FN.adr XDRTOT &> /dev/null
$AUTOSDIR/AUTOSDR.exe < drin
cd $WORKDIR
echo " "
echo "general output written to $FN.drout"
echo "binned DR cross sections written to $FN.bdr"
echo "DR rate coefficients written to $FN.adr"
echo " "

