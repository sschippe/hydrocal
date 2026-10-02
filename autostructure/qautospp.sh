#!/bin/sh
qsub -N p_$1 -j y -cwd -m e -M Stefan.Schippers@uni-giessen.de $HOME/hydrocal/autostructure/autospp.sh $1 $2
