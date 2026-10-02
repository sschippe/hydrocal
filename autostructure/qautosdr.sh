#!/bin/sh
qsub -N a_$1 -j y -cwd -m e -M schippers@jlug.de $HOME/hydrocal/autostructure/autosdr.sh $1 $2
