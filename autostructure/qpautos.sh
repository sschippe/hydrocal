#!/bin/sh
qsub -N a_$1 $HOME/hydrocal/autostructure/pautos.sh $1 $2
