#!/bin/bash
echo " "
echo " Jena Atomic Calculator (https://github.com/OpenJAC/JAC.jl)"
echo " "
if [[ $# -eq 0 ]] ; then
	echo "  Usage: 'JAC <filename>', where <filename>.jl is the JAC input file"
	echo " "
	exit 0
fi
infile=$1
JACcommand="$(echo /usr/local/bin/julia -e "'using JenaAtomicCalculator; include(\""$infile".jl\")'")"
echo $JACcommand
eval $JACcommand
