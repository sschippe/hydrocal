#!/bin/bash
sudo true || exit 1

version=1.10
subversion=9
juliaversion=julia-$version.$subversion
tarfile=$juliaversion-linux-x86_64.tar.gz

tput bold && tput setaf 6; echo "Downloading and installing $juliaversion"
tput sgr0   
sudo wget https://julialang-s3.julialang.org/bin/linux/x64/$version/$tarfile
sudo tar -xvzf $tarfile -C /usr/local
sudo ln -s /usr/local/$juliaversion/bin/julia /usr/local/bin/julia
sudo rm $tarfile
tput bold && tput setaf 6; echo 'Installing jupyter'
tput sgr0   
sudo apt-get install jupyter
tput bold && tput setaf 6; echo 'Installing JAC and its dependencies'
tput sgr0
/usr/local/bin/julia -e 'using Pkg; Pkg.add("PyPlot"); Pkg.add("IJulia"); Pkg.add("Pluto"); Pkg.add("JenaAtomicCalculator");'
ln -s $PWD/JAC.bash ~/bin/JAC
JAC
