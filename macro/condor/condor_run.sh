#!/bin/bash
# Condor job wrapper: clone, build, run.
# Usage: condor_run.sh <executable in build dir> <title> <ClusterId> <ProcId>
set -e
echo GIT CLONE, CMAKE, MAKE, RUN

git clone https://gitlab.cern.ch/mkwon/RooParticle.git
cd RooParticle
mkdir build && cd build
cmake3 ..
make -j4

./"$1" "$2" "$3" "$4"
