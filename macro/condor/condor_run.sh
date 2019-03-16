echo GIT CLONE, CMAKE, MAKE, RUN

git clone https://gitlab.cern.ch/mkwon/RooParticle.git
cd RooParticle
mkdir build && cd build
cmake3 ..
make

$1 $2 $3 $4
