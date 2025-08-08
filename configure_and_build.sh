#!/bin/bash

# On coolmuc2 cluster following modules can be used:
# module load intel-toolkit fftw


# On your own machine (ubuntu)
# install dependencies (as root) by using
# apt-get update && apt-get upgrade --yes && xargs apt-get install --yes < dependencies.txt

# Configure a release build
cmake -S . -B build/ -D CMAKE_BUILD_TYPE=Release

# Build release binaries
cmake --build build/

# create documentation
# cd wo-coulomb
# doxygen docs_config
# rm -f documentation.html && ln -s docs/html/index.html documentation.html
# cd ../

# build a docker container and compile there
# docker build -t wannier-optics .

# run the container
# docker run -ti wannier-optics

# create debian package
# cd build/ && cpack && cd -
