#!/bin/bash

# On your own machine (ubuntu)
# install dependencies (as root) by using
# sudo apt-get update && sudo apt-get upgrade --yes && xargs sudo apt-get install --yes < dependencies.txt

# Configure a release build
cmake -S . -B build/ -D CMAKE_BUILD_TYPE=Release

# Build release binaries
cmake --build build/

# create documentation
#doxygen docs_config
#rm -f documentation.html && ln -s docs/html/index.html documentation.html

# build a docker container and compile there
# docker build -t dev-coulomb .

# run the container
# docker run -ti dev-coulomb

# create debian package
# cd build/ && cpack && cd -
