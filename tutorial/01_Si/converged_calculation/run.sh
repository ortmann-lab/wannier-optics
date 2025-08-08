#!/bin/bash    

mpirun -np 56 ./wo-coulomb.x input.ini
python ./wo-tight-binding.py . --optimize
mpirun -np 28 ./wo-optics.x
./wo-absorption.py output -o spectrum.dat
./compare_to_exp.py
