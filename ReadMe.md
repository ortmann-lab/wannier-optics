# Wannier Optics

Publication: [Merkel, K. & Ortmann, F. Linear scaling approach for optical excitations using maximally localized Wannier functions. Journal of Physics: Materials 7, 015001 (2023).]( https://dx.doi.org/10.1088/2515-7639/ad06cd )


The program implements the linear scaling Wannier optics method as described in the paper.

Starting from two separate Wannier90 calculations for electrons and holes, the code calculates the electron-hole transition matrix elements (optical transition dipoles) and all contributions of the exciton Hamiltonian, i.e. the tight-binding model, the Coulomb interaction and local field effects. These contributions are then used to construct the exciton Hamiltonian as a sparse matrix and to calculate correlation functions and optical absorption/reflection spectra using an efficient Lanczos approach.

All calculations use MPI and can run massively parallelized on high performance clusters.

# Compilation

### Software requirements
- mandatory requirements: FFTW3, MPI, OpenMP, cmake, python
- optional requirements: gtest (for testing), doxygen (to generate documentation)

On Ubuntu/Debian you need to install:
```bash
sudo apt install cmake libfftw3-dev build-essential mpi-default-dev

# for testing (optional):
sudo apt install libgtest-dev

# for documentation (optional):
sudo apt install doxygen graphviz
```
### Compilation
```bash
./configure_and_build.sh

# make documentation:
doxygen doc_config
```


# Usage

**Requirements:**

You need to calculate Wannier functions for the upper valence band and the lower conduction band in two separate calculations using `wannier90`. Both calculations should contain
- Wannier functions as *.xsf files
- Wannier Hamiltonian (wannier90_hr.dat)
- Band structure for verification

### 1. Generate input files and shift Wannier functions to home cell

The first step is to generate the configuration files, import the *.xsf files and shift the Wannier functions to home cell.
```bash
wo-coulomb.x -g input.ini
```
You will be ask to provide the path to the valence and conduction Wannier functions and the corresponding wannier90 seednames.

Output:
```
...
Generate configuration files for you...
[+] Write configuration file input.ini
[+] Write plan file CUSTOM
Path to valence Wannier functions? <path/to/valence/WF>
wannier90 seedname for valence WF? [default=wannier90]
[+] Write vmapping.txt
...
Path to conduction Wannier functions? <path/to/conduction/WF>
wannier90 seedname for conduction WF? [default=wannier90]
[+] Write cmapping.txt
...
[+] Read mappings from file
...
```
This will generate the `input.ini`, `vmapping.txt`, `cmapping.txt` and `POSFILE` in the current directory, which are needed in the following steps.

If the Wannier functions in the *.xsf files are not already centered in the home unit cell, then the program applies a shifts and saves them as `*_shifted.xsf` files in the corresponding Wannier simulation directories.

### 2. Generate tight-binding model

We need to import the single particle Hamiltonians from `wannier90` and apply the same shifts so that they are compatible with the previously calculated Coulomb and LFE matrix elements. This is done by
```bash
wo-tight-binding.py <path/to/coulomb-calculation>
```
This will read the `vmapping.txt`, `cmapping.txt` and `POSFILE` files and generate `TINFILE_c`, `TINFILE_v`, `ONSITE_c` and `ONSITE_v`.

### 3. Calculate transition matrix elements, Coulomb integrals and local field effects

From the (shifted) *.xsf files containing the Wannier functions on a real space grid, we now need to calculate all electron-hole interactions and transition matrix elements. The central input file is `input.ini`, which was created in the previous step. The calculation can be performed by
```bash
mpirun -np 10 wo-coulomb.x input.ini
```

This will generate the `TRANSITION`, `COULOMB` and `LOCALFIELDEFFECTS` files that are needed to construct the exciton Hamiltonian
At this point we have calculated all terms of the exciton Hamiltonian.

### 4. Calculate optical absorption spectrum

In the final step, we need to construct the exciton Hamiltonian and calculate correlation functions and optical spectra using:
```bash
mpirun -np 10 ./wo-optics.x
```
The calculation can be controled using `PARFILE`.

To obtain the optical absorption spectrum please run:
```bash
wo-absorption.py output/
```