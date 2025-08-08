#!/usr/bin/env python
"""
Some utility functions to read and write files.
"""

import os
import re
import numpy as np


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)
    
def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)
    
def read_wsvec(path="wannier90_wsvec.dat"):
    """
    Parses the wannier90_wsvec.dat. This file contains corrections for the interpolation
    that are calculated when use_ws_distance=true is set in wannier90.
    """
    #print("read file ", path)
    with open(path) as file:
        lines = [line.strip() for line in file]

    # determine maximal degeneracy
    ndegenx = 0
    num_wann = 0
    num_distances = 0
    for line in lines:
        if line[0] == '#':  # skip comment line
            continue    

        values = re.split('[ ]+',line)
        if len(values) == 1:
            ndegenx = max([np.int32(values[0]), ndegenx])
        elif len(values) == 5:  # R, i,j
            num_distances += 1
            num_wann = max([num_wann, np.int32(values[3])])

    nrpts = num_distances//(num_wann**2)

    if ndegenx < 1 or num_wann < 1 or num_distances != nrpts*num_wann**2:
        raise NameError(f"Error while parsing {path}. Could not read degeneracies, num_wann or nrpts.")
    # print(f"ndegenx = {ndegenx}")
    # print(f"num_wann = {num_wann}")
    # print(f"nrpts = {nrpts}")

    # create data structures to store information from file
    ws_dist = np.zeros((3,ndegenx,num_wann, num_wann, nrpts), dtype=np.int16)
    num_deg = np.zeros((num_wann, num_wann, nrpts), dtype=np.int16)
    Rvec_dict = {}
    # Rvec = np.zeros((3,nrpts), dtype=np.int16)

    i = 0
    j = 0
    nR = -1
    ndeg = 0
    for line in lines:
        if line[0] == '#':  # skip comment line
            continue 

        values = re.split('[ ]+',line)
        if len(values) == 5:  # R, i,j
            i = np.int32(values[3])-1
            j = np.int32(values[4])-1

            if i==0 and j ==0:
                nR += 1
                # Rvec[:,nR] = np.array(values[:3], dtype=np.int32)
                Rvec_dict[(int(values[0]), int(values[1]), int(values[2]))] = nR
            
        elif len(values) == 1:
            ndeg = 0
            num_deg[i,j,nR] = np.int32(values[0])

        elif len(values)==3:
            ws_dist[:,ndeg,i,j,nR] = np.array(values, dtype=np.int32)
            ndeg += 1

        else:
            print(line)
            print(values)
            raise NameError(f"Error while parsing {path}. Number of elements in line is not like expected.")
        
    # print(ws_dist[:,:,4,5,0].T)
    # print(num_deg[4,5,0])
    #print(Rvec.T)

    return Rvec_dict, num_deg, ws_dist

def write_tinfile(filename, my_model):
    print("write", filename)
    NTi = len(my_model._hoppings)
    with open(filename, 'w') as f:
        f.write(f"{2*NTi}\nl1 \tl2 \tDx \tDy \tDz \tReal(TI) (eV) \t\tImag(TI) (eV)\n")

        for ti, l1, l2, R in my_model._hoppings:
            f.write(f"{l1+1}\t{l2+1}\t{R[0]}\t{R[1]}\t{R[2]}\t{ti.real:.10f}\t {ti.imag:.10f}\n")
            f.write(f"{l2+1}\t{l1+1}\t{-R[0]}\t{-R[1]}\t{-R[2]}\t{ti.real:.10f}\t {-ti.imag:.10f}\n")

            # if ti.imag > 1e-3:
            #     print("WARNING: TI has non-vanishing imaginary part", ti)


def read_posfile(path):
    unitcell = np.zeros((3,3))
    with open(path) as file:
        file.readline()  # first two lines are not interesting
        file.readline()

        for i in range(3):
            unitcell[i,:] = np.array(file.readline().strip().split("\t"),dtype=np.double)

        # print("Unitcell: ")
        # print(unitcell)

        Nval = int(file.readline().strip())
        positions_val = np.zeros((Nval,3))
        modus = file.readline().strip()
        if modus != 'C' and modus != 'c':
            raise NameError("Can only read POSFILES with cartesian coordinates!")
        for i in range(Nval):
            positions_val[i,:] = np.array(file.readline().strip().split("\t"),dtype=np.double)

        # print("Positions of valence WF:")
        # print(positions_val)

        Ncond = int(file.readline().strip())
        positions_cond = np.zeros((Ncond,3))
        modus = file.readline().strip()
        if modus != 'C' and modus != 'c':
            raise NameError("Can only read POSFILES with cartesian coordinates!")

        for i in range(Ncond):
            positions_cond[i,:] = np.array(file.readline().strip().split("\t"),dtype=np.double)

        # print("Positions of conduction WF:")
        # print(positions_cond)
    return unitcell, positions_val, positions_cond


def write_posfile(filename, my_model):
    unitcell = np.array(my_model.get_lat())
    rel_pos = np.array(my_model.get_orb())
    Nsites = my_model._norb
    positions = np.zeros((Nsites, 3))
    for i in np.arange(Nsites):
        positions[i,:] = unitcell.T @ rel_pos[i,:]

    print("write", filename)
    with open(filename, 'w') as f:
        f.write(" 20 20 20 y y y\nmodel\n")
        for v1, v2, v3 in unitcell:
            f.write(f" {v1}\t{v2}\t{v3}\n")
        f.write(f"{Nsites}\nC\n")
        for v1, v2, v3 in positions:
            f.write(f" {v1}\t{v2}\t{v3}\n")

def write_onsite(filename, my_model, shift=0.0):
    if (abs(shift) > 1e-10):
        print("shift onsite energy by", shift)

    print("write", filename)
    Nsites = my_model._norb
    with open(filename, 'w') as f:
        f.write(f"{Nsites}\n")
        for v1 in my_model._site_energies:
            f.write(f"{v1+shift}\n")

def merge_posfiles(posfile_v, posfile_c):

    with open(posfile_v, 'r') as f:
        posfile_v_data = f.readlines()

    with open(posfile_c, 'r') as f:
        posfile_c_data = f.readlines()
    
    if len(posfile_v_data) < 8:
        print("Error: Posfile needs to have at least 8 lines!!!")
        return
    
    if len(posfile_c_data) < 8:
        print("Error: Posfile needs to have at least 8 lines!!!")
        return

    new_posfile = posfile_v_data
    for line in posfile_c_data[5:]:
        new_posfile.append(line)

    print("write POSFILE")
    with open("POSFILE", 'w') as f:
        for line in new_posfile:
            f.write(line)

def write_parfile(filename="PARFILE"):

    parfile = """ [general]
 #------------------------------------------------------------------------------
 #                              TECHNICAL PARAMETERS
 #------------------------------------------------------------------------------
 
 # Number of neighboring connections for every site:
 # values < 0 : use heuristics (might not work or is inefficient)
 NNEIGH       =           -1
 
 # Additional tests and output? (T/F)
 DEBUGGING    =  F
 
 # Enforce Hamiltonian to be hermitian? (T/F)
 FORCE_HERMITIAN    =  F
 
 # Check matrix elements of the Hamiltonian for dublicates? (T/F)
 # (Could be very slow!)
 CHECK_DUBLICATES    =  F
 
 
 #------------------------------------------------------------------------------
 #                                 TIME EVOLUTION
 #------------------------------------------------------------------------------
 
 # Length of time step (in internal units = largest transfer integral):
 T            =   0.20000000000000001     
 
 # Number of time steps:
 nT           =         1000
 
 # Number of Chebyshev polynomials:
 NPOL         =           20
 
 # After which time steps you want to write output?
 TWRITE       =            1
 
 
 #------------------------------------------------------------------------------
 #                              DENSITY OF STATES (DOS)
 #------------------------------------------------------------------------------
 
 # Number of Lanczos recursions
 NRECURS      =         1000
 
 # Lorentzian broadening of the peaks
 EPS          =    9.9999999999999995E-007
 
 # Energy steps and boundarys for DOS output file:
 N_ENERGIES   =        32000
 DOSEMIN      =   -10.000000000000000     
 DOSEMAX      =    10.000000000000000     
 
 # Use projected DOS? (values = none,x,y,z)
 DOS_PROJECTION = none                
 
 
 #------------------------------------------------------------------------------
 #                                      OPTIC
 #------------------------------------------------------------------------------
 
 # Polarization direction (values = x,y,z)
 DIRECTION    = x
 
 
 #------------------------------------------------------------------------------
 #                                     DISORDER
 #------------------------------------------------------------------------------
 
 # Using additional disorder for diagonal parts of the Hamiltonian.
 # This is only experimental! There are no derivations or extensive tests for
 # exciton Hamiltonians!!!!
 # Available disorder models = none, ANDERSON
 DISORDER_MODEL    = none                
 DISORDER_STRENGTH =    0.0000000000000000     """
    
    if not os.path.exists(filename):
        print("write", filename)
        # If the file doesn't exist, write the large string to it
        with open(filename, 'w') as file:
            file.write(parfile)
    else:
        print("skip PARFILE (already exists)")