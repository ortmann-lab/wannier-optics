#!/usr/bin/env python
"""
Postprocessing script that calculates the susceptibility / dielectric function
from a time domain signal (from exciton code) via FFT.

Citation: Merkel, K. & Ortmann, F. Journal of Physics: Materials 7, 015001 (2024)
https://dx.doi.org/10.1088/2515-7639/ad06cd
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
import argparse, os
import configparser
matplotlib.rcParams["savefig.directory"] = os.path.dirname(os.getcwd())

HBAR = 65.82119569  # reduced Planck constant in eV*fs
ARB_TO_EV = 14.39964535 # e^2/(4pi epsilon_0) in eV*Angström

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def read_OUTFILE(OUTFILE):
    """
    reads in all settings from OUTFILE.
    """

    config = configparser.ConfigParser()
    config.read(OUTFILE)

    o = config['output']
    o.getfloat('volume')

    nT = o.getint('nT')
    TWRITE = o.getint('TWRITE')
    Tstep = o.getfloat('T')
    engy_skale = o.getfloat('gamma')
    volume = o.getfloat('volume')
    dipole_norm = o.getfloat('dipole_norm')
    return nT,TWRITE,Tstep,engy_skale, volume, dipole_norm


if __name__=="__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', type=dir_path, help='path to output directory')
    parser.add_argument('--output', '-o', type=str, help='output file for spectrum', default=None)

    parser.add_argument('--broadening', type=float, help="broadening parameter in units of meV (default=80meV)", default=None)

    parser.add_argument('--tau', type=float, help="broadening parameter in units of fs (only if --broadening is not set)", default=None)
    parser.add_argument('--interp', help="interpolates the time signal before taking the FFT", action="store_true", default=False)
    parser.add_argument('--gauss', help="uses gaussian broadening instead of lorzenz broadening", action="store_true", default=False)
    parser.add_argument('--dos', help="plots (joint) density of states", action="store_true", default=False)
    parser.add_argument('--max_nT', type=int, help="restrict evaluation to a maximal number of time steps (default=all)", default=None)
    parser.add_argument('--largerTimeStep', type=int, help="restrict evaluation to every n-th time step (default=1)", default=1)


    args = parser.parse_args()

    if args.broadening is None and args.tau is None:  # nothing specified --> use default values
        args.broadening = 80  # [meV] set default

    if args.broadening is not None:
        args.tau = HBAR/args.broadening*1000
    else:
        args.broadening = HBAR/args.tau/1000  # [meV]


    outfile = os.path.join(args.path,"OUTFILE")
    print("open ", outfile)
    nT,TWRITE,Tstep,engy_skale, Omega, dipole_norm = read_OUTFILE(outfile)
    print("\nParameters from OUTIFLE:")
    print("Energy scale (eV)    = ", engy_skale)
    print("Number of time steps = ", nT+1)
    print("Length of time step  = ", Tstep)
    print("Omega (Ang^3)        = ", Omega)
    print("Norm of dipoles      = ", dipole_norm)
    print("")

    overlap_file = os.path.join(args.path,"T000KOverlap")
    print("open ", overlap_file)
    data = np.loadtxt(overlap_file)
    time = data[:,0] * HBAR/engy_skale  # fs


    print("Actual number of time steps found: ", len(time))

    if args.max_nT is not None and args.max_nT < len(time):
        print("Restrict used time steps to ", args.max_nT)
        time = time[:args.max_nT]
        data = data[:args.max_nT,:]

    if args.largerTimeStep > 1:
        print(f"Only use {args.largerTimeStep} time step in the evaluation.")
        print("Effective time step: ", Tstep*args.largerTimeStep)
        time = time[::args.largerTimeStep]
        data = data[::args.largerTimeStep,:]

    dt = time[1]-time[0]
    print("\nActual parameters used:")
    print(f"Number of time steps  = {len(time)}")
    print(f"Largest time used     = {time[-1]/1000} ps")
    print(f"Time step             = {dt} fs \t= {dt/HBAR*engy_skale} (internal units)")
    print(f"Broadening            = {args.tau} fs \t= {HBAR/args.tau*1000} meV")
    print("")

    # broadening
    if args.gauss:
        broadening = np.exp(-time**2/(2*args.tau**2)) #*np.cos(np.pi/2*it/len(time))
    else:
        broadening = np.exp(-time/ (args.tau))
    signal = data[:,2] * broadening

    # interpolate signal
    if args.interp:
        time_interp = np.linspace(0, np.max(time), len(time)*10)
        dt_interp = time_interp[1] - time_interp[0]
        interp_func = interp1d(time, signal, kind='quadratic')

        signal_interp = interp_func(time_interp)

        # Fourier transform signal
        spectrum_interp = 1- 2*8.*np.pi*ARB_TO_EV/(HBAR*Omega) *dipole_norm *np.fft.rfft(signal_interp) *dt_interp
        energy_interp = np.arange(len(spectrum_interp))* HBAR * 2*np.pi / np.max(time_interp)

    spectrum = 1 - 2*8.*np.pi*ARB_TO_EV/(HBAR*Omega) *dipole_norm *np.fft.rfft(signal) *dt   # additional factor 2 from spin
    energy = np.arange(len(spectrum))* HBAR * 2*np.pi / np.max(time)

    # save to file
    if args.output is not None:
        output_data = np.zeros((len(energy), 3))
        output_data[:,0] = energy[:]
        output_data[:,2] = np.real(spectrum)
        output_data[:,1] = np.imag(spectrum)
        print("Save file: " + args.output)
        np.savetxt(args.output,output_data, header="energy (eV), imag part, real part")

    # create plot
    fig = plt.figure(figsize=(8,8))

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.set_title("Time domain")
    ax1.set_xlabel("Time (fs)")
    ax1.set_ylabel("Amplitude (arb. units)")

    ax2.set_title("Frequency domain")
    ax2.set_xlabel("Energy (eV)")
    ax2.set_ylabel("Absorption (arb. units)")
    ax2.set_xlim([0,min([10., np.max(energy)])])


    ax1.plot(time,data[:,2], label="raw data")
    ax1.plot(time,signal, label="with broadening")
    # ax2.plot(energy, np.imag(spectrum)-np.imag(spectrum[-1]), label="plain vanilla")
    ax2.plot(energy, np.imag(spectrum), label="plain vanilla")


    print("Energy resolution (meV) = ", (energy[1] - energy[0])*1000)

    # print("Integral over epsilon", np.trapz(np.imag(spectrum), energy))
    # print("f-sum rule", np.trapz(energy*np.imag(1./spectrum), energy))
    # print("f-sum rule", 1./(2*np.pi)*np.trapz(energy/spectrum, energy))

    # aB = 0.52917721090380   # bohr radius in angström
    # # physical constants for Si
    # rs = 2.0                # dimensionless electron gas paramter (from yellium model)
    # vF = 20.9               # Fermi velocity (angström/fs)
    # qTF = (12.0/np.pi)**(1./3.) / (np.sqrt(rs) * aB)  # Thomas-Fermi wave vector in 1/angström
    # omegaP = qTF*vF/np.sqrt(3.)                           # plasmon frequency (1/fs)

    # print("plasma-frequency^2", omegaP**2)

    if args.interp:
        ax1.plot(time_interp,signal_interp, label="interpolated")
        ax2.plot(energy_interp, np.imag(spectrum_interp)-np.imag(spectrum_interp[-1]), label="interpolated")


    plt.tight_layout()
    ax1.legend()
    ax2.legend()
    plt.show()

    if args.dos:
        DOS_file = os.path.join(args.path,"T000Kdensity.dat")
        print("open ", DOS_file)
        data = np.loadtxt(DOS_file)

        data[:,0] = data[:,0]*engy_skale

        plt.plot(data[:,0], data[:,1])
        plt.yscale("log")
        plt.show()
