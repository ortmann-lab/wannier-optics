#!/usr/bin/env python
"""
Imports tight-binding model from wannier90 (wannier90_hr.dat) and generates wannier-optics input files (POSFILE, TINFILE, ONSITE_ENERGY).

Citation: Merkel, K. & Ortmann, F. Journal of Physics: Materials 7, 015001 (2024)
https://dx.doi.org/10.1088/2515-7639/ad06cd
"""
from utils import *
from pythtb import * # import TB model class
from copy import deepcopy
import warnings
import argparse, os
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["savefig.directory"] = os.path.dirname(os.getcwd())


def shift_sites(my_model, shift):
    """
    Shifts the sites of the tb-model by a unitcell vector.

    shift is an np.array of dimensions (my_model._norb,3)
    """

    if my_model._norb != shift.shape[0] or 3 != shift.shape[1]:
        raise NameError("Shift array has wrong dimensions.")

    pos = my_model.get_orb()
    # shift postions
    my_model._orb = pos-shift

    # shift hoppings
    for i in range(my_model._norb):
        for h in range(len(my_model._hoppings)):
            if my_model._hoppings[h][1]==i:
                my_model._hoppings[h][3]-=shift[i,:]
            if my_model._hoppings[h][2]==i:
                my_model._hoppings[h][3]+=shift[i,:]


def optimize_tb(my_model, min_hopping, max_hopping, mesh, sparsity_penalty=0.05):

    def prune_model(my_model, min_hopping):
        pruned = deepcopy(my_model)
        new_hoppings = [h for h in pruned._hoppings if np.abs(h[0]) >= min_hopping]
        pruned._hoppings = new_hoppings
        return pruned

    k_vec = my_model.k_uniform_mesh(mesh)
    e1 = my_model.solve_all(k_vec)

    def loss(log_min_hopping):
        min_hopping = np.exp(log_min_hopping)
        model_pruned = prune_model(my_model, min_hopping)
        e2 = model_pruned.solve_all(k_vec)
        return np.mean((e1-e2)**2) + sparsity_penalty*len(model_pruned._hoppings)/2.

    bounds = (np.log(min_hopping), np.log(max_hopping))
    res = minimize_scalar(loss, bounds, bounds)

    if res.success:
        print(f"   --> optimized minimal hopping (cut-off): {np.exp(res.x)}eV")
        best_model = prune_model(my_model, np.exp(res.x))
        return best_model
    else:
        print(res)
        print("Optimization has not converged. Continue with unoptimized model.")
        return my_model


def print_summary(my_model):
    print("   Number of WF:", my_model.get_num_orbitals())
    print("   Number of transfer integrals (pruned):", len(my_model._hoppings))
    print("   Unit cell: ")
    for line in my_model.get_lat():
        print("    ",line)
    print("   Wannier centers: ")
    for line in my_model._orb:
        print("    ",line)
    #print("   Supercell: TODO")


def plot_transfer_integrals(my_model,wannier_model, min_hopping, max_distance, ax):
    dist = np.zeros(len(my_model._hoppings))
    ham = np.zeros(len(my_model._hoppings))
    unitcell = np.array(my_model.get_lat())
    dist_raw,ham_raw=wannier_model.dist_hop()
    for i, (ti, destination, origin, vec) in enumerate(my_model._hoppings):
        ham[i] = np.abs(ti)
        r1 = unitcell.T @ my_model.get_orb()[destination]
        r2 = unitcell.T @ my_model.get_orb()[origin] + unitcell.T @ vec
        dist[i] = np.sqrt(np.dot(r1-r2, r1-r2))

    #plot hopping terms as a function of distance on a log scale
    ax.scatter(dist_raw,np.log10(np.abs(ham_raw)), label="all (Wannier90)")
    ax.hlines(np.log10(min_hopping), 0, np.max(dist_raw))
    if max_distance is not None:
        ax.vlines(max_distance, np.min(np.log10(ham_raw)), np.max(np.log10(ham_raw)))
    ax.scatter(dist,np.log10(np.abs(ham)), color='C1', label="used (TB model)")
    ax.set_xlabel("Distance (A)")
    ax.set_ylabel(r"$\log_{10} ϵ_{ij}$ (eV)")
    ax.legend()
    # fig.tight_layout()
    #plt.show()

def plot_bandstructure(val_model, val_wannier_model, cond_model, cond_wannier_model):
    # plot band structure
    fig2, ax2 = plt.subplots()
    fig2.suptitle("Electronic structure")
    (val_w90_kpt,val_w90_evals)=val_wannier_model.w90_bands_consistency()
    for i in range(val_w90_evals.shape[0]):
        ax2.plot(list(range(val_w90_evals.shape[1])),val_w90_evals[i],"k-",zorder=-100)
    (cond_w90_kpt,cond_w90_evals)=cond_wannier_model.w90_bands_consistency()
    for i in range(cond_w90_evals.shape[0]):
        ax2.plot(list(range(cond_w90_evals.shape[1])),cond_w90_evals[i],"k-",zorder=-100)

    # now interpolate from the model on the same path in k-space
    val_evals=val_model.solve_all(val_w90_kpt)
    cond_evals=cond_model.solve_all(cond_w90_kpt)
    for i in range(val_evals.shape[0]):
        ax2.plot(list(range(val_evals.shape[1])),val_evals[i],"C0-",zorder=-50)

    for i in range(cond_evals.shape[0]):
        ax2.plot(list(range(cond_evals.shape[1])),cond_evals[i],"C1-",zorder=-50)

    # for legend
    ax2.plot([0],[0], color="k", label="Wannier90")
    ax2.plot([0],[0], color="C0", label="valence TB")
    ax2.plot([0],[0], color="C1", label="conduciton TB")

    ax2.set_xlim(0,val_evals.shape[1]-1)
    ax2.set_xlabel("k-path from Wannier90")
    ax2.set_ylabel("Band energy (eV)")
    ax2.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    ax2.legend()
    fig2.tight_layout()

    #plt.show()

def openW90Files(path, seedname, min_hopping, ignorable_imaginary_part, no_ws_dist):
    # Open w90 model
    print(f"   - read {path}")
    wannier_model=w90(path,seedname)
    (dist_raw,ham_raw)=wannier_model.dist_hop()
    print("   - number of transfer integrals (original): ", len(ham_raw))
    print(f"   - enforce minimum hopping of {min_hopping}eV")
    if ignorable_imaginary_part:
        print("   - ignore imaginary part of transfer integrals")

    # get tb model in which some small terms are ignored
    my_model=wannier_model.model(zero_energy=0.0,min_hopping_norm=min_hopping,max_distance=None,ignorable_imaginary_part=ignorable_imaginary_part)

    # ws_dist corrections
    if not no_ws_dist:
        ws_dist_file = os.path.join(path, f'{seedname}_wsvec.dat')
        print("   - use ws_distance corrections from file ", ws_dist_file)
        Rvec_dict, ws_num_deg, ws_dist = read_wsvec(ws_dist_file)
        for i, (ti, destination, origin, vec) in enumerate(my_model._hoppings):
            # TODO implement degenerates using ws_num_deg

            # update TI
            T = ws_dist[:,0,destination,origin, Rvec_dict[(vec[0],vec[1],vec[2])]]
            my_model._hoppings[i] = [ti, destination, origin, vec + T]


    return my_model, wannier_model


def calcShifts(my_model, pos_posfile):
    # go through all positions and determine shifts
    pos_d = my_model.get_orb()  # in direct coordinates
    shift = np.zeros(pos_d.shape, dtype=int)
    unitcell_tb = np.array(my_model.get_lat())
    inv_unitcell = np.linalg.inv(unitcell_tb.T)

    if np.any(pos_posfile.shape != pos_d.shape):
        raise NameError("Shapes of position arrays from TB-Model and POSFILE are not compatible.")

    for i, p in enumerate(pos_d):
        shift[i,:] = np.round(p - inv_unitcell @ pos_posfile[i])

    return shift

def pruneModel(my_model, max_distance, limitR):
    # enforce maximal distance
    if max_distance is not None:
        print("   - enforce maximal distance of transfer integrals")
        acceptedHoppings = np.zeros(len(my_model._hoppings), dtype=np.bool8)
        unitcell = np.array(my_model.get_lat())
        for i, (ti, destination, origin, vec) in enumerate(my_model._hoppings):
            ham = np.abs(ti)
            r1 = unitcell.T @ my_model.get_orb()[destination]
            r2 = unitcell.T @ my_model.get_orb()[origin] + unitcell.T @ vec
            dist = np.sqrt(np.dot(r1-r2, r1-r2))

            if dist <= max_distance:
                acceptedHoppings[i] = True

        my_model._hoppings = np.array(my_model._hoppings)
        my_model._hoppings = my_model._hoppings[acceptedHoppings]
    else:
        print("   - skip maximal distance of transfer integrals")


    # enforce limit R
    if limitR is not None:
        print("   - enforce R limit.")
        acceptedHoppings = np.zeros(len(my_model._hoppings), dtype=np.bool8)
        for i, (ti, destination, origin, vec) in enumerate(my_model._hoppings):
            if abs(vec[0]) <= limitR[0] and abs(vec[1]) <= limitR[1] and abs(vec[2]) <= limitR[2]:
                acceptedHoppings[i] = True

        my_model._hoppings = np.array(my_model._hoppings)
        my_model._hoppings = my_model._hoppings[acceptedHoppings]
    else:
        print("   - skip R limit.")

    return my_model


def cleanUpModel(my_model):
    # remove dublicates / create TB model again from scratch to make sure that everything is perfect
    my_model_ = tb_model(my_model._dim_k, my_model._dim_r,my_model.get_lat(), my_model.get_orb())
    my_model_.set_onsite(my_model._site_energies,mode="reset")
    for h in my_model._hoppings:
        try:
            my_model_.set_hop(h[0],h[1],h[2],h[3], mode="reset", allow_conjugate_pair=False)
        except:
            print(f"Ignoring TI={h} (conjugate pair already set!)")
    return my_model_  # update TB model

def remove_unnecessary_wf(myModel_raw, rm_wf):
    print("You want to remove/rename Wannier functions from the tight-binding model.")
    print("This is an experimental feature! Please make sure that you know what you are doing.")
    print("Especially, make sure that the POSFILE and used *.xsf are in agreement with the new indexes / removal.")
    print("There are no checks or savegards!")

    orbitals_raw = myModel_raw.get_orb()
    N = len(orbitals_raw) - len(rm_wf)

    orbitals = np.zeros((N,3))
    mapping = {}
    j = 0
    for i in np.arange(len(orbitals_raw)):
        if i in rm_wf:
            print(i+1, " (removed)")
            continue
        else:
            print(i+1, "-->", j+1)
            mapping[i] = j
            orbitals[j] = orbitals_raw[i]
            j += 1

    myModel = tb_model(3,3,myModel_raw.get_lat(), orbitals)

    for h in range(len(myModel_raw._hoppings)):
        #print(myModel._hoppings[h])
        if myModel_raw._hoppings[h][1] in rm_wf or myModel_raw._hoppings[h][2] in rm_wf:
            continue
        else:
            value = myModel_raw._hoppings[h][0]
            i = mapping[myModel_raw._hoppings[h][1]]
            j = mapping[myModel_raw._hoppings[h][2]]
            R = myModel_raw._hoppings[h][3]

            # print("add hopping", value, i,j,R)

            myModel.set_hop(value, i, j, R,"set",True)

    onsite_raw = myModel_raw._site_energies
    for i in np.arange(len(orbitals_raw)):
        if i in rm_wf:
            continue
        else:
            myModel.set_onsite(onsite_raw[i],mapping[i], "reset")

    return myModel

if __name__=="__main__":

    parser = argparse.ArgumentParser(description=__doc__, epilog="Have fun :D")
    parser.add_argument('path', type=dir_path, help='path to Coulomb calculation')
    parser.add_argument('--no_plots', help="Turn off plots", action="store_true", default=False)


    # shifts and corrections
    g2 = parser.add_argument_group('Correction options')
    g2.add_argument('--no_ws_dist', help="Turn off ws_dist corrections (using wannier90_wsvec.dat)", action="store_true", default=False)

    # settings for pruning the tight-binding model
    g3 = parser.add_argument_group('Pruning options')
    g3.add_argument('--ignorable_imaginary_part', type=float, help="threshold for minimal imaginary part of transfer integral [eV] (default=None)", default=None)
    g3.add_argument('--min_hopping', type=float, help="threshold for minimal transfer integral [eV] (default=0.001)", default=0.001)
    g3.add_argument('--max_distance', type=float, help="maximal distance for a transfer integral [A] (default=None)", default=None)
    g3.add_argument('--limitR', type=int, nargs=3, help="Restrict the maximal shift R (positive and negative) (default=None)", default=None)

    # use optimizer to obtain min_hopping
    g4 = parser.add_argument_group('Optimizer options (also pruning)')
    g4.add_argument('--optimize', help="Use an optimizer to select cut-off criterion for minimal hopping values", action="store_true", default=False)
    g4.add_argument('--opt_upper_bound', type=float, help="maximimal value of min_hopping for the optimizer (default=0.01)", default=0.01)
    g4.add_argument('--opt_sparsity_penalty', type=float, help="sparsity penalty for the optimizer (default=0.05)", default=0.05)
    g4.add_argument('--opt_grid', type=int, nargs=3, help="k-grid for the optimization (default=[5,5,5])", default=[5,5,5])

    args = parser.parse_args()

    print("Use output directory of Coulomb calculation: ", args.path)
    vmapping_path = os.path.join(args.path, "vmapping.txt")
    cmapping_path = os.path.join(args.path, "cmapping.txt")
    POSFILE_path = os.path.join(args.path, "POSFILE")
    if not os.path.exists(vmapping_path):
        raise FileNotFoundError(vmapping_path)

    if not os.path.exists(cmapping_path):
        raise FileNotFoundError(cmapping_path)

    if not os.path.exists(POSFILE_path):
        raise FileNotFoundError(POSFILE_path)

    print("[+] Read mapping files.")
    print(f"   - read {vmapping_path}")
    with warnings.catch_warnings():  # got some deprecated warnings that are irrelevant
        warnings.simplefilter("ignore")
        warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
        vpaths = np.loadtxt(vmapping_path,str, comments="#", delimiter=">",usecols=1, ndmin=1)

    if len(vpaths) == 0:
        raise NameError("Cannot parse vmappings.")
    valence_dir = os.path.dirname(vpaths[0]).strip()
    valence_seedname = os.path.basename(vpaths[0]).split("_0")[0]
    print(f"   - found w90 calculation for valence states: {valence_dir}, seedname: {valence_seedname}")

    print(f"   - read {cmapping_path}")
    with warnings.catch_warnings():  # got some deprecated warnings that are irrelevant
        warnings.simplefilter("ignore")
        warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
        cpaths = np.loadtxt(cmapping_path,str, comments="#", delimiter=">",usecols=1, ndmin=1)
    if len(cpaths) == 0:
        raise NameError("Cannot parse cmappings.")
    conduction_dir = os.path.dirname(cpaths[0]).strip()
    conduction_seedname = os.path.basename(cpaths[0]).split("_0")[0]
    print(f"   - found w90 calculation for conduction states: {conduction_dir}, seedname: {conduction_seedname}")

    # read in tight-binding models from wannier90 output files
    print("[+] Read Wannier Hamiltonian from w90 calculations.")
    val_model, val_wannier_model = openW90Files(valence_dir, valence_seedname, args.min_hopping, args.ignorable_imaginary_part, args.no_ws_dist)
    cond_model, cond_wannier_model = openW90Files(conduction_dir, conduction_seedname, args.min_hopping, args.ignorable_imaginary_part, args.no_ws_dist)

    if np.any(np.abs(np.array(val_model.get_lat()) - np.array(cond_model.get_lat())) > 1e-5):
        print("valence and conduction Wannier functions have different unit cells!")
        print("valence:")
        print(val_model.get_lat())
        print("conduction:")
        print(val_model.get_lat())
        raise NameError("Wannier90 calculations are not compatible to each other!")

    # get some basic information
    # print("\n--> Obtained tight-binding model")
    # print_summary(my_model)

    # custom changes to the TB-model (remove/rename Wannier functions)
    # rm_wf = np.array([ 8, 12, 13, 14, 15], dtype=np.int16)-1
    # cond_model = remove_unnecessary_wf(cond_model, rm_wf)


    # read POSFILE
    print(f"[+] Read POSFILE from {POSFILE_path}")
    unitcell, positions_val, positions_cond = read_posfile(POSFILE_path)

    print("   - check compatibility with w90 calculations")
    if np.any(np.abs(np.array(val_model.get_lat()) - unitcell) > 1e-5):
        print("POSFILE has different unit cell than Wannier90 calculations!")
        print("Wannier90:")
        print(val_model.get_lat())
        print("POSFILE:")
        print(unitcell)
        raise NameError("POSFILE is not compatible with Wannier90 calculations!")


    # check if number of WF agrees with POSFILE
    if (len(positions_val) != val_model.get_num_orbitals()):
        raise NameError("Different number of valence Wannier functions in POSFILE and Wannier90.")
    if (len(positions_cond) != cond_model.get_num_orbitals()):
        raise NameError("Different number of conduction Wannier functions in POSFILE and Wannier90.")

    unitcell_tb = np.array(val_model.get_lat())
    print(f"[+] Make tight-binding model compatible with POSFILE (compatible with Coulomb calculations)")
    print("   - shift valence Wannier functions")
    val_shift = calcShifts(val_model, positions_val)
    # print(val_shift)
    shift_sites(val_model, val_shift)
    print("   - difference to POSFILE (after shift):")
    pos_d = val_model.get_orb()  # update coordinates
    for i, p in enumerate(pos_d):
        diff = unitcell_tb.T @ p - positions_val[i]
        print(f"    {i}: \t{diff} \tnorm = {np.linalg.norm(diff)}")

    print("   - shift conduction Wannier functions")
    cond_shift = calcShifts(cond_model, positions_cond)
    # print(cond_shift)
    shift_sites(cond_model, cond_shift)
    print("   - difference to POSFILE (after shift):")
    pos_d = cond_model.get_orb()  # update coordinates
    for i, p in enumerate(pos_d):
        diff = unitcell_tb.T @ p - positions_cond[i]
        print(f"    {i}: \t{diff} \tnorm = {np.linalg.norm(diff)}")


    print("[+] Additional pruning of the tight-binding model")
    print("   - valence states")
    val_model = pruneModel(val_model, args.max_distance, args.limitR)
    print("   - conduction states")
    pruneModel(cond_model, args.max_distance, args.limitR)

    print("[+] Clean up TB-model.")
    print("   - valence states")
    val_model = cleanUpModel(val_model)
    print("   - conduction states")
    cond_model = cleanUpModel(cond_model)



    # optimize min_hopping
    if args.optimize:
        print("[+] Optimize min_hopping.")
        print(f"   - try to find optimal cut-off value for min_hopping")
        print(f"   - searching interval (eV) [{args.min_hopping}, {args.opt_upper_bound}]")
        print(f"   - sparsity penalty: ", args.opt_sparsity_penalty)
        print(f"   - k-grid for optimization: ", args.opt_grid)

        print("   - optimize valence states ...")
        val_model = optimize_tb(val_model, args.min_hopping, args.opt_upper_bound, args.opt_grid, args.opt_sparsity_penalty)
        print("   - optimize conduction states ...")
        cond_model = optimize_tb(cond_model, args.min_hopping, args.opt_upper_bound, args.opt_grid, args.opt_sparsity_penalty)


    # print("--> Final tight-binding model:")
    # print_summary(val_model)
    # print_summary(cond_model)


    # plot for control
    if not args.no_plots:
        print("[+] Plot results.")
        fig = plt.figure(figsize=(10,5))
        ax1= fig.add_subplot(121)
        ax2= fig.add_subplot(122)

        ax1.set_title("Valence (hole) transfer integrals")
        plot_transfer_integrals(val_model,val_wannier_model, args.min_hopping, args.max_distance, ax1)

        ax2.set_title("Conduction (electron) transfer integrals")
        plot_transfer_integrals(cond_model,cond_wannier_model, args.min_hopping, args.max_distance, ax2)

        fig.tight_layout()
        plt.show(block=False)

        plot_bandstructure(val_model, val_wannier_model, cond_model, cond_wannier_model)
        plt.show(block=False)

    # output wannier-optics files
    print("[+] Write wannier-optics input files.")
    write_tinfile(os.path.join(args.path,"TINFILE_v"), val_model)
    write_onsite(os.path.join(args.path,"ONSITE_ENERGY_v"), val_model)
    # write_posfile("POSFILE_v", val_model)

    write_tinfile(os.path.join(args.path,"TINFILE_c"), cond_model)
    write_onsite(os.path.join(args.path,"ONSITE_ENERGY_c"), cond_model)
    # write_posfile("POSFILE_c", cond_model)

    write_parfile(os.path.join(args.path,"PARFILE"))

    # print("\nThe usual convention is:")
    # print("\t\t valence states \t==> _v (outer coordinates)")
    # print("\t\t conduction states \t==> _c (inner coordinates)")
    # print("Please notice that the valence (outer) transfer integrals will be multiplied with (-1) in the exciton hamiltonian!")

    print("Finish.")
    if not args.no_plots:
        plt.show()