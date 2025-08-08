#ifndef CHGCAR_H
#define CHGCAR_H

/**
 * @file CHGCAR.h
 * @author Konrad Merkel
 * @brief Data structures and functionalities to deal with CHGCAR data from a vasp simulation.
 */


#include "meshgrid.h"
#include "algebra.h"
#include "wannierfunction.h"

using namespace std;

class CHGCAR
{
private:
    unique_ptr<double[],free_deleter> value{};    //!< data grid that stores the actual function (3D array)
    shared_ptr<RealMeshgrid> meshgrid{};          //!< real space grid where values are stored (in contrast to WannierFunction this instance is unique and not used by other instances)


public:
    CHGCAR(shared_ptr<RealMeshgrid> newMeshgrid, unique_ptr<double[],free_deleter> newValue)
    : value{std::move(newValue)}, meshgrid{newMeshgrid}
    {}

    void reset() {
        value.reset();
        meshgrid.reset();
    }

    // getter and setter
    RealMeshgrid const* getMeshgrid() const {return meshgrid.get(); }
    double getNumElectrons() const {

        double num_electron=0;
        double dV = this->meshgrid->getdV();
        int N = this->meshgrid->getNumDataPoints();

        // integrate over whole density
        for (int i=0; i<N; i++) {
            num_electron += value[i] * dV;
        }
        return num_electron;
    }
    double const* getValue() const { return value.get(); } //!< Returns all data points

    void coarseGrain(uint fx, uint fy, uint fz) {

        if ((fx==0) || (fy==0) || (fz==0))
            throw runtime_error("coarseGrain cannot divide by zero");

        // get new dimensions (number of grid points)
        vector<int> oldDim = meshgrid->getDim();
        vector<int> newDim(3);
        newDim[0] = oldDim[0] / fx;
        newDim[1] = oldDim[1] / fy;
        newDim[2] = oldDim[2] / fz;

        // TODO: check if fx is a divisor of oldDim

        shared_ptr<RealMeshgrid> newMeshgrid = make_shared<RealMeshgrid>(newDim,meshgrid->getLattice(),meshgrid->getOrigin());


        // generate new data array
        double *newValue = (double *)malloc(sizeof(double) * newMeshgrid->getNumDataPoints());

        // go through every grid point
        for (int k = 0; k < newDim[2]; k++) { // z
            int kk = fz*k;  // index in the old (denser) grid
            for (int j = 0; j < newDim[1]; j++) { // y
                int jj = fy*j;
                for (int i = 0; i < newDim[0]; i++) { // x
                    int ii = fx*i;
                    newValue[i + newDim[0] * (j + newDim[1] * k)] = value[ii + oldDim[0] * (jj + oldDim[1] * kk)];
                }
            }
        }


        if (meshgrid->hasMeshgridArrays())
            newMeshgrid->createMeshgridArrays();

        // update meshgrid and values
        meshgrid = newMeshgrid;
        value.reset(newValue);
    }

    /**
     * @brief Expands the charge density in a supercell.
     *
     * The charge density (value) is repeated in every primitive cell troughout the supercell.
     * This is the main difference to a supercell of a Wannier function.
     *
     * @param nx
     * @param ny
     * @param nz
     */
    void createSupercell(int nx, int ny, int nz, vector<double> origin)
    {
        if ((nx <= 1) && (ny <= 1) && (nz <= 1)) { return; }

        // cout << "TEST: create supercell CHG: " << nx << "  " << ny << "  " << nz << "\n";

        vector<vector<double>> unitcell = meshgrid->getLattice();
        // vector<vector<double>> old_lattice = meshgrid->getLattice();

        // get new dimensions (number of grid points)
        vector<int> newDim(3);
        newDim[0] = meshgrid->getDim()[0] * nx;
        newDim[1] = meshgrid->getDim()[1] * ny;
        newDim[2] = meshgrid->getDim()[2] * nz;

        // calculate new (larger) lattice of the supercell (not unitcell!)
        vector<vector<double>> new_lattice(3);
        for (int i=0;i<3; i++)
            new_lattice[i] = vector<double>(3);
        for (int j = 0; j < 3; j++)
        {
            new_lattice[0][j] = unitcell[0][j] *nx;
            new_lattice[1][j] = unitcell[1][j] *ny;
            new_lattice[2][j] = unitcell[2][j] *nz;
        }


        shared_ptr<RealMeshgrid> newMeshgrid = make_shared<RealMeshgrid>(newDim, new_lattice, origin);


        // generate new data array
        double *newValue = (double *)malloc(sizeof(double) * newMeshgrid->getNumDataPoints());

        vector<int> oldDim = meshgrid->getDim();

        // go through every unitcell in the supercell
        for (int Rx=0; Rx<nx; Rx++) {
            for (int Ry=0; Ry<ny; Ry++) {
                for (int Rz=0; Rz<nz; Rz++) {

                    // go through every grid point in the primitive unitcell
                    for (int k = 0; k < oldDim[2]; k++) { // z
                        int kk = k+oldDim[2]*Rz;  // grid point in the new supercell

                        for (int j = 0; j < oldDim[1]; j++) { // y
                            int jj = j+oldDim[1]*Ry;  // grid point in the new supercell

                            for (int i = 0; i < oldDim[0]; i++) { // x
                                int ii = i+oldDim[0]*Rx;  // grid point in the new supercell

                                newValue[ii + newDim[0] * (jj + newDim[1] * kk)] = value[i + oldDim[0] * (j + oldDim[1] * k)];
                            }
                        }
                    }
                }
            }
        }

        if (meshgrid->hasMeshgridArrays())
            newMeshgrid->createMeshgridArrays();

        // update meshgrid and values
        meshgrid = newMeshgrid;
        value.reset(newValue);
    }

    void createSupercell(int nx, int ny, int nz) {
        this->createSupercell(nx,ny,nz,this->meshgrid->getOrigin());
    }

    /**
     * @brief Checks if two Wannier functions are compatible / can be used together
     *
     * Wannier functions are compatible if they are defined on the same real space grid.
     * The RealMeshgrid objects of both Wannier functions needs to be the same or a copy.
     *
     * @param otherWannier
     * @return true     compatible
     * @return false    not compatible
     */
    bool isCompatible(WannierFunction const& otherWannier) const
    {
        if (*(meshgrid.get()) != *(otherWannier.getMeshgrid()))
            return false;

        return true;
    }

    bool makeCompatible(WannierFunction const& wf) {

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // usually the CHGCAR file has double as many points in the SUPERCELL than the Wannier funcitons
        // with respect to the unit cell this might be not an integer multiple
        vector<double> supercell = wf.getLatticeInUnitcellBasis();
        this->createSupercell(round(supercell[0]), round(supercell[1]), round(supercell[2]), wf.getMeshgrid()->getOrigin());

        vector<int> dim_wf = wf.getMeshgrid()->getDim();
        vector<int> dim_chg = this->meshgrid->getDim();

        if (rank==0) {
            cout << "Dim WF: " << dim_wf[0] << "\t" << dim_wf[1] << "\t" << dim_wf[2] << "\n";
            cout << "Dim CHGCAR: " << dim_chg[0] << "\t" << dim_chg[1] << "\t" << dim_chg[2] << "\n";
            cout << "Coarse graining factor: " << dim_chg[0]/dim_wf[0] << "\t" << dim_chg[1]/dim_wf[1] << "\t" << dim_chg[2]/dim_wf[2] << "\n";

        }
        this->coarseGrain(round(dim_chg[0]/dim_wf[0]) ,round(dim_chg[1]/dim_wf[1]), round(dim_chg[2]/dim_wf[2]));


        // cout << "\n\n\n\n CHGCAR:" << endl;
        // this->getMeshgrid()->printSummary();

        // cout << "\n\n\n\n WF:" << endl;
        // wf->getMeshgrid()->printSummary();



        return isCompatible(wf);
    }
};


double calcExpectationDensity(CHGCAR const& chg, WannierFunction const& wf)
{
    if (!chg.isCompatible(wf))
        throw runtime_error("CHGCAR is not compatible with Wannier funciton.");

    int N = chg.getMeshgrid()->getNumDataPoints();
    double dV = chg.getMeshgrid()->getdV();

    const double* n_GS = chg.getValue();
    const double* phi = wf.getValue();

    double meanDens = 0.0;
    for (int i=0; i<N; i++) {
        meanDens += phi[i] * n_GS[i] * phi[i] * dV;
    }
    return meanDens;
}



#endif // CHGCAR_H