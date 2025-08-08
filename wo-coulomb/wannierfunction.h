#ifndef WANNIERFUNCTION_H
#define WANNIERFUNCTION_H

/**
 * @file wannierfunction.h
 * @author Konrad Merkel
 * @brief Data structures and functionalities for single Wannier functions on a real space grid.
 */

#include "meshgrid.h"
#include "algebra.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <omp.h>


using namespace std;

/**
 * @brief Monopole moment of a charge density
 */
struct Monopole
{
    double x,y,z;   //!< cartesian coordinates of charge centre
    double charge;  //!< total charge
};

/**
 * @brief Single Wannier function on a real space grid
 *
 * This is on of the central objects of the code. It stores a single Wannier function \f$w_{m0}(x)\f$
 * within a supercell on a real space grid. For every m (variable 'id' in the code) you would need a single
 * WannierFunction object \f$w_{m0}(x)\f$. Wannier functions for different unit cell vectors \f$R\f$
 * can always be obtained by using the shifting property of Wannier functions, i.e.
 * \f$w_{mR}(x) = w_{m0}(x-R)\f$. This is implemented in the getShiftedValue() method.
 *
 * To read and write Wannier function in the *.xsf file format please use the XSF_controller class.
 *
 */
class WannierFunction
{
private:
    int id;                               //!< index / quantum number
    unique_ptr<double[], free_deleter> value{};         //!< data grid that stores the actual function (3D array)
    shared_ptr<RealMeshgrid> meshgrid{};  //!< real space grid where values are stored
    vector<vector<double>> unitcell{};    //!< primitive unit cell (3x3 matrix)
    vector<double> latticeInUnitcellBasis{};  //!< dimensions of the supercell in terms of unit cells
    vector<double> originInUnitcellBasis{};   //!< origin of the supercell in terms of unit cells

protected:

    /**
     * @brief updates the private variable latticeInUnitcellBasis and originInUnitcellBasis
     *
     * This function is called every time when meshgrid or unitcell is changing.
     *
     */
    void updateLatticeInUnitcellBasis() {
        latticeInUnitcellBasis = vector<double>{0., 0., 0.};
        originInUnitcellBasis  = vector<double>{0., 0., 0.};

        vector<double> origin = meshgrid->getOrigin();
        vector<vector<double>> lattice = meshgrid->getLattice();
        vector<vector<double>> invUnitcell = invMat3x3(unitcell);

        for (int i = 0; i < 3; i++)
        {
            for (int n = 0; n < 3; n++)
            {
                latticeInUnitcellBasis[i] += invUnitcell[n][i] * lattice[i][n];
                originInUnitcellBasis[i]  += invUnitcell[n][i] * origin[n];
            }
        }
    }

public:
    /**
     * @brief Creates an empty Wannier function
     *
     * @param newID index / quantum number
     */
    WannierFunction(int newID=0) : id(newID), value{}, meshgrid{}, latticeInUnitcellBasis{}, originInUnitcellBasis{}
    {}

    /**
     * @brief Creates a Wannier function
     *
     * @param newID       index / quantum number
     * @param newMeshgrid real space grid where values are stored
     * @param newValue    data grid that stores the actual function (3D array)
     * @param newUnitcell primitive unit cell (3x3 matrix)
     */
    WannierFunction(int newID, shared_ptr<RealMeshgrid> newMeshgrid, unique_ptr<double[], free_deleter> newValue, vector<vector<double>> const& newUnitcell)
        : id{newID}, value{std::move(newValue)}, meshgrid{newMeshgrid}, unitcell{newUnitcell}
    {
        if (unitcell.size() != 3)
            throw runtime_error("Unitcell needs to be a 3x3 matrix!");

        for (int i = 0; i < 3; i++)
        {
            if (unitcell[i].size() != 3)
                throw runtime_error("Unitcell needs to be a 3x3 matrix!");
        }
        updateLatticeInUnitcellBasis();
    }


    // getter and setter
    RealMeshgrid* getMeshgrid() const {return meshgrid.get(); }
    shared_ptr<RealMeshgrid> getSharedMeshgridPtr() const {return meshgrid; }
    vector<vector<double>> getUnitcell() const { return unitcell; }
    int getId() const { return this->id; }
    double getVunitcell() const { return det3x3(unitcell); } //!< Volume of the supercell

    /**
     * @brief Calculates the norm of the Wannier function
     */
    double getNorm() const
    {
        double norm = 0.0;
        double dV = this->meshgrid->getdV();
        for (int i = 0; i < this->meshgrid->getNumDataPoints(); i++)
        {
            norm += value[i] * value[i] * dV;
        }
        return norm;
    }

    double const* getValue() const { return value.get(); } //!< Returns all data points  // TODO replace by operator(size_t i)
    double& operator()(int i) { return value[i]; }
    double& operator()(int i, int j, int k) { return value[meshgrid->getGlobId(i,j,k)]; }

    /**
     * @brief Get dimensions of the lattice (supercell) in terms of unit cells
     */
    vector<double> getLatticeInUnitcellBasis() const { return this->latticeInUnitcellBasis;}

    /**
     * @brief Get origin of the lattice (supercell) in terms of unit cells
     */
    vector<double> getOriginInUnitcellBasis() const { return this->originInUnitcellBasis; }

    /**
     * @brief Returns the value of a shifted Wannier function
     *
     * Returns the value n of a shifted Wannier function. R is the shift
     * vector in units of unitcell vectors.
     * Please note that this is not the value at 'n+R' but rather 'n-R'.
     *
     * DO NOT USE THIS ROUTINE WITHIN A LOOP BECAUSE indexOffset WILL THEN BE CALCULATED
     * EVERY TIME. INSTEAD USE
     *      RealMeshgrid::getIndexOffset()
     * OUTSIDE THE LOOP AND
     *      RealMeshgrid::getShiftedGlobalIndex()
     * INSIDE THE LOOP FOR BETTER PERFORMANCE.
     *
     * @param n       global index of the data array
     * @param R       shift vector in units of lattice vectors
     * @return double value of that grid point
     */
    double getShiftedValue(const int n, const vector<int>& R) const
    {
        vector<int> indexOffset = meshgrid->getIndexOffset(R, this->getLatticeInUnitcellBasis());
        int globId = meshgrid->getShiftedGlobalIndex(n,indexOffset);

        if (globId < 0) return 0.0;
        else return value[globId];
    }



    // vector<double> getDipol() const {
    //   vector<double> dipol{0.,0.,0.};
    //   Monopole mono = this->getMonopole();

    //   int Npoints = meshgrid->getNumDataPoints();
    //   double dV = meshgrid->getdV();
    //   const double* dens = this->getDensity();

    //   // calculate charge center
    //   double x,y,z;
    //   for (int j=0; j<Npoints; j++){
    //       meshgrid->xyz(j,x,y,z);
    //       dipol[0] += dens[j] * (x-mono.x) * dV / mono.charge;
    //       dipol[1] += dens[j] * (y-mono.y) * dV / mono.charge;
    //       dipol[2] += dens[j] * (z-mono.z) * dV / mono.charge;
    //   }

    //   return dipol;
    // }

    void setValue(double *newValue) { value.reset(newValue); }   //!< Assigns new data points (on the same real space grid) to WF object
    void setMeshgrid(shared_ptr<RealMeshgrid> newMeshgrid) {  //!< Assigns new real space mesh grid
        this->meshgrid = newMeshgrid;
        updateLatticeInUnitcellBasis();
    }
    void setUnitcell(vector<vector<double>> const& newUnitcell) {  //!< Set new primitive unit cell
        this->unitcell = newUnitcell;
        updateLatticeInUnitcellBasis();
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

        if (*(meshgrid) != *(otherWannier.getMeshgrid()))
            return false;

        // check if unit cell is the same
        vector< vector<double> > otherLattice = otherWannier.getUnitcell();
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                if (abs(otherLattice[i][j] - this->unitcell[i][j])> 1e-8) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Prints a small summary in std output
     *
     */
    void printSummary() const
    {
        cout << fixed;
        cout << setprecision(12);
        cout << "Wannierfunction : " << id << endl;

        cout << "Meshgrid:" << endl;
        meshgrid->printSummary();

        cout << endl
             << "Unit cell: " << endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                cout << "\t" << unitcell[i][j];
            }
            cout << endl;
        }

        cout << "Supercell\t: ";
        for (int i = 0; i < 3; i++)
        {
            cout << latticeInUnitcellBasis[i] << "  ";
        }
        cout << endl;
        cout << "Origin (in terms of unit cells)\t: ";
        for (int i = 0; i < 3; i++)
        {
            cout << originInUnitcellBasis[i] << "  ";
        }
        cout << endl;

        cout << "Vunitcell\t: " << getVunitcell() << endl;
        cout << "Norm\t\t: " << getNorm() << endl;
        cout << "Created mesh arrays\t\t: " << meshgrid->hasMeshgridArrays() << endl;
        // Monopole mono = getMonopole();
        // cout << "Charge center\t\t: " << mono.x << "  " << mono.y << "  " << mono.z << "  "  << endl;
        // cout << "Support radius\t\t: " << getSupportRadius();

        cout << "\nHead value (first 10 values): " << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << value[i] << "\t";
        }
        cout << "..." << endl;
    }

    /**
     * @brief Normalizes the Wannier function to one.
     */
    void normalize(double charge=1.0)
    {
        double scale = sqrt(getNorm());

        if (abs(scale - sqrt(abs(charge))) < 1e-10) // already normalized
            return;

        int N = meshgrid->getNumDataPoints();
        for (int i = 0; i < N; i++)
        {
            value[i] = value[i] / scale * sqrt(abs(charge));
        }
    }
};


/**
 * @brief Checks if all WannierFunctions in the maps are compatible to each other
 * and makes them use the same RealMeshgrid instance.
 *
 * This is useful to save memory since Meshgrid information are only stored once.
 *
 * @param vWannMap WannierFunctions of the valence bands
 * @param cWannMap WannierFunctions of the conduction bands
 * @return true   WannierFunctions are compatible and share the same RealMeshgrid
 * @return false  WannierFuncitons are not compatible. Nothing has been changed.
 */
bool tryToShareMeshgrid( map<int, WannierFunction>& vWannMap, map<int, WannierFunction>& cWannMap)
{
    WannierFunction const& firstWann = vWannMap.begin()->second;
    shared_ptr<RealMeshgrid> mesh = firstWann.getSharedMeshgridPtr();

    // check if Wannier functions and their underlying meshgrids are compatible to each other
    // if not there is nothing we can do!
    auto isCompatible = [&firstWann](std::pair<const int, WannierFunction> const& p){
        return p.second.isCompatible(firstWann);
    };
    if (! all_of(vWannMap.begin(), vWannMap.end(), isCompatible)) return false;
    if (! all_of(cWannMap.begin(), cWannMap.end(), isCompatible)) return false;

    // set the same meshgrid for every Wannier function
    auto setMesh = [&mesh](std::pair<const int, WannierFunction>& p){ p.second.setMeshgrid(mesh); };
    for_each(vWannMap.begin(), vWannMap.end(), setMesh);
    for_each(cWannMap.begin(), cWannMap.end(), setMesh);

    return true;
}

/**
 * @brief Creates a WannierFunction object in a larger supercell
 *
 * An existing Wannier function is placed in a larger supercell. New values
 * are filled with zeros (padding). This can be important for FFTs of Wannier
 * functions.
 *
 * @param wann    Wannier function in a smaller supercell. This will be overwritten.
 * @param supercellDim  Dimensions of the new supercell in units of the old supercell
 * @return void
 */
void createLargerSupercell(WannierFunction& wann, vector<int> supercellDim)
{

    if ((supercellDim[0] <= 1) && (supercellDim[1] <= 1) && (supercellDim[2] <= 1))
    {
        return;
    }

    const RealMeshgrid* oldMeshgrid = wann.getMeshgrid();
    // vector<vector<double>> unitcell = wann.getUnitcell();
    vector<vector<double>> old_lattice = oldMeshgrid->getLattice();
    vector<int> oldDim = oldMeshgrid->getDim();

    // get new dimensions (number of grid points)
    vector<int> newDim(3);
    newDim[0] = oldMeshgrid->getDim()[0] * supercellDim[0];
    newDim[1] = oldMeshgrid->getDim()[1] * supercellDim[1];
    newDim[2] = oldMeshgrid->getDim()[2] * supercellDim[2];

    // calculate new (larger) lattice of the supercell (not unitcell!)
    vector<vector<double>> new_lattice(3);
    for (int i = 0; i < 3; i++)
    {
        new_lattice[i] = vector<double>(3);
        for (int j = 0; j < 3; j++)
        {
            new_lattice[i][j] = old_lattice[i][j] / (oldDim[i]) * (newDim[i]);
        }
    }

    shared_ptr<RealMeshgrid> newMeshgrid = make_shared<RealMeshgrid>(newDim, new_lattice, wann.getMeshgrid()->getOrigin());

    // generate new data array
    double *data = (double *)malloc(sizeof(double) * newMeshgrid->getNumDataPoints());
    double const* olddata = wann.getValue();
    printf("Create supercell with %d threads\n", omp_get_max_threads());
    #pragma omp parallel for shared(data, olddata, newDim, oldDim)
    for (int k = 0; k < newDim[2]; k++)
    { // z
        for (int j = 0; j < newDim[1]; j++)
        { // y
            for (int i = 0; i < newDim[0]; i++)
            { // x

                if ((i < oldDim[0]) && (j < oldDim[1]) && (k < oldDim[2]))
                {
                    data[i + newDim[0] * (j + newDim[1] * k)] = olddata[i + oldDim[0] * (j + oldDim[1] * k)];
                }
                else
                {
                    data[i + newDim[0] * (j + newDim[1] * k)] = 0;
                }
            }
        }
    }

    // update meshgrid and data arrays
    wann.setMeshgrid(newMeshgrid);
    wann.setValue(data);
    return;
}

/**
 * @brief Creates a copy of a Wannier function that is shifted in units of unit cells
 *
 * Data points are wrapped around the supercell such that no data is lost. This is different
 * than the usual WannierFunction::getShiftedValue() method or joinedDensity().
 *
 * @param wann              Wannier function that is shifted
 * @param shift             Shift vector in terms of unit cells
 * @return WannierFunction New copy of Wannier function
 */
void rotateWannierFunction(WannierFunction &wann, vector<int> shift)
{
    int N = wann.getMeshgrid()->getNumDataPoints();
    const RealMeshgrid* mesh = wann.getMeshgrid();
    vector<int> dim = mesh->getDim();

    vector<double> R_cart{0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            R_cart[i] += -wann.getUnitcell()[k][i] * shift[k];
            // minus because we want to get the values at n for the shifted WF
            // and not the value at n+R
        }
    }

    double *value = (double *)malloc(N * sizeof(double)); // new data array for shifted WF
    double const* data = wann.getValue();                      // old (unshifted) values
    int x, y, z;
    for (int i = 0; i < N; i++)
    {

        mesh->getShiftedIndexTipple(i, R_cart, x, y, z);

        // wrap around the unit cell
        x = x % dim[0];
        y = y % dim[1];
        z = z % dim[2];

        if (x < 0)
            x = dim[0] + x;
        if (y < 0)
            y = dim[1] + y;
        if (z < 0)
            z = dim[2] + z;

        int globId = mesh->getGlobId(x, y, z);

        if ((globId >= N) || (globId < 0))
        {
            cerr << x << " " << y << " " << z << endl;
            cerr << globId << endl;
            throw runtime_error("Problems while wraping around the supercell in rotateWannierFunction() (Invalid global ID!)");
        }

        value[i] = data[globId];
    }

    wann.setValue(value);

    // return new WannierFunction(wann->getId(), wann->getMeshgrid(), value, wann->getUnitcell());
}



#endif  // WANNIERFUNCTION_H