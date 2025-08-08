#ifndef MESHGRID_H
#define MESHGRID_H

/**
 * @file meshgrid.h
 * @author Konrad Merkel
 * @brief Meshgrids and coordinate systems
 *
 */

#include "algebra.h"

#include <iostream>
//#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <math.h>


using namespace std;

/**
 * @brief General grid object
 *
 * All grid objects store all necessary information about the supercell, discretization
 * and provide methods to calculate coordinates.
 *
 */
class Meshgrid
{
protected:
    vector<vector<double>> lattice; //!< supercell lattice
    vector<double> origin;          //!< origin of the supercell
    vector<int> dim;                //!< number of data points in each direction
    unique_ptr<double[], free_deleter>XX{};              //!< meshgrid array for x
    unique_ptr<double[], free_deleter>YY{};              //!< meshgrid array for y
    unique_ptr<double[], free_deleter>ZZ{};              //!< meshgrid array for z

public:

    Meshgrid(vector<int> newDimensions, vector<vector<double>> newLattice, vector<double> newOrigin)
     : lattice{ newLattice }, origin{ newOrigin }, dim{ newDimensions }
    {
        if (newDimensions.size() != 3)
            throw runtime_error("Number of dimensions needs to be 3 (in Meshgrid constructor)!");

        if (newLattice.size() != 3)
            throw runtime_error("Lattice needs to be a 3x3 matrix (in Meshgrid constructor)!");
        for (int i = 0; i < 3; i++)
        {
            if (newLattice[i].size() != 3)
                throw runtime_error("Lattice needs to be a 3x3 matrix (in Meshgrid constructor)!");
        }

        if (newOrigin.size() != 3)
            throw runtime_error("Origin on data field needs to have 3 dimensions (in Meshgrid constructor)!");

        //int N = newDimensions[0] * newDimensions[1] * newDimensions[2];
    }

    virtual ~Meshgrid() =default;
    /**
     * @brief Create meshgrid arrays that can be used instead of calling xyz(...) for every grid point
     *
     * The routine allocates memory and fills the arrays only if they are not already filled.
     */
    void createMeshgridArrays()
    {
        if (this->hasMeshgridArrays()) return; // arrays are already filled

        //cout << "create meshgrid arrays (increases memory consumption)\n";
        int N = this->getNumDataPoints();

        XX.reset((double *)malloc(sizeof(double) * N));
        YY.reset((double *)malloc(sizeof(double) * N));
        ZZ.reset((double *)malloc(sizeof(double) * N));

        // fill arrays
        for (int i = 0; i < N; i++)
        {
            this->xyz(i, XX[i], YY[i], ZZ[i]);
        }
    }

    bool hasMeshgridArrays() const {return (XX.get() != nullptr) && (YY.get() != nullptr) && (ZZ.get() != nullptr);}

    const double *getXX() const { return XX.get(); }
    const double *getYY() const { return YY.get(); }
    const double *getZZ() const { return ZZ.get(); }

    double getVgrid() const { return det3x3(lattice); }                                                    //!< calculates the volume of the entire supercell
    double getdV() const { return this->getVgrid() / ((this->dim[0]) * (this->dim[1]) * (this->dim[2])); } //!< discretized volume element
    vector<int> getDim() const { return this->dim; }
    vector<double> getOrigin() const { return this->origin; }
    vector<vector<double>> getLattice() const { return this->lattice; }

    /**
     * @brief Calculates cartesian coordinates from global array index
     *
     * @param n global index of a data array
     * @param x,y,z  returned cartesian coordinates
     */
    virtual void xyz(int n, double &x, double &y, double &z) const = 0;

    int inline getNumDataPoints() const { return this->dim[0] * this->dim[1] * this->dim[2]; } //!< total number of grid points

    /**
     * @brief Checks if two Meshgrid objects describe the same supercell
     */
    bool operator==(const Meshgrid &other) const
    {
        // dimensions
        vector<int> otherDim = other.getDim();
        if ((this->dim[0] != otherDim[0]) || (this->dim[1] != otherDim[1]) || (this->dim[2] != otherDim[2]))
            return false;

        // origin
        vector<double> otherOrig = other.getOrigin(); // *.xsf files have only 7 decimals
        if ((abs(this->origin[0] - otherOrig[0]) > 1e-6) || (abs(this->origin[1] - otherOrig[1]) > 1e-6) || (abs(this->origin[2] - otherOrig[2]) > 1e-6))
            return false;

        // lattice
        vector<vector<double>> otherLattice = other.getLattice();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (abs(otherLattice[i][j] - this->lattice[i][j]) > 5e-6)  // *.xsf files have only 7 decimals
                {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Checks if two Meshgrid objects do not describe the same supercell
     */
    bool operator!=(const Meshgrid &other) const
    {
        return !this->operator==(other);
    }

    /**
     * @brief Prints a small summary in std output
     *
     */
    void printSummary() const
    {
        cout << fixed;
        cout << setprecision(12);
        cout << "Cell of Meshgrid: " << endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                cout << "\t" << this->lattice[i][j];
            }
            cout << endl;
        }

        cout << "Meshgrid dimensions\t: ";
        for (int i = 0; i < 3; i++)
        {
            cout << this->dim[i] << "  ";
        }
        cout << endl;

        cout << "Meshgrid origin\t: ( ";
        for (int i = 0; i < 3; i++)
        {
            cout << this->origin[i] << " ";
        }
        cout << ")" << endl;

        cout << "Vgrid\t\t: " << this->getVgrid() << endl;
        cout << "dV\t\t: " << this->getdV() << endl;
    }

    /**
     * @brief Calculates i,j,k indexes from global array index
     *
     * i,j,k are the indexes along the different directions
     * in a data array. This is the inverse operation to getGlobId()
     *
     * The global index (data index) can be obtained from: n = i + n0*( j + n1*k )
     *
     * @param n global index (of a data array)
     * @return i,j,k array index in every direction
     * @return j array index in 2nd direction
     * @return k array index in 3rd direction
     */
    void inline getIndexTipple(const int n, int &i, int &j, int &k) const
    {
        i = (n) % (this->dim[0]);
        j = ((n - i) / this->dim[0]) % (this->dim[1]);
        k = ((n - i) - j * this->dim[0]) / (this->dim[1] * this->dim[0]);
    }

    /**
     * @brief Calculates the global index from i,j,k indexes
     *
     * i,j,k are the indexes along the different directions
     * in a data array. This is the inverse operation to getIndexTipple(..))
     *
     * @param i array index in 1st direction
     * @param j array index in 2nd direction
     * @param k array index in 3rd direction
     * @return int global index (of a data array)
     */
    int inline getGlobId(const int i, const int j, const int k) const
    {
        return i + dim[0] * (j + dim[1] * k);
    }
};

/**
 * @brief Real space grid
 *
 */
class RealMeshgrid : public Meshgrid
{
private:
    vector<vector<double>> invLattice; //!< inverse cell

public:
    RealMeshgrid(vector<int> const& newDimensions, vector<vector<double>> const& newLattice, vector<double> const& newOrigin)
        : Meshgrid(newDimensions, newLattice, newOrigin), invLattice{invMat3x3(transpose3x3(newLattice))}
        // cache for later usage in getShiftedIndexTipple(...)
        // Please note that we transpose the matrix so that we can later use the matVecMul3x3(...)
        // directly. This is different to the definition of the reciprocal lattice below, which
        // is not transposed.
    {}
    ~RealMeshgrid() =default;

    void inline xyz(const int n, double &x, double &y, double &z) const override
    {
        int i, j, k;
        getIndexTipple(n, i, j, k);
        // lattice is the same as in the POSFILE/POSCAR
        // it needs to be transposed, which is done implicitly here
        x = i * lattice[0][0] / (dim[0]) + j * lattice[1][0] / (dim[1]) + k * lattice[2][0] / (dim[2]) + origin[0];
        y = i * lattice[0][1] / (dim[0]) + j * lattice[1][1] / (dim[1]) + k * lattice[2][1] / (dim[2]) + origin[1];
        z = i * lattice[0][2] / (dim[0]) + j * lattice[1][2] / (dim[1]) + k * lattice[2][2] / (dim[2]) + origin[2];
    }

    void inline xyz(const int i, const int j, const int k, double &x, double &y, double &z) const
    {
        x = i * lattice[0][0] / (dim[0]) + j * lattice[1][0] / (dim[1]) + k * lattice[2][0] / (dim[2]) + origin[0];
        y = i * lattice[0][1] / (dim[0]) + j * lattice[1][1] / (dim[1]) + k * lattice[2][1] / (dim[2]) + origin[1];
        z = i * lattice[0][2] / (dim[0]) + j * lattice[1][2] / (dim[1]) + k * lattice[2][2] / (dim[2]) + origin[2];
    }

    /**
     * @brief Calculates the indexes of the meshgrid that belong to a certain cartesian position (Reverse operation to xyz(...))
     *
     * @param x,y,z  cartesian coordinates
     * @param i,j,k  returned indexes
     * @return true  found correct index
     * @return false position is out of range / not contained in the meshgrid
     */
    bool estimateIndexFromPosition(double x, double y, double z, int& i, int&j, int& k) const
    {
        vector<double> id = matVecMul3x3(invLattice, vector<double>{x - origin[0], y - origin[1], z - origin[2]});

        i = round(id[0]*dim[0]);
        j = round(id[1]*dim[1]);
        k = round(id[2]*dim[2]);

        // check boundaries
        if ((i<0) || (i>=dim[0])) return false;
        if ((j<0) || (j>=dim[1])) return false;
        if ((k<0) || (k>=dim[2])) return false;

        return true;
    }

    /**
     * @brief Gets i,j,k coordinates of a shifted data point (old and slow)
     *
     * THIS ROUTINE IS SLOW! ONLY USE IT FOR SINGLE POINTS, NOT INSIDE A LOOP!!!
     * See also RealMeshgrid::getIndexOffset() and Real Meshgrid::getShiftedGlobalIndex()
     * For more optimized workflows.
     *
     * Calculates the index triple (x,y,z) for a lattice point n that is shifted by +R
     * i.e. (x,y,z) <=> 'n+R'. Depending on R and the used supercell it is not guaranteed
     * that (x,y,z) is a valid data point.
     *
     * @param n global index (data index)
     * @param R shift vector
     * @param x,y,z index for each direction
     * @return true if (x,y,z) is within the data grid and false otherwise.
     */
    bool getShiftedIndexTipple(int n, vector<double> R, int &x, int &y, int &z) const
    {
        this->getIndexTipple(n, x, y, z);
        // vector< vector<double> > invLattice = invMat3x3(transpose3x3(lattice));
        vector<double> R_direct = matVecMul3x3(invLattice, R); // invLattice was already transposed in the constructor

        x += round(R_direct[0] * dim[0]);
        y += round(R_direct[1] * dim[1]);
        z += round(R_direct[2] * dim[2]);

        // check if index is in data range
        if ((x >= 0) && (x < dim[0]) && (y >= 0) && (y < dim[1]) && (z >= 0) && (z < dim[2]))
            return true;
        else
            return false;
    }

    /**
     * @brief Shift the index of a grid point by a given offset in x,y,z-direction
     *
     * The indexOffset can be obtained by RealMeshgrid::getIndexOffset()
     *
     * @param n     global index that should be shifted
     * @param indexOffset   x-,y-,z- offset in terms of indexes
     * @return int  resulting global index if shift leads to a valid grid point or -1 if new point is outside the grid
     */
    int getShiftedGlobalIndex(const int n, const vector<int>& indexOffset) const
    {
        int x,y,z;
        this->getIndexTipple(n, x, y, z);

        x += indexOffset[0];
        y += indexOffset[1];
        z += indexOffset[2];

        // check if index is in data range
        if ((x >= 0) && (x < dim[0]) && (y >= 0) && (y < dim[1]) && (z >= 0) && (z < dim[2]))
            return this->getGlobId(x, y, z);
        else
            return -1;  // value not contained in the array
    }

    /**
     * @brief Get the index offset for a shift in terms of unitcells
     *
     * @param R  shift in terms of lattice vectors
     * @param supercell  how many unit cells are contained in the supercell for each direction
     * can be obtained by WannierFunction::getLatticeInUnitcellBasis()
     * @return vector<int>
     */
    vector<int> getIndexOffset(const vector<int>& R, vector<double> const& supercell) const
    {
        vector<int> ijk(3);
        for (int i=0; i<3; i++){
            ijk[i] = -R[i] * round(dim[i] / supercell[i]);
        }

        return ijk;
    }
};

/**
 * @brief Meshgrid in reciprocal space
 *
 */
class ReciprocalMeshgrid : public Meshgrid
{
public:
    explicit ReciprocalMeshgrid(const RealMeshgrid *realMesh)
        : Meshgrid(realMesh->getDim(), realMesh->getLattice(), vector<double>{0.,0.,0.})
    {
        // We define the reciprocal lattice as the inverse of the lattice as given from
        // POSFILE/POSCAR/*.xsf file (not the transpose).
        // The order of transposition and inversion does not matter: (A^(-1))^T = (A^T)^(-1)

        vector<vector<double>> recLattice = invMat3x3(realMesh->getLattice());
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                recLattice[i][j] = 2 * M_PI * recLattice[i][j];
            }
        }

        this->lattice = recLattice; // update lattice to reciprocal lattice
    }
    ~ReciprocalMeshgrid() = default;

    void inline xyz(int n, double &x, double &y, double &z) const override
    {
        int i, j, k;
        getIndexTipple(n, i, j, k); // get indexes in x,y and z direction
        xyz(i, j, k, x, y, z);      // get reciprocal coordinates for that index
    }

    void inline xyz(int i, int j, int k, double &x, double &y, double &z) const
    {

        // shift reciprocal lattice such that we have positive and negative wave vectors
        if (i > dim[0] / 2) i = i - dim[0];
        if (j > dim[1] / 2) j = j - dim[1];
        if (k > dim[2] / 2) k = k - dim[2];

        // multiply with reciprocal lattice vectors
        x = lattice[0][0] * i + lattice[0][1] * j + lattice[0][2] * k;
        y = lattice[1][0] * i + lattice[1][1] * j + lattice[1][2] * k;
        z = lattice[2][0] * i + lattice[2][1] * j + lattice[2][2] * k;

        // for test (does not give the correct results)
        // x = lattice[0][0] * i + lattice[1][0] * j + lattice[2][0] * k;
        // y = lattice[0][1] * i + lattice[1][1] * j + lattice[2][1] * k;
        // z = lattice[0][2] * i + lattice[1][2] * j + lattice[2][2] * k;
    }
};

#endif