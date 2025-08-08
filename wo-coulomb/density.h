#ifndef DENSITY_H
#define DENSITY_H

/**
 * @file density.h
 * @author Konrad Merkel
 * @copyright Copyright (c) 2024
 *
 * @brief Contains all data structures and functions that deal with (overlap) densities
 *
 * A charge density is defined as
 *
 *      rho(x) = w1(x) * w2(x-R),
 *
 * where w1 and w2 are Wannier functions and R is the relative shift vector in terms of
 * unit cells.
 *
 * This file contains parallel and non-parallel functions to
 *      - calculate the density on a meshgrid (as flattened array)
 *      - calculate indicator data, such as charge center, absolute charge, extend, ...
 *        which can be used to estimate values of the integrals to evaluate if they
 *        need to be calculated in full detail (e.g. MBIE-0 procedure).
 *
 */

#include "wannierfunction.h"
#include "mpi_tags.h"
#include <mpi.h>
#include <fftw3.h>
#include <limits.h>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cassert>

/**
 * Data structure for labels of a density (descriptors)
 *
 * This only stores the labels in terms of Wannier function ids and relative shift vector.
 * It does not store the actual meshgrid data (array).
 */
struct Density_descr
{
    int id1;            //!< Wannier function 1 (usually conduction WF)
    int id2;            //!< Wannier function 2 (usually valence WF)
    vector<int> R;      //!< relative shift vector between WFs

    Density_descr()
     : id1{INT_MAX}, id2{INT_MAX}, R{ vector<int>{INT_MAX,INT_MAX,INT_MAX} }
    {}

    Density_descr(int i, int j, vector<int> const& nR)
     : id1{i}, id2{j}, R{ nR }
    {}

    string toString() const {
        ostringstream out;
        out << id1 << "\t" << id2 << "\t" << R[0] << "\t" << R[1] << "\t" << R[2];
        return std::move(out).str();
    }

    vector<int> toVector() const {
        return vector<int>{id1, id2, R[0], R[1], R[2]};
    }

    bool operator==(const Density_descr& b) const {
        return (b.id1 == id1) && (b.id2 == id2) && (b.R[0] == R[0]) && (b.R[1] == R[1]) && (b.R[2] == R[2]);
    }

    bool operator!=(const Density_descr& b) const {
        return ! operator==(b);
    }

    bool operator<(const Density_descr& b) const {  // needed to use it as key in a map
        return this->toVector() < b.toVector();
    }

    bool operator<=(const Density_descr& b) const {  // needed for Hermicity
        return this->toVector() <= b.toVector();
    }

    Density_descr conjugated() const {
        return Density_descr(id2,id1,vector<int>{-R[0],-R[1],-R[2]});
    }

    bool isOverlapDensity() const {
        return (id1!=id2) || (R[0]!=0) || (R[1]!=0) || (R[2]!=0);
    }

    bool isClassicalDensity() const {
        return this->isOverlapDensity() == false;
    }
};

/**
 * @brief Data structure to store necessary information about a (overlap) charge density.
 *
 * This information can be used to estimate integrals.
 *
 */
struct Density_indicator
{
    double x,y,z;       //!< charge center
    double absCharge;   //!< absolute charge = integral dr |rho(r)|
    double extend;      //!< extend of the charge density in Angstrom

    Density_indicator()
     : x(0), y(0), z(0), absCharge(0), extend(0) {}

    Density_indicator(double x, double y, double z, double absCharge, double extend)
     : x(x), y(y), z(z), absCharge(absCharge), extend(extend) {}

    vector<double> center() const { return vector<double>{x,y,z}; }

    vector<double> toVector() const { return vector<double>{x,y,z,absCharge,extend}; }
};

/**
 * @brief Calculates the (overlap) density of two Wannier functions
 *
 * Calculates the (overlap) density of two Wannier functions and a shift vector R (in units of unit cells)
 * The shift vector acts on wann2
 *
 * rho(x) = w1(x) * w2(x-R)
 *
 * If R=0 this gives a classical charge density.
 * This is the fftw version. Maybe in the future we can only have one version
 *
 * @param wann1
 * @param wann2
 * @param R
 * @param density  output density (already allocated array)
 */
void joinedDensity(WannierFunction const& wann1, WannierFunction const& wann2, vector<int> const& R, fftw_complex* density)
{
    if (density == nullptr)
        throw runtime_error("fftw_complex array for density has to be allocated!");

    if (!wann1.isCompatible(wann2)) {
        throw runtime_error("Cannot create joined density because Wannier functions are not compatible");
    }

    int N  = wann1.getMeshgrid()->getNumDataPoints();
    double const* value1 = wann1.getValue();
    double const* value2 = wann2.getValue();

    #pragma omp parallel shared(density, value1, value2, wann1, wann2) firstprivate(N,R)
    {
        if ((R[0] ==0 ) && (R[1] ==0 ) && (R[2] ==0 )) {

            #pragma omp parallel for
            for (int i=0; i<N;  i++) {
                density[i][0] = value1[i]* value2[i];  // charge density
                density[i][1] = 0.0;
            }

        } else {

            vector<double> sc = wann1.getLatticeInUnitcellBasis();

            if ((round(sc[0])<=abs(R[0])) || (round(sc[1])<=abs(R[1])) || (round(sc[2])<=abs(R[2]))) {
                cout << "Shift vector for joint density is larger than supercell\n";
                cout << "R = " << R[0] << " " << R[1] << " " << R[2] << "\n";
                cout << "supercell = " << sc[0] << " " << sc[1] << " " << sc[2] << "\n";
                #pragma omp for
                for (int i=0; i<N;  i++) {
                    density[i][0] = 0.0;
                    density[i][1] = 0.0;
                }
            } else {

                // calculate index offset
                const RealMeshgrid* mesh = wann2.getMeshgrid();
                vector<int> indexOffset = mesh->getIndexOffset(R, wann2.getLatticeInUnitcellBasis());
                vector<int> dim = mesh->getDim();

                // initialize all values with zeros
                #pragma omp for
                for (int i=0; i<N;  i++) {
                    density[i][0] = 0.0;
                    density[i][1] = 0.0;
                }

                // set boundaries for loops
                int boundi, boundj, boundk;
                if (indexOffset[0] < 0) boundi = dim[0];
                else boundi = dim[0] - indexOffset[0];

                if (indexOffset[1] < 0) boundj = dim[1];
                else boundj = dim[1] - indexOffset[1];

                if (indexOffset[2] < 0) boundk = dim[2];
                else boundk = dim[2] - indexOffset[2];

                int globId;
                #pragma omp for private(globId)
                for (int k=max(0,-indexOffset[2]); k<boundk;  k++) {
                    for (int j=max(0,-indexOffset[1]); j<boundj;  j++) {
                        for (int i=max(0,-indexOffset[0]); i<boundi;  i++) {

                            // int globId = mesh->getGlobId(i,j,k);
                            // int globId2 = mesh->getGlobId(i+ indexOffset[0],j+ indexOffset[1],k+ indexOffset[2]); //mesh->getGlobId(x,y,z);
                            // density[globId][0] = value1[globId]*value2[globId2] ;  // overlap density

                            globId = i + dim[0] * (j + dim[1] * k);
                            density[globId][0] = value1[globId]*value2[(i + indexOffset[0]) + dim[0] * ((j + indexOffset[1]) + dim[1]*(k + indexOffset[2]))] ;  // overlap density
                        }
                    }
                }
            }
        }
    }
}


/**
 * @brief Calculates the (overlap) density of two Wannier functions
 *
 * Calculates the (overlap) density of two Wannier functions and a shift vector R (in units of unit cells)
 * The shift vector acts on wann2
 * rho(x) = w1(x) * w2(x-R)
 *
 * If R=0 this gives the classical charge density.
 *
 * @param wann1
 * @param wann2
 * @param R
 * @return double*
 */
unique_ptr<double[], free_deleter> joinedDensity(WannierFunction const& wann1, WannierFunction const& wann2, vector<int> const& R)
{

    if (!wann1.isCompatible(wann2)) {
        throw runtime_error("Cannot create joined density because Wannier functions are not compatible");
    }

    int N  = wann1.getMeshgrid()->getNumDataPoints();
    unique_ptr<double[], free_deleter> density{ (double*) malloc(sizeof(double)* N) };
    double const* value1 = wann1.getValue();
    double const* value2 = wann2.getValue();

    #pragma omp parallel shared(density, value1, value2, wann1, wann2) firstprivate(N,R)
    {
        if ((R[0] ==0 ) && (R[1] ==0 ) && (R[2] ==0 )) {

            #pragma omp for
            for (int i=0; i<N;  i++) {
                density[i] = value1[i]* value2[i];  // charge density
            }

        } else {

            vector<double> sc = wann1.getLatticeInUnitcellBasis();

            if ((round(sc[0])<=abs(R[0])) || (round(sc[1])<=abs(R[1])) || (round(sc[2])<=abs(R[2]))) {
                cout << "Shift vector for joint density is larger than supercell\n";
                cout << "R = " << R[0] << " " << R[1] << " " << R[2] << "\n";
                cout << "supercell = " << sc[0] << " " << sc[1] << " " << sc[2] << "\n";
                #pragma omp for
                for (int i=0; i<N;  i++) {
                    density[i] = 0.0;
                }
            } else {
                const RealMeshgrid* mesh = wann2.getMeshgrid();
                vector<int> indexOffset = mesh->getIndexOffset(R, wann2.getLatticeInUnitcellBasis());
                vector<int> dim = mesh->getDim();

                // initialize all values with zeros
                #pragma omp for
                for (int i=0; i<N;  i++) {
                    density[i] = 0.0;
                }

                // set boundaries for loops
                int boundi, boundj, boundk;
                if (indexOffset[0] < 0) boundi = dim[0];
                else boundi = dim[0] - indexOffset[0];

                if (indexOffset[1] < 0) boundj = dim[1];
                else boundj = dim[1] - indexOffset[1];

                if (indexOffset[2] < 0) boundk = dim[2];
                else boundk = dim[2] - indexOffset[2];

                int globId;
                #pragma omp for private(globId)
                for (int k=max(0,-indexOffset[2]); k<boundk;  k++) {
                    for (int j=max(0,-indexOffset[1]); j<boundj;  j++) {
                        for (int i=max(0,-indexOffset[0]); i<boundi;  i++) {

                            // int globId = mesh->getGlobId(i,j,k);
                            // int globId2 = mesh->getGlobId(i+ indexOffset[0],j+ indexOffset[1],k+ indexOffset[2]); //mesh->getGlobId(x,y,z);
                            // density[globId][0] = value1[globId]*value2[globId2] ;  // overlap density

                            globId = i + dim[0] * (j + dim[1] * k);
                            density[globId] = value1[globId]*value2[(i + indexOffset[0]) + dim[0] * ((j + indexOffset[1]) + dim[1]*(k + indexOffset[2]))] ;  // overlap density
                        }
                    }
                }
            }
        }
    }
    return density;
}


/**
 * @brief Checks if an array only contains zeros.
 *
 */
bool isZero(const double* density, int N) {
    bool isZero = true;
    for (int j=0; j<N; j++) {
        if (abs(density[j]) > 1e-10) {
            isZero = false;
            break;
        }
    }
    return isZero;
}


/**
 * @brief Checks if an array only contains zeros.
 *
 */
bool isZero(const fftw_complex* density, int N) {
    bool isZero = true;
    for (int j=0; j<N; j++) {
        if ((abs(density[j][0]) > 1e-10) || (abs(density[j][1]) > 1e-10)) {
            isZero = false;
            break;
        }
    }
    return isZero;
}


/**
 * @brief Get the Monopole object from a density array
 *
 * @param density (overlap) density
 * @param mesh real space grid, where the density is defined
 * @return Monopole
 */
Monopole getMonopole(const double* density, const RealMeshgrid* mesh)
{
    // shared variables
    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double charge = 0.0;
    double norm =0.0;
    double xm=0.0;
    double ym=0.0;
    double zm=0.0;
    #pragma omp parallel shared(density, charge, xm,ym,zm) firstprivate(dV, Npoints)
    {
        #pragma omp for reduction(+:charge) reduction(+:norm)
        for (int j=0; j<Npoints; j++){
            charge += density[j] * dV;
            norm += abs(density[j]) * dV;
        }

        if (abs(charge) > 1e-10) {  // to make sure we do not divide by zero

            // calculate charge center

            if (mesh->hasMeshgridArrays()){
                const double* XX = mesh->getXX();
                const double* YY = mesh->getYY();
                const double* ZZ = mesh->getZZ();

                # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
                for (int j=0; j<Npoints; j++){
                    xm += abs(density[j]) * XX[j] * dV / norm;  // overlap densities can get negative!!!
                    ym += abs(density[j]) * YY[j] * dV / norm;
                    zm += abs(density[j]) * ZZ[j] * dV / norm;
                }
            }else{
                double x,y,z;
                # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
                for (int j=0; j<Npoints; j++){
                    mesh->xyz(j,x,y,z);
                    xm += abs(density[j]) * x * dV / norm;  // overlap densities can get negative!!!
                    ym += abs(density[j]) * y * dV / norm;
                    zm += abs(density[j]) * z * dV / norm;
                }
            }
        }
    }

    Monopole mono;
    mono.x = xm;
    mono.y = ym;
    mono.z = zm;
    mono.charge=charge;

    return mono;
}


/**
 * @brief Get the Monopole object from a density array
 *
 * @param density (overlap) density
 * @param mesh real space grid, where the density is defined
 * @return Monopole
 */
Monopole getMonopole(const fftw_complex* density, const RealMeshgrid* mesh)
{
    // shared variables
    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double charge = 0.0;
    double norm =0.0;
    double xm=0.0;
    double ym=0.0;
    double zm=0.0;
    #pragma omp parallel shared(density, charge, xm,ym,zm) firstprivate(dV, Npoints)
    {
        #pragma omp for reduction(+:charge) reduction(+:norm)
        for (int j=0; j<Npoints; j++){
            charge += density[j][0] * dV;
            norm += abs(density[j][0]) * dV;
        }

        if (abs(charge) > 1e-10) {  // to make sure we do not divide by zero

            if (mesh->hasMeshgridArrays()){
                const double* XX = mesh->getXX();
                const double* YY = mesh->getYY();
                const double* ZZ = mesh->getZZ();

                # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
                for (int j=0; j<Npoints; j++){
                    xm += abs(density[j][0]) * XX[j] * dV / norm;  // overlap densities can get negative!!!
                    ym += abs(density[j][0]) * YY[j] * dV / norm;
                    zm += abs(density[j][0]) * ZZ[j] * dV / norm;
                }
            }else{
                double x,y,z;
                # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
                for (int j=0; j<Npoints; j++){
                    mesh->xyz(j,x,y,z);
                    xm += abs(density[j][0]) * x * dV / norm;  // overlap densities can get negative!!!
                    ym += abs(density[j][0]) * y * dV / norm;
                    zm += abs(density[j][0]) * z * dV / norm;
                }
            }
        }
    }

    Monopole mono;
    mono.x = xm;
    mono.y = ym;
    mono.z = zm;
    mono.charge=charge;

    return mono;
}


/**
 * @brief Calculates covariance matrix for a given charge density
 *
 */
vector< vector<double>> getCovariance(const double* density, const RealMeshgrid* mesh, vector<double> chargeCenter)
{
    vector< vector<double>> var(3);
    for (int i=0;i<3;i++) {
        var[i] = vector<double>{0,0,0};
    }

    vector<double> xyz{0,0,0};
    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();

    for (int n=0; n<Npoints; n++){
        mesh->xyz(n,xyz[0],xyz[1],xyz[2]);
        for (int i=0;i<3;i++) {
            for (int j=0;j<3; j++) {
                var[i][j] = density[n] * dV * (xyz[i] - chargeCenter[i])*(xyz[j] - chargeCenter[j]);
            }
        }
    }

    return var;
}


/**
 * @brief Calculates mean of the eigen-values of the covariance matrix from the trace of the matrix.
 *
 */
double getMeanVariance(const fftw_complex* density, const RealMeshgrid* mesh, vector<double> chargeCenter)
{
    double var=0.0;
    vector<double> xyz{0,0,0};
    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();

    # pragma omp parallel for shared(mesh) firstprivate(xyz, dV, chargeCenter) reduction(+:var)
    for (int n=0; n<Npoints; n++){
        mesh->xyz(n,xyz[0],xyz[1],xyz[2]);
        for (int i=0;i<3;i++) {
            var += density[n][0] * dV * pow(xyz[i] - chargeCenter[i],2);  // trace of the covariance matrix
        }
    }

    return abs(var)/3.;  // for a general overlap density the variance could be negative
}


/**
 * @brief Scalar-product between two Wannier functions
 *
 * @param wann1,wann2 Wannier functions
 * @param R  relative shift in units of primitive lattice vectors
 */
double scalarProduct(WannierFunction const& wann1, WannierFunction const& wann2, vector<int> const& R = vector<int>{0,0,0}) {

  if (! wann1.isCompatible(wann2)) {
    throw runtime_error("Cannot calcualte scalar product between two WF that are not compatible.");
  }

  const RealMeshgrid* mesh = wann1.getMeshgrid();
  double dV = mesh->getdV();
  int N = mesh->getNumDataPoints();

//   double* value1 = wann1->getValue();
//   double* value2 = wann2->getValue();

  unique_ptr<double[], free_deleter>value{ joinedDensity(wann1, wann2, R) };

  double ret = 0.0;
  for (int i=0; i<N; i++) {
    // ret += value1[i]*value2[i] *dV;
    ret += value[i] *dV;
  }
  return ret;
}


/**
 * @brief Calculates the support radius of a Wannier function
 *
 * Data points with larger distance to the Wannier centre than the support radius
 * have zero values (within the threshold).
 *
 * @param wann
 * @param threshold
 * @param mono centre of the Wannier function
 * @return double
 */
double getSupportRadius(WannierFunction const& wann, double threshold, Monopole* mono=0)
{
    const RealMeshgrid* meshgrid = wann.getMeshgrid();
    if (mono==0) {
        unique_ptr<double[], free_deleter>density{ joinedDensity(wann, wann, vector<int>{0,0,0}) };
        Monopole m = getMonopole(density.get(),meshgrid);
        mono = &m;
    }

    double const* values = wann.getValue();

    int Npoints = meshgrid->getNumDataPoints();
    double supp_radius = 0;
    double x,y,z,r;
    for (int i=0; i<Npoints; i++) {
        if (abs(values[i]) > threshold) {
            meshgrid->xyz(i,x,y,z);
            r = sqrt(pow(mono->x - x, 2) + pow(mono->y - y, 2) + pow(mono->z - z, 2));
            if (r>supp_radius) supp_radius=r;
        }
    }
    return supp_radius;
}


/**
 * @brief Calculates the absolute monopole of an overlap density
 *
 * The absolute monopole is defined as
 * abs_charge    = integral dr   |rho(r)|
 * charge center = integral dr r*|rho(r)| / abs_charge
 *
 * @param density
 * @param mesh
 * @return Monopole
 */
Monopole getAbsMonopole(const double* density, const RealMeshgrid* mesh)
{

    if (! mesh->hasMeshgridArrays()) {
        throw runtime_error("Meshgrid has no MeshgridArrays!");
    }

    if (density == nullptr) {
        Monopole mono;
        mono.x = 0.0;
        mono.y = 0.0;
        mono.z = 0.0;
        mono.charge=0.0;

        return mono;
    }

    // shared variables
    const double* XX = mesh->getXX();
    const double* YY = mesh->getYY();
    const double* ZZ = mesh->getZZ();

    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double norm =0.0;
    double xm=0.0;
    double ym=0.0;
    double zm=0.0;

    #pragma omp parallel shared(density, norm, xm,ym,zm, XX, YY, ZZ) firstprivate(dV, Npoints)
    {
        #pragma omp for reduction(+:norm)
        for (int j=0; j<Npoints; j++){
            norm += abs(density[j]) * dV;
        }

        if (abs(norm) > 1e-10) {  // to make sure we do not divide by zero

            // calculate charge center
            # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
            for (int j=0; j<Npoints; j++){
                xm += abs(density[j]) * XX[j] * dV / norm;  // overlap densities can get negative!!!
                ym += abs(density[j]) * YY[j] * dV / norm;
                zm += abs(density[j]) * ZZ[j] * dV / norm;
            }
        }
    }

    Monopole mono;
    mono.x = xm;
    mono.y = ym;
    mono.z = zm;
    mono.charge=norm;

    return mono;
}


/**
 * @brief Calculates the absolute monopole of an overlap density
 *
 * The absolute monopole is defined as
 * abs_charge    = integral dr   |rho(r)|
 * charge center = integral dr r*|rho(r)| / abs_charge
 *
 * This is the FFTW version.
 *
 * @param density
 * @param mesh
 * @return Monopole
 */
Monopole getAbsMonopole(const fftw_complex* density, RealMeshgrid* mesh)
{

    if (! mesh->hasMeshgridArrays()) {
        throw runtime_error("Meshgrid has no MeshgridArrays!");
    }

    if (density == nullptr) {
        Monopole mono;
        mono.x = 0.0;
        mono.y = 0.0;
        mono.z = 0.0;
        mono.charge=0.0;

        return mono;
    }

    // shared variables
    const double* XX = mesh->getXX();
    const double* YY = mesh->getYY();
    const double* ZZ = mesh->getZZ();

    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double norm =0.0;
    double xm=0.0;
    double ym=0.0;
    double zm=0.0;

    #pragma omp parallel shared(density, norm, xm,ym,zm, XX, YY, ZZ) firstprivate(dV, Npoints)
    {
        #pragma omp for reduction(+:norm)
        for (int j=0; j<Npoints; j++){
            norm += abs(density[j][0]) * dV;
        }

        if (abs(norm) > 1e-10) {  // to make sure we do not divide by zero

            // calculate charge center
            # pragma omp for reduction(+:xm) reduction(+:ym) reduction(+:zm)
            for (int j=0; j<Npoints; j++){
                xm += abs(density[j][0]) * XX[j] * dV / norm;  // overlap densities can get negative!!!
                ym += abs(density[j][0]) * YY[j] * dV / norm;
                zm += abs(density[j][0]) * ZZ[j] * dV / norm;
            }
        }
    }

    Monopole mono;
    mono.x = xm;
    mono.y = ym;
    mono.z = zm;
    mono.charge=norm;

    return mono;
}


vector<double> getAbsVariance(const double* density, const RealMeshgrid* mesh, vector<double> chargeCenter)
{

    if (! mesh->hasMeshgridArrays()) {
        throw runtime_error("Meshgrid has no MeshgridArrays!");
    }

    // shared variables
    const double* XX = mesh->getXX();
    const double* YY = mesh->getYY();
    const double* ZZ = mesh->getZZ();


    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double vx=0.0, vy=0.0, vz=0.0;

    # pragma omp parallel for shared(XX,YY,ZZ) firstprivate(dV, chargeCenter) reduction(+:vx) reduction(+:vy) reduction(+:vz)
    for (int n=0; n<Npoints; n++){
        // diagonal of the covariance matrix
        vx += abs(density[n]) * dV * pow(XX[n] - chargeCenter[0],2);
        vy += abs(density[n]) * dV * pow(YY[n] - chargeCenter[1],2);
        vz += abs(density[n]) * dV * pow(ZZ[n] - chargeCenter[2],2);
    }

    return vector<double>{vx, vy, vz};;
}


vector<double> getAbsVariance(const fftw_complex* density, const RealMeshgrid* mesh, vector<double> chargeCenter)
{

    if (! mesh->hasMeshgridArrays()) {
        throw runtime_error("Meshgrid has no MeshgridArrays!");
    }

    // shared variables
    const double* XX = mesh->getXX();
    const double* YY = mesh->getYY();
    const double* ZZ = mesh->getZZ();


    const int Npoints = mesh->getNumDataPoints();
    double dV = mesh->getdV();
    double vx=0.0, vy=0.0, vz=0.0;

    # pragma omp parallel for shared(XX,YY,ZZ) firstprivate(dV, chargeCenter) reduction(+:vx) reduction(+:vy) reduction(+:vz)
    for (int n=0; n<Npoints; n++){
        // diagonal of the covariance matrix
        vx += abs(density[n][0]) * dV * pow(XX[n] - chargeCenter[0],2);
        vy += abs(density[n][0]) * dV * pow(YY[n] - chargeCenter[1],2);
        vz += abs(density[n][0]) * dV * pow(ZZ[n] - chargeCenter[2],2);
    }

    return vector<double>{vx, vy, vz};;
}


double getExtend(const double* density, const RealMeshgrid* mesh, Monopole const& mono, double criterion=1e-4)
{
    vector<double> variance = getAbsVariance(density, mesh, vector<double>{mono.x,mono.y,mono.z});

    auto max_iter = max_element(variance.begin(), variance.end());
    // cout << "variance = " << variance[0] << " " << variance[1] << " " << variance[2] << endl;
    //int max_id = distance(variance.begin(), max_iter);
    //cout << "max id: " <<  max_id << ": " << *max_iter << " = " << variance[max_id] << endl;

    double alpha = 1./(2* (*max_iter));
    double extend = sqrt(-log(criterion)/(2*alpha));
    int direction = distance(variance.begin(), max_iter);
    // cout << "direction = " << direction << endl;
    // cout << "extend (initial guess) = " << extend << endl;

    // get initial interval
    auto lattice = mesh->getLattice();
    double dx = sqrt(pow(lattice[direction][0],2) + pow(lattice[direction][1],2) + pow(lattice[direction][2],2)) / mesh->getDim()[0];
    int delta = ceil(extend / dx) ;

    vector<int> id{0,0,0};
    mesh->estimateIndexFromPosition(mono.x, mono.y, mono.z, id[0], id[1], id[2]);  // TODO: we could also use minus extend

    // copy index
    vector<int> ida(3), idb(3), idc(3);
    for (int l=0; l<3; l++) {
        ida[l] = id[l];
        idb[l] = id[l];
        idc[l] = id[l];
    }

    // TODO: what is with id-delta ? that would also be a valid choice and maybe leads to a larger extend
    ida[direction] = id[direction]+delta*0.4;
    idb[direction] = min(mesh->getDim()[direction]-1, int(id[direction] + 3*delta/2));
    double fa = abs(density[mesh->getGlobId(ida[0], ida[1], ida[2])]);
    double fb = abs(density[mesh->getGlobId(idb[0], idb[1], idb[2])]);

    // check intervall and extend it if needed
    if ((fa-criterion) < 0.0) {
        ida[direction] = id[direction];  // set to center of charges
        fa = abs(density[mesh->getGlobId(ida[0], ida[1], ida[2])]);
    }
    while((fb-criterion) > 0.0)
    {
        if (idb[direction]>=mesh->getDim()[direction]-1) {
            throw runtime_error("Unable to find intervall for bisection algorithm.");
        }
        // cout << "[find intervall] " << ida[direction] << " " << idb[direction]  << " fb=" << fb-criterion << endl;

        ida[direction] = idb[direction];

        idb[direction] = ceil((idb[direction] + mesh->getDim()[direction]-1 )/2);
        fb = abs(density[mesh->getGlobId(idb[0], idb[1], idb[2])]);
    }

    if (sgn<double>((fa-criterion)) == sgn<double>((fb-criterion))) // interval is not valid
    {
        //cerr << "[WARNING] cannot find interval for bisection method. Just use initial guess from variance." << endl;
        return extend;
    }

    if (ida[direction] >= idb[direction]) {
        cout << "a=" << ida[direction] << " b=" << idb[direction] << " i=" << id[direction] << " delta=" << delta << endl;
        throw runtime_error("Problem in bisection");
        return extend;
    }

    // for (int l=a; l<=b; l++) {
    //     int globId2 = mesh->getGlobId(l,j,k);
    //     cout << l << " : fx=" << abs(density[globId2]) << endl;
    // }


    // cout << "fa = " << fa << " fb = " << fb << endl;

    // cout << "(initial guess) : a=" << ida[direction] << " b=" << idb[direction] << " extend=" << extend << endl;
    int globId;
    for (int n=0; n<100; n++) {

        idc[direction] = (ida[direction] + idb[direction]) / 2;  // middle of the interval

        globId = mesh->getGlobId(idc[0], idc[1], idc[2]);
        double fc = abs(density[globId]);

        // cout << n << " : a=" << ida[direction] << " b=" << idb[direction] <<  " f(x)=" << fc-criterion << endl;
        if ((abs(fc-criterion)< 1e-5) || ((idb[direction]-ida[direction]) <= 1)) { // found result
            // cout << "found solution : " << idc[direction] << " f(x)=" << fc-criterion << endl;
            break;
        }

        int globId2 = mesh->getGlobId(ida[0], ida[1], ida[2]);
        double fa = abs(density[globId2]);

        if (sgn<double>((fc-criterion)) == sgn<double>((fa-criterion))) {
            ida[direction] = idc[direction];
        } else {
            idb[direction] = idc[direction];
        }
    }

    double posX, posY, posZ;
    mesh->xyz(globId,posX, posY, posZ);
    extend = sqrt( pow(mono.x - posX, 2) + pow(mono.y - posY, 2) + pow(mono.z - posZ, 2) );
    // cout << "extend = " << extend << endl;

    return extend;
}


double getExtend(fftw_complex* density, const RealMeshgrid* mesh, Monopole mono, double criterion=1e-4)
{
    vector<double> variance = getAbsVariance(density, mesh, vector<double>{mono.x,mono.y,mono.z});

    auto max_iter = max_element(variance.begin(), variance.end());
    //int max_id = distance(variance.begin(), max_iter);
    //cout << "max id: " <<  max_id << ": " << *max_iter << " = " << variance[max_id] << endl;

    double alpha = 1./(2* (*max_iter));
    double extend = sqrt(-log(criterion)/(2*alpha));
    int direction = distance(variance.begin(), max_iter);
    // cout << "direction = " << direction << endl;
    // cout << "extend (initial guess) = " << extend << endl;

    // get initial interval
    auto lattice = mesh->getLattice();
    double dx = sqrt(pow(lattice[direction][0],2) + pow(lattice[direction][1],2) + pow(lattice[direction][2],2)) / mesh->getDim()[0];
    int delta = ceil(extend / dx) ;

    vector<int> id{0,0,0};
    mesh->estimateIndexFromPosition(mono.x, mono.y, mono.z, id[0], id[1], id[2]);  // TODO: we could also use minus extend

    // copy index
    vector<int> ida(3), idb(3), idc(3);
    for (int l=0; l<3; l++) {
        ida[l] = id[l];
        idb[l] = id[l];
        idc[l] = id[l];
    }


    // TODO: what is with id-delta ? that would also be a valid choice and maybe leads to a larger extend
    ida[direction] = id[direction]+delta*0.4;
    idb[direction] = min(mesh->getDim()[direction]-1, int(id[direction] + 3*delta/2));
    double fa = abs(density[mesh->getGlobId(ida[0], ida[1], ida[2])][0]);
    double fb = abs(density[mesh->getGlobId(idb[0], idb[1], idb[2])][0]);

    // check intervall and extend it if needed
    if ((fa-criterion) < 0.0) {
        ida[direction] = id[direction];  // set to center of charges
        fa = abs(density[mesh->getGlobId(ida[0], ida[1], ida[2])][0]);
    }
    while((fb-criterion) > 0.0)
    {
        if (idb[direction]>=mesh->getDim()[direction]-1) {
            throw runtime_error("Unable to find intervall for bisection algorithm.");
        }
        // cout << "[find intervall] " << ida[direction] << " " << idb[direction]  << " fb=" << fb-criterion << endl;

        ida[direction] = idb[direction];

        idb[direction] = ceil((idb[direction] + mesh->getDim()[direction]-1 )/2);
        fb = abs(density[mesh->getGlobId(idb[0], idb[1], idb[2])][0]);
    }

    if (sgn<double>((fa-criterion)) == sgn<double>((fb-criterion))) // interval is not valid
    {
        cerr << "[WARNING] cannot find interval for bisection method. Just use initial guess from variance." << endl;
        return extend;
    }

    if (ida[direction] >= idb[direction]) {
        cout << "a=" << ida[direction] << " b=" << idb[direction] << " i=" << id[direction] << " delta=" << delta << endl;
        cerr << "Problem in bisection\n";
        return extend;
    }


    // cout << "(initial guess) : a=" << ida[direction] << " b=" << idb[direction] << " extend=" << extend << endl;
    int globId;
    for (int n=0; n<100; n++) {

        idc[direction] = (ida[direction] + idb[direction]) / 2;  // middle of the interval

        globId = mesh->getGlobId(idc[0], idc[1], idc[2]);
        double fc = abs(density[globId][0]);

        // cout << n << " : a=" << ida[direction] << " b=" << idb[direction] <<  " f(x)=" << fc-criterion << endl;
        if ((abs(fc-criterion)< 1e-5) || ((idb[direction]-ida[direction]) <= 1)) { // found result
            // cout << "found solution : " << idc[direction] << " f(x)=" << fc-criterion << endl;
            break;
        }

        int globId2 = mesh->getGlobId(ida[0], ida[1], ida[2]);
        double fa = abs(density[globId2][0]);

        if (sgn<double>((fc-criterion)) == sgn<double>((fa-criterion))) {
            ida[direction] = idc[direction];
        } else {
            idb[direction] = idc[direction];
        }
    }

    double posX, posY, posZ;
    mesh->xyz(globId,posX, posY, posZ);
    extend = sqrt( pow(mono.x - posX, 2) + pow(mono.y - posY, 2) + pow(mono.z - posZ, 2) );
    // cout << "extend = " << extend << endl;

    return extend;
}


/**
 * @brief Calculates the indicator data for a charge density explicitly given by the Wannier functions and the shift vector.
 *
 * @param wann1
 * @param wann2
 * @param R
 * @param criterion_extend
 * @return Density_indicator
 */
Density_indicator calc_indicator(WannierFunction const& wann1, WannierFunction const& wann2, vector<int> const& R, const double criterion_extend) {

    Density_indicator data;

    unique_ptr<double[], free_deleter>density{ joinedDensity(wann1, wann2, R)};

    Monopole mono = getAbsMonopole(density.get(), wann1.getMeshgrid());
    double extend = getExtend(density.get(),  wann1.getMeshgrid(), mono, criterion_extend);

    data.x = mono.x;
    data.y = mono.y;
    data.z = mono.z;
    data.absCharge  = mono.charge;
    data.extend     = extend;

    return data;
}


map<Density_descr,Density_indicator> calcIndicator_estimates(const map< int,WannierFunction >& wannMap, const vector<vector<vector<int>>>& shells, double criterion_extend)
{
    map<Density_descr, Density_indicator> indicators{};

    for (auto itr1 = wannMap.begin(); itr1 != wannMap.end(); itr1++) {
        size_t v1 = itr1->first;
        WannierFunction const& wann1 = itr1->second;

        for (auto itr2 = wannMap.begin(); itr2 != wannMap.end(); itr2++) {
            size_t v2 = itr2->first;

            if (v2>v1) break;  // make sure we do not calculate hermitian conjugated

            WannierFunction const& wann2 = itr2->second;
            double maxMonopole = 1.0;

            // cout << "id1 = " << v1 << ", id2 = " << v2 << endl;

            for (auto const& shell : shells) {
                // cout << "new Shell: " << shell.size() << " maxMonopole: " << maxMonopole << endl;

                if (maxMonopole < 0.1) break;
                maxMonopole = 0.0;

                for (const vector<int> & R : shell) {

                    Density_descr ids = Density_descr(v1,v2,R);
                    auto data = calc_indicator(wann1, wann2, R, criterion_extend);
                    indicators.insert({ids, data});


                    if (v1 != v2) {  // add hermitian conjugated
                        indicators.insert({ids.conjugated(),data});
                    }

                    maxMonopole = max(maxMonopole, data.absCharge);
                }
            }
        }
    }
    return indicators;
}

/**
 * @brief Calculates the indicator data for all important (depending on thresholds) densities.
 *
 * Evaluates the charge density for either valence or conduction overlap densities depending on
 * wannMap.
 * It does not calculate overlap densities between valence and conduction band. To do this pleas use
 * calcLFE_estimates_parallel().
 *
 * @param wannMap               all Wannier functions for valence or conduction band (map: id->WF)
 * @param shells                possible shift vectors (Rv or Rc) in terms of shells (see createShells())
 * @param criterion_extend      threshold for extend
 * @param ABSCHARGE_THRESHOLD   threshold for absolute charge
 * @return map<Density_descr,Density_indicator>*
 */
map<Density_descr,Density_indicator> calcIndicator_estimates_parallel(map< int,WannierFunction > const& wannMap, vector<vector<vector<int>>> const& shells, const double criterion_extend, const double ABSCHARGE_THRESHOLD){

    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    // cout << "calc indicators rank = " << rank << " / " << num_worker << endl;

    map<Density_descr, Density_indicator> indicators{};
    int counter = -1;

    for (auto itr1 = wannMap.begin(); itr1 != wannMap.end(); itr1++) {
        size_t v1 = itr1->first;
        WannierFunction const& wann1 = itr1->second;

        for (auto itr2 = wannMap.begin(); itr2 != wannMap.end(); itr2++) {
            size_t v2 = itr2->first;
            //if (v2>v1) break;  // TODO: BUG make sure we do not calculate hermitian conjugated

            counter++;
            if (counter >= num_worker) counter=0;
            if (counter != rank) continue;

            WannierFunction const& wann2 = itr2->second;

            // cout << "rank = " << rank << " id1 = " << v1 << ", id2 = " << v2 << endl;

            double maxMonopole = 1.0;
            for (auto const& shell : shells) {

                // cout << "rank = " << rank << " new Shell: " << shell.size() << " maxMonopole: " << maxMonopole << endl;

                if (maxMonopole < ABSCHARGE_THRESHOLD) break;

                maxMonopole = 0.0;
                for (const vector<int> & R : shell) {

                    Density_descr ids = Density_descr(v1,v2,R);
                    auto indicator = calc_indicator(wann1, wann2, R, criterion_extend);
                    indicators.insert({ids, indicator});

                    maxMonopole = max(maxMonopole, indicator.absCharge);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);


    // collect all results  //TODO: what happens if there are more workers than combinations?
    if (rank!=0) {  // worker
        //serialize data to send via MPI

        int N = indicators.size();

        if (N==0) {  // worker has no data
            MPI_Send(&N, 1, MPI_INT, 0, TAG_INDICATORS_NUM, MPI_COMM_WORLD);  // inform master process
            return map<Density_descr,Density_indicator>{};  // return empty object
        }

        double* data = (double*)malloc(N*5*sizeof(double));
        int* indexes = (int*)malloc(N*5*sizeof(int));

        int i = 0;
        for (auto itr = indicators.begin(); itr != indicators.end(); itr++) {
            indexes[5*i+0] = itr->first.id1;
            indexes[5*i+1] = itr->first.id2;
            indexes[5*i+2] = itr->first.R[0];
            indexes[5*i+3] = itr->first.R[1];
            indexes[5*i+4] = itr->first.R[2];

            vector<double> indicator_data = itr->second.toVector();

            for (int k=0; k<5; k++) data[5*i+k] = indicator_data[k];
            i++;
        }

        MPI_Send(&N, 1, MPI_INT, 0, TAG_INDICATORS_NUM, MPI_COMM_WORLD);
        MPI_Send(indexes, N*5, MPI_INT, 0, TAG_INDICATORS_INDEXES, MPI_COMM_WORLD);
        MPI_Send(data, N*5, MPI_DOUBLE, 0, TAG_INDICATORS_DATA, MPI_COMM_WORLD);

        free(data);
        free(indexes);

        return map<Density_descr,Density_indicator>{};
    } else {  // master process

        // cout << "Collect data from workers\n";
        int* indexes = nullptr;
        double* data = nullptr;

        for (int worker = 1; worker < num_worker; worker++) {
            // cout << "worker " << worker << " of " << num_worker-1 << endl;
            int N = 0;
            MPI_Recv(&N, 1, MPI_INT, worker, TAG_INDICATORS_NUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // cout << "rank=" << rank <<" worker " << worker << " N=" << N << endl;

            if (N==0) continue; // worker has no data

            indexes=(int*) malloc(sizeof(int)*N*5);
            MPI_Recv(indexes, N*5, MPI_INT, worker, TAG_INDICATORS_INDEXES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            data=(double*) malloc(sizeof(double)*N*5);
            MPI_Recv(data, N*5, MPI_DOUBLE, worker, TAG_INDICATORS_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            // add recived data to map
            for (int i=0; i<N; i++) {
                Density_descr ids = Density_descr(indexes[5*i+0], indexes[5*i+1], vector<int>{indexes[5*i+2],indexes[5*i+3],indexes[5*i+4]});
                // cout << "density: " << ids.toString() << endl;
                Density_indicator indicator(data[5*i+0],data[5*i+1],data[5*i+2],data[5*i+3],data[5*i+4]);
                indicators.insert({ids,indicator});

                // add hermitian conjugate TODO: here is a bug!!!!
                // if (ids.id1 != ids.id2) {

                //     // TODO: the position should be shifted: r = r+R
                //     vector<double> indicator_conj = vector<double>{
                //         data[5*i+0],data[5*i+1],data[5*i+2],data[5*i+3],data[5*i+4]
                //     };

                //     indicators.insert({ids.conjugated(),indicator_conj});
                // }
            }

            // free memory
            free(data);
            free(indexes);
            data = nullptr;
            indexes = nullptr;
        }

        return indicators;
    }
}


map<Density_descr,Density_indicator> calcLFE_estimates(const map< int,WannierFunction >& cWannMap, const map< int,WannierFunction >& vWannMap, vector<vector<vector<int>>> const& shells, const double absCharge_threshold)
{
    const RealMeshgrid* mesh = cWannMap.begin()->second.getMeshgrid();
    map<Density_descr, Density_indicator> indicators{};

    for (auto itr1 = cWannMap.begin(); itr1 != cWannMap.end(); itr1++) {
        size_t c1 = itr1->first;
        WannierFunction const& cWann = itr1->second;

        for (auto itr2 = vWannMap.begin(); itr2 != vWannMap.end(); itr2++) {
            size_t v1 = itr2->first;

            WannierFunction const& vWann = itr2->second;
            double maxMonopole = 1.0;

            // cout << "id1 = " << v1 << ", id2 = " << v2 << endl;

            for (auto const& shell : shells) {
                // cout << "new Shell: " << shell.size() << " maxMonopole: " << maxMonopole << endl;

                if (maxMonopole < absCharge_threshold) break;
                maxMonopole = 0.0;

                for (const vector<int> & R : shell) {

                    // Here we calculate the overlap densities rho(x) = w_c1,0(x) w_v1,R(x) (not w_c1,0(x) w_v1,-S1(x) !!!)
                    // The Indicator-map has (c1,v1,R) as keys. This is the same nomenclature as for all other overlap densities (c1,c2,Rc)
                    // It is not the same as in the LocalFieldEffectsSolver, where (c1,v1,S1) with S1=-R is used!

                    Density_descr ids = Density_descr(c1,v1,R);
                    Density_indicator data;

                    unique_ptr<double[], free_deleter>density{ joinedDensity(cWann, vWann, R) };
                    cout << c1 << " " << v1 << " " << R[0] << " " << R[1] << " " << R[2] << endl;

                    Monopole mono = getAbsMonopole(density.get(), mesh);

                    indicators.insert({ids,Density_indicator(
                        mono.x, mono.y, mono.z, mono.charge, -1  // not set
                    )});

                    maxMonopole = max(maxMonopole, mono.charge);

                }
            }
        }
    }
    return indicators;
}

/**
 * @brief Calculates the indicator data for overlap densities between valence and conduction band
 *
 * This is only needed for local field effects (LFE) calculation.
 *
 * @param cWannMap      all Wannier functions for conduction band (map: id->WF)
 * @param vWannMap      all Wannier functions for valence band (map: id->WF)
 * @param shells        possible shift vectors (S, S') in terms of shells (see createShells())
 * @param ABSCHARGE_THRESHOLD
 * @return map<Density_descr,Density_indicator>*
 */
map<Density_descr,Density_indicator> calcLFE_estimates_parallel(map< int,WannierFunction > const& cWannMap, map< int,WannierFunction > const& vWannMap, vector<vector<vector<int>>> const& shells, const double ABSCHARGE_THRESHOLD){

    int rank, num_worker;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);

    const RealMeshgrid* mesh = cWannMap.begin()->second.getMeshgrid();
    map<Density_descr, Density_indicator> indicators{};
    int counter = -1;

    for (auto itr1 = cWannMap.begin(); itr1 != cWannMap.end(); itr1++) {
        size_t c1 = itr1->first;
        WannierFunction const& cWann = itr1->second;

        for (auto itr2 = vWannMap.begin(); itr2 != vWannMap.end(); itr2++) {
            size_t v1 = itr2->first;

            counter++;
            if (counter >= num_worker) counter=0;
            if (counter != rank) continue;

            WannierFunction const& vWann = itr2->second;

            // cout << "rank = " << rank << " c1 = " << c1 << ", v2 = " << v1 << endl;

            double maxMonopole = 1.0;
            for (auto const& shell : shells) {
                // cout << "rank = " << rank << " new Shell: " << shell.size() << " maxMonopole: " << maxMonopole << endl;

                if (maxMonopole < ABSCHARGE_THRESHOLD) break;

                maxMonopole = 0.0;
                for (const vector<int> & R : shell) {

                    // Here we calculate the overlap densities rho(x) = w_c1,0(x) w_v1,R(x) (not w_c1,0(x) w_v1,-S1(x) !!!)
                    // The Indicator-map has (c1,v1,R) as keys. This is the same nomenclature as for all other overlap densities (c1,c2,Rc)
                    // It is not the same as in the LocalFieldEffectsSolver, where (c1,v1,S1) with S1=-R is used!

                    Density_descr ids = Density_descr(c1,v1,R);
                    unique_ptr<double[], free_deleter>density{ joinedDensity(cWann, vWann, R) };
                    // cout << c1 << " " << v1 << " " << R[0] << " " << R[1] << " " << R[2] << endl;

                    Monopole mono = getAbsMonopole(density.get(), mesh);
                    indicators.insert({ids,Density_indicator(mono.x,mono.y,mono.z, mono.charge, -1)});

                    maxMonopole = max(maxMonopole, mono.charge);
                }
            }
        }
    }

    // cout << "rank " << rank << " at barrier\n";

    MPI_Barrier(MPI_COMM_WORLD);


    // collect all results  //TODO: what happens if there are more workers than combinations?
    if (rank!=0) {  // worker
        //serialize data to send via MPI

        int N = indicators.size();

        if (N==0) {  // worker has no data
            MPI_Send(&N, 1, MPI_INT, 0, TAG_INDICATORS_NUM, MPI_COMM_WORLD);
            return indicators;  // empty map
        }

        double* data = (double*)malloc(N*5*sizeof(double));
        int* indexes = (int*)malloc(N*5*sizeof(int));

        int i = 0;
        for (auto itr = indicators.begin(); itr != indicators.end(); itr++) {
            indexes[5*i+0] = itr->first.id1;
            indexes[5*i+1] = itr->first.id2;
            indexes[5*i+2] = itr->first.R[0];
            indexes[5*i+3] = itr->first.R[1];
            indexes[5*i+4] = itr->first.R[2];

            vector<double> indicator_data = itr->second.toVector();

            for (int k=0; k<5; k++) data[5*i+k] = indicator_data[k];
            i++;
        }

        MPI_Send(&N, 1, MPI_INT, 0, TAG_INDICATORS_NUM, MPI_COMM_WORLD);
        MPI_Send(indexes, N*5, MPI_INT, 0, TAG_INDICATORS_INDEXES, MPI_COMM_WORLD);
        MPI_Send(data, N*5, MPI_DOUBLE, 0, TAG_INDICATORS_DATA, MPI_COMM_WORLD);

        free(data);
        free(indexes);

        return map<Density_descr,Density_indicator>{};
    } else {  // master process

        // cout << "Collect data from workers\n";
        int* indexes = nullptr;
        double* data = nullptr;

        for (int worker = 1; worker < num_worker; worker++) {
            // cout << "worker " << worker << " of " << num_worker-1 << endl;
            int N = 0;
            MPI_Recv(&N, 1, MPI_INT, worker, TAG_INDICATORS_NUM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // cout << "rank=" << rank <<" worker " << worker << " N=" << N << endl;

            if (N==0) continue; // worker has no data

            indexes=(int*) malloc(sizeof(int)*N*5);
            MPI_Recv(indexes, N*5, MPI_INT, worker, TAG_INDICATORS_INDEXES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            data=(double*) malloc(sizeof(double)*N*5);
            MPI_Recv(data, N*5, MPI_DOUBLE, worker, TAG_INDICATORS_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            // add recived data to map
            for (int i=0; i<N; i++) {
                Density_descr ids = Density_descr(indexes[5*i+0], indexes[5*i+1], vector<int>{indexes[5*i+2],indexes[5*i+3],indexes[5*i+4]});
                // cout << "density: " << ids.toString() << endl;
                Density_indicator indicator(data[5*i+0],data[5*i+1],data[5*i+2],data[5*i+3],data[5*i+4]);
                indicators.insert({ids,indicator});
            }

            // free memory
            free(data);
            free(indexes);
            data = nullptr;
            indexes = nullptr;
        }

        return indicators;
    }
}

/**
 * @brief Checks if the indicators are compatible with the known Wannier functions.
 *
 * When we read indicator data file we need to check if they are compatible with the loaded Wannier functions.
 * I.e. all Wannier functions need to be contained and all contained WF need to be known.
 *
 *
 * @param indicators
 * @param wannMap
 */
void checkCompatibilityIndicator(map<Density_descr,Density_indicator> const& indicators, map<int, WannierFunction> const& wannMap) {

    // check if all Wannier functions are contained in the indicator data
    for (auto itr = wannMap.begin(); itr != wannMap.end(); ++itr) {
        bool found = false;
        for (auto itr2 = indicators.begin(); itr2 != indicators.end(); ++itr2) {
            Density_descr dens = itr2->first;

            if ((itr->first == dens.id1) || (itr->first == dens.id2)) {
                found = true;
                break;
            }
        }
        if (!found) {
            cerr << "Wannier funciton: " << itr->first << endl;
            throw runtime_error("Some Wannier funcitons are not contained in the indicator data.");
        }
    }

    // check if indicator data contains unknown Wannier functions
    for (auto itr = indicators.begin(); itr != indicators.end(); ++itr) {
        Density_descr dens = itr->first;

        auto itr2 = wannMap.find(dens.id1);
        if (itr2 == wannMap.end()) {  // not found
            cerr << "Density with unknown Wannier function: " <<dens.toString() << endl;
            throw runtime_error("indicator parameters contain densities with Wannier functions that are not known.");
        }
        itr2 = wannMap.find(dens.id2);
        if (itr2 == wannMap.end()) {  // not found
            cerr << "Density with unknown Wannier function: " <<dens.toString() << endl;
            throw runtime_error("indicator parameters contain densities with Wannier functions that are not known.");
        }
    }
}

/**
 * @brief Same as checkCompatibilityIndicator() but for local field effects (LFE)
 *
 * @param indicators
 * @param vWannMap
 * @param cWannMap
 */
void checkCompatibilityIndicator_lfe(map<Density_descr,Density_indicator> const& indicators, map<int, WannierFunction> const& vWannMap, map<int, WannierFunction> const& cWannMap) {

    // check if all valence Wannier functions are contained in the indicator data
    for (auto itr = vWannMap.begin(); itr != vWannMap.end(); ++itr) {
        bool found = false;
        for (auto itr2 = indicators.begin(); itr2 != indicators.end(); ++itr2) {
            Density_descr dens = itr2->first;

            if (itr->first == dens.id2) {
                found = true;
                break;
            }
        }
        if (!found) {
            cerr << "Wannier funciton: " << itr->first << endl;
            throw runtime_error("Some Wannier funcitons are not contained in the indicator data.");
        }
    }

    // check if all conduction Wannier functions are contained in the indicator data
    for (auto itr = cWannMap.begin(); itr != cWannMap.end(); ++itr) {
        bool found = false;
        for (auto itr2 = indicators.begin(); itr2 != indicators.end(); ++itr2) {
            Density_descr dens = itr2->first;

            if (itr->first == dens.id1) {
                found = true;
                break;
            }
        }
        if (!found) {
            cerr << "Wannier funciton: " << itr->first << endl;
            throw runtime_error("Some Wannier funcitons are not contained in the indicator data.");
        }
    }

    for (auto itr = indicators.begin(); itr != indicators.end(); ++itr) {
        Density_descr dens = itr->first;

        auto itr2 = cWannMap.find(dens.id1);
        if (itr2 == cWannMap.end()) {  // not found
            cerr << "LFE-Density with unknown Wannier function: " <<dens.toString() << endl;
            throw runtime_error("LFE-indicator parameters contain densities with unknown conduction Wannier function.");
        }
        itr2 = vWannMap.find(dens.id2);
        if (itr2 == vWannMap.end()) {  // not found
            cerr << "LFE-Density with unknown Wannier function: " <<dens.toString() << endl;
            throw runtime_error("LFE-indicator parameters contain densities with unknown valence Wannier function.");
        }
    }
}

/**
 * @brief Calculates the center of all classical charge densities
 *
 * This is used for Density-density integrals only. Here we do not need to calculate the
 * entire indicator with all possible overlaps.
 *
 * @param WannMap       all Wannier functions for valence or conduction band (map: id->WF)
 * @param id            maps array-index to wannier-id
 * @return vector<vector<double>>  shape: [Number of WF, 3 (x,y,z)]
 */
vector<vector<double>> calcWannierCenter(map< int,WannierFunction > const& WannMap, vector<int> const& id )
{
    const RealMeshgrid* mesh = WannMap.begin()->second.getMeshgrid();

    // create vectors for Wannier centers
    vector<vector<double>> pos(id.size());

    // calculate charge centers
    for (size_t i=0; i<id.size(); i++) {

        auto itr = WannMap.find(id[i]);
        if (itr ==  WannMap.end()) {
            throw runtime_error("id and WannMap are not compatible in calcWannierCenter");
        }

        unique_ptr<double[], free_deleter>density{ joinedDensity(itr->second, itr->second, vector<int>{0,0,0}) };
        Monopole m = getMonopole(density.get(), mesh);
        pos[i] = vector<double>{m.x,m.y,m.z};
    }

    return pos;
}


#endif // DENSITY_H