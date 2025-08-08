#ifndef SOLVER_H
#define SOLVER_H

/**
 * @file solver.h
 * @author Konrad Merkel
 * @brief Actual implementation of the Coulomb integrals and local field effects
 */


#include "wannierfunction.h"
#include "meshgrid.h"
#include "potential.h"
#include "algebra.h"
#include "coulombIntegral.h"

#include <complex>
#include <vector>
#include <map>
#include <fftw3.h>
#include <float.h>  // for definition of DBL_MAX
#include <math.h>
#include <mpi.h>
#include <omp.h>  // for parallelization

using namespace std;


/**
 * @brief Multiplication of two complex numbers
 */
void inline multiply(const fftw_complex &a, const fftw_complex &b, fftw_complex &result) {
    result[0] = a[0] * b[0] -a[1] * b[1] ;
    result[1] = a[0] * b[1] + a[1] * b[0];
}


/**
 * @brief Actual calculation of Coulomb and LFE integrals.
 *
 * This is the interface that any algorithm that calculates integrals needs to implement.
 * Every worker has one instance of this class which it uses to calculate a list of
 * integrals.
 * The solver does not store known or calculated integrals (this is the job of the Scheduler
 * which only lives at the master process). It only updates Integral::value.
 * However, some solvers might cache overlap densities for better performance.
 *
 */
class Solver
{
private:
    string name; //!< name of the implementation. This will be used in the output.

protected:
    bool verbose;                           //!< enable additional output
    map< int,WannierFunction > const& vWannMap;  //!< valence Wannier functions
    map< int,WannierFunction > const& cWannMap;  //!< conduction Wannier functions
    vector<vector<double>> unitcell_T;      //!< transpose of the unit cell

    void msg(string const& s) const {
        if (verbose) {
            #pragma omp critical
            {
                cout << "("<< omp_get_thread_num() << "):  " << s << endl;
            }
        }
    }

    void setName(const string& newName) {name = newName;}


public:

    /**
     * @brief Construct a new Solver object
     *
     * @param name_         Name of the scheduler that is shown in standard output or in files
     * @param newVWannMap   Wannier functions for the valence bands
     * @param newCWannMap   Wannier functions for the conduction bands
     * @param verbose_      enable more output
     */
    Solver(string const& name_, map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap, bool verbose_=false)
    : name(name_), verbose(verbose_), vWannMap(newVWannMap), cWannMap(newCWannMap)
    {
        this->unitcell_T = transpose3x3(newCWannMap.begin()->second.getUnitcell());

        // check consistency
        WannierFunction const& wann = vWannMap.begin()->second;
        for (auto const& [key, val] : vWannMap) {
            if (! val.isCompatible(wann)) throw runtime_error("Solver cannot be initialized: Wannier functions are not compatible!");
        }
        for (auto const& [key, val] : cWannMap) {
            if (! val.isCompatible(wann)) throw runtime_error("Solver cannot be initialized: Wannier functions are not compatible!");
        }
    }
    virtual ~Solver() {}
    string getName() const { return this->name; }

    void setVerbose(bool verbose_) {this->verbose=verbose_;}
    bool getVerbose() const {return this->verbose;}

    /**
     * @brief Actual calculation of integrals
     *
     * The calculation is performed in small batches (batch_size = integrals.size()), because it might be faster to calculate
     * multiple integrals in the same run and maybe cache intermediate results. However, it depends on the actual algorithm and
     * implementation if the solver can utilize such things.
     *
     * The value of each integral is stored in integrals.value.
     *
     * @param integrals         array of integrals that needs to be calculated
     * @param verbose           enables more user output if true
     * @param numOuterThreads
     * @param numInnerThreads
     */
    virtual void calculate(vector<Integral>& integrals, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) =0;

    /**
     * @brief Calculation of a single integral
     *
     * This is just a shortcut for a single integral. Please see Solver::calculate()
     *
     */
    virtual void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) {
        auto tmp = vector<Integral>{integral};
        this->calculate(tmp,verbose,numOuterThreads,numInnerThreads);
        integral.value = tmp[0].value;  // copy result
    }

};

/**
 * @brief This is just for test purposes.
 *
 * The values of the integrals are just dummy values, i.e. they are wrong!
 *
 */
class DummySolver : public Solver
{
private:
    map< int,WannierFunction > tmp{};
    /* data */
public:
    DummySolver() : Solver("Dummy",tmp,tmp) {}

    void calculate(vector<Integral>& integrals, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (size_t i=0; i<integrals.size(); i++) {
            integrals[i].value = rank  + 0.001 * i;
        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};


/**
 * @brief Solves the Coulomb integral in real space
 *
 * This is only for Coulomb integrals. It cannot solve LFE integrals!
 *
 * THIS SOLVER SHOULD ONLY BE USED FOR TEST PURPOSES! DO NOT USE IT IN PRODUCTION!
 * It directly solves the double integral (int dx int dx' ...) by splitting it
 * into small intervals and a subsequent summation. This is very slow!!!
 *
 *
 * THIS SOLVER IS NUMERICALLY UNSTABLE WHEN YOU USE A COULOMB POTENTIAL!
 * The value of integrals might diverge (become infinite) if you use a Coulomb potential
 * together with overlapping electron and hole densities. This is because the Coulomb potential
 * has a 1/|x-x'| divergence which is numerically unstable. However this is not a problem
 * if electron and hole densities don't have a shared support or the Yukawa or Gauss potential
 * (which do not diverge) are used.
 *
 */
class RealSpaceSolver : public Solver
{
private:
    int Npoints;        //!< number of points in the meshgrid (size of vDensity and cDensity array)
    double dV;          //!< volume of small interval (Anstrom^3)
    const RealMeshgrid* mesh; //!< meshgrid of the Wannier functions
    const Potential* potential = nullptr;      //!< interaction potential


    Density_descr last_vDens, last_cDens;  //!< indexes of the densities that are currently stored (for caching)
    unique_ptr<double[], free_deleter> vDensity{};  //!< valence density (cached)
    unique_ptr<double[], free_deleter> cDensity{};  //!< conduction density (cached)
public:
    RealSpaceSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap, const Potential* pot)
    : Solver("RealSpace", newVWannMap, newCWannMap), mesh{newVWannMap.begin()->second.getMeshgrid()}, potential(pot), vDensity{}, cDensity{}
    {
        Npoints = mesh->getNumDataPoints();
        dV = mesh->getdV();
        last_vDens = Density_descr();
        last_cDens = Density_descr();
    }

    void calculate(vector<Integral>& integrals, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        if (integrals.size()==0) return;

        omp_set_num_threads(numInnerThreads);
        // solving Coulomb integrals:
        for (size_t i=0; i<integrals.size(); i++){

            if (integrals[i].isEmpty()) { // no integral
                continue;
            }

            vector<int> RD{integrals[i].indexes[4],integrals[i].indexes[5],integrals[i].indexes[6]};
            vector<int> Rc{integrals[i].indexes[7],integrals[i].indexes[8],integrals[i].indexes[9]};
            vector<int> Rv{integrals[i].indexes[10],integrals[i].indexes[11],integrals[i].indexes[12]};

            Density_descr cDens(integrals[i].indexes[0], integrals[i].indexes[1], Rc);
            Density_descr vDens(integrals[i].indexes[2], integrals[i].indexes[3], Rv);


            // update valence density
            if (vDens == this->last_vDens) {
                msg("Reuse valence density!");
            } else {
                msg("Update valence-WF");
                msg(" create joined density");
                vDensity = joinedDensity(vWannMap.find(vDens.id1)->second , vWannMap.find(vDens.id2)->second ,Rv);
                this->last_vDens = vDens;
            }

            if (isZero(vDensity.get(), Npoints)) {
                integrals[i].value = 0.0;
                continue;
            }

            // update conduction density
            if (cDens == this->last_cDens) {
                msg("Reuse conduction density!");
            } else {
                msg("Update conduction-WF");
                msg(" create joined density");
                cDensity = joinedDensity(cWannMap.find(cDens.id1)->second , cWannMap.find(cDens.id2)->second ,Rc);
                this->last_cDens = cDens;
            }

            if (isZero(cDensity.get(), Npoints)) {
                integrals[i].value = 0.0;
                continue;
            }

            // actual calculation
            vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);
            double xj,yj,zj, xl,yl,zl, W=0.0;
            #pragma omp parallel for private(xj,yj,zj, xl,yl,zl) shared(cDensity, vDensity) reduction(+:W)
            for (int j=0; j<Npoints; j++){
                if (abs(cDensity[j]) < 1e-10)
                    continue;
                mesh->xyz(j,xj,yj,zj);
                for (int l=0;l<Npoints; l++){
                    if (abs(vDensity[l]) < 1e-10)
                        continue;
                    mesh->xyz(l,xl,yl,zl);
                    W += cDensity[j] * vDensity[l] * this->potential->realCart(
                        xj- xl - vec_shift[0],
                        yj- yl - vec_shift[1],
                        zj- zl - vec_shift[2]) * dV * dV;
                }
            }
            integrals[i].value = W;
        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};

/**
 * @brief Calculation of the monopole-monopole approximation of an integral
 *
 * THIS IS AN APPROXIMATION!
 * This is only for Coulomb integrals. It cannot solve LFE integrals!
 *
 * The monopole moments for valence and conduction densities are calculated.
 * They are stored (cached) for later usage. The integral value is obtained
 * from the monopole-monopole interaction. Higher multipoles are neglected.
 *
 * Only density-density integral (2-center) have a monopole-monopole interaction.
 * Therefore you should only calculate this kind of integrals.
 *
 */
class MonoMonoSolver : public Solver
{
private:
    map< Density_descr,Monopole > vMonoMap{};  //!< store calculated monopoles for valence densities
    map< Density_descr,Monopole > cMonoMap{};  //!< store calculated monopoles for conduction densities
    const CoulombPotential potential; //!< interaction potential

public:
    MonoMonoSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap)
    : Solver("MonoMono", newVWannMap, newCWannMap)
    { }

    void calculate(vector<Integral>& integrals, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        if (integrals.size()==0) return;

        omp_set_num_threads(numInnerThreads);
        // solving Coulomb integrals:
        for (size_t i=0; i<integrals.size(); i++){

            if (integrals[i].isEmpty()) { // no integral
                continue;
            }

            vector<int> RD{integrals[i].indexes[4],integrals[i].indexes[5],integrals[i].indexes[6]};
            vector<int> Rc{integrals[i].indexes[7],integrals[i].indexes[8],integrals[i].indexes[9]};
            vector<int> Rv{integrals[i].indexes[10],integrals[i].indexes[11],integrals[i].indexes[12]};

            Density_descr cDens(integrals[i].indexes[0], integrals[i].indexes[1], Rc);
            Density_descr vDens(integrals[i].indexes[2], integrals[i].indexes[3], Rv);

            Monopole vMono, cMono;

            auto itr = vMonoMap.find(vDens);
            if (itr != vMonoMap.end()) {  // density is known
                vMono = itr->second;
            } else {  // calculate density and monopole
                msg("calculate valence density");
                WannierFunction const& wann1 = vWannMap.find(vDens.id1)->second;
                WannierFunction const& wann2 = vWannMap.find(vDens.id2)->second;
                unique_ptr<double[], free_deleter>density{ joinedDensity(wann1, wann2, Rv) };
                vMono = getMonopole(density.get(), wann1.getMeshgrid());
                vMonoMap.insert({vDens, vMono});
            }

            if (abs(vMono.charge)< 1e-10)  {
                integrals[i].value = 0.0;
                continue;
            }

            itr = cMonoMap.find(cDens);
            if (itr != cMonoMap.end()) {  // density is known
                cMono = itr->second;
            } else {  // calculate density and monopole
                msg("calculate conduction density");
                WannierFunction const& wann1 = cWannMap.find(cDens.id1)->second;
                WannierFunction const& wann2 = cWannMap.find(cDens.id2)->second;
                unique_ptr<double[], free_deleter>density{ joinedDensity(wann1, wann2, Rc) };
                cMono = getMonopole(density.get(), wann1.getMeshgrid());
                cMonoMap.insert({cDens, cMono});
            }

            if (abs(cMono.charge)< 1e-10)  {
                integrals[i].value = 0.0;
                continue;
            }

            // calculate monopol-monopol interaction
            vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);
            integrals[i].value = vMono.charge * cMono.charge * this->potential.realCart(cMono.x-vMono.x-vec_shift[0], cMono.y-vMono.y-vec_shift[1], cMono.z-vMono.z-vec_shift[2]);
        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};

/**
 * @brief Solver of the local field effects integral (LFE)
 *
 * This is only for LFE integrals. It cannot solve Coulomb integrals!
 *
 * The calculation is done in Fourier space with a short range Coulomb potential.
 * You can only use the Coulomb potential. It automatically than assumes the short
 * range version. Other potentials are not supported.
 *
 */
class LocalFieldEffectsSolver : public Solver
{
private:
    Density_descr last_Dens1, last_Dens2;  //!< labels of the last densities (used for caching)
    ReciprocalMeshgrid recMesh;
    const CoulombPotential potential;      //!< interaction potential
    double dV;
    double V_unitcell;
    int N;
    vector<int> dim, supercell;

    fftw_plan p_f1, p_f2;
    fftw_complex* f1;           //!< fourier transformed densities (between electron and hole), cached
    fftw_complex* f2;           //!< fourier transformed densities (between electron and hole), cached

public:
    LocalFieldEffectsSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap)
    : Solver("LocalFieldEffects", newVWannMap, newCWannMap), recMesh{newVWannMap.begin()->second.getMeshgrid()},
      supercell{vector<int>(3)}
    {
        // check supercell
        vector<double> supercell_d = newVWannMap.begin()->second.getLatticeInUnitcellBasis();
        for (int j=0;j<3;j++) {
            supercell[j] = round(supercell_d[j]);   // will be important to have an integer

            if (abs(supercell_d[j]-supercell[j]) > 1e-5) {
                throw runtime_error("Supercell is not compatible with unitcell. Cannot calculate local field effects.");
            }
        }

        // reset memory in case it was initialized before
        this->last_Dens1 = Density_descr();
        this->last_Dens2 = Density_descr();

        // get information about meshgrid und cell
        RealMeshgrid const* realMesh{newVWannMap.begin()->second.getMeshgrid()};
        dV = realMesh->getdV();
        V_unitcell = newVWannMap.begin()->second.getVunitcell();
        N = realMesh->getNumDataPoints();
        dim = recMesh.getDim();

        // allocate arrays for densities and FFT
        // use same arrays for real and reciprocal space (inplace transformation)
        f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        // plan needs to be created before initializing inputs !!!
        // dimensions are reversed because our datagrids are in Fortran (column-major) format!
        // planing is not thread safe unless void fftw_make_planner_thread_safe(void); is invoked first!!!
        #pragma omp critical
        {
            p_f1 = fftw_plan_dft_3d(dim[2], dim[1], dim[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE); // f1(+q) FORWARD corresponds to exp(-iqx)
            p_f2 = fftw_plan_dft_3d(dim[2], dim[1], dim[0], f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);  // f2(-q)  FORWARD corresponds to exp(+iqx)
        }
    }
    ~LocalFieldEffectsSolver() {
        fftw_destroy_plan(p_f1);
        fftw_destroy_plan(p_f2);
        if (f1!=nullptr) fftw_free(f1);
        if (f2!=nullptr) fftw_free(f2);
    }

    void calculate(vector<Integral>& integrals, const bool verbose=true,
        const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        if (integrals.size()==0) return;

        omp_set_num_threads(numInnerThreads);
        // solving Coulomb integrals:
        for (size_t i=0; i<integrals.size(); i++){

            if (integrals[i].isEmpty()) { // no integral
                continue;
            }

            vector<int> RD{integrals[i].indexes[4],integrals[i].indexes[5],integrals[i].indexes[6]};
            if ((RD[0]!=0) || (RD[1]!=0) || (RD[2]!=0)) {
                integrals[i].setFailed("RD has to be zero when calculating local field effects!");
                // throw runtime_error("RD has to be zero when calculating local field effects!");
                continue;
            }

            // Nomenclature:
            // We use Rc=S1, Rv=S2, RD=0 in the integrals object. With this definition we make sure that the obtained
            // LFE file has S1 and S2 as columns (instead of -S1 and -S2).
            // For calculating the overlap densities using the joinedDensity function this means that we need to use
            // the negative shift vector
            // rho(x) = w_c1,0(x) w_v1,R1(x) = w_c1,0(x) w_v1,-S1(x)  ==> R1 = -S1 and S1=Rc
            vector<int> R1{ -integrals[i].indexes[7], -integrals[i].indexes[8], -integrals[i].indexes[9]}; // =-S1
            vector<int> R2{ -integrals[i].indexes[10],-integrals[i].indexes[11],-integrals[i].indexes[12]}; // =-S2

            // overlap densities between valence and conduction WF!
            Density_descr Dens1(integrals[i].indexes[0], integrals[i].indexes[2], R1);
            Density_descr Dens2(integrals[i].indexes[1], integrals[i].indexes[3], R2);

            // update density1 and corresponding Fourier quantities
            if (Dens1 == this->last_Dens1) {
                msg("Reuse 1st joined density (c1,v1,S1))!");
                if (isZero(f1, N)) {
                    msg(" joined density is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }
            } else {
                msg("Update 1st joined density (c1,v1,S1)");
                msg(" create joined density");
                WannierFunction const& c1 = cWannMap.find(Dens1.id1)->second;
                WannierFunction const& v1 = vWannMap.find(Dens1.id2)->second;

                joinedDensity(c1,v1,Dens1.R,f1);
                this->last_Dens1 = Dens1;

                if (isZero(f1, N)) {
                    msg(" joined density is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }else{
                    // Fourier transform wave functions
                    msg(" perform Fourier transform");
                    fftw_execute(p_f1);
                }
            }

            // update density2 and corresponding Fourier quantities
            if (Dens2 == this->last_Dens2) {
                msg("Reuse 2nd joined density (c2,v2,S2)!");
                if (isZero(f2, N)) {
                    msg(" joined density is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }
            } else {
                msg("Update 2nd joined density (c2,v2,S2)");
                msg(" create joined density");
                WannierFunction const& c2 = cWannMap.find(Dens2.id1)->second;
                WannierFunction const& v2 = vWannMap.find(Dens2.id2)->second;

                joinedDensity(c2,v2,Dens2.R,f2);
                this->last_Dens2 = Dens2;

                if (isZero(f2, N)) {
                    msg(" joined density is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }else{
                    // Fourier transform wave functions
                    msg(" perform Fourier transform");
                    fftw_execute(p_f2);
                }
            }


            //
            // perform fourier sums
            //
            double qx,qy,qz;
            fftw_complex tmp, mat_rec_ana;
            mat_rec_ana[0] = 0.0;
            mat_rec_ana[1] = 0.0;
            msg(" - Perform Fourier sums ...");
            int l;
            for (int i_=0; i_<dim[0]; i_+=supercell[0]) {
                for (int j=0; j<dim[1]; j+=supercell[1]) {
                    for (int k=0; k<dim[2]; k+=supercell[2]) {

                        if ((i_==0) && (j==0) && (k==0)) continue;  // short range Coulomb potential v(G=0)=0

                        recMesh.xyz(i_,j,k,qx,qy,qz);
                        l = recMesh.getGlobId(i_,j,k);
                        multiply(f1[l], f2[l], tmp);

                        mat_rec_ana[0] += tmp[0] * potential.fourierCart(qx, qy, qz) * dV * dV / pow(2*M_PI, 3);
                        mat_rec_ana[1] += tmp[1] * potential.fourierCart(qx, qy, qz) * dV * dV / pow(2*M_PI, 3);
                    }
                }
            }

            integrals[i].value = mat_rec_ana[0] / V_unitcell;
            // Reason for the V_unitcell prefactor:
            //   Short range Coulomb potental uses continues Fourier transform with normalization Omega = Volume of supercell
            //   Wannier functions use discreate Fourier transform with normalization N_Omega = Number of prim. unit cells in supercell
            //   --> prefactor: N_Omega / Omega = 1./V_unitcell

            // the integral should be real, because the Wannier functions are real.
            if(abs(mat_rec_ana[1]) > 1e-8 ) {
                cerr << "Integral value = " << mat_rec_ana[0] << " +i* " << mat_rec_ana[1] << endl;
                throw runtime_error("LFE integral has non-zero imaginary part!");
            }

        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};


/**
 * @brief Coulomb integral solver in Fourier space
 *
 * This is only for Coulomb integrals. It cannot solve LFE integrals!
 *
 * This is the simplest implementation of a solver in Fourier space. It is
 * similar to the RealSpaceSolver but much faster.
 *
 * THIS SOLVER IS NUMERICALLY UNSTABLE WHEN YOU USE A COULOMB POTENTIAL!
 * This is because the Coulomb potential is has a 1/q^2 divergence in Fourier
 * space. This solver is not able to deal with this correctly!
 * However, this is not a problem with other (non-divergent) potentials.
 *
 * THIS SOLVER SHOULD ONLY BE USED WITH YUKAWA OR GAUSS POTENTIALS!
 *
 */
class YukawaSolver : public Solver
{
private:
    ReciprocalMeshgrid recMesh;
    double dV;
    int N;
    vector<double> origin;
    vector<double> supercell;
    const Potential* potential{};      //!< interaction potential
    map<int, double> const vMeanDensity;
    map<int, double> const cMeanDensity;
    double RELATIVE_PERMITTIVITY, SCREENING_ALPHA;

    fftw_plan p_f1, p_f2;
    fftw_complex* f1{};       //!< fourier transformed densities (cached)
    fftw_complex* f2{};       //!< fourier transformed densities (cached)

    Density_descr last_vDens, last_cDens;  //!< labels of the last densities (used for caching)

protected:

    /**
     * @brief Checks for aliasing when doing a Fourier method
     *
     * To protect against aliasing you can use padding as implemented in createLargerSupercell()
     *
     * @param wann
     * @param R unit cell vector of the Wannier function w_mR(x)
     */
    bool checkAliasing(const vector<int>& R) {
        for (int i=0; i<3; i++) {
            if (supercell[i] < abs(R[i])/2) {
                return false;
                //throw runtime_error("The supercell is not large enough to protect against aliasing!");
            }
        }
        return true;
    }

    bool checkMaps(map<int, WannierFunction> const& WannMap, map<int, double> const& MeanDensity) {
        // Check if any of the keys in vWannMap are missing from MeanDensity
        for (const auto& entry : WannMap) {
            int key = entry.first;

            // If the key is not found in vMeanDensity, return false
            if (MeanDensity.find(key) == MeanDensity.end()) {
                return false;
            }
        }

        // All keys in WannMap are found in MeanDensity
        return true;
    }

    YukawaSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap, const Potential* pot,
         map<int, double> const& newVMeanDensity, map<int, double> const& newCMeanDensity,
        double RELATIVE_PERMITTIVITY_, double SCREENING_ALPHA_)
    : Solver("Fourier",newVWannMap, newCWannMap), recMesh{newVWannMap.begin()->second.getMeshgrid()}, potential{pot},
      vMeanDensity{newVMeanDensity}, cMeanDensity{newCMeanDensity},
      RELATIVE_PERMITTIVITY(RELATIVE_PERMITTIVITY_), SCREENING_ALPHA(SCREENING_ALPHA_)
    {
        this->last_vDens = Density_descr();
        this->last_cDens = Density_descr();

        RealMeshgrid* realMesh = newVWannMap.begin()->second.getMeshgrid();
        this->supercell = newVWannMap.begin()->second.getLatticeInUnitcellBasis();

        // get information about meshgrid
        dV = realMesh->getdV();
        N = realMesh->getNumDataPoints();
        origin = realMesh->getOrigin();

        // allocate arrays for densities and FFT
        // use same arrays for real and reciprocal space (inplace transformation)
        f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        // create plans
        #pragma omp critical
        {
            p_f1 = fftw_plan_dft_3d(recMesh.getDim()[2], recMesh.getDim()[1], recMesh.getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE); // f1(+q) FORWARD corresponds to exp(-iqx)
            p_f2 = fftw_plan_dft_3d(recMesh.getDim()[2], recMesh.getDim()[1], recMesh.getDim()[0], f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);  // f2(-q)  FORWARD corresponds to exp(+iqx)
        }
    }

public:

    // use generalized model screening that depends on the Wannier functions
    YukawaSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap,
        map<int, double> const& newVMeanDensity, map<int, double> const& newCMeanDensity,
        double RELATIVE_PERMITTIVITY_, double SCREENING_ALPHA_)
    : YukawaSolver(newVWannMap, newCWannMap, nullptr, newVMeanDensity, newCMeanDensity, RELATIVE_PERMITTIVITY_, SCREENING_ALPHA_)
    {
        // check if MeanDensities are compatible with WannMap
        if (((!checkMaps(vWannMap, vMeanDensity)) || (!checkMaps(cWannMap, cMeanDensity)))) {
            throw runtime_error("MeanDensity maps are not compatible with Wannier maps.");
        }

    }

    // use a predefined potential (either yukawa potential or gauss potential for testing)
    YukawaSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap, const Potential* pot)
    : YukawaSolver(newVWannMap, newCWannMap, pot, {}, {}, 0.0, 0.0)
    {
        if (pot->getName() == "Coulomb3D") throw runtime_error("YukawaSolver cannot work with Coulomb potentials. Use the CoulombSolver instead!");
    }

    ~YukawaSolver() {
        // clean up
        fftw_destroy_plan(p_f1); fftw_destroy_plan(p_f2);
        fftw_free(f1); fftw_free(f2);
    }

    void calculate(vector<Integral>& integrals, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        if (integrals.size()==0) return;

        omp_set_num_threads(numInnerThreads);
        // solving Coulomb integrals:
        for (size_t i=0; i<integrals.size(); i++){

            if (integrals[i].isEmpty()) { // no integral
                continue;
            }

            vector<int> RD{integrals[i].indexes[4],integrals[i].indexes[5],integrals[i].indexes[6]};
            vector<int> Rc{integrals[i].indexes[7],integrals[i].indexes[8],integrals[i].indexes[9]};
            vector<int> Rv{integrals[i].indexes[10],integrals[i].indexes[11],integrals[i].indexes[12]};

            Density_descr cDens(integrals[i].indexes[0], integrals[i].indexes[1], Rc);
            Density_descr vDens(integrals[i].indexes[2], integrals[i].indexes[3], Rv);


            // prepare shift vectors
            if (!checkAliasing(RD)) {
                integrals[i].setFailed("The supercell is not large enough to protect against aliasing!");
                continue;
            }


            // update valence density and corresponding Fourier quantities
            if (vDens == this->last_vDens)
            {
                msg("Reuse valence density");
                if (isZero(f2, N)) {
                    msg(" joined density (valence) is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }
            } else {
                msg("Update valence-WF");
                msg(" create joined density");
                joinedDensity(vWannMap.find(vDens.id1)->second , vWannMap.find(vDens.id2)->second ,Rv, f2);
                this->last_vDens = vDens;

                if (isZero(f2, N)) {
                    msg(" joined density (valence) is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }else{
                    // Fourier transform wave functions
                    msg(" perform Fourier transform");
                    fftw_execute(p_f2);
                }
            }

            // update conduction density and corresponding Fourier quantities
            if (cDens == this->last_cDens)
            {
                msg("Reuse conducton density and auxiliary functions!");
                if (isZero(f1, N)) {
                    msg(" joined density (conduction) is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }
            } else {
                msg("Update conduction-WF");
                msg(" create joined density");
                joinedDensity(cWannMap.find(cDens.id1)->second , cWannMap.find(cDens.id2)->second,Rc,f1);
                this->last_cDens = cDens;

                if (isZero(f1, N)) {
                    msg(" joined density (conduction) is zero.");
                    integrals[i].value = 0.0;
                    continue;
                }else{

                    // Fourier transform wave functions
                    msg(" perform Fourier transform");
                    fftw_execute(p_f1);
                }
            }

            // set potential
            const Potential* pot = nullptr;
            if (this->potential != nullptr){
                pot=this->potential;
            }else{  // use model screening where the charge density depends on the Wannier functions
                double nv1 = vMeanDensity.at(vDens.id1);
                double nv2 = vMeanDensity.at(vDens.id2);
                double nc1 = cMeanDensity.at(cDens.id1);
                double nc2 = cMeanDensity.at(cDens.id2);

                double mean_density = sqrt(sqrt(nv1*nv2) * sqrt(nc1*nc2));
                double yuk_exp_factor = YukawaPotential::calc_yukawa_screening_factor(mean_density, RELATIVE_PERMITTIVITY, SCREENING_ALPHA);
                // cout << "mean density: " << mean_density << "   yuk_exp_factor: " << yuk_exp_factor << endl;
                pot = new YukawaPotential(yuk_exp_factor);
            }

            // perform fourier sums
            msg("Summation in Fourier space ...");
            vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);
            fftw_complex tmp, tmp2, phase;
            double qx,qy,qz;
            fftw_complex mat_rec_ana;
            mat_rec_ana[0] = 0.;
            mat_rec_ana[1] = 0.;
            for (int k=0; k<N; k++){  // TODO: OMP parallelization
                recMesh.xyz(k,qx,qy,qz);
                phase[0] = cos(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);
                phase[1] = sin(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);

                multiply(f1[k], f2[k], tmp2);
                multiply(tmp2, phase, tmp);
                mat_rec_ana[0] += tmp[0] * pot->fourierCart(qx, qy, qz) *dV /N;
                mat_rec_ana[1] += tmp[1] * pot->fourierCart(qx, qy, qz) *dV /N;
            }

            if (abs(mat_rec_ana[1]) > 1e-8) {
                cerr << "Integral = " << mat_rec_ana[0] << " +i* " << mat_rec_ana[1] << endl;
                throw runtime_error("Imaginary part of integral should be zero!");
            }

            integrals[i].value = mat_rec_ana[0];

            if (this->potential != pot){ delete pot; }
        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};


/**
 * @brief Coulomb integral solver in Fourier space using auxiliary Gauss functions
 *
 * This is only for Coulomb integrals. It cannot solve LFE integrals!
 *
 * This is a numerical stable implementation for the Coulomb potential. It cannot be used
 * for other potentials.
 *
 * For more information see:
 *      Konrad Merkel and Frank Ortmann 2024 J. Phys. Mater. 7 015001
 *      DOI 10.1088/2515-7639/ad06cd
 *
 */
class CoulombSolver : public Solver
{
private:
    // parameter for auxiliary Gaussians
    const double STD_PER_CELL       = 11;    //!< number of standard deviations that have to be contained along the smallest lattice vector
    const double POINTS_PER_STD     = 2;     //!< number of mesh points per standard deviation of the gaussian
    const double THRESHOLD_DENSITY  = 1e-7;  //!< threshold for f(q=0) to throw an error
    const bool wrapAux;                      //!< wrapping around auxiliary function

    const RealMeshgrid* realMesh{};
    ReciprocalMeshgrid recMesh;
    double dV, std_min, std_max;
    int N;
    vector<double> origin;
    vector<double> supercell;
    const CoulombPotential potential;      //!< interaction potential


    fftw_plan p_f1, p_f2;
    fftw_complex* f1{};                       //!< fourier transformed densities (cached)
    fftw_complex* f2{};                       //!< fourier transformed densities (cached)

    Monopole m1{}, m2{};
    double alpha{}, gamma{};

    Density_descr last_vDens, last_cDens;   //!< labels of the last densities (used for caching)

protected:

    double inline gaussian3D(double x, double y, double z, double x0, double y0, double z0, double alpha, double norm=1.) const {
        return norm * pow(alpha / M_PI, 3./2.) * exp(-alpha * (  pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2) ) );
    }

    void inline recGaussian3D(double qx, double qy, double qz, double x0, double y0, double z0, double alpha, double norm, fftw_complex& result) const {
        result[0] = norm * exp( -(pow(qx, 2) + pow(qy, 2) + pow(qz, 2))/(4.*alpha) ) *cos(qx * x0 + qy * y0 + qz * z0);
        result[1] = norm * exp( -(pow(qx, 2) + pow(qy, 2) + pow(qz, 2))/(4.*alpha) ) *-sin(qx * x0 + qy * y0 + qz * z0); // minus because of FORWARD transform --> exp(-iqx)
    }

    /**
     * @brief Checks for aliasing when doing a Fourier method
     *
     * To protect against aliasing you can use padding as implemented in createLargerSupercell()
     *
     * @param wann
     * @param R unit cell vector of the Wannier function w_mR(x)
     */
    bool checkAliasing(const vector<int>& R) {
        for (int i=0; i<3; i++) {
            if (supercell[i] < abs(R[i])/2) {
                return false;
                //throw runtime_error("The supercell is not large enough to protect against aliasing!");
            }
        }
        return true;
    }

    void add_auxiliary(RealMeshgrid const* realMesh, fftw_complex* data, Monopole mono, double alpha ) const
    {
        double x,y,z;
        double norm = mono.charge * pow(alpha / M_PI, 3./2.);
        //int N = realMesh->getNumDataPoints();

        if (wrapAux) {
            vector<vector<double>> shifts(27);
            int l=0;

            // all possible shift vectors
            vector< vector<double> > lattice_T = transpose3x3(realMesh->getLattice());
            for (int i=-1; i<=1; i++) {
                for (int j=-1; j<=1; j++) {
                    for (int k=-1; k<=1; k++) {
                        shifts[l] = matVecMul3x3(lattice_T,vector<int>{i,j,k});
                        l++;
                    }
                }
            }

            double param;
            #pragma omp parallel for shared(data,realMesh, shifts) firstprivate(alpha, N, mono, norm) private(x,y,z, param)
            for (int i=0; i<N; i++){
                realMesh->xyz(i,x,y,z);  // get coordinates
                param = DBL_MAX;

                for (int l=0; l<27; l++){
                    param = min(param,  pow(x-mono.x-shifts[l][0], 2) + pow(y-mono.y-shifts[l][1], 2) + pow(z-mono.z-shifts[l][2], 2) );
                }

                // real part
                data[i][0] -= norm *  exp(-alpha * param);
            }
        } else {

            #pragma omp parallel for shared(data,realMesh) firstprivate(alpha, N, mono, norm) private(x,y,z)
            for (int i=0; i<N; i++){
                realMesh->xyz(i,x,y,z);  // get coordinates

                // real part
                data[i][0] -= norm *  exp(-alpha *  (pow(x-mono.x, 2) + pow(y-mono.y, 2) + pow(z-mono.z, 2)));

            }
        }
    }


    bool prepareDensity2() {

        msg(" create auxiliary gaussian");
        m2 = getMonopole(f2,realMesh);
        double std2 = sqrt(getMeanVariance(f2, realMesh, vector<double>{m2.x,m2.y,m2.z}));

        // enforce boundaries for standard deviation
        if (std2<std_min) std2=std_min;
        else if (std2>std_max) std2=std_max;
        gamma = 1./(2*pow(std2,2));
        add_auxiliary(realMesh, f2, m2, gamma );

        // Fourier transform wave functions
        msg(" perform Fourier transform");
        fftw_execute(p_f2);

        // check if overall charge of density-gaussian is actually zero
        if (abs(f2[0][0]*dV) > THRESHOLD_DENSITY) {
            cerr << "f2(q=0) is not zero! f2(q=0) = " << f2[0][0]*dV << endl;
            // throw runtime_error("Solver failed!");
            return false;
        }
        return true;
    }

    bool prepareDensity1() {
        msg(" create auxiliary gaussian");
        m1 = getMonopole(f1,realMesh);
        double std1 = sqrt(getMeanVariance(f1, realMesh, vector<double>{m1.x,m1.y,m1.z}));

        // enforce boundaries for standard deviation
        if (std1<std_min) std1=std_min;
        else if (std1>std_max) std1=std_max;
        alpha = 1./(2*pow(std1,2));


        // TODO Improve: check if m1 and m2 are separated by a unit cell vector.
        // if so then add 0.01 (this should be valid for all choices of RD)

        // distance between auxiliary Gaussians
        //double r = sqrt(( pow(m1.x-m2.x-vec_shift[0],2) + pow(m1.y-m2.y-vec_shift[1],2) + pow(m1.z-m2.z-vec_shift[2],2) ));
        // double r = sqrt(( pow(m1.x-m2.x,2) + pow(m1.y-m2.y,2) + pow(m1.z-m2.z,2) ));
        // if (r<1e-10) {  // shift to make sure we do not divide by zero later
        //     m1.x = m1.x + 0.01;
        //     m1.y = m1.y + 0.01;
        //     m1.z = m1.z + 0.01;
        // }


        add_auxiliary(realMesh, f1, m1, alpha );

        // Fourier transform wave functions
        msg(" perform Fourier transform");
        fftw_execute(p_f1);

        // check if overall charge of density-gaussian is actually zero
        if (abs(f1[0][0]*dV) > THRESHOLD_DENSITY) {
            cerr << "f1(q=0) is not zero! f1(q=0) = " << f1[0][0]*dV <<endl;
            //throw runtime_error("Solver failed!");
            return false;
        }
        return true;
    }

    complex<double> fourier_sum(vector<int> const& RD) {
        // check if overall charge of density-gaussian is actually zero
        if (abs(f1[0][0]*dV) > THRESHOLD_DENSITY) {
            // TODO: failed
            cerr << "[FAIL] abs(f1[0][0]*dV) > THRESHOLD_DENSITY" << endl;
        }

        vector<double> vec_shift = matVecMul3x3(unitcell_T, RD);

        //
        // perform fourier sums
        //
        msg("Summation in Fourier space ...");
        fftw_complex mat_rec_ana;
        mat_rec_ana[0] = 0.;
        mat_rec_ana[1] = 0.;
        double x,y,z, V_q;
        fftw_complex tmp2, tmp, phase, gauss;

        #pragma omp parallel for private(x,y,z, phase, tmp, tmp2, gauss, V_q) shared(f1, f2) firstprivate(N, alpha, gamma, m1, m2, vec_shift, dV, origin) reduction(+:mat_rec_ana)
        for (int k=1; k<N; k++){
            recMesh.xyz(k,x,y,z);                   // get reciprocal coordinates
            V_q = potential.fourierCart(x,y,z);    // potential in reciprocal space
            phase[0] = cos(x * vec_shift[0] + y * vec_shift[1] + z * vec_shift[2]);  // phase due to shift of unit cells
            phase[1] = sin(x * vec_shift[0] + y * vec_shift[1] + z * vec_shift[2]);

            // 1st term (density1-density2)
            multiply(f1[k], f2[k], tmp2);
            multiply(tmp2, phase, tmp);
            mat_rec_ana[0] += tmp[0] * V_q *dV /N;
            mat_rec_ana[1] += tmp[1] * V_q *dV /N;

            // 2nd term (density1-gauss2)
            recGaussian3D(-x,-y,-z, m2.x-origin[0], m2.y-origin[1], m2.z-origin[2], gamma, m2.charge, gauss);
            // shift by origin is necessary because f1[k] and f2[k] dont know the origin.
            multiply(f1[k], gauss, tmp2);
            multiply(tmp2, phase, tmp);
            mat_rec_ana[0] += tmp[0] * V_q /N;
            mat_rec_ana[1] += tmp[1] * V_q /N;

            // // 3rd term (gauss1-density2)
            recGaussian3D(x,y,z, m1.x-origin[0], m1.y-origin[1], m1.z-origin[2], alpha, m1.charge, gauss);
            multiply(gauss, f2[k], tmp2);
            multiply(tmp2, phase, tmp);
            mat_rec_ana[0] += tmp[0] * V_q /N;
            mat_rec_ana[1] += tmp[1] * V_q /N;
        }

        // distance between auxiliary functions
        double r = sqrt(( pow(m1.x-m2.x-vec_shift[0],2) + pow(m1.y-m2.y-vec_shift[1],2) + pow(m1.z-m2.z-vec_shift[2],2) ));

        // 4th term (gauss1 - gauss2) --> analytic
        double analytic_part = ARB_TO_EV * m1.charge* m2.charge* erf(sqrt(alpha*gamma/(alpha+gamma)) * r) / r;
        mat_rec_ana[0] += analytic_part;

        return complex<double>(mat_rec_ana[0], mat_rec_ana[1]);
    }



public:
    CoulombSolver(map< int,WannierFunction > const& newVWannMap, map< int,WannierFunction > const& newCWannMap, bool wrapAux_=true)
    : Solver("FourierGauss",newVWannMap, newCWannMap), wrapAux(wrapAux_), realMesh{newVWannMap.begin()->second.getMeshgrid()},
      recMesh{ newVWannMap.begin()->second.getMeshgrid() }
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // reset memory in case it was initialized before
        this->last_vDens = Density_descr();
        this->last_cDens = Density_descr();
        this->supercell = newVWannMap.begin()->second.getLatticeInUnitcellBasis();

        // get information about meshgrid
        dV = realMesh->getdV();
        N = realMesh->getNumDataPoints();
        origin = realMesh->getOrigin();

        // determine smallest lattice vector and corresponding discretization length
        vector<double> uvec = matVecMul3x3(realMesh->getLattice(), vector<int>{1,0,0});  // x
        double length = sqrt(pow(uvec[0],2) + pow(uvec[1],2) + pow(uvec[2],2));
        double min_length = length;
        double discret_length = length/realMesh->getDim()[0];

        uvec = matVecMul3x3(realMesh->getLattice(), vector<int>{0,1,0});  // y
        length = sqrt(pow(uvec[0],2) + pow(uvec[1],2) + pow(uvec[2],2));
        discret_length = min(discret_length, length/realMesh->getDim()[1]);
        min_length = min(min_length, length);

        uvec = matVecMul3x3(realMesh->getLattice(), vector<int>{0,0,1});  // z
        length = sqrt(pow(uvec[0],2) + pow(uvec[1],2) + pow(uvec[2],2));
        discret_length = min(discret_length, length/realMesh->getDim()[2]);
        min_length = min(min_length, length);

        // determine smallest and largest possible standard deviations for auxiliary Gaussians
        // that meets our restrictions
        this->std_min = POINTS_PER_STD * discret_length;
        this->std_max = min_length / STD_PER_CELL;

        if (rank==0) {
            cout << "\nNumerical parameter of the solver:\n";
            cout << "POINTS_PER_STD    = " << POINTS_PER_STD << endl;
            cout << "STD_PER_CELL      = " << STD_PER_CELL << endl;
            cout << "THRESHOLD_DENSITY = " << THRESHOLD_DENSITY << endl;
            cout << "wrap aux. func.   = " << wrapAux << endl;
            cout << "\nFor the given meshgrid we get:\n";
            cout << "length of smallest lattice vector = " << min_length << endl;
            cout << "length per mesh point             = " << discret_length << endl;
            cout << "std_min                           = " << std_min << endl;
            cout << "std_max                           = " << std_max << endl;
        }

        if (std_max<=std_min) {  // check consistency
            throw runtime_error("std_min > std_max : Not possible to construct any auxiliary function!");
        }

        // allocate arrays for densities and FFT
        // use same arrays for real and reciprocal space (inplace transformation)
        f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        // create plans
        #pragma omp critical
        {
            p_f1 = fftw_plan_dft_3d(recMesh.getDim()[2], recMesh.getDim()[1], recMesh.getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE); // f1(+q) FORWARD corresponds to exp(-iqx)
            p_f2 = fftw_plan_dft_3d(recMesh.getDim()[2], recMesh.getDim()[1], recMesh.getDim()[0], f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);  // f2(-q)  FORWARD corresponds to exp(+iqx)
        }
    }

    ~CoulombSolver() {
        // clean up
        //fftw_cleanup_threads();
        fftw_destroy_plan(p_f1); fftw_destroy_plan(p_f2);
        fftw_free(f1); fftw_free(f2);
    }

    void calculate(vector<Integral>& integrals, const bool verbose=false, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        if (integrals.size()==0) return;

        omp_set_num_threads(numInnerThreads);
        // solving Coulomb integrals:
        for (size_t i=0; i<integrals.size(); i++){

            if (integrals[i].isEmpty()) { // no integral
                continue;
            }

            vector<int> RD{integrals[i].indexes[4],integrals[i].indexes[5],integrals[i].indexes[6]};
            vector<int> Rc{integrals[i].indexes[7],integrals[i].indexes[8],integrals[i].indexes[9]};
            vector<int> Rv{integrals[i].indexes[10],integrals[i].indexes[11],integrals[i].indexes[12]};

            Density_descr cDens(integrals[i].indexes[0], integrals[i].indexes[1], Rc);
            Density_descr vDens(integrals[i].indexes[2], integrals[i].indexes[3], Rv);


            // prepare shift vectors
            if (!checkAliasing(RD)) {
                integrals[i].setFailed("The supercell is not large enough to protect against aliasing!");
                continue;
            }


            //
            // update valence density and corresponding Fourier quantities
            //
            if (vDens == this->last_vDens)
            {
                msg("Reuse valence density and auxiliary functions!");
            } else {
                msg("Update valence-WF");
                msg(" create joined density");
                joinedDensity(vWannMap.find(vDens.id1)->second , vWannMap.find(vDens.id2)->second ,Rv, f2);

                // add gauss and perform FFT
                if (!prepareDensity2()) {
                    integrals[i].setFailed("Solver failed! f2(q=0) is not zero!");
                    this->last_vDens = Density_descr();  // reset memory
                    continue;
                }
                this->last_vDens = vDens;
            }

            //
            // update conduction density and corresponding Fourier quantities
            //
            if (cDens == this->last_cDens)
            {
                msg("Reuse conduction density and auxiliary functions!");
            } else {
                msg("Update conduction-WF");
                msg(" create joined density");
                joinedDensity(cWannMap.find(cDens.id1)->second , cWannMap.find(cDens.id2)->second,Rc,f1);

                // add gauss and perform FFT
                if (!prepareDensity1(/*RD*/)) {
                    integrals[i].setFailed("Solver failed! f1(q=0) is not zero!");
                    this->last_cDens = Density_descr();  // reset memory
                    continue;
                }
                this->last_cDens = cDens;
            }


            // Fourier sum
            complex<double> W = fourier_sum(RD);
            integrals[i].value = W.real();

            if(abs(W.imag()) > 1e-8 ) {
                cerr << "Integral value = " << W.real() << " +i* " << W.imag() << endl;
                throw runtime_error("Coulomb integral has non-zero imaginary part!");
            }
        }
    }

    void calculate(Integral& integral, const bool verbose=true, const unsigned int numOuterThreads=1, const unsigned int numInnerThreads=1) override
    {
        Solver::calculate(integral, verbose, numOuterThreads, numInnerThreads);
    }
};


#endif  //SOLVER_H