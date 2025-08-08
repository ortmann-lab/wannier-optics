#ifndef POTENTIAL_H
#define POTENTIAL_H

/**
 * @file potential.h
 * @author Konrad Merkel
 * @brief Potentials that can be used within Coulomb integrals.
 * 
 */

#include <cmath>
#include <string>
using namespace std;

#define ARB_TO_EV 14.39964535

/**
 * @brief Potential that can be used in the Coulomb integrals to calculate the interaction.
 * 
 * All potentials depend on the distance (not on actual coordinates!!!). Potentials could also
 * depend on parameters.
 * 
 */
class Potential
{
private:
    /**
     * @brief Name of the potential (e.g. Coulomb, Yukawa,...)
     * 
     * This name will be used in all output files and for std-output.
     * The name may contain used parameters of the potential.
     */
    string name;

public:
    explicit Potential(string const& name_) : name(name_) {}
    virtual ~Potential() {}

    /**
     * @brief Value of the potential function given a real space distance
     * 
     * @param r     real space distance
     */
    virtual double inline real(double r) const = 0;
    double inline realCart(double x, double y, double z) const { return this->real( sqrt(x*x+y*y+z*z) ); }

    /**
     * @brief Value of the potential function given a reciprocal space distance
     * 
     * @param q     reciprocal space distance
     */
    virtual double inline fourier(double q) const = 0;
    double inline fourierCart(double qx, double qy, double qz) const { return this->fourier( sqrt(qx*qx+qy*qy+qz*qz) ); }

    string inline getName() const { return this->name; }
};


/**
 * @brief 3D unscreened Coulomb potential
 * 
 */
class CoulombPotential: public Potential
{
public:
    CoulombPotential() :Potential("Coulomb3D") {}
    double inline real(double r) const override { return ARB_TO_EV /  abs(r);}
    double inline fourier(double q) const override { return 4.*M_PI/pow(q,2.) * ARB_TO_EV; }

};

/**
 * @brief 3D Yukawa potential
 * 
 */
class YukawaPotential : public Potential
{
private:
    double alpha;  //!< screening parameter
public:
    explicit YukawaPotential(double alpha_) : Potential("Yukawa3D("  + to_string(alpha_) + ")"), alpha(alpha_) {}

    double inline real(double r) const override { return exp(-alpha * abs(r)) /  abs(r) * ARB_TO_EV; }
    double inline fourier(double q) const override { return 4*M_PI/(q*q + alpha*alpha)  * ARB_TO_EV;}

    double inline getAlpha() const { return alpha; }

    static inline double calc_yukawa_screening_factor(double mean_density, double epsilon, double alpha) {

        // general physical constants
        //const double HBAR = 65.82119569;      // reduced Planck constant in eV*fs
        const double aB = 0.52917721090380;   // Bohr radius in Angström

        // dimensionless electron gas paramter (from yellium model; see Bechstedt 2015)
        //double rs = (3./(4*np.pi*density*aB**3))**(1./3.)
        double rs = pow((3.0 / (4.0 * M_PI * mean_density * pow(aB, 3))), (1.0 / 3.0));
        // if (rank==0) cout << "Dimensionless electron gas paramter rs = "<< rs << endl;


        // Thomas-Fermi wave vector in 1/angström
        // double qTF = (12.0/np.pi)**(1./3.) / (np.sqrt(rs) * aB)
        double qTF = pow(12.0/M_PI, 1.0/3.0) / (sqrt(rs) * aB);
        // if (rank==0) cout << "Thomas-Fermi wave vector qTF = " << qTF << " Ang^-1" << endl;

        if (epsilon<=1) throw runtime_error("Cannot calculate Yukawa-Screening for epsilon <= 1");


        // Yukawa screening
        double pot_screening = qTF / sqrt((1-1./epsilon)*alpha);
        // if (rank==0) cout << "Yukawa screening factor pot_screening = " << pot_screening << " Ang^-1" << endl;

        return pot_screening;
    }

    static inline double calc_yukawa_screening_factor(double N_electron, double volume, double epsilon, double alpha) {
        // electron density
        double mean_density = N_electron / volume;
        // if (rank==0) cout << "Mean density of valence electrons = " << mean_density << " Ang^-3" << endl;

        return calc_yukawa_screening_factor(mean_density, epsilon, alpha);
    }
};

/**
 * @brief 3D isotropic Gauss function as potential
 * 
 * This is mainly for debugging purposes.
 */
class GaussPotential : public Potential
{
private:
    double a;  //!< parameter in the gauss function
public:
    explicit GaussPotential(double a_) : Potential("Gauss3D("  + to_string(a_) + ")"), a(a_) {}

    // double real(double r) const { return exp( -pow(r, 2) / (2*pow(sigma, 2))) /  pow(sqrt(2*M_PI)*sigma,3); }
    // double fourier(double q) const { return exp(-0.5*pow(q*sigma, 2)); }

    double inline real(double r) const override { return exp( -pow(a*r, 2)) * ARB_TO_EV; }
    double inline fourier(double q) const override { return pow(M_PI, 3./2) / pow(a,3) * exp(-pow(q, 2)/(4*pow(a,2))) * ARB_TO_EV; }

    double inline getA() const { return a; }
};


#endif