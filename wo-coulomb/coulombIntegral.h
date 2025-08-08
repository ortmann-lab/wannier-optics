#ifndef COULOMBINTEGRAL_H
#define COULOMBINTEGRAL_H

/**
 * @file coulombIntegral.h
 * @author Konrad Merkel
 * @brief Data structures for single Coulomb integrals and local field effects integrals
 *
 */

#include "density.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cassert>
#include <limits.h>
#include <cmath>
using namespace std;


enum class IntegralStatus { normal, failed };


/**
 * @brief Data structure for a single Coulomb / LFE integral
 *
 * Every integral can be described by its Wannier function ids and relative shift
 * vectors (in terms of unit cells). We use the same data structure for Coulomb and
 * LFE integrals. However the role of the indexes is slightly different.
 *
 * For Coulomb integrals: (see J. Phys. Mater. 7 015001 (2024))
 * c1,c2  ... conduction Wannier function ids
 * v1,v2  ... valence Wannier function ids
 * RD     ... shift between valence and conduction densities (in terms of primitive unit cells)
 * Rc     ... shift between conduction densities (in terms of primitive unit cells)
 * Rv     ... shift between valence densities (in terms of primitive unit cells)
 *
 * For LFE integrals: (see J. Phys. Mater. 7 015001 (2024))
 * c1,c2  ... conduction Wannier function ids
 * v1,v2  ... valence Wannier function ids
 * RD=0   ... has no meaning, but needs to be zero
 * Rc=S1  ... relative shift between WF of density 1 (c1,v1,S1)
 * Rv=S2  ... relative shift between WF of density 2 (c2,v2,S2)
 *
 * The value of the integral is stored in Integral::value. It is NAN if it was
 * not calculated.
 *
 * Every integral has a status and error_msg that is used when the solver fails
 * to calculate the integral.
 *
 * All integrals with c1<0 (or indexes[0] < 0) are not considered as usual integrals
 * but as placeholders, i.e. empty integral or NO-integrals. Such empty integrals will
 * be skipped by the solver. They are useful for the Scheduler when a fixed number of
 * integrals is expected but all relevant integrals are already calculated.
 *
 *
 */
struct Integral
{
    vector<int> indexes;            //!< Wannier function labels and relative shift vectors (size=13)
    double value;                   //!< value of the integral (will be set by the Solver)
    IntegralStatus currentStatus;   //!< determines if the evaluation of the integral has failed
    string error_msg;               //!< error message in terms of failure

    Integral(int c1,int c2,int v1,int v2,int  RDx,int RDy,int RDz,int  Rcx,int Rcy,int Rcz,int  Rvx,int Rvy,int Rvz) {
        this->indexes = vector<int>{c1,c2,v1,v2, RDx,RDy,RDz, Rcx,Rcy,Rcz, Rvx,Rvy,Rvz};
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }
    Integral(const Density_descr& cDens, const Density_descr& vDens, const vector<int>& RD) {
        assert(RD.size() == 3);
        this->indexes = vector<int>{cDens.id1,cDens.id2,vDens.id1,vDens.id2, RD[0],RD[1],RD[2], cDens.R[0],cDens.R[1],cDens.R[2], vDens.R[0],vDens.R[1],vDens.R[2]};
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }
    explicit Integral(vector<int> &ids) {
        if (ids.size() != 13) throw runtime_error("Integral needs index of length 13");
        this->indexes = ids;
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }
    Integral(int c1,int c2,int v1,int v2, const vector<int>& RD, const vector<int>& Rc, const vector<int>& Rv) {
        assert(RD.size() == 3);
        assert(Rc.size() == 3);
        assert(Rv.size() == 3);
        this->indexes = vector<int>{c1,c2,v1,v2, RD[0],RD[1],RD[2], Rc[0],Rc[1],Rc[2], Rv[0],Rv[1],Rv[2]};
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }
    Integral() {
        this->indexes = vector<int>(13);
        this->indexes[0] = 0; // make sure it is not a NO-integral
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }

    bool isEmpty() const {
        return this->indexes[0]<0;  // No integral;
    }
    bool isFailed() const {return currentStatus==IntegralStatus::failed;}

    vector<int> getRD() const {return vector<int>{this->indexes[4],this->indexes[5],this->indexes[6]};}
    vector<int> getRc() const {return vector<int>{this->indexes[7],this->indexes[8],this->indexes[9]};}
    vector<int> getRv() const {return vector<int>{this->indexes[10],this->indexes[11],this->indexes[12]};}

    /**
     * @brief Returns the labels of the valence density
     *
     * ONLY FOR COULOMB INTEGRALS
     */
    Density_descr getValDensity() const {
        return Density_descr(this->indexes[2], this->indexes[3], this->getRv());
    }

    /**
     * @brief Returns the labels of the conduction density
     *
     * ONLY FOR COULOMB INTEGRALS
     */
    Density_descr getConDensity() const {
        return Density_descr(this->indexes[0], this->indexes[1], this->getRc());
    }

    /**
     * @brief Create an empty-integral (placeholder)
     *
     * Empty integrals (c1<0) are not evaluated by the Solver. They are only
     * placeholders for padding.
     *
     * @return Integral
     */
    static Integral CreateEmptyIntegral() {
        Integral i = Integral();
        i.setEmpty();
        return i;
    }

    // Integral getHermitianConj() {
    //     // TODO
    // }

    /**
     * @brief Convert an integral to an empty-integral (placeholder)
     *
     * empty integrals as specified by c1<0 (indexes[0]<0). Additionally all other indexes and the value is reset.
     *
     */
    void setEmpty() {  // no indexes and values, this is just a placeholder object
        this->indexes = vector<int>(13);
        this->indexes[0] = INT_MIN;  // No integral
        this->value = NAN;
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }

    /**
     * @brief Changes state to failed.
     *
     * @param error_msg
     */
    void setFailed(string const& error_msg="") {
        this->value = NAN;
        this->currentStatus = IntegralStatus::failed;
        this->error_msg = error_msg;
    }

    /**
     * @brief Changes state to normal (is default).
     *
     */
    void setNormal() {
        this->currentStatus = IntegralStatus::normal;
        this->error_msg = "";
    }

    void print() const {
        cout << this->toString() << endl;
    }

    const string toString() const {

        if (this->isEmpty()) {
            return "<< Empty Integral >>";
        }

        ostringstream out;
        out << fixed;
        out << setprecision(12);

        for (size_t j=0; j<indexes.size(); j++) {
            out << to_string(indexes[j]) << "\t";
        }
        if (currentStatus==IntegralStatus::normal) {
            if (isnan(value)) {
                out << "(planned)";
            }else{
                out << value;
            }
        }else if (currentStatus==IntegralStatus::failed) {
            out << "FAILED! Error Message: " << error_msg;
        }
        return std::move(out).str();
    }

    bool operator==(const Integral& other) const {
        return this->indexes == other.indexes && this->value == other.value;
    }

    bool operator!=(const Integral& other) const {
        return !(*this == other);
    }
};


#endif // COULOMBINTEGRAL_H