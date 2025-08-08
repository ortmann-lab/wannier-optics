#ifndef OPTICALDIPOLE_H
#define OPTICALDIPOLE_H

/**
 * @file opticalDipole.h
 * @author Konrad Merkel
 * @brief 
 */

#include <vector>
#include <sstream>
#include <iomanip>


enum class Status { planned, calculated, failed };

/**
 * @brief Data structure to store and calculate optical transition matrix elements between conduction and valence WF
 * 
 * In contrast to Coulomb integrals we do not separate the implementation and storage, since the calculation of
 * transition dipoles is much easier and can be done easily in real space.
 * 
 */
struct OpticalDipole
{
    vector<int> indexes;          //!< ids of WF and relative shifts
    vector<double> dipole;        //!< dipole moment in cartesian coordinates
    Status currentStatus;         //!< is values already calculated?


    OpticalDipole(vector<int> indexes_, vector<double> dipole_)
     : indexes(indexes_), dipole(dipole_), currentStatus(Status::calculated)
    {
        assert(indexes_.size() == 5);
        assert(dipole_.size() == 3);
    }

    explicit OpticalDipole(vector<int> indexes_)
     : indexes(indexes_), currentStatus(Status::planned)
    {
        assert(indexes_.size() == 5);
        dipole = vector<double>{0.,0.,0.};
    }

    OpticalDipole(int idc, int idv, vector<int> shift)
    {
        assert(shift.size() == 3);
        currentStatus = Status::planned;
        dipole = vector<double>{0.,0.,0.};
        indexes = vector<int>{idc,idv,shift[0],shift[1],shift[2]};
    }


    vector<int> getShift() const {
        return vector<int>{indexes[2],indexes[3],indexes[4]};
    }
    int getIdc() const { return indexes[0]; }
    int getIdv() const { return indexes[1]; }

    bool isCalculated() const { return currentStatus == Status::calculated; }
    Status getStatus() const {return this->currentStatus;}

    vector<double> getDipole() const { return dipole; }
    void setDipole(vector<double> const& new_dipole) {
        this->dipole=new_dipole;
        this->currentStatus = Status::calculated;
    }


    string toString() const {
        ostringstream out;
        out << fixed;
        out << setprecision(12);

        for (int j=0; j<5; j++) {
            out << to_string(indexes[j]) << "\t";
        }
        
        if (currentStatus == Status::calculated) {
            for (int j=0; j<3; j++) {
                out << this->dipole[j] << "\t";
            }
        }
        else if (currentStatus == Status::failed) {
            out << "FAILED! ";
        } else if (currentStatus == Status::planned) { 
            out << "PLANNED";
        }
        
        return std::move(out).str();
    }

    bool operator<(const OpticalDipole& b) const {  // needed for ordering
        return this->indexes < b.indexes;
    }

};

#endif // OPTICALDIPOLE_H
