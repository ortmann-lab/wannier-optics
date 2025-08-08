

#include <wannierfunction.h>
#include <opticalDipole.h>
#include <generator.cpp>
#include <gtest/gtest.h>
#include <iostream>
using namespace std;


// TEST(OpticalDipoleTest, calcDipole) {  // TODO this test can be improved. Now its basicly shown that we get a number
//     /**
//      * TODO
//      **/

//     double L_unitcell = 3.1;
//     double L = 5;

//     vector< vector<double> > unitcell(3);
//     unitcell[0] = vector<double>{L_unitcell,0,0};
//     unitcell[1] = vector<double>{0,L_unitcell,0};
//     unitcell[2] = vector<double>{0,0,L_unitcell};

//     vector<int> dim{50,51,41};
//     vector<int> R{1,0,0};
//     vector<double> origin{0,0,0.};
//     vector< vector<double> > lattice(3);
//     lattice[0] = vector<double>{L,0,0};
//     lattice[1] = vector<double>{0,L,0};
//     lattice[2] = vector<double>{0,0,L_unitcell};
//     RealMeshgrid meshgrid(dim, lattice, origin);

//     WannierFunction* wann0 = generateCube(0, &meshgrid, lattice, 1.8, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0);
//     WannierFunction* wann1 = generateCube(0, &meshgrid, lattice, 1.8, L_unitcell/2+L_unitcell/4., L_unitcell/2, L_unitcell/2, 1.0);
//     //WannierFunction* wann1 = generateSphere(1, &meshgrid, lattice, 1.0, 4.0, L_unitcell/2, L_unitcell/2, 1.0);
//     // WannierFunction* wann1_shifted_ = generateSphere(1, &meshgrid, lattice, 0.95, 0.95+L_unitcell, L_unitcell/2, L_unitcell/2, 1.0);

//     wann0->setUnitcell(unitcell);
//     //wann1->setUnitcell(unitcell);
//     // wann1_shifted_->setUnitcell(unitcell);


//     OpticalDipole* M = new OpticalDipole(wann0, wann1, vector<int>{0,0,0});
//     vector<double> dipole = M->getDipole();

//     for (int j=0; j<3; j++) {
//         cout << dipole[j] << "\t";
//     }
// }