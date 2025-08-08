// tests.cpp
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <gtest/gtest.h>

// all tests:
#include "MeshgridTest.cpp"
#include "SchedulerTest.cpp"
#include "WannierfunctionTest.cpp"
#include "FourierShiftTest.cpp"
#include "XSF_controllerTest.cpp"
#include "ImplementationCoulombSolverTest.cpp"
#include "ImplementationTest.cpp"
#include "LocalFieldEffectsTest.cpp"
#include "OpticalDipoleTest.cpp"
#include "DensityTest.cpp"
#include "CHGCARTest.cpp"


using namespace std;


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    cout << fixed;
    cout << setprecision(12);
    testing::InitGoogleTest(&argc, argv);
    auto ret = RUN_ALL_TESTS();

    MPI_Finalize();
    return ret;
}