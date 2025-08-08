#include <gtest/gtest.h>
#include <exception>
#include "filehandler.h"
#include <generator.cpp>


TEST(XSF_controllerTest, FileNotFound) {
    ASSERT_THROW(XSF_controller::read_file("NON.xsf"), std::runtime_error);
}

TEST(XSF_controllerTest, WriteFileNotPossible)
{
    vector<int> R{0,0,0};
    vector<int> dim{31,31,31};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,0,0};
    lattice[1] = vector<double>{0,2.,0};
    lattice[2] = vector<double>{0,0,2.};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);
    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, 1, 1, 1, 1.0)};

    ASSERT_FALSE(XSF_controller::save_file("",wann0));
    ASSERT_FALSE(XSF_controller::save_file(".",wann0));
}

TEST(XSF_controllerTest, ReadWriteFile)
{
    vector<int> R{0,0,0};
    vector<int> dim{31,31,31};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,0,0};
    lattice[1] = vector<double>{0,2.,0};
    lattice[2] = vector<double>{0,0,2.};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, 1, 1, 1, 1.0)};
    ASSERT_TRUE(XSF_controller::save_file("testfile.xsf",wann0));

    WannierFunction wann1{ XSF_controller::read_file("testfile.xsf")};
    ASSERT_NE(&wann0, &wann1);  // not the very same object in memory
    ASSERT_NE(wann0.getMeshgrid(), wann1.getMeshgrid()); // not the very same object in memory
    ASSERT_TRUE(wann0.isCompatible(wann1));

    int N = wann0.getMeshgrid()->getNumDataPoints();
    double const* data0 = wann0.getValue();
    double const* data1 = wann1.getValue();
    ASSERT_NE(data0, data1); // not the very same object in memory
    for (int i=0; i<N; i++) {
        ASSERT_NEAR(data0[i], data1[i], 1e-5);  // only 6 decimal places are written
    }
}

TEST(XSF_controllerTest, ReadWriteFile2)
{
    vector<int> R{0,0,0};
    vector<int> dim{31,36,47};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{4.,1,0};
    lattice[1] = vector<double>{0,2.5,0};
    lattice[2] = vector<double>{0.5,0,3.};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);


    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 0.8, 1.2, 1.5, 1, 1.0)};

    ASSERT_TRUE(XSF_controller::save_file("testfile.xsf",wann0));

    WannierFunction wann1{ XSF_controller::read_file("testfile.xsf")};
    ASSERT_NE(&wann0, &wann1);  // not the very same object in memory
    ASSERT_NE(wann0.getMeshgrid(), wann1.getMeshgrid()); // not the very same object in memory
    ASSERT_TRUE(wann0.isCompatible(wann1));

    int N = wann0.getMeshgrid()->getNumDataPoints();
    double const* data0 = wann0.getValue();
    double const* data1 = wann1.getValue();
    ASSERT_NE(data0, data1); // not the very same object in memory
    for (int i=0; i<N; i++) {
        ASSERT_NEAR(data0[i], data1[i], 1e-5);  // only 6 decimal places are written
    }
}