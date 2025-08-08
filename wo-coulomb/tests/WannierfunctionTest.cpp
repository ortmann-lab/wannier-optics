// tests.cpp
#include <wannierfunction.h>
#include <generator.cpp>
#include <vector>
#include <gtest/gtest.h>

TEST(WannierfunctionTest, LatticeSupercell) {

    double L_unitcell = 4.;
    double L = 3*L_unitcell;

    vector<int> dim{100,100,100};
    vector<int> R1{1,0,0};
    vector<double> origin{-L_unitcell,-L_unitcell,-L_unitcell};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{L_unitcell,0,0};
    unitcell[1] = vector<double>{0,L_unitcell,0};
    unitcell[2] = vector<double>{0,0,L_unitcell};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L_unitcell};
    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{generateCube(0, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    wann0.setUnitcell(unitcell);

    vector<double> shouldBe{3,3,1};
    vector<double> is = wann0.getLatticeInUnitcellBasis();

    for (int i=0; i<3; i++) {
        EXPECT_NEAR(shouldBe[i], is[i], 1e-10);
    }

}


TEST(WannierfunctionTest, LatticeSupercell2) {

    double L_unitcell = 4.;

    vector<int> dim{100,61,8};
    vector<double> origin{-16.577331, -28.712780, -0.185185};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{16.485744,0,0};
    unitcell[1] = vector<double>{-8.242872,14.277073,0};
    unitcell[2] = vector<double>{0,0,20.0};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{82.42872, 0, 0};
    lattice[1] = vector<double>{-24.728616,42.831219,0};
    lattice[2] = vector<double>{0,0,30.0};
    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    wann0.setUnitcell(unitcell);

    vector<double> shouldBe{5,3,1.5};
    vector<double> is = wann0.getLatticeInUnitcellBasis();

    for (int i=0; i<3; i++) {
        EXPECT_NEAR(shouldBe[i], is[i], 1e-10);
    }

}




TEST(WannierfunctionTest, CreateSupercell) {

    vector<int> dim{100,61,82};
    vector<double> origin{-16.577331, -28.712780, 0.0};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{16.485744,0,0};
    unitcell[1] = vector<double>{-8.242872,14.277073,0};
    unitcell[2] = vector<double>{0,0,20.0};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{49.457232, 0, 0};
    lattice[1] = vector<double>{-24.728616,42.831219,0};
    lattice[2] = vector<double>{0,0,20.0};
    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);
    shared_ptr<RealMeshgrid> meshgrid2 = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{generateCube(0, meshgrid1, lattice, 5.0, -10., -10, 10., 1.0)};
    wann0.setUnitcell(unitcell);

    WannierFunction wann_supercell{generateCube(0, meshgrid2, lattice, 5.0, -10., -10, 10., 1.0)};
    wann_supercell.setUnitcell(unitcell);

    ASSERT_NE(wann0.getValue(), wann_supercell.getValue());  // check that its not the identical object
    createLargerSupercell(wann_supercell, vector<int>{2, 2, 2});

    // XSF_controller::save_file("wann0.xsf",wann0);
    // XSF_controller::save_file("wann_supercell.xsf",wann_supercell);

    vector<double> wann0Latt = wann0.getLatticeInUnitcellBasis();
    vector<double> wann_supercellLatt = wann_supercell.getLatticeInUnitcellBasis();

    EXPECT_NEAR(3.0, wann0Latt[0], 1e-10);
    EXPECT_NEAR(3.0, wann0Latt[1], 1e-10);
    EXPECT_NEAR(1.0, wann0Latt[2], 1e-10);
    EXPECT_NEAR(wann_supercellLatt[0],6.0, 1e-10);
    EXPECT_NEAR(wann_supercellLatt[1],6.0, 1e-10);
    EXPECT_NEAR(wann_supercellLatt[2],2.0, 1e-10);

    // test meshgrid and unitcell
    for (int i=0; i<3; i++) {
        ASSERT_EQ(wann0.getMeshgrid()->getDim()[i]*2, wann_supercell.getMeshgrid()->getDim()[i]);
        ASSERT_EQ(wann0.getMeshgrid()->getOrigin()[i], wann_supercell.getMeshgrid()->getOrigin()[i]);
        for (int j=0; j<3; j++) {
            ASSERT_EQ(wann0.getUnitcell()[i][j], wann_supercell.getUnitcell()[i][j]);
        }
    }
    EXPECT_NEAR(wann0.getMeshgrid()->getdV(), wann_supercell.getMeshgrid()->getdV(), 1e-10);
    EXPECT_NEAR(wann0.getMeshgrid()->getVgrid()*8, wann_supercell.getMeshgrid()->getVgrid(), 1e-10);

    //wann_supercell->getMeshgrid()->setupCoordinates();
    double x,y,z, x2,y2,z2;
    for (int k=0; k<wann_supercell.getMeshgrid()->getDim()[2]; k++) {  // z
        for (int j=0; j<wann_supercell.getMeshgrid()->getDim()[1]; j++) {  // y
            for (int i=0; i<wann_supercell.getMeshgrid()->getDim()[0]; i++) {  // x

                if ((i<wann0.getMeshgrid()->getDim()[0]) && (j<wann0.getMeshgrid()->getDim()[1]) && (k<wann0.getMeshgrid()->getDim()[2])) {
                    wann0.getMeshgrid()->xyz(i + dim[0]*( j + dim[1]*k ),x,y,z);
                    wann_supercell.getMeshgrid()->xyz(i + wann_supercell.getMeshgrid()->getDim()[0]*( j + wann_supercell.getMeshgrid()->getDim()[1]*k ), x2,y2,z2);
                    ASSERT_NEAR(x, x2, 1e-10);
                    ASSERT_NEAR(y, y2, 1e-10);
                    ASSERT_NEAR(z, z2, 1e-10);
                }
            }
        }
    }


    // test data
    EXPECT_NEAR(wann0.getNorm(), wann_supercell.getNorm(), 1e-10);
    // TODO more!

}


// TEST(WannierfunctionTest, Dipol1) {

//     double L = 5.;

//     vector<int> dim{101,101,101};
//     vector<double> origin{0,0,0};

//     vector< vector<double> > lattice(3);
//     lattice[0] = vector<double>{L,0,0};
//     lattice[1] = vector<double>{0,L,0};
//     lattice[2] = vector<double>{0,0,L};
//     RealMeshgrid meshgrid(dim, lattice, origin);

//     WannierFunction* wann = generateSphere(0, &meshgrid, lattice, L/2, L/2, L/2, L/2, -1.3);

//     Monopole mono = wann->getMonopole();
//     vector<double> dipol = wann->getDipol();

//     EXPECT_NEAR(0.0, dipol[0], 1e-5);
//     EXPECT_NEAR(0.0, dipol[1], 1e-5);
//     EXPECT_NEAR(0.0, dipol[2], 1e-5);
// }


TEST(WannierfunctionTest, getShiftedValue) {

    vector<int> dim{25,30, 21};
    vector<double> origin{0, 0, 0.0};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{2,0,0};
    unitcell[1] = vector<double>{-1,3,0};
    unitcell[2] = vector<double>{0,0,20.0};

    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, unitcell, origin);
    unique_ptr<double[], free_deleter> values{ (double*) malloc(sizeof(double)*meshgrid1->getNumDataPoints()) };
    unique_ptr<double[], free_deleter> values2{ (double*) malloc(sizeof(double)*meshgrid1->getNumDataPoints()) };
    for (int i=0; i<meshgrid1->getNumDataPoints(); i++) {
        values[i] = 1.0;
        values2[i] = 1.0;
    }
    WannierFunction wann0(0,meshgrid1, std::move(values), unitcell);
    WannierFunction wann_supercell(0,meshgrid1,std::move(values2),unitcell);
    ASSERT_NE(wann0.getValue(), wann_supercell.getValue());  // check that its not the identical object
    createLargerSupercell(wann_supercell, vector<int>{2, 2, 2});

    // XSF_controller::save_file("test.xsf", wann_supercell);

    vector<double> wann0Latt = wann0.getLatticeInUnitcellBasis();
    // vector<double> wann_supercellLatt = wann_supercell->getLatticeInUnitcellBasis();

    EXPECT_NEAR(1.0, wann0Latt[0], 1e-10);
    EXPECT_NEAR(1.0, wann0Latt[1], 1e-10);
    EXPECT_NEAR(1.0, wann0Latt[2], 1e-10);

    RealMeshgrid* mesh_sc = wann_supercell.getMeshgrid();

    for (int i=-1; i<2; i++) {
        for (int j=-1; j<2; j++) {
            for (int k=0; k<2; k++) {
                for (int x=0; x <dim[0]; x++){
                    for (int y=0; y <dim[1]; y++){
                        for (int z=0; z<dim[2]; z++){
                            if (i==0 && j==0 && k==0) {
                                ASSERT_NEAR(1.0, wann_supercell.getShiftedValue(mesh_sc->getGlobId(x,y,z),vector<int>{i,j,k}), 1e-10);
                            } else
                                ASSERT_NEAR(0.0, wann_supercell.getShiftedValue(mesh_sc->getGlobId(x,y,z),vector<int>{i,j,k}), 1e-10);
                        }
                    }
                }
            }
        }
    }
}

// double inline gaussian3D(double x, double y, double z, double x0, double y0, double z0, double alpha, double norm=1.) {
//     return norm * pow(alpha / M_PI, 3./2.) * exp(-alpha * (  pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2) ) );
// }

TEST(WannierfunctionTest, getShiftedValue2) {

    vector<int> dim{20,20,20};
    vector<double> origin{0, 0, 0};


    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{10,0,0};
    lattice[1] = vector<double>{0,10,0};
    lattice[2] = vector<double>{0,0,10};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{2,0,0};
    unitcell[1] = vector<double>{0,2,0};
    unitcell[2] = vector<double>{0,0,2};

    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);
    double x,y,z;
    unique_ptr<double[], free_deleter> values{ (double*) malloc(sizeof(double)*meshgrid1->getNumDataPoints()) };
    for (int i=0; i<meshgrid1->getNumDataPoints(); i++) {
        meshgrid1->xyz(i,x,y,z);
        // if ((x<2) && (y<2) && (z<2))
        //     values[i] = 1.;//gaussian3D(x,y,z,1,1,1,1);
        // else
        //     values[i] = 0.0;
        values[i] = cube(x,y,z,1,1,1,1);
    }

    unique_ptr<double[], free_deleter> values_shifted{ (double*) malloc(sizeof(double)*meshgrid1->getNumDataPoints()) };
    for (int i=0; i<meshgrid1->getNumDataPoints(); i++) {
        meshgrid1->xyz(i,x,y,z);

        // if ((x<4) && (x>=2) && (y<6) && (y>=4) && (z<2))
        //     values_shifted[i] = 1.;//gaussian3D(x,y,z,3,5,1,1);
        // else
        //     values_shifted[i] = 0.0;
        values_shifted[i] = cube(x,y,z,1+2,1+4,1,1);
    }


    WannierFunction wann(0,meshgrid1,std::move(values),unitcell);
    WannierFunction wann_shifted(0,meshgrid1,std::move(values_shifted),unitcell);

    ASSERT_NEAR(wann.getNorm(), wann_shifted.getNorm(), 1e-10);


    // double* values_shifted2 = (double*) malloc(sizeof(double)*meshgrid1->getNumDataPoints());
    // for (int i=0; i<meshgrid1->getNumDataPoints(); i++) {
    //     values_shifted2[i] = wann->getShiftedValue(i,vector<int>{1,2,0});  // TODO warum negativ???
    // }
    // WannierFunction* wann_shifted2 = new WannierFunction(0,meshgrid1,values_shifted2,unitcell);


    // XSF_controller::save_file("wann.xsf", wann);
    // XSF_controller::save_file("wann_s.xsf", wann_shifted);
    // XSF_controller::save_file("wann_s2.xsf", wann_shifted2);
    int i,j,k;
    int N = meshgrid1->getNumDataPoints();
    for (int n=0; n<N; n++) {
        // cout << n << endl;
        meshgrid1->getIndexTipple(n,i,j,k);
        // cout << i << " " << j << " " << k << "\n";
        ASSERT_NEAR(wann_shifted.getValue()[n], wann.getShiftedValue(n,vector<int>{1,2,0}), 1e-10);

        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{5,0,0}), 1e-10);
        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{0,5,0}), 1e-10);
        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{0,0,5}), 1e-10);

        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{-1,0,0}), 1e-10);
        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{0,-1,0}), 1e-10);
        ASSERT_NEAR(0.0, wann.getShiftedValue(n,vector<int>{0,0,-1}), 1e-10);
    }

}
