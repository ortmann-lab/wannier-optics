#include <gtest/gtest.h>
#include <CHGCAR.h>
#include <filehandler.h>

TEST(CHGCARTest, NumElectrons) {

    vector<int> dim{100,100,100};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{4.,0,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0,0,4.0};
    shared_ptr<RealMeshgrid> mesh = make_shared<RealMeshgrid>(dim, lattice, origin);

    int N = mesh->getNumDataPoints();
    double NELECT = 4;
    double x,y,z;
    int N_inside = 0;

    unique_ptr<double[],free_deleter> density{ (double*)malloc(sizeof(double)* N) };
    for (int i=0; i<N; i++) {
        density[i] = 0.0;

        mesh->xyz(i,x,y,z);

        if ((x>=0) && (x<2.0)) {
            if ((y>=0) && (y<3.0)) {
                if ((z>=0) && (z<2.0)) {
                    //cout << x << "\t" << y << "\t" << z << endl;
                    density[i] = NELECT/(2.*3.*2.0);
                    N_inside++;
                }
            }
        }
    }

    CHGCAR c(mesh, std::move(density));

    // XSF_controller::save_file("test1.xsf",mesh,density, lattice);

    // cout << "NELECT: " << NELECT << endl;
    // cout << "NumElectron: " << c.getNumElectrons() << endl;
    // cout << "volume: " << 2.*3.*2.0 << endl;
    // cout << "approx volume: " << N_inside*1.0/N * mesh->getVgrid() << endl;

    ASSERT_NEAR(c.getNumElectrons(), NELECT, 1e-8);

    c.createSupercell(2,3,1);
    // XSF_controller::save_file("test2.xsf",c.getMeshgrid(),c.getValue(), lattice);

    ASSERT_EQ(c.getMeshgrid()->getDim()[0], dim[0]*2);
    ASSERT_EQ(c.getMeshgrid()->getDim()[1], dim[1]*3);
    ASSERT_EQ(c.getMeshgrid()->getDim()[2], dim[2]*1);
    ASSERT_NEAR(c.getNumElectrons(), NELECT*2*3*1, 1e-8);

}

TEST(CHGCARTest, coaseGrain) {

    vector<int> dim{200,200,200};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{4.,0,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0,0,4.0};
    shared_ptr<RealMeshgrid> mesh = make_shared<RealMeshgrid>(dim, lattice, origin);

    int N = mesh->getNumDataPoints();
    double NELECT = 4;
    double x,y,z;
    int N_inside = 0;

    unique_ptr<double[],free_deleter> density{ (double*)malloc(sizeof(double)* N) };
    for (int i=0; i<N; i++) {
        density[i] = 0.0;

        mesh->xyz(i,x,y,z);

        if ((x>=0) && (x<2.0)) {
            if ((y>=0) && (y<3.0)) {
                if ((z>=0) && (z<2.0)) {
                    //cout << x << "\t" << y << "\t" << z << endl;
                    density[i] = NELECT/(2.*3.*2.0);
                    N_inside++;
                }
            }
        }
    }

    CHGCAR c(mesh, std::move(density));

    // XSF_controller::save_file("test1.xsf",mesh,density, lattice);

    // cout << "NELECT: " << NELECT << endl;
    // cout << "NumElectron: " << c.getNumElectrons() << endl;
    // cout << "volume: " << 2.*3.*2.0 << endl;
    // cout << "approx volume: " << N_inside*1.0/N * mesh->getVgrid() << endl;

    ASSERT_NEAR(c.getNumElectrons(), NELECT, 1e-8);

    c.coarseGrain(2,2,2);
    // XSF_controller::save_file("test2.xsf",c.getMeshgrid(),c.getValue(), lattice);

    ASSERT_EQ(c.getMeshgrid()->getDim()[0], dim[0]/2);
    ASSERT_EQ(c.getMeshgrid()->getDim()[1], dim[1]/2);
    ASSERT_EQ(c.getMeshgrid()->getDim()[2], dim[2]/2);
    ASSERT_NEAR(c.getNumElectrons(), NELECT, 1e-8);

}


TEST(CHGCARTest, read_write) {
    GTEST_SKIP();

    CHGCAR chg {read_CHGCAR("/home/kmerkel/wo-test_data/test_chgcar/CHGCAR")};

    const RealMeshgrid* mesh = chg.getMeshgrid();

    ASSERT_EQ(mesh->getDim()[0], 48);
    ASSERT_EQ(mesh->getDim()[1], 48);
    ASSERT_EQ(mesh->getDim()[2], 48);


    // cout << "NELECT: " << chg->getNumElectrons() << endl;

    ASSERT_NEAR(chg.getNumElectrons(), 8, 1e-6);

    XSF_controller::save_file("test3.xsf", chg.getMeshgrid(), chg.getValue());

    chg.createSupercell(2,2,2);
    XSF_controller::save_file("test4.xsf", chg.getMeshgrid(), chg.getValue());
}


TEST(CHGCARTest, expectationValue) {
    GTEST_SKIP();

    CHGCAR chg{read_CHGCAR("/home/kmerkel/wo-test_data/test_chgcar/CHGCAR")};
    WannierFunction wf{ XSF_controller::read_file("/home/kmerkel/Desktop/Si/valence_11x11x11/wannier90_00001.xsf") };

    wf.normalize();

    EXPECT_THROW(calcExpectationDensity(chg, wf), runtime_error);

    chg.makeCompatible(wf);

    // cout << "\n\n\n\n CHGCAR:" << endl;
    // chg->getMeshgrid()->printSummary();

    // cout << "\n\n\n\n WF:" << endl;
    // wf->getMeshgrid()->printSummary();

    double value = calcExpectationDensity(chg, wf);
    // cout << "expectation value: " << value << endl;
    ASSERT_NEAR(value, 0.2773445207917277, 1e-8);
}