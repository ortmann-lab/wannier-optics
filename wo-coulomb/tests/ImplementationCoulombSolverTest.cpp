

#include <coulombIntegral.h>
#include <wannierfunction.h>
#include <potential.h>
#include <solver.h>
#include <generator.cpp>
#include <mpi.h>
#include <gtest/gtest.h>
#include <iostream>
using namespace std;


WannierFunction generateGauss2(int id, shared_ptr<RealMeshgrid> meshgrid, vector< vector<double> > const& lattice, double alpha, double x0, double y0, double z0, double norm) {

    // fill arrays (initialize inputs)
    // column major order (Fortan style !!!)
    unique_ptr<double[], free_deleter> data{ (double*) malloc(sizeof(double)*meshgrid->getNumDataPoints()) };
    double x,y,z;
    for (int i=0; i < meshgrid->getNumDataPoints(); i++){
        meshgrid->xyz(i,x,y,z);
        data[i] = sqrt(gaussian3D(x,y,z, x0, y0, z0, alpha));
    }

    //normalizeWann(wannier, norm);
    return WannierFunction(id, meshgrid, std::move(data), lattice);
}

TEST(ImplementationFourierGaussTest, GaussFunction)  // TODO: maybe remove
/*
    Tests the norm and fourier transform of the gauss function
*/
{

    GTEST_SKIP();

    double L = 10.;
    vector<int> R{0,0,0};
    vector<int> dim{151,151,151};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);
    int N = meshgrid->getNumDataPoints();

    double x0 = L/2-1.6;
    double y0 = x0;
    double z0 = x0;
    double alpha = 1.5;

    // real space gaussian
    double x,y,z;
    double* values = (double*) malloc(sizeof(double)*N);
    for (int i=0; i<N; i++) {
        meshgrid->xyz(i,x,y,z);
        values[i] = gaussian3D(x,y,z, x0, y0, z0, alpha, 1.0);
    }
    // integrate
    double norm = 0;
    for (int i=0; i<N; i++) {
        norm += values[i] * meshgrid->getdV();
    }

    EXPECT_NEAR(norm, 1.0, 1e-6);


    // reciprocal space gaussian
    ReciprocalMeshgrid recMesh(meshgrid.get());
    //recMesh.setupCoordinates();
    // const double* qx = recMesh.getX();
    // const double* qy = recMesh.getY();
    // const double* qz = recMesh.getZ();
    for (int i=0; i<N; i++) {
        values[i] = 1.0; //gaussian3D(qx[i], qy[i], qz[i], 0.0, 0.0, 0.0, 1./(4.*alpha), 1.0);  //pow( M_PI / alpha, 3./2.)
    }
    // integrate
    norm = 0;
    for (int i=0; i<N; i++) {
        norm += values[i] * recMesh.getdV();
    }
    cout << "\n\n\nREC:\n" << recMesh.getdV() << " ?= " << 2*M_PI/meshgrid->getdV() << endl;
    EXPECT_NEAR(norm, 1.0, 1e-6);

}


TEST(ImplementationFourierGaussTest, FourierRealGaussCubicCell)
{
    double L = 15.;
    vector<int> R{0,0,0};
    vector<int> dim{100,100,109};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double a = 1.5;

    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1+origin[0], L/2-a/2+origin[1], L/2-a/2+origin[2], 1.0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, generateCube(1, meshgrid, lattice, a, L/2+a/2+origin[0], L/2+a/2+origin[1], L/2+a/2+origin[2], 1.0)});


    Potential* pot = new CoulombPotential();
    auto WFouier = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
    auto WRealSpace = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    CoulombSolver implFourier(vWannMap, cWannMap);
    RealSpaceSolver implRealSpace(vWannMap, cWannMap, pot);

    implFourier.calculate(WFouier, false,1,10);
    implRealSpace.calculate(WRealSpace, false,1,10);

    cout << WFouier.value << " ?= " << WRealSpace.value << endl;
    EXPECT_NEAR(WFouier.value, WRealSpace.value, 1e-4);

}

TEST(ImplementationFourierGaussTest, FourierRealGaussCubicCellOrigin)
{
    double L = 15.;
    vector<int> R{0,0,0};
    vector<int> dim{100,100,109};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};

    vector<double> origin{0,0,0};
    double a = 1.5;

    Potential* pot = new CoulombPotential();
    map< int,WannierFunction > cWannMap1{};
    map< int,WannierFunction > vWannMap1{};
    Integral WFouier, WRealSpace;

    // zero origin
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);


    cWannMap1.insert({0, generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1+origin[0], L/2-a/2+origin[1], L/2-a/2+origin[2], 1.0)});
    vWannMap1.insert({1, generateCube(1, meshgrid, lattice, a, L/2+a/2+origin[0], L/2+a/2+origin[1], L/2+a/2+origin[2], 1.0)});

    WFouier = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
    WRealSpace = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    CoulombSolver implFourier1(vWannMap1, cWannMap1);
    RealSpaceSolver implRealSpace1(vWannMap1, cWannMap1, pot);

    implFourier1.calculate(WFouier, false,1,10);
    implRealSpace1.calculate(WRealSpace, false,1,10);

    cout << WFouier.value << " ?= " << WRealSpace.value << endl;
    EXPECT_NEAR(WFouier.value, WRealSpace.value, 1e-4);


    for (double i=-3.32; i<9.; i+=5.1) {
        for (double j=-2.7; j<10.5; j+=6.2) {
            for (double k=-1.9; k<9.; k+=7.5) {
                origin[0] = i;
                origin[1] = j;
                origin[2] = k;

                cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
                meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

                // create mappings
                map< int,WannierFunction > cWannMap{};
                map< int,WannierFunction > vWannMap{};
                cWannMap.insert({0, generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1+origin[0], L/2-a/2+origin[1], L/2-a/2+origin[2], 1.0)});
                vWannMap.insert({1, generateCube(1, meshgrid, lattice, a, L/2+a/2+origin[0], L/2+a/2+origin[1], L/2+a/2+origin[2], 1.0)});

                // set mappings in solver

                CoulombSolver implFourier(vWannMap, cWannMap);
                RealSpaceSolver implRealSpace(vWannMap, cWannMap, pot);


                // define Integrals
                WFouier = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
                WRealSpace = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

                // solve integrals
                implFourier.calculate(WFouier, false,1,10);
                implRealSpace.calculate(WRealSpace, false,1,10);


                cout << WFouier.value << " ?= " << WRealSpace.value << endl;
                ASSERT_NEAR(WFouier.value, WRealSpace.value, 1e-4);

            }
        }
    }
}


TEST(ImplementationFourierGaussTest, FourierRealGaussCubicCell2) {

    //GTEST_SKIP(); // takes too long


    double L = 15.;
    vector<int> R{0,1,0};
    vector<int> dim{101,101,101};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);
    Potential* pot = new CoulombPotential();
    double a = 3;

    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0, generateCube(0, meshgrid, lattice, a, L/2, L/2, L/2, 1.0)});
    vWannMap.insert({1, generateCube(0, meshgrid, lattice, a, L/2, L/2, L/2, 1.0)});

    WannierFunction wann0_sc{ generateCube(0, meshgrid, lattice, a, L/2, L/2, L/2, 1.0)};
    createLargerSupercell(wann0_sc, vector<int>{2,2,2});
    WannierFunction wann0_sc2{ generateCube(0, meshgrid, lattice, a, L/2, L/2, L/2, 1.0)};
    createLargerSupercell(wann0_sc2, vector<int>{2,2,2});

    map< int,WannierFunction > cWannMap_sc{};
    map< int,WannierFunction > vWannMap_sc{};
    cWannMap_sc.insert({0, std::move(wann0_sc)});
    vWannMap_sc.insert({1, std::move(wann0_sc2)});

    CoulombSolver implFourier(vWannMap_sc, cWannMap_sc);
    RealSpaceSolver implRealSpace(vWannMap, cWannMap, pot);

    auto WFouier = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
    auto WRealSpace = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    implFourier.calculate(WFouier, false,1,10);
    implRealSpace.calculate(WRealSpace, false,1,10);

    cout << WFouier.value << " ?= " << WRealSpace.value << endl;
    EXPECT_NEAR(WFouier.value, WRealSpace.value, 1e-3);
}

TEST(ImplementationFourierGaussTest, FourierRealGaussCubicCellFiniteR)
{
    double L = 15.;
    vector<int> dim{101,101,101};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double a = 1.4;
    vector<int> R{1,0,0};

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1, L/2-a/2-0.1, L/2-a/2-0.1, 1.0)};
    WannierFunction wann1{ generateCube(1, meshgrid, lattice, a, L/2+a/2, L/2+a/2, L/2+a/2, 1.0)};

    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0, std::move(wann0)});
    vWannMap.insert({1, std::move(wann1)});


    WannierFunction wann0_sc{ generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1, L/2-a/2-0.1, L/2-a/2-0.1, 1.0)};
    WannierFunction wann1_sc{ generateCube(1, meshgrid, lattice, a, L/2+a/2, L/2+a/2, L/2+a/2, 1.0)};
    createLargerSupercell(wann0_sc,vector<int>{2,2,2});
    createLargerSupercell(wann1_sc,vector<int>{2,2,2});
    // wann0_sc->getMeshgrid()->setupCoordinates();
    // wann1_sc->getMeshgrid()->setupCoordinates();

    map< int,WannierFunction > cWannMap_sc{};
    map< int,WannierFunction > vWannMap_sc{};
    cWannMap_sc.insert({0, std::move(wann0_sc)});
    vWannMap_sc.insert({1, std::move(wann1_sc)});


    Potential* pot = new CoulombPotential();
    auto WFouier = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
    auto WRealSpace = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    CoulombSolver implFourier(vWannMap_sc, cWannMap_sc);
    RealSpaceSolver implRealSpace(vWannMap, cWannMap, pot);

    implFourier.calculate(WFouier, false,1,10);
    implRealSpace.calculate(WRealSpace, false,1,10);

    cout << WFouier.value << " ?= " << WRealSpace.value << endl;
    EXPECT_NEAR(WFouier.value, WRealSpace.value, 1e-4);
}


TEST(ImplementationFourierGaussTest, CubicAnalytic)
{
    double L = 25.;
    vector<int> R{0,0,0};
    vector<int> dim{40,40,40};
    vector<double> origin{100,-52.,7};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double alpha = 0.3;
    double gamma = 0.4;

    double x0 = L/2 + origin[0];
    double y0 = L/2 + origin[1];
    double z0 = L/2 + origin[2];

    double deltax = 1.0;
    double deltay = -2.11;
    double deltaz = 0.7;

    double r = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
    double shouldbe = ARB_TO_EV * erf(sqrt(alpha*gamma/(alpha+gamma)) * r) / r;

    WannierFunction wann0_sc{ generateGauss2(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0) };
    WannierFunction wann1_sc{ generateGauss2(1, meshgrid, lattice, gamma, x0+deltax,y0+deltay, z0+deltaz, 1.0) };

    EXPECT_NEAR(wann0_sc.getNorm(), 1.0, 1e-10);
    EXPECT_NEAR(wann1_sc.getNorm(), 1.0, 1e-10);

    // supercell
    createLargerSupercell(wann0_sc,vector<int>{3,3,3});
    createLargerSupercell(wann1_sc,vector<int>{3,3,3});

    EXPECT_NEAR(wann0_sc.getNorm(), 1.0, 1e-10);
    EXPECT_NEAR(wann1_sc.getNorm(), 1.0, 1e-10);


    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, std::move(wann0_sc)});
    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, std::move(wann1_sc)});

    auto WFouier_sc = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    CoulombSolver implFourier(vWannMap, cWannMap);

    implFourier.calculate(WFouier_sc, false,1,10);

    cout << WFouier_sc.value << " ?= " << shouldbe << endl;

    EXPECT_NEAR(WFouier_sc.value, shouldbe, 1e-4);

}


TEST(ImplementationFourierGaussTest, HexagonalAnalytic)
{
    double L = 25.;
    vector<int> R{0,0,0};
    vector<int> dim{40,40,40};
    vector<double> origin{-2.3,-5.2,3.3333};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{-L/2.,L*sqrt(3.)/2.,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double alpha = 0.2;
    double gamma = 0.3;

    double x0 = L/4. + origin[0];
    double y0 = L*sqrt(3.)/4. + origin[1];
    double z0 = L/2. + origin[2];

    double deltax = 1.0;
    double deltay = -2.11;
    double deltaz = 0.7;

    double r = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
    double shouldbe = ARB_TO_EV * erf(sqrt(alpha*gamma/(alpha+gamma)) * r) / r;

    WannierFunction wann0{ generateGauss2(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0)};
    WannierFunction wann1{ generateGauss2(1, meshgrid, lattice, gamma, x0+deltax,y0+deltay, z0+deltaz, 1.0)};

    // XSF_controller::save_file("test.xsf", wann0);

    EXPECT_NEAR(wann0.getNorm(), 1.0, 1e-10);
    EXPECT_NEAR(wann1.getNorm(), 1.0, 1e-10);

    // CoulombIntegral* WFouier = new CoulombIntegral(implFourier, wann0, wann0, wann1, wann1, R, pot);
    // cout << WFouier->getValue().real() << " ?= " << shouldbe << endl;

    // ASSERT_LT(abs(WFouier->getValue().imag()), 1e-10);
    // EXPECT_NEAR(WFouier->getValue().real(), shouldbe, 1e-3);


    // supercell
    createLargerSupercell(wann0,vector<int>{2,2,2}); // set to 5
    createLargerSupercell(wann1,vector<int>{2,2,2});

    EXPECT_NEAR(wann0.getNorm(), 1.0, 1e-10);
    EXPECT_NEAR(wann1.getNorm(), 1.0, 1e-10);


    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, std::move(wann0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, std::move(wann1)});

    auto WFouier_sc = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

    CoulombSolver implFourier(vWannMap, cWannMap);
    implFourier.calculate(WFouier_sc, false,1,10);

    cout << WFouier_sc.value << " ?= " << shouldbe << endl;
    EXPECT_NEAR(WFouier_sc.value, shouldbe, 1e-4);
}




// TEST(ImplementationFourierGaussTest, HexagonalAnalyticSupercellSeries) {
//     /*
//         These tests are for convergence plots. They are not software tests in
//         the closer sense.
//     */

//     GTEST_SKIP();

//     Solver* implFourier = new CoulombSolver();
//     Potential* pot = new CoulombPotential();

//     double L = 25.;
//     vector<int> R{0,0,0};
//     vector<int> dim{60,60,60};
//     vector<double> origin{-2.3,-5.2,3.3333};
//     vector< vector<double> > lattice(3);
//     lattice[0] = vector<double>{L,0,0};
//     lattice[1] = vector<double>{-L/2.,L*sqrt(3.)/2.,0};
//     lattice[2] = vector<double>{0,0,L};
//     shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

//     double alpha = 0.2;
//     double gamma = 0.3;

//     double x0 = L/4. + origin[0];
//     double y0 = L*sqrt(3.)/4. + origin[1];
//     double z0 = L/2. + origin[2];

//     double deltax = 1.0;
//     double deltay = -2.11;
//     double deltaz = 0.7;

//     double r = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
//     double shouldbe = ARB_TO_EV * erf(sqrt(alpha*gamma/(alpha+gamma)) * r) / r;
//     cout << "Analytic result (eV): " << shouldbe << endl;
//     cout << "i\trelative error(%)\n";

//     WannierFunction* wann0 = generateGauss2(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0);
//     WannierFunction* wann1 = generateGauss2(1, meshgrid, lattice, gamma, x0+deltax,y0+deltay, z0+deltaz, 1.0);

//     EXPECT_NEAR(wann0->getNorm(), 1.0, 1e-10);
//     EXPECT_NEAR(wann1->getNorm(), 1.0, 1e-10);


//     CoulombIntegral* WFouier = new CoulombIntegral(wann0, wann0, wann1, wann1, R, pot);
//     implFourier->calculate(WFouier, false,1,10);
//     cout << 1 << "\t" << (WFouier->getValue().real()-shouldbe)/shouldbe*100 << endl;

//     ASSERT_LT(abs(WFouier->getValue().imag()), 1e-10);
//     // EXPECT_NEAR(WFouier->getValue().real(), shouldbe, 1e-5);


//     // supercell
//     WannierFunction* wann0_sc=nullptr;
//     WannierFunction* wann1_sc=nullptr;
//     CoulombIntegral* WFouier_sc=nullptr;
//     for (int i=2; i<10; i++) {
//         wann0_sc = createLargerSupercell(wann0,vector<int>{i,i,i});
//         wann1_sc = createLargerSupercell(wann1,vector<int>{i,i,i});

//         EXPECT_NEAR(wann0_sc->getNorm(), 1.0, 1e-10);
//         EXPECT_NEAR(wann1_sc->getNorm(), 1.0, 1e-10);

//         WFouier_sc = new CoulombIntegral(wann0_sc, wann0_sc, wann1_sc, wann1_sc, R, pot);
//         implFourier->calculate(WFouier_sc, false,1,10);
//         cout << i << "\t" << (WFouier_sc->getValue().real()-shouldbe)/shouldbe*100 << endl;

//         ASSERT_LT(abs(WFouier_sc->getValue().imag()), 1e-10);
//         // EXPECT_NEAR(WFouier_sc->getValue().real(), shouldbe, 1e-5);

//         delete wann0_sc, wann1_sc, WFouier_sc;
//     }

// }



// TEST(ImplementationFourierGaussTest, HexagonalAnalyticDimSeries) {
    /*
        These tests are for convergence plots. They are not software tests in
        the closer sense.
    */

    //GTEST_SKIP();

    // Solver* implFourier = new CoulombSolver();
    // Potential* pot = new CoulombPotential();

    // double L = 25.;
    // vector<int> R{0,0,0};
    // vector<double> origin{-2.3,-5.2,3.3333};
    // vector< vector<double> > lattice(3);
    // lattice[0] = vector<double>{L,0,0};
    // lattice[1] = vector<double>{-L/2.,L*sqrt(3.)/2.,0};
    // lattice[2] = vector<double>{0,0,L};

    // double alpha = 0.2;
    // double gamma = 0.3;

    // double x0 = L/4. + origin[0];
    // double y0 = L*sqrt(3.)/4. + origin[1];
    // double z0 = L/2. + origin[2];

    // double deltax = 1.0;
    // double deltay = -2.11;
    // double deltaz = 0.7;

    // double r = sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
    // double shouldbe = ARB_TO_EV * erf(sqrt(alpha*gamma/(alpha+gamma)) * r) / r;

    // cout << "Analytic result (eV): " << shouldbe << endl;
    // cout << "i\trelative error(%)\n";


    // // supercell
    // vector<int> dim(3);
    // RealMeshgrid* meshgrid = nullptr;
    // WannierFunction* wann0 = nullptr;
    // WannierFunction* wann1 = nullptr;
    // WannierFunction* wann0_sc=nullptr;
    // WannierFunction* wann1_sc=nullptr;
    // CoulombIntegral* WFouier_sc=nullptr;
    // for (int i=15; i<100; i+=5) {

    //     dim[0] = i;
    //     dim[1] = i;
    //     dim[2] = i;

    //     meshgrid = new RealMeshgrid(dim, lattice, origin);
    //     wann0 = generateGauss2(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0);
    //     wann1 = generateGauss2(1, meshgrid, lattice, gamma, x0+deltax,y0+deltay, z0+deltaz, 1.0);
    //     // EXPECT_NEAR(wann0->getNorm(), 1.0, 1e-8);
    //     // EXPECT_NEAR(wann1->getNorm(), 1.0, 1e-8);

    //     // create supercell
    //     wann0_sc = createLargerSupercell(wann0,vector<int>{2,2,2});
    //     wann1_sc = createLargerSupercell(wann1,vector<int>{2,2,2});
    //     // EXPECT_NEAR(wann0_sc->getNorm(), 1.0, 1e-8);
    //     // EXPECT_NEAR(wann1_sc->getNorm(), 1.0, 1e-8);

    //     WFouier_sc = new CoulombIntegral(wann0_sc, wann0_sc, wann1_sc, wann1_sc, R, pot);
    //     implFourier->calculate(WFouier_sc, false,1,10);

    //     cout << i << "\t" << (WFouier_sc->getValue().real()-shouldbe)/shouldbe*100 << endl;

    //     // ASSERT_LT(abs(WFouier_sc->getValue().imag()), 1e-10);
    //     // EXPECT_NEAR(WFouier_sc->getValue().real(), shouldbe, 1e-5);

    //     delete wann0, wann1, meshgrid, wann0_sc, wann1_sc, WFouier_sc;
    // }

// }