

#include <coulombIntegral.h>
#include <wannierfunction.h>
#include <potential.h>
#include <solver.h>
#include <generator.cpp>
#include <gtest/gtest.h>
#include <iostream>
using namespace std;


TEST(ImplementationTest, Consistency1)
{
    Potential* pot = new CoulombPotential();

    double L = 15.;
    vector<int> dim{100,100,109};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double a = 1.5;

    map< int,WannierFunction > cWannMap{};
    cWannMap.insert( {0, generateCube(0, meshgrid, lattice, a, 4, 10, 5, 1.0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert( {1, generateCube(1, meshgrid, lattice, a, 4.4, 10, 5.3, 1.0)});

    lattice[0][1] = 1.0;
    shared_ptr<RealMeshgrid> meshgrid2 = make_shared<RealMeshgrid>(dim, lattice, origin);
    vWannMap.insert({2, generateCube(2, meshgrid2, lattice, a, L/2+a/2+origin[0], L/2+a/2+origin[1], L/2+a/2+origin[2], 1.0)});

    ASSERT_THROW(RealSpaceSolver(vWannMap, cWannMap, pot), runtime_error);

}

TEST(ImplementationTest, Consistency2)
{
    Potential* pot = new CoulombPotential();

    double L = 15.;
    vector<int> dim{100,100,109};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    map< int,WannierFunction > cWannMap{};
    //cWannMap->insert({0, generateCube(0, meshgrid, lattice, a, 4, 10, 5, 1.0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, generateCube(1, meshgrid, lattice, 1.0, 4.4, 10, 5.3, 1.0)});

    lattice[2][0] = 1.0;
    cWannMap.insert({2, generateCube(2, meshgrid, lattice, 1.0, 4.4, 10, 5.3, 1.0)});

    ASSERT_THROW(RealSpaceSolver(vWannMap, cWannMap, pot), runtime_error);
}

TEST(ImplementationTest, Consistency3)
{
    Potential* pot = new CoulombPotential();

    double L = 15.;
    vector<int> dim{100,100,109};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,1.,0};
    lattice[1] = vector<double>{0,L,0.2};
    lattice[2] = vector<double>{0.4,1.,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    dim[1] = 80;
    shared_ptr<RealMeshgrid> meshgrid2 = make_shared<RealMeshgrid>(dim, lattice, origin);

    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, generateCube(0, meshgrid, lattice, 1.0, 4, 10, 5, 1.0)});
    cWannMap.insert({2, generateCube(2, meshgrid, lattice, 1.0, 6.4, 9, 4.3, 1.0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, generateCube(1, meshgrid2, lattice, 1.0, 4.4, 10, 5.3, 1.0)});


    ASSERT_THROW(RealSpaceSolver(vWannMap, cWannMap, pot), runtime_error);
}

TEST(ImplementationTest, Consistency4)
{
    Potential* pot = new CoulombPotential();

    double L = 15.;
    vector<int> dim{100,100,109};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,1.,0};
    lattice[1] = vector<double>{0,L,0.2};
    lattice[2] = vector<double>{0.4,1.,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    origin[2] = 3.;
    shared_ptr<RealMeshgrid> meshgrid2 = make_shared<RealMeshgrid>(dim, lattice, origin);

    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, generateCube(0, meshgrid, lattice, 1.0, 4, 10, 5, 1.0)});
    cWannMap.insert({2, generateCube(2, meshgrid, lattice, 1.0, 6.4, 9, 4.3, 1.0)});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({1, generateCube(1, meshgrid2, lattice, 1.0, 4.4, 10, 5.3, 1.0)});


    ASSERT_THROW(RealSpaceSolver(vWannMap, cWannMap, pot), runtime_error);
}


TEST(ImplementationTest, Monopole1) {

    double L = 5.;

    vector<int> dim{51,101,51};
    vector<double> origin{0,0,0};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann{ generateCube(0, meshgrid, lattice, L/2-0.2, L/2+0.1, L/2-0.2, L/2, 341.3)};

    unique_ptr<double[], free_deleter> density{ joinedDensity(wann,wann,vector<int>{0,0,0}) };
    Monopole mono = getMonopole(density.get(),wann.getMeshgrid());

    EXPECT_NEAR(L/2+0.1, mono.x, L/(dim[0]-1)/2);
    EXPECT_NEAR(L/2-0.2, mono.y, L/(dim[1]-1)/2);
    EXPECT_NEAR(L/2, mono.z, L/(dim[2]-1)/2);
    EXPECT_NEAR(341.3, mono.charge, 1e-8);
}

TEST(ImplementationTest, Monopole2) {

    double L = 5.;

    vector<int> dim{102,101,105};
    vector<double> origin{0,0,0};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,1,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,2,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann{ generateSphere(0, meshgrid, lattice, 1., L/2-0.2, L/2+0.1, L/2+0.1, 22.)};

    XSF_controller::save_file("wann.xsf", wann);

    unique_ptr<double[], free_deleter> density{ joinedDensity(wann,wann,vector<int>{0,0,0}) };
    Monopole mono = getMonopole(density.get(),wann.getMeshgrid());    EXPECT_NEAR(L/2-0.2, mono.x, L/(dim[0]-1)/2);
    EXPECT_NEAR(L/2+0.1, mono.y, L/(dim[1]-1)/2);
    EXPECT_NEAR(L/2+0.1, mono.z, L/(dim[2]-1)/2);
    EXPECT_NEAR(22., mono.charge, 1e-8);
}

TEST(ImplementationTest, Monopole3) {

    double L = 5.;
    vector<int> dim{101,101,101};
    vector<double> origin{0,0,0};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann{ generateSphere(0, meshgrid, lattice, L/2, L/2, L/2, L/2, -1.3)};

    // XSF_controller::save_file("wann.xsf", wann);

    unique_ptr<double[], free_deleter> density{ joinedDensity(wann,wann,vector<int>{0,0,0}) };
    Monopole mono = getMonopole(density.get(),wann.getMeshgrid());
    EXPECT_NEAR(L/2, mono.x, L/(dim[0]-1)/2);
    EXPECT_NEAR(L/2, mono.y, L/(dim[1]-1)/2);
    EXPECT_NEAR(L/2, mono.z, L/(dim[2]-1)/2);
    EXPECT_NEAR(1.3, mono.charge, 1e-8);  // charge is always positive
}

TEST(ImplementationTest, Monopole4)
{
    vector<int> R{0,0,0};
    vector<int> dim{191,191,191};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{8.2,0,0};
    lattice[1] = vector<double>{-4.1,7.14,0};
    lattice[2] = vector<double>{-5.,0.0,10.0};

    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateSphere(0, meshgrid, lattice, 1., 1.0, 1.00, 1, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.,2.5, 3.7, 4.1, 1.0)};

    // wann1->printSummary();
    unique_ptr<double[], free_deleter> density0{ joinedDensity(wann0,wann0,vector<int>{0,0,0}) };
    unique_ptr<double[], free_deleter> density1{ joinedDensity(wann1,wann1,vector<int>{0,0,0}) };

    Monopole mono0 = getMonopole(density0.get(), wann0.getMeshgrid());
    Monopole mono1 = getMonopole(density1.get(), wann1.getMeshgrid());

    // cout << mono0.charge << endl;
    // cout << mono0.x << " " << mono0.y << " " << mono0.z << endl;
    // cout << mono1.charge << endl;
    // cout << mono1.x << " " << mono1.y << " " << mono1.z << endl;

    ASSERT_NEAR(mono0.charge, 1.0, 1e-8);
    ASSERT_NEAR(mono0.x, 1.0, 1e-3);
    ASSERT_NEAR(mono0.y, 1.0, 1e-3);
    ASSERT_NEAR(mono0.z, 1.0, 1e-3);

    ASSERT_NEAR(mono1.charge, 1.0, 1e-8);
    ASSERT_NEAR(mono1.x, 2.5, 1e-3);
    ASSERT_NEAR(mono1.y, 3.7, 1e-3);
    ASSERT_NEAR(mono1.z, 4.1, 1e-3);


    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    Potential* pot = new CoulombPotential();
    MonoMonoSolver impl1(vWannMap, cWannMap);

    auto W = Integral(0, 0, 1, 1, vector<int>{0,0,0}, vector<int>{0,0,0}, vector<int>{0,0,0});

    impl1.calculate(W, false,1,10);
    ASSERT_NEAR(W.value, pot->realCart(mono0.x-mono1.x, mono0.y-mono1.y, mono0.z-mono1.z), 1e-6);

    for (int i=-4; i<=4; i++) {

        if (i==0) continue;

        vector<Integral> W_all = vector<Integral>(5);
        W_all[0] = Integral(0, 0, 1, 1, vector<int>{i,0,0}, vector<int>{0,0,0}, vector<int>{0,0,0});
        W_all[1] = Integral(0, 0, 1, 1, vector<int>{0,i,0}, vector<int>{0,0,0}, vector<int>{0,0,0});
        W_all[2] = Integral(0, 0, 1, 1, vector<int>{i,i,0}, vector<int>{0,0,0}, vector<int>{0,0,0});
        W_all[3] = Integral(0, 0, 1, 1, vector<int>{0,0,i}, vector<int>{0,0,0}, vector<int>{0,0,0});
        W_all[4] = Integral(0, 0, 1, 1, vector<int>{i,i,i}, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W_all, false,1,10);

        EXPECT_NEAR(W_all[0].value, pot->realCart(mono0.x-(mono1.x+8.2*i), mono0.y-mono1.y, mono0.z-mono1.z), 1e-6);
        EXPECT_NEAR(W_all[1].value, pot->realCart(mono0.x-(mono1.x-4.1*i), mono0.y-(mono1.y+7.14*i), mono0.z-(mono1.z-0.)), 1e-6);
        EXPECT_NEAR(W_all[2].value, pot->realCart(mono0.x-(mono1.x+8.2*i-4.1*i), mono0.y-(mono1.y+7.14*i), mono0.z-(mono1.z-0.)), 1e-6);
        EXPECT_NEAR(W_all[3].value, pot->realCart(mono0.x-(mono1.x-5.0*i), mono0.y-(mono1.y), mono0.z-(mono1.z+10.*i)), 1e-6);
        EXPECT_NEAR(W_all[4].value, pot->realCart(mono0.x-(mono1.x+8.2*i-4.1*i-5*i), mono0.y-(mono1.y+7.14*i), mono0.z-(mono1.z+10.*i)), 1e-6);
    }


    delete pot;
}


TEST(ImplementationTest, FourierRealGaussCubicCell) {


    double L = 10.;
    vector<int> R{0,0,0};
    vector<int> dim{61,61,61};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.5, L/2+1.6, L/2+1.7, L/2, 1.0)};
    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    for (int i=0; i<5; i++) {

        pot = new GaussPotential(0.8 + i * 0.1);
        // set mappings and potential in solver
        YukawaSolver impl1(vWannMap, cWannMap, pot);
        RealSpaceSolver impl2(vWannMap, cWannMap, pot);

        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        EXPECT_NEAR(W1.value, W2.value, 1e-6);

        delete pot;
    }
}

TEST(ImplementationTest, FourierRealGaussArbCell)
{
    double L = 10.;
    vector<int> R{0,0,0};
    vector<int> dim{51,51,51};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,1,0};
    lattice[1] = vector<double>{2,L,0};
    lattice[2] = vector<double>{0,0,0.9*L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};

    // XSF_controller::save_file("wann0.xsf", wann0);
    // XSF_controller::save_file("wann1.xsf", wann1);
    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    for (int i=0; i<5; i++) {

        pot = new GaussPotential(1.0 + i * 0.1);
        YukawaSolver impl1(vWannMap, cWannMap, pot);
        RealSpaceSolver impl2(vWannMap, cWannMap, pot);


        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});


        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        EXPECT_NEAR(W1.value, W2.value, 1e-6);

        delete pot;
    }
}


TEST(ImplementationTest, FourierRealGaussUnitcell)
{
    double L = 10.;
    vector<int> R{0,0,0};
    vector<int> dim{51,51,51};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{L/2,0,0};
    unitcell[1] = vector<double>{0,L/2,0};
    unitcell[2] = vector<double>{0,1.,L};

    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2-1.1, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.5, L/3, L/2, L/2+1.1, 1.0)};

    wann0.setUnitcell(unitcell);
    wann1.setUnitcell(unitcell);
    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    for (int i=0; i<5; i++) {

        pot = new GaussPotential(0.6 + i * 0.1);
        // set mappings and potential in solver
        YukawaSolver impl1(vWannMap, cWannMap, pot);
        RealSpaceSolver impl2(vWannMap, cWannMap, pot);

        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        EXPECT_NEAR(W1.value, W2.value, 1e-6);

        delete pot;
    }
}


TEST(ImplementationTest, FourierGaussDiffR1)
{
    double L_unitcell = 4.;
    double L = 4*L_unitcell;

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
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    WannierFunction wann1_shifted{ generateSphere(1, meshgrid1, lattice, 1.5, L_unitcell/2+L_unitcell, L_unitcell/2, L_unitcell/2, 1.0)};

    wann0.setUnitcell(unitcell);
    wann1.setUnitcell(unitcell);
    wann1_shifted.setUnitcell(unitcell);

    // XSF_controller::save_file("wann0.xsf", wann0);
    // XSF_controller::save_file("wann1.xsf", wann1);
    // XSF_controller::save_file("wann1_shifted.xsf", wann1_shifted);
    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});
    vWannMap.insert({2,std::move(wann1_shifted)});

    for (int i=0; i<1; i++) {

        pot = new GaussPotential(0.8 + i * 0.1);
        YukawaSolver impl1(vWannMap, cWannMap, pot);

        auto W1 = Integral(0, 0, 1, 1, R1, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 2, 2, vector<int>{0,0,0}, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl1.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        EXPECT_NEAR(W1.value, W2.value, 1e-6);

        delete pot;
    }
}

TEST(ImplementationTest, RealGaussDiffR)
{
    double L_unitcell = 4.;
    double L = 3*L_unitcell;

    vector<int> dim{81,81,81};
    vector<int> R1{1,1,0};
    vector<double> origin{-L_unitcell,-L_unitcell,-L_unitcell};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{L_unitcell,0,0};
    unitcell[1] = vector<double>{0,L_unitcell,0};
    unitcell[2] = vector<double>{0,0,L_unitcell};

    vector<double> vec_shift = matVecMul3x3(unitcell, R1);

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid1 = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0{ generateCube(0, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid1, lattice, 1.5, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0)};
    WannierFunction wann1_shifted{ generateSphere(1, meshgrid1, lattice, 1.5, L_unitcell/2+vec_shift[0], L_unitcell/2+vec_shift[1], L_unitcell/2+vec_shift[2], 1.0)};

    wann0.setUnitcell(unitcell);
    wann1.setUnitcell(unitcell);
    wann1_shifted.setUnitcell(unitcell);

    // XSF_controller::save_file("wann0.xsf", wann0);
    // XSF_controller::save_file("wann1.xsf", wann1);
    // XSF_controller::save_file("wann1_shifted.xsf", wann1_shifted);

    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});
    vWannMap.insert({2,std::move(wann1_shifted)});

    for (int i=0; i<1; i++) {

        pot = new GaussPotential(0.8 + i * 0.1);

        // set mappings and potential in solver
        RealSpaceSolver impl1(vWannMap, cWannMap, pot);


        auto W1 = Integral(0, 0, 1, 1, R1, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 2, 2, vector<int>{0,0,0}, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl1.calculate(W2, false,1,10);


        // cout << W1.value << " ?= " << W2.value << endl;
        ASSERT_NEAR(W1.value, W2.value, 1e-5);

        delete pot;
    }
}


// TEST(ImplementationTest, FourierRealGaussDiffR) {
//     /**
//      * Test Fourier method for a specific shift when lattice of the
//      * supercell is not compatible with unitcell.
//      *
//      * Larger super cell is created.
//      **/
//     Solver* impl1 = new YukawaSolver();
//     Solver* impl2 = new RealSpaceSolver();

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
//     shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

//     WannierFunction* wann0_ = generateCube(0, meshgrid, lattice, 1.8, L_unitcell/2, L_unitcell/2, L_unitcell/2, 1.0);
//     WannierFunction* wann1_ = generateSphere(1, meshgrid, lattice, 0.95, 0.95, L_unitcell/2, L_unitcell/2, 1.0);
//     WannierFunction* wann1_shifted_ = generateSphere(1, meshgrid, lattice, 0.95, 0.95+L_unitcell, L_unitcell/2, L_unitcell/2, 1.0);

//     wann0_->setUnitcell(unitcell);
//     wann1_->setUnitcell(unitcell);
//     wann1_shifted_->setUnitcell(unitcell);

//     // vector<double> tmp = wann0_->getLatticeInUnitcellBasis();
//     // cout << tmp[0] << "\t" << tmp[1] << "\t" << tmp[2] << endl;

//     WannierFunction* wann0 = createLargerSupercell(wann0_, vector<int>{3,3,3});
//     WannierFunction* wann1 = createLargerSupercell(wann1_, vector<int>{3,3,3});
//     WannierFunction* wann1_shifted = createLargerSupercell(wann1_shifted_, vector<int>{3,3,3});

//     // XSF_controller::save_file("wann0.xsf", wann0);
//     // XSF_controller::save_file("wann1.xsf", wann1);
//     // XSF_controller::save_file("wann1_shifted.xsf", wann1_shifted);


//     CoulombIntegral* W1 = 0;
//     CoulombIntegral* W2 = 0;
//     CoulombIntegral* W_fail = 0;
//     Potential* pot = 0;

//     for (int i=0; i<1; i++) {

//         pot = new GaussPotential(1.0 + i * 0.1);

//         W_fail = new CoulombIntegral(wann0_, wann0_, wann1_, wann1_, vector<int>{3,0,0}, pot);  // should fail because supercell is not large enough

//         impl1->calculate(W_fail, false,1,10);
//         ASSERT_EQ(W_fail->getStatus(), Status::failed);
//         ASSERT_EQ(W_fail->getErrorMsg(), "The supercell is not large enough to protect against aliasing!");

//         //W1 = new CoulombIntegral(impl1, wann0, wann0, wann1, wann1, R, pot);
//         W1 = new CoulombIntegral(wann0, wann0, wann1_shifted, wann1_shifted, vector<int>{0,0,0}, pot);
//         W2 = new CoulombIntegral(wann0_, wann0_, wann1_, wann1_, R, pot);
//         //W2 = new CoulombIntegral(impl2, wann0, wann0, wann1, wann1, R, pot);
//         //W2 = new CoulombIntegral(impl2, wann0_, wann0_, wann1_shifted_, wann1_shifted_, vector<int>{0,0,0}, pot);

//         impl1->calculate(W1);
//         impl2->calculate(W2);

//         ASSERT_LT(abs(W1->getValue().imag()), 1e-10);
//         ASSERT_LT(abs(W2->getValue().imag()), 1e-10);
//         cout << W1->getValue().real() << " ?= " << W2->getValue().real() << endl;
//         ASSERT_NEAR(W1->getValue().real(), W2->getValue().real(), 1e-4);

//         delete W1, W2, pot;
//     }
// }

TEST(ImplementationTest, FourierGaussDiffR3)
{
    double L = 5.;
    vector<int> dim{50,52,51};
    vector<int> R{1,0,0};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0_{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1_{ generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    createLargerSupercell(wann0, vector<int>{4,3,3});
    createLargerSupercell(wann1, vector<int>{4,3,3});

    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap_{};
    map< int,WannierFunction > vWannMap_{};
    cWannMap_.insert({0,std::move(wann0_)});
    vWannMap_.insert({1,std::move(wann1_)});

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    for (int i=0; i<1; i++) {

        pot = new GaussPotential(2/L + i * 0.1);
        // set mappings and potential in solver
        YukawaSolver impl1(vWannMap, cWannMap, pot); // use supercell
        RealSpaceSolver impl2(vWannMap_, cWannMap_, pot);

        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        EXPECT_NEAR(W1.value, W2.value, 1e-6);

        delete pot;
    }
}


TEST(ImplementationTest, FourierRealYukawaDiffR)
{
    double L = 5.;
    vector<int> dim{51,51,51};
    vector<int> R{0,1,0};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    WannierFunction wann0_{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1_{ generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann1{ generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    createLargerSupercell(wann0, vector<int>{4,4,4});
    createLargerSupercell(wann1, vector<int>{4,4,4});

    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap_{};
    map< int,WannierFunction > vWannMap_{};
    cWannMap_.insert({0,std::move(wann0_)});
    vWannMap_.insert({1,std::move(wann1_)});

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0,std::move(wann0)});
    vWannMap.insert({1,std::move(wann1)});

    for (int i=0; i<1; i++) {

        //cout << "i = " << i << endl;

        pot = new YukawaPotential(0.8 + i * 0.1);
        // set mappings and potential in solver
        YukawaSolver impl1(vWannMap, cWannMap, pot);  // use supercell
        RealSpaceSolver impl2(vWannMap_, cWannMap_, pot);  // normal cell

        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);


        cout << W1.value << " ?= " << W2.value << endl;
        ASSERT_NEAR(W1.value, W2.value, 1e-5);

        delete pot;

    }
}


// TEST(ImplementationTest, MonoGaussDiffR)
// {
//     double L = 7.;
//     vector<int> dim{51,51,51};
//     vector<int> R{0,-1,1};
//     vector<double> origin{0,0,0};
//     vector< vector<double> > lattice(3);
//     lattice[0] = vector<double>{L,0,0};
//     lattice[1] = vector<double>{0,L,0};
//     lattice[2] = vector<double>{0,0,L};
//     shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

//     WannierFunction* wann0 = generateSphere(0, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0);
//     WannierFunction* wann1 = generateSphere(1, meshgrid, lattice, 1.5, L/2, L/2, L/2, 1.0);

//     Potential* pot = 0;

//     // create mappings
//     map< int,WannierFunction* >* cWannMap = new map< int,WannierFunction* >;
//     map< int,WannierFunction* >* vWannMap = new map< int,WannierFunction* >;
//     cWannMap->insert({0, wann0});
//     vWannMap->insert({1, wann1});

//     for (int i=0; i<3; i++) {

//         pot = new GaussPotential(0.8 + i * 0.1);

//         // set mappings and potential in solver
//         Solver* impl1 = new MonoMonoSolver(vWannMap, cWannMap, pot);  // MonoMono is now restricted to CoulombPotential
//         Solver* impl2 = new RealSpaceSolver(vWannMap, cWannMap, pot);

//         auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
//         auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

//         impl1->calculate(W1, false,1,10);
//         impl2->calculate(W2, false,1,10);

//         ASSERT_NEAR(W1.value, W2.value, 1e-8);

//         delete pot;
//         delete impl1;
//         delete impl2;
//     }
// }

TEST(ImplementationTest, FourierYukawa)
{
    double L = 10.;
    vector<int> R{0,0,0};
    vector<int> dim{61,61,61};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);
    Potential* pot = nullptr;

    // create mappings
    map< int,WannierFunction > cWannMap{};
    map< int,WannierFunction > vWannMap{};
    cWannMap.insert({0, generateCube(0, meshgrid, lattice, 1., 3.9, L/2, L/2, 1.0)});
    vWannMap.insert({1, generateSphere(1, meshgrid, lattice, 1., 6.1, L/2, L/2, 1.0)});

    for (int i=0; i<5; i++) {

        //cout << "i=" << i << endl;

        pot = new YukawaPotential(1.5 + i * 0.1);
        // set mappings and potential in solver
        YukawaSolver impl1(vWannMap, cWannMap, pot);
        RealSpaceSolver impl2(vWannMap, cWannMap, pot);

        auto W1 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});
        auto W2 = Integral(0, 0, 1, 1, R, vector<int>{0,0,0}, vector<int>{0,0,0});

        impl1.calculate(W1, false,1,10);
        impl2.calculate(W2, false,1,10);

        cout << W1.value << " ?= " << W2.value << endl;
        ASSERT_NEAR(W1.value, W2.value, 1e-4);

        delete pot;
    }
}

TEST(YukawaPotential, calc_yukawa_screening_factor)
{
    double mean_density = 0.037516174403821975;
    double epsilon = 2.0;
    double alpha = 1.563;

    ASSERT_NEAR( YukawaPotential::calc_yukawa_screening_factor(mean_density,epsilon, alpha), 1.7856510467611366, 1e-10);
}