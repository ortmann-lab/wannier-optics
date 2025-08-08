

#include <coulombIntegral.h>
#include <wannierfunction.h>
#include <potential.h>
#include <solver.h>
#include <generator.cpp>
#include <gtest/gtest.h>
#include <iostream>
using namespace std;


WannierFunction generateGauss3(int id, shared_ptr<RealMeshgrid> meshgrid, vector< vector<double> > lattice, double alpha, double x0, double y0, double z0, double norm) {

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

TEST(ImplementationLocalFieldEffectsTest, CubicAnalytic)
{
    double L = 11;
    vector<int> dim{40,40,40};
    vector<double> origin{100,-52.,7};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double alpha = 1.0;

    double x0 = L/2 + origin[0];
    double y0 = L/2 + origin[1];
    double z0 = L/2 + origin[2];


    double shouldbe = 42.08457689412335 / (L*L*L);  // calculated with python script ('test_case_local_field_effects.py')

    WannierFunction wann0{ generateGauss3(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0) };
    EXPECT_NEAR(wann0.getNorm(), 1.0, 1e-10);

    // supercell
    WannierFunction wann0_sc{ generateGauss3(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0) };
    createLargerSupercell(wann0_sc,vector<int>{2,3,2});
    EXPECT_NEAR(wann0_sc.getNorm(), 1.0, 1e-10);

    // create mappings
    map< int,WannierFunction > cWannMap{};
    cWannMap.insert({0, std::move(wann0_sc)});

    // create the same object twice (currently there is no copy operator)
    WannierFunction wann0_sc2{ generateGauss3(0, meshgrid, lattice, alpha, x0,y0,z0, 1.0) };
    createLargerSupercell(wann0_sc2,vector<int>{2,3,2});

    map< int,WannierFunction > vWannMap{};
    vWannMap.insert({0, std::move(wann0_sc2)});

    // set mappings in solver
    LocalFieldEffectsSolver implFourier(vWannMap, cWannMap);

    // define Integrals
    auto WFouier = Integral(0, 0, 0, 0, vector<int>{0,0,0}, vector<int>{0,0,0}, vector<int>{0,0,0});
    implFourier.calculate(WFouier, false,1,10);

    cout << WFouier.value << " ?= " << shouldbe << endl;

    EXPECT_NEAR(WFouier.value, shouldbe, 1e-8);
}