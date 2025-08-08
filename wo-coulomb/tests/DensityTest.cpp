// tests.cpp
#include <density.h>
#include <meshgrid.h>
#include <wannierfunction.h>
#include <solver.h>
#include <vector>
#include <fftw3.h>
#include <gtest/gtest.h>


TEST(DensityTest, getExtend_cube)
{
    double L = 15.;
    vector<int> R{0,0,0};
    vector<int> dim{100,100,100};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double a = 1.5;

    WannierFunction wann0{ generateCube(0, meshgrid, lattice, a, L/2-a/2-0.1+origin[0], L/2-a/2+origin[1], L/2-a/2+origin[2], 1.0)};
    RealMeshgrid* mesh = wann0.getMeshgrid();
    mesh->createMeshgridArrays();

    unique_ptr<double[], free_deleter>density{ joinedDensity(wann0, wann0, vector<int>{0,0,0}) };
    Monopole mono = getMonopole(density.get(), mesh);
    cout << "center = " << mono.x << " " << mono.y << " " << mono.z << endl;


    double extend = getExtend(density.get(), mesh, mono);

    ASSERT_LE(extend, 0.8);
    ASSERT_GE(extend, 0.6);
}


TEST(DensityTest, getExtend_sphere)
{
    double L = 15.;
    vector<int> R{0,0,0};
    vector<int> dim{100,100,100};
    vector<double> origin{-3,5.,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid = make_shared<RealMeshgrid>(dim, lattice, origin);

    double a = 1.5;

    WannierFunction wann0{ generateSphere(0, meshgrid, lattice, a, L/2-a/2-0.1+origin[0], L/2-a/2+origin[1], L/2-a/2+origin[2], 1.0)};
    RealMeshgrid* mesh = wann0.getMeshgrid();
    mesh->createMeshgridArrays();

    unique_ptr<double[], free_deleter>density{ joinedDensity(wann0, wann0, vector<int>{0,0,0}) };
    Monopole mono = getMonopole(density.get(), mesh);
    cout << "center = " << mono.x << " " << mono.y << " " << mono.z << endl;

    double extend = getExtend(density.get(), mesh, mono);

    ASSERT_LE(extend, 1.6);
    ASSERT_GE(extend, 1.4);
}