#ifndef GENERATOR_CPP
#define GENERATOR_CPP

/**
 * @file generator.cpp
 * @author Konrad Merkel
 * @brief Generate test (Wannier) function.
 *
 * This is only for test purposes.
 *
 */

#include <fstream>
#include <iostream>
#include <cmath>
#include "wannierfunction.h"
#include "filehandler.h"
using namespace std;

double sphere(double x, double y, double z, double x0, double y0, double z0, double radius=1.)
{
    if (pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2) <= pow(radius,2))
        return 1.0;
    else
        return 0.0;
}

double gauss(double x, double y, double z, double x0, double y0, double z0, double sigma)
{
    return exp(- (pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2))/(2* sigma*sigma)) /
        (pow(2*M_PI, 3./2)*pow(sigma, 3));
}

double cube(double x, double y, double z, double x0, double y0, double z0, double A=1.)
{
    if ((abs(x-x0) < A/2) && (abs(y-y0) < A/2) && (abs(z-z0) < A/2))
        return 1.;
    else
        return 0.0;
}


WannierFunction generateSphere(int id, shared_ptr<RealMeshgrid> meshgrid, vector< vector<double> > lattice, double radius, double x0, double y0, double z0, double charge) {

    // fill arrays (initialize inputs)
    // column major order (Fortan style !!!)
    unique_ptr<double[], free_deleter> data{ (double*) malloc(sizeof(double)*meshgrid->getNumDataPoints()) };
    double x,y,z;
    for (int i=0; i < meshgrid->getNumDataPoints(); i++){
        meshgrid->xyz(i,x,y,z);
        data[i] = sphere(x,y,z, x0, y0, z0, radius);
    }

    WannierFunction wannier{ WannierFunction(id, meshgrid, std::move(data), lattice) };
    wannier.normalize(charge);
    return wannier;
}

WannierFunction generateGauss(int id, shared_ptr<RealMeshgrid> meshgrid, vector< vector<double> > lattice, double sigma, double x0, double y0, double z0, double charge) {

    // fill arrays (initialize inputs)
    // column major order (Fortan style !!!)
    unique_ptr<double[], free_deleter> data{ (double*) malloc(sizeof(double)*meshgrid->getNumDataPoints()) };
    double x,y,z;
    for (int i=0; i < meshgrid->getNumDataPoints(); i++){
        meshgrid->xyz(i,x,y,z);
        data[i] = gauss(x,y,z, x0, y0, z0, sigma);
    }

    WannierFunction wannier{ WannierFunction(id, meshgrid, std::move(data), lattice) };
    wannier.normalize(charge);
    // normalizeWann(wannier, charge);
    return wannier;
}


WannierFunction generateCube(int id, shared_ptr<RealMeshgrid> meshgrid, vector< vector<double> > lattice, double A, double x0, double y0, double z0, double charge) {

    // fill arrays (initialize inputs)
    // column major order (Fortan style !!!)
    unique_ptr<double[], free_deleter> data{ (double*) malloc(sizeof(double)*meshgrid->getNumDataPoints()) };
    double x,y,z;
    for (int i=0; i < meshgrid->getNumDataPoints(); i++){
        meshgrid->xyz(i,x,y,z);
        data[i] = cube(x,y,z, x0, y0, z0, A);
    }

    WannierFunction wannier{ WannierFunction(id, meshgrid, std::move(data), lattice) };
    wannier.normalize(charge);
    // normalizeWann(wannier, charge);
    return wannier;
}


#endif  // GENERATOR_CPP