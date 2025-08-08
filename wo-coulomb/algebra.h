#ifndef ALGEBRA_H
#define ALGEBRA_H

/**
 * @file algebra.h
 * @author Konrad Merkel
 * @brief Some algebra routines to deal with 3x3 matrixes
 *
 */

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

using namespace std;

struct free_deleter{
    template <typename T>
    void operator()(T *p) const {
        std::free(const_cast<std::remove_const_t<T>*>(p));
    }
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double det3x3(vector< vector<double> > const& mat){
    return mat[0][0]*mat[1][1]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1]
        - mat[2][0]*mat[1][1]*mat[0][2] - mat[1][0]*mat[0][1]*mat[2][2] - mat[0][0]*mat[2][1]*mat[1][2];

}


double det3x3(double const mat[3][3]){
    return mat[0][0]*mat[1][1]*mat[2][2] + mat[0][1]*mat[1][2]*mat[2][0] + mat[0][2]*mat[1][0]*mat[2][1]
        - mat[2][0]*mat[1][1]*mat[0][2] - mat[1][0]*mat[0][1]*mat[2][2] - mat[0][0]*mat[2][1]*mat[1][2];

}

template< typename T>
T crossProduct(T const& vect_A, T const& vect_B)
{
    T cross_P(3);
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    return cross_P;
}

double dotProduct(double const vect_A[], double const vect_B[])
{
    double product = 0;

    // Loop for calculate cot product
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}

double dotProduct(vector<double> const& vect_A, vector<double> const& vect_B)
{
    double product = 0;
    if (vect_A.size() != vect_B.size()) {
        throw runtime_error("Vectors need to be of same dimension.");
    }

    // Loop for calculate cot product
    for (size_t i = 0; i < vect_A.size(); i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}


vector< vector<double> > invMat3x3(vector< vector<double> > const& mat){

    // check if matrix is a 3x3 matrix
    assert(mat.size() == 3);
    assert(mat[0].size() == 3);
    assert(mat[1].size() == 3);
    assert(mat[2].size() == 3);

    vector< vector<double> > inv;
    inv.resize(3, vector<double>(3,0.0));

    double det = det3x3(mat);
    inv[0][0] = (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) / det;
    inv[0][1] = (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2]) / det;
    inv[0][2] = (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]) / det;
    inv[1][0] = (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]) / det;
    inv[1][1] = (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0]) / det;
    inv[1][2] = (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2]) / det;
    inv[2][0] = (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]) / det;
    inv[2][1] = (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1]) / det;
    inv[2][2] = (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]) / det;

    return inv;
}

// void invMat3x3(double mat[3][3], double inv[3][3]){
//     double det = det3x3(mat);
//     inv[0][0] = (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) / det;
//     inv[0][1] = (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2]) / det;
//     inv[0][2] = (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]) / det;
//     inv[1][0] = (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]) / det;
//     inv[1][1] = (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0]) / det;
//     inv[1][2] = (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2]) / det;
//     inv[2][0] = (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]) / det;
//     inv[2][1] = (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1]) / det;
//     inv[2][2] = (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]) / det;
// }

vector< vector<double> > matMul3x3(vector< vector<double> > const& A, vector< vector<double> > const& B){
    vector< vector<double> > C(3);
    for (int i=0; i<3; i++) {
        C[i] = vector<double>{0,0,0};
    }

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<3; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

void matMul3x3(double const A[3][3], double const B[3][3], double C[3][3]){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            C[i][j] = 0.0;
        }
    }
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<3; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// void matVecMul3x3(double A[3][3], double v[], double res[]){
//     for (int i=0; i<3; i++){
//         res[i] = 0.0;
//     }
//     for (int i=0; i<3; i++){
//         for (int k=0; k<3; k++){
//             res[i] += A[i][k] * v[k];
//         }
//     }
// }

vector<double> inline matVecMul3x3(vector< vector<double> > const& A, vector<double> const& v){
    vector<double> res{0,0,0};
    for (int i=0; i<3; i++){
        for (int k=0; k<3; k++){
            res[i] += A[i][k] * v[k];
        }
    }
    return res;
}

vector<double> matVecMul3x3(vector< vector<double> > const& A, vector<int> const& v){
    vector<double> res{0,0,0};

    for (int i=0; i<3; i++){
        for (int k=0; k<3; k++){
            res[i] += A[i][k] * v[k];
        }
    }
    return res;
}

vector< vector<double> > transpose3x3(const vector< vector<double> >& A){
    vector< vector<double> > res(3);
    for (int i=0; i<3; i++){
        res[i] = vector<double>{0,0,0};
    }
    for (int i=0; i<3; i++){
        for (int k=0; k<3; k++){
            res[i][k] += A[k][i];
        }
    }
    return res;
}

void printMat3x3(double const A[3][3]){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            cout << "\t" << A[i][j];
        }
        cout << endl;
    }
}


void printMat3x3(vector<vector<double>> const& A){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            cout << "\t" << A[i][j];
        }
        cout << endl;
    }
}

double L2Norm(vector<double> const& v) {
    double norm = 0;
    for (auto c:v) norm = norm + c*c;
    return sqrt(norm);
}

#endif