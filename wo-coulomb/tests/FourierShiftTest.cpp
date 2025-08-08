// tests.cpp
#include <wannierfunction.h>
#include <generator.cpp>
#include <meshgrid.h>
#include <vector>
#include <gtest/gtest.h>
#include <fftw3.h>


void multiply_(fftw_complex &a, fftw_complex &b, fftw_complex &result) {
    result[0] = a[0] * b[0] -a[1] * b[1] ;
    result[1] = a[0] * b[1] + a[1] * b[0];
}

TEST(FourierShiftTest, ShiftTheorem) {

    double L = 5.;
    vector<int> dim{50,81,43};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid_real = make_shared<RealMeshgrid>(dim, lattice, origin);
    ReciprocalMeshgrid meshgrid_rec(meshgrid_real.get());


    // double dx = L/dim[0];
    double dy = L/dim[1];
    // double dz = L/dim[2];
    double x0 = 1.;
    double y0 = 1.;
    double z0 = 1;
    vector<double> vec_shift{0.,10*dy,0.};  // shift by exactly 10 data points


    WannierFunction wann{ generateCube(0, meshgrid_real, lattice, 2, x0, y0, z0, 1.0)};
    WannierFunction wann_shifted{ generateCube(0, meshgrid_real, lattice, 2, x0+vec_shift[0], y0+vec_shift[1], z0+vec_shift[2], 1.0)};
    double N = meshgrid_rec.getNumDataPoints();

    fftw_complex* f1_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f2_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f1_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f2_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* phase = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // plan needs to be created before initializing inputs !!!
    // dimensions are reversed because our datagrids are in Fortran (column-major) format!
    fftw_plan p_f1 = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f1_r, f1_q, FFTW_FORWARD, FFTW_ESTIMATE); // BACKWARD because we need f1(-q) !!!!
    fftw_plan p_f2 = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f2_r, f2_q, FFTW_FORWARD, FFTW_ESTIMATE);  // f2(+q)

    vector< vector<double> > recLattice = meshgrid_rec.getLattice();

    // cout << N << endl;
    // cout << x0 << " " <<  y0 << " " << z0 << " " << endl;
    // cout << vec_shift[0] << " " <<  vec_shift[1] << " " << vec_shift[2] << " " << endl;
    // cout << "recLattice: \n";
    // printMat3x3(recLattice);

    // fill arrays (initialize inputs)
    double qx,qy,qz;
    for (int i=0; i<N; i++){
        meshgrid_rec.xyz(i,qx,qy,qz);
        // real part
        f1_r[i][0] = wann.getValue()[i];
        f2_r[i][0] = wann_shifted.getValue()[i];
        phase[i][0] = cos(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);


        // imaginary part
        f1_r[i][1] = 0.0;
        f2_r[i][1] = 0.0;
        phase[i][1] = -sin(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);  // minus because of FFTW_FORWARD

    }

    // XSF_controller::save_file("phase_r.xsf", &meshgrid, phase_r);
    // XSF_controller::save_file("phase_i.xsf", &meshgrid, phase_i);

    // Fourier transform wave functions
    fftw_execute(p_f1);
    fftw_execute(p_f2);

    fftw_complex value;
    for (int i=0; i<N; i++) {

        multiply_(f1_q[i], phase[i], value);

        ASSERT_NEAR(f2_q[i][0], value[0], 1e-10);
        ASSERT_NEAR(f2_q[i][1], value[1], 1e-10);
    }

}




TEST(FourierShiftTest, ShiftInReciprocalSpace) {

    double L = 5.;
    vector<int> dim{50,50,50}; // does not work for 51,51,51 because mesh and shift are then not compatible
    vector<double> vec_shift{0.5,0,0};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid_real = make_shared<RealMeshgrid>(dim, lattice, origin);
    ReciprocalMeshgrid meshgrid_rec(meshgrid_real.get());

    WannierFunction wann{ generateCube(0, meshgrid_real, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann_fourier_shifted{ generateCube(0, meshgrid_real, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    WannierFunction wann_shifted{ generateCube(0, meshgrid_real, lattice, 1.5, L/2+vec_shift[0], L/2+vec_shift[1], L/2+vec_shift[2], 1.0)};
    double N = meshgrid_rec.getNumDataPoints();

    fftw_complex* f1_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f1_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f_shift_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f_shift_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* phase = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // plan needs to be created before initializing inputs !!!
    // dimensions are reversed because our datagrids are in Fortran (column-major) format!
    fftw_plan p_f1 = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f1_r, f1_q, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan back = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f_shift_q, f_shift_r, FFTW_BACKWARD, FFTW_ESTIMATE);

    // fill arrays (initialize inputs)
    double qx,qy,qz;
    for (int i=0; i<N; i++){
        meshgrid_rec.xyz(i,qx,qy,qz);
        // real part
        f1_r[i][0] = wann.getValue()[i];
        phase[i][0] = cos(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);

        // imaginary part
        f1_r[i][1] = 0.0;
        phase[i][1] = -sin(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);

    }

    // Fourier transform wave functions
    fftw_execute(p_f1);

    // shift in Fourier space by multiplying with a phase
    for (int i=0; i<N; i++) {
        multiply_(f1_q[i], phase[i], f_shift_q[i]);
    }

    fftw_execute(back);  // back trafo
    double* data = (double*) malloc(sizeof(double) * N);
    double const* shouldbe = wann_shifted.getValue();
    for (int i=0; i<N; i++){

        // imaginary part
        ASSERT_LT(abs(f_shift_r[i][1]), 1e-8);

        // real part
        data[i] = f_shift_r[i][0] / N;
        ASSERT_NEAR(data[i], shouldbe[i], 1e-8);
    }
    wann_fourier_shifted.setValue(data);

    // XSF_controller::save_file("wann_fourier_shifted.xsf", wann_fourier_shifted);


}




TEST(FourierShiftTest, CellShiftInReciprocalSpace) {
    double L = 5.;
    vector<int> dim{60,54,46};
    vector<double> vec_shift{L,0,0};  // this is a workaround. Is this really correct?

    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    shared_ptr<RealMeshgrid> meshgrid_real = make_shared<RealMeshgrid>(dim, lattice, origin);
    ReciprocalMeshgrid meshgrid_rec(meshgrid_real.get());

    WannierFunction wann{ generateCube(0, meshgrid_real, lattice, 1.5, L/2, L/2, L/2, 1.0)};
    double N = meshgrid_rec.getNumDataPoints();

    fftw_complex* f1_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f1_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f_shift_r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* f_shift_q = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* phase = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // plan needs to be created before initializing inputs !!!
    // dimensions are reversed because our datagrids are in Fortran (column-major) format!
    fftw_plan p_f1 = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f1_r, f1_q, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan back = fftw_plan_dft_3d(meshgrid_real->getDim()[2], meshgrid_real->getDim()[1], meshgrid_real->getDim()[0], f_shift_q, f_shift_r, FFTW_BACKWARD, FFTW_ESTIMATE);

    // fill arrays (initialize inputs)
    double qx,qy,qz;
    for (int i=0; i<N; i++){
        meshgrid_rec.xyz(i,qx,qy,qz);
        // real part
        f1_r[i][0] = wann.getValue()[i];
        phase[i][0] = cos(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);

        // imaginary part
        f1_r[i][1] = 0.0;
        phase[i][1] = -sin(qx * vec_shift[0] + qy * vec_shift[1] + qz * vec_shift[2]);

    }

    // Fourier transform wave functions
    fftw_execute(p_f1);

    // shift in Fourier space by multiplying with a phase
    for (int i=0; i<N; i++) {
        multiply_(f1_q[i], phase[i], f_shift_q[i]);
    }

    fftw_execute(back);  // back trafo
    double* data = (double*) malloc(sizeof(double) * N);
    double const* shouldbe = wann.getValue();
    for (int i=0; i<N; i++){

        // imaginary part
        ASSERT_LT(abs(f_shift_r[i][1]), 1e-8);

        // real part
        data[i] = f_shift_r[i][0] / N;
        ASSERT_NEAR(data[i], shouldbe[i], 1e-8);
    }
}