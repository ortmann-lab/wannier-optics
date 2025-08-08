// tests.cpp
#include <meshgrid.h>
#include <solver.h>
#include <vector>
#include <fftw3.h>
#include <gtest/gtest.h>
 
double inline gaussian3D(double x, double y, double z, double x0, double y0, double z0, double alpha, double norm=1.) {
    return norm * pow(alpha / M_PI, 3./2.) * exp(-alpha * (  pow(x-x0, 2) + pow(y-y0, 2) + pow(z-z0, 2) ) );
}

void inline recGaussian3D(double qx, double qy, double qz, double x0, double y0, double z0, double alpha, double norm, fftw_complex& result) {
    result[0] = norm * exp( -(pow(qx, 2) + pow(qy, 2) + pow(qz, 2))/(4.*alpha) ) *cos(qx * x0 + qy * y0 + qz * z0);
    result[1] = norm * exp( -(pow(qx, 2) + pow(qy, 2) + pow(qz, 2))/(4.*alpha) ) *-sin(qx * x0 + qy * y0 + qz * z0); // minus because of FORWARD transform --> exp(-iqx)
}


TEST(MeshgridTest, EqualOperator) {

    vector<int> dim{81,61,101};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,1.,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0,0.5,5.};
    RealMeshgrid meshgrid1(dim, lattice, origin);

    lattice[0][1] = 2.35;
    RealMeshgrid meshgrid2(dim, lattice, origin);

    ASSERT_FALSE(meshgrid1 == meshgrid2);
    ASSERT_TRUE(meshgrid1 != meshgrid2);

    origin[1] = 1.;
    RealMeshgrid meshgrid3(dim, lattice, origin);
    ASSERT_FALSE(meshgrid3 == meshgrid2);

}


TEST(MeshgridTest, getIndexTipple) {

    vector<int> dim{71,241,10};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,1.,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0,0.5,5.};
    RealMeshgrid meshgrid(dim, lattice, origin);

    int N = meshgrid.getNumDataPoints();
    ASSERT_EQ(N, 71*241*10);

    vector<int> dim2 = meshgrid.getDim();
    ASSERT_EQ(dim2[0], 71);
    ASSERT_EQ(dim2[1], 241);
    ASSERT_EQ(dim2[2], 10);


    int i,j,k;
    for (int l=0; l<N; l++) {
        meshgrid.getIndexTipple(l,i,j,k);
        //cout << ijk[0] << " " << ijk[1] << " " << ijk[2] << endl;
        ASSERT_EQ(i + dim[0]*( j + dim[1]*k ), l);
    }

}


TEST(MeshgridTest, getIndexTipple2) {

    vector<int> dim{51,50,50};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{7.,1.,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0.1,0.5,5.};
    RealMeshgrid meshgrid(dim, lattice, origin);

    int N = meshgrid.getNumDataPoints();


    int i,j,k;
    for (int l=0; l<N; l++) {
        meshgrid.getIndexTipple(l,i,j,k);
        //cout << ijk[0] << " " << ijk[1] << " " << ijk[2] << endl;
        ASSERT_EQ(i + dim[0]*( j + dim[1]*k ), l);
    }

}

TEST(MeshgridTest, ReciprocalMesh) {

    vector<int> dim{80,60,100};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,0.,0};
    lattice[1] = vector<double>{0,3.,0};
    lattice[2] = vector<double>{0,0,5.};
    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);
    //rec_mesh.setupCoordinates();
    int N  = real_mesh.getNumDataPoints();

    //ASSERT_EQ(N, 81*61*101);
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    double  val = 0.;
    vector<double> vec{0,0,0};
    vector<double> rec_vec{0,0,0};
    double x,y,z;
    for (int i=0; i< real_mesh.getDim()[0]; i++) {
        rec_mesh.xyz(i,x,y,z);
        for (int xj=-10; xj < 10; xj++) {
            for (int yj=-10; yj < 10; yj++) {
                for (int zj=-10; zj < 10; zj++) {
                    vec = matVecMul3x3(transpose3x3(lattice), vector<double>{double(xj),double(yj),double(zj)});
                    rec_vec[0] = x;
                    rec_vec[1] = y;
                    rec_vec[2] = z;

                    val = dotProduct(vec, rec_vec)  / (2*M_PI);
                    ASSERT_NEAR( val, round(val), 1e-10);  // check if integer
                }
            }
        }

    }

}


TEST(MeshgridTest, ReciprocalMesh2) {

    vector<int> dim{85,80,40};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,0.,2};
    lattice[1] = vector<double>{1,3.,0};
    lattice[2] = vector<double>{0.3,0,5.};
    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);
    //rec_mesh.setupCoordinates();
    int N  = real_mesh.getNumDataPoints();
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    double  val = 0.;
    vector<double> vec{0,0,0};
    vector<double> rec_vec{0,0,0};
    double x,y,z;
    for (int i=0; i< real_mesh.getDim()[0]; i++) {
        rec_mesh.xyz(i,x,y,z);
        for (int xj=-10; xj < 10; xj++) {
            for (int yj=-10; yj < 10; yj++) {
                for (int zj=-10; zj < 10; zj++) {
                    vec = matVecMul3x3(transpose3x3(lattice), vector<double>{double(xj),double(yj),double(zj)});
                    rec_vec[0] = x;
                    rec_vec[1] = y;
                    rec_vec[2] = z;

                    val = dotProduct(vec, rec_vec)  / (2*M_PI);
                    ASSERT_NEAR( val, round(val), 1e-10);  // check if integer
                }
            }
        }

    }

}


TEST(MeshgridTest, ReciprocalLattice) {

    vector<int> dim{21,5,3};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,-1.,0};
    lattice[1] = vector<double>{1.,3.,0};
    lattice[2] = vector<double>{0,0.3,5.};

    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);

    int N  = real_mesh.getNumDataPoints();
    ASSERT_EQ(N, 21*5*3);
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    vector< vector<double>> real_lattice = real_mesh.getLattice();
    vector< vector<double>> rec_lattice = rec_mesh.getLattice();

    vector< vector<double>> product = matMul3x3(real_lattice,rec_lattice);

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // ASSERT_NEAR(dotProduct(real_lattice[i], rec_lattice[j]), int(i==j)*2*M_PI, 1e-10);
            ASSERT_NEAR(product[i][j], int(i==j)*2*M_PI, 1e-10);
        }
    }

    product = matMul3x3(rec_lattice, real_lattice);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // ASSERT_NEAR(dotProduct(real_lattice[i], rec_lattice[j]), int(i==j)*2*M_PI, 1e-10);
            ASSERT_NEAR(product[i][j], int(i==j)*2*M_PI, 1e-10);
        }
    }
}


TEST(MeshgridTest, ReciprocalLattice2) {

    vector<int> dim{21,5,3};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    double L = 2.3;

    lattice[0] = vector<double>{L,0.,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0.,L};

    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);

    int N  = real_mesh.getNumDataPoints();
    ASSERT_EQ(N, 21*5*3);
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    vector< vector<double>> real_lattice = real_mesh.getLattice();
    vector< vector<double>> rec_lattice = rec_mesh.getLattice();

    printMat3x3(rec_lattice);

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            if (i==j)
                ASSERT_NEAR(rec_lattice[i][j], 2*M_PI/L, 1e-10);
            else
                ASSERT_NEAR(rec_lattice[i][j], 0.0, 1e-10);
        }
    }
}


TEST(MeshgridTest, ReciprocalLattice3) {

    vector<int> dim{21,5,3};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    double L = 2.3;

    lattice[0] = vector<double>{L,0.,0};
    lattice[1] = vector<double>{-L/2.,L*sqrt(3.0)/2.,0};
    lattice[2] = vector<double>{0,0.,L};

    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);

    int N  = real_mesh.getNumDataPoints();
    ASSERT_EQ(N, 21*5*3);
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    vector< vector<double>> real_lattice = real_mesh.getLattice();
    vector< vector<double>> rec_lattice = rec_mesh.getLattice();


    vector< vector<double> > shouldbe(3);
    shouldbe[0] = vector<double>{2*M_PI/L,0.,0};
    shouldbe[1] = vector<double>{2*M_PI/(sqrt(3.)*L),4*M_PI/(sqrt(3.)*L),0};
    shouldbe[2] = vector<double>{0,0.,2*M_PI/L};

    printMat3x3(shouldbe);
    printMat3x3(rec_lattice);

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            ASSERT_NEAR(rec_lattice[i][j], shouldbe[i][j], 1e-10);
        }
    }
}


TEST(MeshgridTest, ReciprocalLatticeOrder) {

    vector<int> dim{20,20,20};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,0.,0};
    lattice[1] = vector<double>{0.,3.,0};
    lattice[2] = vector<double>{0,0.0,5.};

    // lattice[0] = vector<double>{16.4857,0.,0};
    // lattice[1] = vector<double>{-8.24287,14.2771,0};
    // lattice[2] = vector<double>{0,0.0,20.0};

    RealMeshgrid real_mesh(dim, lattice, origin);
    ReciprocalMeshgrid rec_mesh(&real_mesh);

    int N  = real_mesh.getNumDataPoints();
    ASSERT_EQ(N, rec_mesh.getNumDataPoints());

    double x0,y0,z0;
    rec_mesh.xyz(0,x0,y0,z0);
    ASSERT_NEAR(x0,0,1e-10);
    ASSERT_NEAR(y0,0,1e-10);
    ASSERT_NEAR(z0,0,1e-10);

    double x,y,z;
    int i,j,k;
    for (int n=1; n<N; n++) {
        rec_mesh.getIndexTipple(n, i,j,k);
        rec_mesh.xyz(n,x,y,z);

        // cout << i << " " << j << " " << k << endl;
        // cout << x << " " << y << " " << z << endl;

        if (i != rec_mesh.getDim()[0]/2+1) {
            ASSERT_GT(x-x0, -1e-8);
        }

        if (j != rec_mesh.getDim()[1]/2+1) {
            ASSERT_GT(y-y0, -1e-8);
        }
        
        if (k != rec_mesh.getDim()[1]/2+1) {
            ASSERT_GT(z-z0, -1e-8);
        }

        x0 = x;
        y0 = y;
        z0 = z;
    }
}




TEST(MeshgridTest, Fourier) {

    double L =10., alpha=1.1;
    vector<int> dim{61,61,61};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0.,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    RealMeshgrid* real_mesh = new RealMeshgrid(dim, lattice, origin);
    ReciprocalMeshgrid* rec_mesh = new ReciprocalMeshgrid(real_mesh);
    //rec_mesh->setupCoordinates();

    int N  = real_mesh->getNumDataPoints();
    ASSERT_EQ(N, 61*61*61);
    ASSERT_EQ(N, rec_mesh->getNumDataPoints());

    fftw_complex* f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p_f1 = fftw_plan_dft_3d(real_mesh->getDim()[2], real_mesh->getDim()[1], real_mesh->getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE);


    double x,y,z;
    for (int i=0; i<N; i++){
        real_mesh->xyz(i,x,y,z);  // get coordinates

        // real part
        f1[i][0] = gaussian3D(x,y,z, L/2, L/2, L/2, alpha, 1.0);

        // imaginary part
        f1[i][1] = 0.0;
        }

    fftw_execute(p_f1);
    fftw_complex gauss;
    double dV = real_mesh->getdV();
    for (int i=0; i<N; i++){
        rec_mesh->xyz(i,x,y,z);
        recGaussian3D(x,y,z, L/2, L/2, L/2, alpha, 1.0, gauss);
        ASSERT_NEAR(f1[i][0], gauss[0]/dV, 1e-8);
        ASSERT_NEAR(f1[i][1], gauss[1]/dV, 1e-8);
    }

}

TEST(MeshgridTest, FourierEven) {

    double L =10., alpha=1.1;
    vector<int> dim{60,60,60};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{L,0.,0};
    lattice[1] = vector<double>{0,L,0};
    lattice[2] = vector<double>{0,0,L};
    RealMeshgrid* real_mesh = new RealMeshgrid(dim, lattice, origin);
    ReciprocalMeshgrid* rec_mesh = new ReciprocalMeshgrid(real_mesh);
    //rec_mesh->setupCoordinates();

    int N  = real_mesh->getNumDataPoints();
    ASSERT_EQ(N, 60*60*60);
    ASSERT_EQ(N, rec_mesh->getNumDataPoints());

    fftw_complex* f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p_f1 = fftw_plan_dft_3d(real_mesh->getDim()[2], real_mesh->getDim()[1], real_mesh->getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE);


    double x,y,z;
    for (int i=0; i<N; i++){
        real_mesh->xyz(i,x,y,z);  // get coordinates

        // real part
        f1[i][0] = gaussian3D(x,y,z, L/2, L/2, L/2, alpha, 1.0);

        // imaginary part
        f1[i][1] = 0.0;
        }

    fftw_execute(p_f1);
    fftw_complex gauss;
    double dV = real_mesh->getdV();
    for (int i=0; i<N; i++){
        rec_mesh->xyz(i,x,y,z);
        recGaussian3D(x,y,z, L/2, L/2, L/2, alpha, 1.0, gauss);
        ASSERT_NEAR(f1[i][0], gauss[0]/dV, 1e-8);
        ASSERT_NEAR(f1[i][1], gauss[1]/dV, 1e-8);
    }

}


TEST(MeshgridTest, Fourier2) {

    double alpha=1.0;
    vector<int> dim{81,81,81};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{16.4,0.,0};
    lattice[1] = vector<double>{-8.3,10.4,0};
    lattice[2] = vector<double>{0,0,10};
    RealMeshgrid* real_mesh = new RealMeshgrid(dim, lattice, origin);
    ReciprocalMeshgrid* rec_mesh = new ReciprocalMeshgrid(real_mesh);
    //rec_mesh->setupCoordinates();

    int N  = real_mesh->getNumDataPoints();
    ASSERT_EQ(N, 81*81*81);
    ASSERT_EQ(N, rec_mesh->getNumDataPoints());

    fftw_complex* f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p_f1 = fftw_plan_dft_3d(real_mesh->getDim()[2], real_mesh->getDim()[1], real_mesh->getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE);


    double x,y,z, val=0.0;
    for (int i=0; i<N; i++){
        real_mesh->xyz(i,x,y,z);  // get coordinates

        // real part
        f1[i][0] = gaussian3D(x,y,z, 5.5, 4.5, 5., alpha, 1.0);
        val += f1[i][0];

        // imaginary part
        f1[i][1] = 0.0;
    }
    fftw_execute(p_f1);

    // XSF_controller::save_file("test2.xsf",rec_mesh,f1,rec_mesh->getLattice());
    // fftw_complex* f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex gauss;
    double dV = real_mesh->getdV();
    for (int i=0; i<N; i++){
        rec_mesh->xyz(i,x,y,z);
        recGaussian3D(x,y,z, 5.5, 4.5, 5., alpha, 1.0, gauss);
        // cout << i << endl;
        ASSERT_NEAR(f1[i][0], gauss[0]/dV, 1e-5);
        ASSERT_NEAR(f1[i][1], gauss[1]/dV, 1e-5);
        // f2[i][0] = gauss[0];
        // f2[i][1] = gauss[1];
    }

    // XSF_controller::save_file("test3.xsf",rec_mesh,f1,rec_mesh->getLattice());
}



TEST(MeshgridTest, Fourier2Even) {

    double alpha=1.0;
    vector<int> dim{80,80,80};
    vector<double> origin{0,0,0};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{16.4,0.,0};
    lattice[1] = vector<double>{8.3,10.4,0};
    lattice[2] = vector<double>{0,0,10};
    RealMeshgrid* real_mesh = new RealMeshgrid(dim, lattice, origin);
    ReciprocalMeshgrid* rec_mesh = new ReciprocalMeshgrid(real_mesh);
    //rec_mesh->setupCoordinates();

    int N  = real_mesh->getNumDataPoints();
    ASSERT_EQ(N, 80*80*80);
    ASSERT_EQ(N, rec_mesh->getNumDataPoints());

    fftw_complex* f1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p_f1 = fftw_plan_dft_3d(real_mesh->getDim()[2], real_mesh->getDim()[1], real_mesh->getDim()[0], f1, f1, FFTW_FORWARD, FFTW_ESTIMATE);


    double x,y,z;
    for (int i=0; i<N; i++){
        real_mesh->xyz(i,x,y,z);  // get coordinates

        // real part
        f1[i][0] = gaussian3D(x,y,z, 11., 4, 5., alpha, 1.0);

        // imaginary part
        f1[i][1] = 0.0;
    }

    // XSF_controller::save_file("test.xsf",real_mesh,f1,real_mesh->getLattice());

    fftw_execute(p_f1);

    // XSF_controller::save_file("test2.xsf",rec_mesh,f1,rec_mesh->getLattice());
    // fftw_complex* f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex gauss;
    double dV = real_mesh->getdV();
    for (int i=0; i<N; i++){
        rec_mesh->xyz(i,x,y,z);
        recGaussian3D(x,y,z, 11., 4., 5., alpha, 1.0, gauss);
        ASSERT_NEAR(f1[i][0], gauss[0]/dV, 1e-5);
        ASSERT_NEAR(f1[i][1], gauss[1]/dV, 1e-5);
        // f2[i][0] = gauss[0];
        // f2[i][1] = gauss[1];
    }

    // XSF_controller::save_file("test3.xsf",rec_mesh,f1,rec_mesh->getLattice());
}



TEST(MeshgridTest, getShiftedTiple) {  // shifts by lattice vectors

    vector<int> dim{19,5,14};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,-1.,0};
    lattice[1] = vector<double>{1.,3.,0};
    lattice[2] = vector<double>{0,0.3,5.};

    RealMeshgrid real_mesh(dim, lattice, origin);
    vector<double> R_cart{0,0,0};
    vector<double> R_test{0,0,0};
    int x1,y1,z1,x2,y2,z2;
    bool ret;
    int N = real_mesh.getNumDataPoints();

    for (int i=-3; i<=3; i++) {
        for (int j=-3; j<=3; j++) {
            for (int k=-3; k<=3; k++) {
                R_cart = matVecMul3x3(transpose3x3(lattice), vector<int>{i,j,k});
                // cout << "R_cart = "<< R_cart[0] << " " << R_cart[1] << " " << R_cart[2] << endl;

                for (int n=0; n<N; n++) {
                    ret = real_mesh.getShiftedIndexTipple(n, R_cart, x1,y1,z1);
                    real_mesh.getIndexTipple(n,x2,y2,z2);

                    ASSERT_EQ(ret, (x1>=0) && (x1<dim[0]) && (y1>=0) && (y1<dim[1]) && (z1>=0) && (z1<dim[2]));

                    // cout << "n = " << n << " x1 = "<< x1 << " " << x1 << " " << x1 << endl;
                    // cout << "n = " << n << " x2 = "<< x2 << " " << x2 << " " << x2 << endl;

                    real_mesh.xyz(x1-x2, y1-y2, z1-z2, R_test[0], R_test[1],R_test[2]);

                    // cout << "R_test = "<< R_test[0] << " " << R_test[1] << " " << R_test[2] << endl;

                    for (int l=0; l<3; l++) {
                        ASSERT_NEAR(R_test[l]-origin[l], R_cart[l], 1e-10);
                    }


                }
            }
        }
    }
}



TEST(MeshgridTest, getShiftedTiple2) {  // shift by arbitary vectors

    vector<int> dim{6,12,14};
    vector<double> origin{-1,2,4};
    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{2.,-1.,0};
    lattice[1] = vector<double>{0,3.,1.};
    lattice[2] = vector<double>{0,0.3,5.};

    RealMeshgrid real_mesh(dim, lattice, origin);
    vector<double> R_cart{0,0,0};
    vector<double> R_test{0,0,0};
    int x1,y1,z1,x2,y2,z2;
    bool ret;
    int N = real_mesh.getNumDataPoints();


    real_mesh.xyz(N/2+1, R_cart[0], R_cart[1], R_cart[2]);  // arbitrary vector
    R_cart[0] -= origin[0];
    R_cart[1] -= origin[1];
    R_cart[2] -= origin[2];
    // cout << "R_cart = "<< R_cart[0] << " " << R_cart[1] << " " << R_cart[2] << endl;

    for (int n=0; n<N; n++) {
        ret = real_mesh.getShiftedIndexTipple(n, R_cart, x1,y1,z1);
        real_mesh.getIndexTipple(n,x2,y2,z2);

        ASSERT_EQ(ret, (x1>=0) && (x1<dim[0]) && (y1>=0) && (y1<dim[1]) && (z1>=0) && (z1<dim[2]));

        // cout << "n = " << n << " x1 = "<< x1 << " " << x1 << " " << x1 << endl;
        // cout << "n = " << n << " x2 = "<< x2 << " " << x2 << " " << x2 << endl;

        real_mesh.xyz(x1-x2, y1-y2, z1-z2, R_test[0], R_test[1],R_test[2]);

        // cout << "R_test = "<< R_test[0] << " " << R_test[1] << " " << R_test[2] << endl;

        for (int l=0; l<3; l++) {
            ASSERT_NEAR(R_test[l]-origin[l], R_cart[l], 1e-10);
        }


    }

}


int getShiftedGlobalIndex_old(const int n, const vector<int> R, RealMeshgrid* meshgrid, vector<vector<double>> unitcell)  // old routine only for test purposes
    {

        vector<double> R_cart{0, 0, 0};
        // vector<vector<double>> A = transpose3x3(unitcell);

        for (int i = 0; i < 3; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                R_cart[i] += -unitcell[k][i] * R[k];
                // minus because we want to get the values at n for the shifted WF
                // and not the value at n+R
            }
        }

        vector<int> ijk{0, 0, 0};
        if (meshgrid->getShiftedIndexTipple(n, R_cart, ijk[0], ijk[1], ijk[2]))
        {
            return meshgrid->getGlobId(ijk[0], ijk[1], ijk[2]);
        }
        else
        {
            return -1; // requested values outside the data grid
        }
    }

TEST(MeshgridTest, getShiftedGlobIndex) {

    // geometry of Si
    vector<int> dim{88, 88, 88};
    vector<double> origin{-27.376250000000, -27.376250000000, -27.376250000000};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{29.865000000,   29.865000000,   0.0000000000};
    lattice[1] = vector<double>{0.0000000000,   29.865000000,   29.865000000};
    lattice[2] = vector<double>{29.865000000,   0.0000000000,   29.865000000};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{2.715000000000,   2.715000000000,   0.000000000000};
    unitcell[1] = vector<double>{0.000000000000,   2.715000000000,   2.715000000000};
    unitcell[2] = vector<double>{2.715000000000,   0.000000000000,   2.715000000000};

    RealMeshgrid* meshgrid1 = new RealMeshgrid(dim, lattice, origin);
    int N = meshgrid1->getNumDataPoints();

    vector<double> supercell{11.0,11.0,11.0};

    int newGlobId, oldGlobId, m,n,l;
    double x,y,z;

    vector<int> R{5,0,10};
    vector<int> indexOffset = meshgrid1->getIndexOffset(R, supercell);
;
    //cout << "offset: " << indexOffset[0] << " " << indexOffset[1] << " " << indexOffset[2] << "\n";

    for (int i=0; i<N; i++) {
        newGlobId = meshgrid1->getShiftedGlobalIndex(i, indexOffset);
        oldGlobId = getShiftedGlobalIndex_old(i, R, meshgrid1, unitcell);
        

        if (newGlobId != oldGlobId) {
            meshgrid1->getIndexTipple(oldGlobId, m,n,l);
            cout << i << " offset (old): " << m << " " << n << " " << l << "\n";
            meshgrid1->getIndexTipple(newGlobId, m,n,l);
            cout << i << " offset (new): " << m << " " << n << " " << l << "\n";
            meshgrid1->xyz(oldGlobId,x,y,z);
            cout << "position (old): " << x << " " << y << " " << z << "\n";
            meshgrid1->xyz(newGlobId,x,y,z);
            cout << "position (new): " << x << " " << y << " " << z << "\n";
            cout << "\n";
        }

        ASSERT_EQ(newGlobId, oldGlobId);
    }

    R = vector<int>{0,-3,-4};
    indexOffset = meshgrid1->getIndexOffset(R, supercell);
    //cout << "offset: " << indexOffset[0] << " " << indexOffset[1] << " " << indexOffset[2] << "\n";

    for (int i=0; i<N; i++) {
        newGlobId = meshgrid1->getShiftedGlobalIndex(i, indexOffset);
        oldGlobId = getShiftedGlobalIndex_old(i, R, meshgrid1, unitcell);
        

        if (newGlobId != oldGlobId) {
            meshgrid1->getIndexTipple(oldGlobId, m,n,l);
            cout << i << " offset (old): " << m << " " << n << " " << l << "\n";
            meshgrid1->getIndexTipple(newGlobId, m,n,l);
            cout << i << " offset (new): " << m << " " << n << " " << l << "\n";
            meshgrid1->xyz(oldGlobId,x,y,z);
            cout << "position (old): " << x << " " << y << " " << z << "\n";
            meshgrid1->xyz(newGlobId,x,y,z);
            cout << "position (new): " << x << " " << y << " " << z << "\n";
            cout << "\n";
        }

        ASSERT_EQ(newGlobId, oldGlobId);
    }
}


TEST(MeshgridTest, getShiftedGlobIndex2) {

    // geometry of stacked COF-CN-1Ph
    vector<int> dim{108 ,108, 56};
    vector<double> origin{-7.580743000000, -17.784735000000, -10.124349000000};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{130.716095465181, 0.000000000000, 0.000000000000};
	lattice[1] = vector<double>{-65.364281047354, 113.199859754875, 0.000000000000};
	lattice[2] = vector<double>{0.597506669856, 11.774209966507, 23.363883186603};

    vector< vector<double> > unitcell(3);
    unitcell[0] = vector<double>{14.524010600000, 0.000000000000, 0.000000000000};
	unitcell[1] = vector<double>{-7.262697900000, 12.577762200000, 0.000000000000};
	unitcell[2] = vector<double>{0.085358100000, 1.682030000000, 3.337697600000};

    RealMeshgrid* meshgrid1 = new RealMeshgrid(dim, lattice, origin);
    int N = meshgrid1->getNumDataPoints();

    vector<double> supercell{9.0,9.0,7.0};
    int newGlobId, oldGlobId, m,n,l;
    double x,y,z;


    vector<int> R{-2,5,1};
    vector<int> indexOffset = meshgrid1->getIndexOffset(R, supercell);
    //cout << "offset: " << indexOffset[0] << " " << indexOffset[1] << " " << indexOffset[2] << "\n";

    for (int i=0; i<N; i++) {
        newGlobId = meshgrid1->getShiftedGlobalIndex(i, indexOffset);
        oldGlobId = getShiftedGlobalIndex_old(i, R, meshgrid1, unitcell);
        

        if (newGlobId != oldGlobId) {
            meshgrid1->getIndexTipple(oldGlobId, m,n,l);
            cout << i << " offset (old): " << m << " " << n << " " << l << "\n";
            meshgrid1->getIndexTipple(newGlobId, m,n,l);
            cout << i << " offset (new): " << m << " " << n << " " << l << "\n";
            meshgrid1->xyz(oldGlobId,x,y,z);
            cout << "position (old): " << x << " " << y << " " << z << "\n";
            meshgrid1->xyz(newGlobId,x,y,z);
            cout << "position (new): " << x << " " << y << " " << z << "\n";
            cout << "\n";
        }

        ASSERT_EQ(newGlobId, oldGlobId);
    }


    R = vector<int>{3,0,-2};
    indexOffset = meshgrid1->getIndexOffset(R, supercell);
    //cout << "offset: " << indexOffset[0] << " " << indexOffset[1] << " " << indexOffset[2] << "\n";

    for (int i=0; i<N; i++) {
        newGlobId = meshgrid1->getShiftedGlobalIndex(i, indexOffset);
        oldGlobId = getShiftedGlobalIndex_old(i, R, meshgrid1, unitcell);
        

        if (newGlobId != oldGlobId) {
            meshgrid1->getIndexTipple(oldGlobId, m,n,l);
            cout << i << " offset (old): " << m << " " << n << " " << l << "\n";
            meshgrid1->getIndexTipple(newGlobId, m,n,l);
            cout << i << " offset (new): " << m << " " << n << " " << l << "\n";
            meshgrid1->xyz(oldGlobId,x,y,z);
            cout << "position (old): " << x << " " << y << " " << z << "\n";
            meshgrid1->xyz(newGlobId,x,y,z);
            cout << "position (new): " << x << " " << y << " " << z << "\n";
            cout << "\n";
        }

        ASSERT_EQ(newGlobId, oldGlobId);
    }
}

TEST(MeshgridTest, estimateIndexFromPosition)
{

        // geometry of stacked COF-CN-1Ph
    vector<int> dim{18 ,27, 16};
    vector<double> origin{-7.580743000000, -17.784735000000, -10.124349000000};

    vector< vector<double> > lattice(3);
    lattice[0] = vector<double>{130.716095465181, 0.000000000000, 0.000000000000};
	lattice[1] = vector<double>{-65.364281047354, 113.199859754875, 0.000000000000};
	lattice[2] = vector<double>{0.597506669856, 11.774209966507, 23.363883186603};

    RealMeshgrid* mesh = new RealMeshgrid(dim, lattice, origin);

    int N = mesh->getNumDataPoints();
    int io,jo,ko, i,j,k;
    double x,y,z;
    bool ret;

    for (int l=0; l<N; l++) {
        mesh->xyz(l, x,y,z);
        mesh->getIndexTipple(l, io,jo,ko);
        ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);

        ASSERT_TRUE(ret);
        ASSERT_EQ(i, io);
        ASSERT_EQ(j, jo);
        ASSERT_EQ(k, ko);
    }

    mesh->xyz(-1,0,0, x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);

    mesh->xyz(0,-1,0, x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);

    mesh->xyz(0,0,-1, x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);


    mesh->xyz(dim[0],0,0, x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);

    mesh->xyz(0,dim[1],0, x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);

    mesh->xyz(0,0,dim[2], x,y,z);
    ret = mesh->estimateIndexFromPosition(x,y,z, i,j,k);
    ASSERT_FALSE(ret);

}