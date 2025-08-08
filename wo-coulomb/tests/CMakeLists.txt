cmake_minimum_required (VERSION 3.16)

# parallelization
find_package(OpenMP 4.0 REQUIRED)
find_package(MPI REQUIRED)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests.x tests.cpp)
target_include_directories(runTests.x PUBLIC ..)
target_link_libraries(runTests.x PUBLIC ${GTEST_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX ${FFTW_THREADS_LIB} ${FFTW_LIB})
target_link_libraries(runTests.x )

# libaries
find_library(FFTW_LIB fftw3 HINTS /dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/haswell/fftw/3.3.10-intel-pbsjdkm/lib/)
find_library(FFTW_THREADS_LIB fftw3_threads HINTS /dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/haswell/fftw/3.3.10-intel-pbsjdkm/lib/)
target_include_directories(${PROJECT_NAME}.x PUBLIC /dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/haswell/fftw/3.3.10-intel-pbsjdkm/include)

add_test(
    NAME runTests.x
    COMMAND runTests.x
)
