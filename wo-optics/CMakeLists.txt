cmake_minimum_required(VERSION 3.16)

project(wo-optics Fortran)

set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED True)

# Set the Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -O3 -cpp -Wunused")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
#add_compile_options(-fallow-argument-mismatch)

# parallelization
find_package(MPI REQUIRED)
find_package(LAPACK REQUIRED)

# Define the source files
set(FORTRAN_SOURCES
    external/cfgio/src/string_conv_mod.f90
    external/cfgio/src/cfgio_mod.f90
    seed_data.f90
    geometry.f90
    test_geometry_mod.f90
    parfile.f90
    potentials.f90
    hamiltonian.f90
    sparseCSR.f90
    algorithms.f90
)
# Gather all Fortran source files in the src directory
# file(GLOB_RECURSE FORTRAN_SOURCES src/*.f90)

# Create the main executable
add_executable(${PROJECT_NAME}.x ${FORTRAN_SOURCES} main.f90 )
target_link_libraries(${PROJECT_NAME}.x PUBLIC ${MPI_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME}.x PRIVATE ${LAPACK_LIBRARIES})
target_include_directories(${PROJECT_NAME}.x PRIVATE ${MPI_Fortran_INCLUDE_DIRS})

# add_executable(test_parfile.x ${FORTRAN_SOURCES} test_parfile.f90 )
# target_link_libraries(test_parfile.x PUBLIC ${MPI_Fortran_LIBRARIES})
# target_include_directories(test_parfile.x PRIVATE ${MPI_Fortran_INCLUDE_DIRS})

# add_executable(test_geometry.x ${FORTRAN_SOURCES} test_geometry.f90 )
# target_link_libraries(test_geometry.x PUBLIC ${MPI_Fortran_LIBRARIES})
# target_include_directories(test_geometry.x PRIVATE ${MPI_Fortran_INCLUDE_DIRS})



# Specify installation rules
install(TARGETS ${PROJECT_NAME}.x RUNTIME DESTINATION bin)
#install(DIRECTORY seed DESTINATION share/wannier-optics/)
