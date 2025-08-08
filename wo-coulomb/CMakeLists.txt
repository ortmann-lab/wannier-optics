cmake_minimum_required(VERSION 3.16)

project (wo-coulomb LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-Wall -O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

#project(CoulombIntegral)
add_executable(${PROJECT_NAME}.x main.cpp)

# parallelization
find_package(OpenMP 4.0 REQUIRED)
find_package(MPI REQUIRED)

# libaries and targets
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW IMPORTED_TARGET REQUIRED fftw3)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

target_include_directories(${PROJECT_NAME}.x PUBLIC PkgConfig::FFTW)
target_link_libraries(${PROJECT_NAME}.x PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::FFTW)

# tell the program its git commit id
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

# Try to find GoogleTest
find_package(GTest QUIET)

# Check if GoogleTest was found
if (GTest_FOUND)
    message(STATUS "GoogleTest found. Enabling tests.")
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
else()
    message(WARNING "GoogleTest not found. Tests will not be built.")
endif()


install(TARGETS ${PROJECT_NAME}.x RUNTIME DESTINATION bin)
