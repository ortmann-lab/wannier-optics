# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.16)

# Define the main project
project(wannier-optics)

# Add subdirectories for each subproject
add_subdirectory(wo-coulomb)
add_subdirectory(wo-optics)


# Specify installation rules
# install(TARGETS coulomb-integral/CoulombIntegral.x coulomb-integral/OpticalDipole.x optics/WannierOptics.x RUNTIME DESTINATION bin)
# install(TARGETS CoulombIntegral.x RUNTIME DESTINATION bin)


# Specify package information
set(CPACK_PACKAGE_NAME "wannier-optics")
set(CPACK_PACKAGE_VERSION "0.0.1")
set(CPACK_PACKAGE_CONTACT "konrad.merkel@tum.de")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Konrad Merkel <${CPACK_PACKAGE_CONTACT}>")

set(CPACK_PACKAGE_DESCRIPTION "Calculation of optical absorption spectra using the linear scaling Wannier optics approach.")
set(CPACK_GENERATOR "DEB")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.md")
#set(CPACK_DEBIAN_PACKAGE_COPYRIGHT "${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT")
#set(CPACK_DEBIAN_PACKAGE_CHANGELOG "${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG")
set(CPACK_STRIP_FILES YES)
set(
    CPACK_INSTALL_DEFAULT_DIRECTORY_PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)


# Set package dependencies
set(CPACK_DEBIAN_PACKAGE_DEPENDS "mpi-default-dev (>= 1.13), mpi-default-bin (>= 1.13), libfftw3-dev (>=3.3.8), libfftw3-bin (>=3.3.8), libc6 (>= 2.29), libgcc-s1 (>= 3.0), libgomp1 (>= 4.9), libstdc++6 (>= 9), libgfortran5 (>= 8)")
#set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS YES)

# Include CPack for packaging
include(CPack)