include(CMakePrintHelpers)

cmake_minimum_required(VERSION 3.15.0)
project (optimizer LANGUAGES Fortran C)

# set paths to configure projects
set(MAINLIB periodicIO)
set(fingerprintlib fingerprint)

# Define flags for intel and gnu compilers. If the compiler is not known, CMake will choose flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -no-wrap-margin")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -check bounds -fpe0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -ip -xHOST")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -fcheck=all -Wall -ffpe-trap=zero,overflow,underflow,invalid -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native -mtune=native")
endif()

# Configure target runnerlib. This is the main library of the runner project.
add_library(${MAINLIB})
add_subdirectory("dev/src/read")
add_subdirectory("dev/src/write")
add_subdirectory("dev/src/periodic_operations")
target_sources(${MAINLIB} PRIVATE ${SOURCES})

add_subdirectory("dev/src/OM-fingerprint")
add_subdirectory("dev/src/transform")
add_subdirectory("dev/src/fingerprint-programs")

find_package(OpenMP)
  if(NOT OpenMP_Fortran_FOUND)
  message("Your Fortran compiler does not support openmp.")
  endif()
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${MAINLIB} PUBLIC OpenMP::OpenMP_Fortran)
endif()
  
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(${fingerprintlib} PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Set Default install directory to $HOME/bin and $HOME/lib
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX $ENV{HOME}
           CACHE PATH "default install path" FORCE)
endif()
install(TARGETS ${MAINLIB} DESTINATION lib)