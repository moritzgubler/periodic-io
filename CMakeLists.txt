
include(CMakePrintHelpers)

cmake_minimum_required(VERSION 3.15.0)
project (optimizer LANGUAGES Fortran C)

# set paths to configure projects
set(EXEC optimizer.x)
set(MAINLIB optimizer)
set(MAINFILE src/main/main.f90)

# find all source files and clean them up
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/fortran/*.f90)

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
target_sources(${MAINLIB} PRIVATE ${SOURCES})
set_target_properties (${MAINLIB} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(${MAINLIB} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)


find_package(OpenMP)
  if(NOT OpenMP_Fortran_FOUND)
  message("Your Fortran compiler does not support openmp.")
  endif()
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${MAINLIB} PUBLIC OpenMP::OpenMP_Fortran)
endif()

  
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(${MAINLIB} PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# if(NOT DEFINED USE_MPI)
#   set(USE_MPI true)
# endif()
# if(USE_MPI)
#   find_package(MPI REQUIRED COMPONENTS Fortran)
#   target_link_libraries(${MAINLIB} PUBLIC MPI::MPI_Fortran)
# endif()

add_executable(${EXEC} ${MAINFILE})
target_link_libraries(${EXEC} ${MAINLIB})

# Set Default install directory to $HOME/bin and $HOME/lib
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX $ENV{HOME}
           CACHE PATH "default install path" FORCE)
endif()
install(TARGETS ${EXEC} DESTINATION bin)
install(TARGETS ${MAINLIB} DESTINATION lib)