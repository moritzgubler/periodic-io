################################################################################
# Makefile of molecularIO library. Feel free configure the variables below
# to your needs.
# Use "make debug" to debug this makefile
################################################################################

# compiler settings
fc := ifort
linking_flags = -g -qopenmp -O2 -no-wrap-margin
compile_flags = $(linking_flags)
lapack = -mkl

# Name of library file
lib_name := libperiodicIO.a

# install directory for executables:
install_exec_dir := $(HOME)/bin

# library install path:
install_lib_dir := $(HOME)/.local/lib

# library settings
lib_command := ar cr
