################################################################################
# Makefile of molecularIO library. Feel free configure the variables below
# to your needs.
# Use "make debug" to debug this makefile
################################################################################

# compiler settings
fc := gfortran
linking_flags = -g -fcheck=all -Wall
compile_flags = $(linking_flags)

# Name of library file
lib_name := libmolecularIO.a

# install directory for executables:
install_exec_dir := $(HOME)/bin

# library install path:
install_lib_dir := $(HOME)/.local/lib

# library settings
lib_command := ar cr

################################################################################
# Don't change anything below this line unless you know what you are doing.
# You should be able to adjust enough by tweaking the paths above.
################################################################################

# path definitions
# location of source files:
src_read_dir := dev/src/read
src_write_dir := dev/src/write
src_transform_dir := dev/src/transform
src_periodic_dir := dev/src/periodic_operations

# location of binaries:
bin_dir := dev/bin
prog_dir := executables


# find all relevant source files.
src_read := $(wildcard $(src_read_dir)/*.f90)
src_write := $(wildcard $(src_write_dir)/*.f90)
src_transform := $(wildcard $(src_transform_dir)/*.f90)
src_periodic := $(wildcard $(src_periodic_dir)/*.f90)

# set location of object files:
obj_read := $(patsubst $(src_read_dir)/%.f90, $(bin_dir)/%.o, $(src_read))
obj_write := $(patsubst $(src_write_dir)/%.f90, $(bin_dir)/%.o, $(src_write))
obj_transform := $(patsubst $(src_transform_dir)/%.f90, $(bin_dir)/%.o, $(src_transform))
obj_periodic := $(patsubst $(src_periodic_dir)/%.f90, $(bin_dir)/%.o, $(src_periodic))

lib_file := $(bin_dir)/$(lib_name)

# set program names (all files from $(src_transform_dir) should be executable)
programs := $(patsubst $(src_transform_dir)/%.f90, $(prog_dir)/%, $(src_transform))

# set phony targets
.PHONY: default debug clean install_bin install_lib install test doc

default: $(programs) $(lib_file)

doc:
	(cd dev && ford doc.md)

# call make debug to debug the following variables.
debug:
	@echo "src_read = $(src_read)"
	@echo "src_write = $(src_write)"
	@echo "src_transform = $(src_transform)"
	@echo "src_periodic = $(src_periodic)"
	@echo "obj_read = $(obj_read)"
	@echo "obj_write = $(obj_write)"
	@echo "obj_transform = $(obj_transform)"
	@echo "obj_periodic = $(obj_periodic)"
	@echo "lib_file = $(lib_file)"
	@echo "programs = $(programs)"

test: $(programs) $(lib_file)
	${MAKE} -C dev/test/test_src_fingerprint
	(cd dev/test && ./test_all.sh)

clean:
	rm -rf $(bin_dir) $(prog_dir)
	$(MAKE) -C dev/test/test_src_fingerprint clean

install: install_lib install_bin

install_bin: $(programs) | $(install_exec_dir)
	cp $(prog_dir)/* $(install_exec_dir)/.

install_lib: $(lib_file) | $(install_lib_dir)
	cp $(lib_file) $(install_lib_dir)/.

#create library:
$(lib_file): $(obj_periodic) $(obj_read) $(obj_write)
	$(lib_command) $(lib_file) $(obj_periodic) $(obj_read) $(obj_write)

# compile src_read_dir
$(obj_read): $(bin_dir)/%.o: $(src_read_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

# compile src_write_dir
$(obj_write): $(bin_dir)/%.o: $(src_write_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

# compile src_transform_dir
$(obj_transform): $(bin_dir)/%.o: $(src_transform_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

# compile src_periodic_dir
$(obj_periodic): $(bin_dir)/%.o: $(src_periodic_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

#link all transformer programs:
$(patsubst $(src_transform_dir)/%.f90, $(prog_dir)/%, $(src_transform)): $(obj_read) $(obj_write) $(obj_periodic) $(obj_transform) | $(prog_dir)
	$(fc) -o $@ $(bin_dir)/$(notdir $@).o $(obj_read) $(obj_write) $(obj_periodic) $(linking_flags)

# create directiories
$(bin_dir):
	mkdir -p $(bin_dir)

$(prog_dir):
	mkdir -p $(prog_dir)

$(install_exec_dir):
	mkdir -p $(install_exec_dir)

$(install_lib_dir):
	mkdir -p $(install_lib_dir)
