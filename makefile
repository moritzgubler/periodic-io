ifneq ("$(wildcard make.inc)","")
    include make.inc
else
    include make-gnu.inc
endif


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
src_om_fp_dir := dev/src/OM-fingerprint
src_fingerprint_programs_dir := dev/src/fingerprint-programs

# location of binaries:
bin_dir := dev/bin
prog_dir := executables


# find all relevant source files.
src_read := $(wildcard $(src_read_dir)/*.f90)
src_write := $(wildcard $(src_write_dir)/*.f90)
src_transform := $(wildcard $(src_transform_dir)/*.f90)
src_periodic := $(wildcard $(src_periodic_dir)/*.f90)
src_om_fp := $(wildcard $(src_om_fp_dir)/*.f90)
src_fingerprint_programs := $(wildcard $(src_fingerprint_programs_dir)/*.f90)

# set location of object files:
obj_read := $(patsubst $(src_read_dir)/%.f90, $(bin_dir)/%.o, $(src_read))
obj_write := $(patsubst $(src_write_dir)/%.f90, $(bin_dir)/%.o, $(src_write))
obj_transform := $(patsubst $(src_transform_dir)/%.f90, $(bin_dir)/%.o, $(src_transform))
obj_periodic := $(patsubst $(src_periodic_dir)/%.f90, $(bin_dir)/%.o, $(src_periodic))
obj_om_fp := $(patsubst $(src_om_fp_dir)/%.f90, $(bin_dir)/%.o, $(src_om_fp))
obj_fingerprint_programs := $(patsubst $(src_fingerprint_programs_dir)/%.f90, $(bin_dir)/%.o, $(src_fingerprint_programs))

lib_file := $(bin_dir)/$(lib_name)

# set program names (all files from $(src_transform_dir) should be executable)
programs := $(patsubst $(src_transform_dir)/%.f90, $(prog_dir)/%, $(src_transform))
fp_programs := $(patsubst $(src_fingerprint_programs_dir)/%.f90, $(prog_dir)/%, $(src_fingerprint_programs))

# set phony targets
.PHONY: default debug clean install_bin install_lib install test doc uninstall

default: $(programs) $(lib_file) $(fp_programs)

doc:
	(cd dev && ford doc.md)

# call make debug to debug the following variables.
debug:
	@echo "src_read = $(src_read)"
	@echo "src_write = $(src_write)"
	@echo "src_transform = $(src_transform)"
	@echo "src_periodic = $(src_periodic)"
	@echo "src_om_fp = $(src_om_fp)"
	@echo "src_fingerprint_programs = $(src_fingerprint_programs)"
	@echo "obj_read = $(obj_read)"
	@echo "obj_write = $(obj_write)"
	@echo "obj_transform = $(obj_transform)"
	@echo "obj_periodic = $(obj_periodic)"
	@echo "obj_om_fp = $(obj_om_fp)"
	@echo "obj_fingerprint_programs = $(obj_fingerprint_programs)"
	@echo "lib_file = $(lib_file)"
	@echo "programs = $(programs)"
	@echo "fp_programs = $(fp_programs)"

test: default
	(cd dev/test && ./test_all.sh)

clean:
	rm -rf $(bin_dir) $(prog_dir)

install: install_lib install_bin

uninstall:
	rm -f $(install_lib_dir)/$(lib_name)
	rm -f $(patsubst $(prog_dir)/%, $(install_exec_dir)/%, $(programs))

install_bin: default | $(install_exec_dir)
	cp $(prog_dir)/* $(install_exec_dir)/.

install_lib: default $(lib_file) | $(install_lib_dir)
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

# compile overlap matrix dir
$(obj_om_fp): $(bin_dir)/%.o: $(src_om_fp_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

# compile src_fingerprint_programs_dir
$(obj_fingerprint_programs): $(bin_dir)/%.o: $(src_fingerprint_programs_dir)/%.f90 | $(bin_dir)
	$(fc) -c $< -o $@ $(compile_flags)

#link all transformer programs:
$(patsubst $(src_transform_dir)/%.f90, $(prog_dir)/%, $(src_transform)): $(obj_read) $(obj_write) $(obj_periodic) $(obj_transform) | $(prog_dir)
	$(fc) -o $@ $(bin_dir)/$(notdir $@).o $(obj_read) $(obj_write) $(obj_periodic) $(linking_flags)

#link all fingerprint programs:
$(patsubst $(src_fingerprint_programs_dir)/%.f90, $(prog_dir)/%, $(src_fingerprint_programs)): $(lib_file) $(obj_om_fp) $(obj_fingerprint_programs) | $(prog_dir)
	$(fc) -o $@ $(bin_dir)/$(notdir $@).o $(lib_file) $(obj_om_fp) $(linking_flags) $(lapack)

# create directiories
$(bin_dir):
	mkdir -p $(bin_dir)

$(prog_dir):
	mkdir -p $(prog_dir)

$(install_exec_dir):
	mkdir -p $(install_exec_dir)

$(install_lib_dir):
	mkdir -p $(install_lib_dir)
