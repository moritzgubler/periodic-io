SRC := .
OBJ := bin

fc := gfortran

linking_flags = -Ofast -llapack -lblas -fopenmp
compile_flags = $(linking_flags)

molecular_io_lib=../../bin/libmolecularIO.a

SOURCES := $(wildcard $(SRC)/*.f90)
OBJECTS := $(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(SOURCES))

default: fingerprint_distance

fingerprint_distance: $(OBJ)/hung.o $(OBJ)/fpd.o $(OBJ)/main.o $(OBJ)/sym2rcov.o
	$(fc) $^ ${molecular_io_lib} -o $@ $(linking_flags) -llapack -lblas

$(OBJ)/%.o: $(SRC)/%.f90 | bin
	$(fc) -c $< -o $@ $(compile_flags)

bin:
	mkdir -p bin

clean:
	rm -rf fingerprint_distance bin


.PHONY: clean default
