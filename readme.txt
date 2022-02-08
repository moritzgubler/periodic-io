This project creates a library that will be able to read and write most of the
atomic file formats. It also contains conversion software that allows you to
convert the atomic easily to other formats.

The cif reader does not support the full cif syntax. Nevertheless it should
be able to read most of the files. It is designed to read cifs from
materialsproject.org

Supported file formats:
- .ascii
- .in (FHI-AIMS)
- .xyz
- .cif
- .gen (DFTB)
- .qe (QUANTUM-ESPRESSO) (not tested, use with care)

Installation:
Set your favorite fortran compiler in the makefile and the linking and
compiling flags (no flags required).
"make install" will compile everything in the src directory and install all the
conversion tools to $HOME/bin. If this is not on path add
if [ -d "$HOME/bin" ] ; then
    export PATH="$HOME/bin:$PATH"
fi
The library will be installed to
$HOME/.local/lib/libmolecularIO.a.
This is not on the $LIBRARY_PATH by default. Add

# set PATH so it includes user's private lib if it exists
if [ -d "$HOME/.local/lib" ] ; then
    export LIBRARY_PATH="$HOME/.local/lib:$LIBRARY_PATH"
fi

to your ~/.profile file. Now reboot or "source ~/.profile".

You will now be able to convert files with
"ascii2in file1.ascii new_file.in" etc.. If the ascii2in is not found try
"source ~/.profile". Now the console should be able to find the programs.
Compiling with the library works like that:
gfortran *.f90 -lmolecularIO

If you have any questions or you want to contribute feel free to ask me.

Adding your reader/writer subroutines:
place the subroutines that read files in dev/src/read and the ones that write files
in dev/src/write. The make file will find them and include them in the library when
the code is recompiled with make.
Put transform programs in dev/src/transform. Once again the makefile will find
them on it's own and compile and/or install them. The filename whithout the .f90 ending will be the name of
the program.