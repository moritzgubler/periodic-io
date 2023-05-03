set(SOURCES
    read_ascii.f90
    read_cif.f90
    read_dftb.f90
    read_extxyz.f90
    read_in.f90
    read_periodic.f90
    read_quantum_espresso.f90
    read_vasp.f90
    read_xyz.f90
)

target_sources(${MAINLIB} PRIVATE ${SOURCES})