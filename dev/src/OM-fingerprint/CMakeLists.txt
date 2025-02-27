add_library(${fingerprintlib} OM.f90 hung.f90)
install(TARGETS ${fingerprintlib} DESTINATION lib)