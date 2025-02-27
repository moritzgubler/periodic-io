set(execname fingerprint-distance)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${fingerprintlib} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname env-diff-max)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${fingerprintlib} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname env-diff-min)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${fingerprintlib} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)
