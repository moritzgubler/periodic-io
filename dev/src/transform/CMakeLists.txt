set(execname back2cell)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname back2cell-wildcard)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname periodic2qe)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname reshape-cell-wildcard)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname reshape-cell)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname set-comment-periodic)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname supercell)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname transform-periodic-wildcard)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)

set(execname transform-periodic)
add_executable(${execname} ${execname}.f90)
target_link_libraries(${execname} ${MAINLIB})
install(TARGETS ${execname} DESTINATION bin)