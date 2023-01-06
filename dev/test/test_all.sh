#!/bin/bash

tempdir=tempdir
mkdir -p $tempdir

tester=../../executables/transform-periodic

tol=0.0001

echo "Self constistent reading / writing check:"
for file in test_dat/*.ascii; do
  basefile="$(basename $file .ascii)"
  ciffile=$tempdir/${basefile}.cif
  infile=$tempdir/${basefile}.in
  genfile=$tempdir/${basefile}.gen
  vaspfile=$tempdir/${basefile}.vasp
  poscar=$tempdir/POSCAR
  extxyz=$tempdir/${basefile}.extxyz
  comparefile=$tempdir/${basefile}.ascii
  rm -f $comparefile
  ./$tester $file $ciffile
  ./$tester $ciffile $poscar
  ./$tester $poscar $infile
  ./$tester $infile $vaspfile
  ./$tester $vaspfile $extxyz
  ./$tester $extxyz $genfile
  ./$tester $genfile $comparefile
  res=$(./../../executables/fingerprint-distance $file $comparefile)
  if python -c "exit(0 if $res < $tol else 1)"; then
    pass="passed"
  else
    pass="failed"
  fi
  echo "Test: $basefile read_ascii, write_cif, read_cif, write_in, read_in, write_gen, read_gen, write_ascii -> $pass with fp_distance: $res"
done

./../../executables/fingerprint-distance test_free_boundary_conditions/CS2H12-1.xyz test_free_boundary_conditions/CS2H12-2.xyz
