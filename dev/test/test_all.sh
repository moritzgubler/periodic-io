#!/bin/bash

tempdir=tempdir
mkdir -p $tempdir

tester=../../executables/transform_periodic

tol=0.0001

echo "Self constistent reading / writing check:"
for file in test_dat/*.ascii; do
  basefile="$(basename $file .ascii)"
  ciffile=$tempdir/${basefile}.cif
  infile=$tempdir/${basefile}.in
  genfile=$tempdir/${basefile}.gen
  vaspfile=$tempdir/${basefile}.vasp
  poscar=$tempdir/POSCAR
  comparefile=$tempdir/${basefile}.ascii
  rm -f $comparefile
  ./$tester $file $ciffile
  ./$tester $ciffile $poscar
  ./$tester $poscar $infile
  ./$tester $infile $vaspfile
  ./$tester $vaspfile $genfile
  ./$tester $genfile $comparefile
  res=$(./test_src_fingerprint/fingerprint_distance $file $comparefile)
  if python -c "exit(0 if $res < $tol else 1)"; then
    pass="passed"
  else
    pass="failed"
  fi
  echo "Test: $basefile read_ascii, write_cif, read_cif, write_in, read_in, write_gen, read_gen, write_ascii -> $pass with fp_distance: $res"
done
