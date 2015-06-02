#!/bin/bash
#
gcc -c -g -I/$HOME/include tetrahedron_integrals_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_integrals_prb.c"
  exit
fi
rm compiler.txt
#
gcc tetrahedron_integrals_prb.o /$HOME/libc/$ARCH/tetrahedron_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_integrals_prb.o."
  exit
fi
#
rm tetrahedron_integrals_prb.o
#
mv a.out tetrahedron_integrals_prb
./tetrahedron_integrals_prb > tetrahedron_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tetrahedron_integrals_prb."
  exit
fi
rm tetrahedron_integrals_prb
#
echo "Program output written to tetrahedron_integrals_prb_output.txt"
