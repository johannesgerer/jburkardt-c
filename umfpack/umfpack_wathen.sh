#!/bin/bash
#
#  Compile
#
gcc -c -I/usr/local/include umfpack_wathen.c
if [ $? -ne 0 ]; then
  echo "Errors compiling umfpack_wathen.c"
  exit
fi
#
#  Link and load
#
gfortran umfpack_wathen.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
  -lsuitesparseconfig -lblas
#gcc umfpack_wathen.o -L/usr/local/lib -L/$HOME/lib/$ARCH -lumfpack -lamd -lcholmod -lcolamd -lm \
#  -lsuitesparseconfig -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading umfpack_wathen.o"
  exit
fi
rm umfpack_wathen.o
mv a.out umfpack_wathen
#
#  Run
#
./umfpack_wathen > umfpack_wathen_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running umfpack_wathen"
  exit
fi
rm umfpack_wathen
#
#  Terminate.
#
echo "Program output written to umfpack_wathen_output.txt"
