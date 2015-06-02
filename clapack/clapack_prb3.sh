#!/bin/bash
#
#  The -I switch allows us to access the include file clapack.h.
#
#myINCLUDE=/usr/local/include/clapack
#myLIB=/usr/local/lib
#gcc -c -I$myINCLUDE clapack_prb.c
#
gcc -c -I/$HOME/include clapack_prb3.c
if [ $? -ne 0 ]; then
  echo "Errors compiling clapack_prb3.c"
  exit
fi
#
#  The -L switch allows us to access 4 libraries associated with CLAPACK.
#
gcc clapack_prb3.o -L$HOME/lib/$ARCH -lclapack_lapack -lclapack_blas -lf2c -lclapack_tmglib -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clapack_prb3.o."
  exit
fi
rm clapack_prb3.o
#
mv a.out clapack_prb3
./clapack_prb3 > clapack_prb3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running clapack_prb3"
  exit
fi
#
rm clapack_prb3
#
echo "clapack_prb3 output written to clapack_prb3_output.txt"

