#!/bin/bash
#
gcc -c -g -I /$HOME/include super_lu_s3.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling super_lu_s3.c"
  exit
fi
rm compiler.txt
#
gcc super_lu_s3.o -L/$HOME/libc/$ARCH -lsuper_lu -lcxml -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading super_lu_s3.o"
  exit
fi
#
rm super_lu_s3.o
#
mv a.out super_lu_s3
./super_lu_s3 < g10_rua.txt > super_lu_s3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running super_lu_s3"
  exit
fi
rm super_lu_s3
#
echo "Program output written to super_lu_s3_output.txt"
