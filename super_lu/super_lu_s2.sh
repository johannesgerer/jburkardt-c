#!/bin/bash
#
gcc -c -g -I /$HOME/include super_lu_s2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling super_lu_s2.c"
  exit
fi
rm compiler.txt
#
gcc super_lu_s2.o -L/$HOME/libc/$ARCH -lsuper_lu -lcxml -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading super_lu_s2.o"
  exit
fi
#
rm super_lu_s2.o
#
mv a.out super_lu_s2
./super_lu_s2 < g10_rua.txt > super_lu_s2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running super_lu_s2"
  exit
fi
rm super_lu_s2
#
echo "Program output written to super_lu_s2_output.txt"
