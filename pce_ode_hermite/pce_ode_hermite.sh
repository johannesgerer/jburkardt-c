#!/bin/bash
#
cp pce_ode_hermite.h /$HOME/include
#
gcc -c -g -I /$HOME/include pce_ode_hermite.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_ode_hermite.c"
  exit
fi
rm compiler.txt
#
mv pce_ode_hermite.o ~/libc/$ARCH/pce_ode_hermite.o
#
echo "Library installed as ~/libc/$ARCH/pce_ode_hermite.o"
