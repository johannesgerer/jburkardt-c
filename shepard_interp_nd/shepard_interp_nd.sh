#!/bin/bash
#
cp shepard_interp_nd.h /$HOME/include
#
gcc -c -g -I/$HOME/include shepard_interp_nd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_nd.c"
  exit
fi
rm compiler.txt
#
mv shepard_interp_nd.o ~/libc/$ARCH/shepard_interp_nd.o
#
echo "Library installed as ~/libc/$ARCH/shepard_interp_nd.o"
