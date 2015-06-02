#!/bin/bash
#
cp fem2d_bvp_serene.h /$HOME/include
#
gcc -c -I /$HOME/include fem2d_bvp_serene.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_serene.c"
  exit
fi
#
mv fem2d_bvp_serene.o ~/libc/$ARCH/fem2d_bvp_serene.o
#
echo "Library installed as ~/libc/$ARCH/fem2d_bvp_serene.o"
