#!/bin/bash
#
cp pwl_interp_2d_scattered.h /$HOME/include
#
gcc -c -g -I /$HOME/include pwl_interp_2d_scattered.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_scattered.c"
  exit
fi
rm compiler.txt
#
mv pwl_interp_2d_scattered.o ~/libc/$ARCH/pwl_interp_2d_scattered.o
#
echo "Library installed as ~/libc/$ARCH/pwl_interp_2d_scattered.o"
