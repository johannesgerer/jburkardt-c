#!/bin/bash
#
cp cube_exactness.h /$HOME/include
#
gcc -c -I/$HOME/include cube_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_exactness.c"
  exit
fi
#
mv cube_exactness.o ~/libc/$ARCH/cube_exactness.o
#
echo "Library installed as ~/libc/$ARCH/cube_exactness.o"
