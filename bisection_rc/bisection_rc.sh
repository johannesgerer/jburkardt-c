#!/bin/bash
#
cp bisection_rc.h /$HOME/include
#
gcc -c -I/$HOME/include bisection_rc.c
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_rc.c"
  exit
fi
#
mv bisection_rc.o ~/libc/$ARCH/bisection_rc.o
#
echo "Library installed as ~/libc/$ARCH/bisection_rc.o"
