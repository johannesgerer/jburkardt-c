#!/bin/bash
#
cp quadrule.h /$HOME/include
#
gcc -c quadrule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule.c."
  exit
fi
#
mv quadrule.o ~/libc/$ARCH/quadrule.o
#
echo "Library installed as ~/libc/$ARCH/quadrule.o"
