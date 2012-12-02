#!/bin/bash
#
cp ode.h /$HOME/include
#
gcc -c -g -I /$HOME/include ode.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode.c"
  exit
fi
rm compiler.txt
#
mv ode.o ~/libc/$ARCH/ode.o
#
echo "Library installed as ~/libc/$ARCH/ode.o"
