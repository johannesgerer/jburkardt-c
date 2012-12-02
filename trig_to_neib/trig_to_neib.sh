#!/bin/bash
#
gcc -c -g trig_to_neib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling trig_to_neib.c"
  exit
fi
#
gcc trig_to_neib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors loading trig_to_neib.c"
  exit
fi
#
rm trig_to_neib.o
mv a.out ~/binc/$ARCH/trig_to_neib
#
echo "Executable installed as ~/binc/$ARCH/trig_to_neib."
