#!/bin/bash
#
gcc -c pce_burgers.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_burgers.c"
  exit
fi
rm compiler.txt
#
gcc pce_burgers.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pce_burgers.o"
  exit
fi
rm pce_burgers.o
#
mv a.out ~/binc/$ARCH/pce_burgers
#
echo "Executable installed as ~/binc/$ARCH/pce_burgers"
