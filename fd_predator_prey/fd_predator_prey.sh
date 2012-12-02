#!/bin/bash
#
gcc -c -g fd_predator_prey.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd_predator_prey.C"
  exit
fi
rm compiler.txt
#
gcc fd_predator_prey.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd_predator_prey.o"
  exit
fi
rm fd_predator_prey.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fd_predator_prey
#
echo "Executable installed as ~/binc/$ARCH/fd_predator_prey"
