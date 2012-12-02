#!/bin/bash
#
gcc -c spring_ode2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode2.c"
  exit
fi
rm compiler.txt
#
gcc spring_ode2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode2.o"
  exit
fi
rm spring_ode2.o
#
mv a.out ~/binc/$ARCH/spring_ode2
#
echo "Executable installed as ~/binc/$ARCH/spring_ode2"
