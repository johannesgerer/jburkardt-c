#!/bin/bash
#
gcc -c spring_ode.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode.c"
  exit
fi
rm compiler.txt
#
gcc spring_ode.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode.o"
  exit
fi
rm spring_ode.o
#
mv a.out ~/binc/$ARCH/spring_ode
#
echo "Executable installed as ~/binc/$ARCH/spring_ode"
