#!/bin/bash
#
gcc -c edgelist.c
if [ $? -= 0 ]; then
  echo "Errors compiling edgelist.c"
  exit
fi
gcc -c geometry.c
if [ $? -= 0 ]; then
  echo "Errors compiling geometry.c"
  exit
fi
gcc -c heap.c
if [ $? -= 0 ]; then
  echo "Errors compiling heap.c"
  exit
fi
gcc -c main.c
if [ $? -= 0 ]; then
  echo "Errors compiling main.c"
  exit
fi
gcc -c memory.c
if [ $? -= 0 ]; then
  echo "Errors compiling memory.c"
  exit
fi
gcc -c output.c
if [ $? -= 0 ]; then
  echo "Errors compiling output.c"
  exit
fi
gcc -c voronoi.c
if [ $? -= 0 ]; then
  echo "Errors compiling voronoi.c"
  exit
fi
#
gcc *.o -lm
if [ $? -= 0 ]; then
  echo "Errors linking and loading the object files"
  exit
fi
#
rm *.o
mv a.out ~/binc/$ARCH/sweep2
#
echo "Executable installed as ~/binc/$ARCH/sweep2"
