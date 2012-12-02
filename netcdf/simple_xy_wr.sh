#!/bin/bash
#
cp ~/include/netcdf.h .
gcc -c -g simple_xy_wr.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_xy_wr.c."
  exit
fi
rm compiler.txt
rm netcdf.h
#
gcc simple_xy_wr.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_xy_wr.o."
  exit
fi
#
rm simple_xy_wr.o
#
mv a.out simple_xy_wr
./simple_xy_wr > simple_xy_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_xy_wr."
  exit
fi
rm simple_xy_wr
#
echo "Program output written to simple_xy_wr_output.txt"
