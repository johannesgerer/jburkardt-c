#!/bin/bash
#
cp ~/include/netcdf.h .
gcc -c -g pres_temp_4D_wr.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pres_temp_4D_wr.c."
  exit
fi
rm compiler.txt
rm netcdf.h
#
gcc pres_temp_4D_wr.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pres_temp_4D_wr.o."
  exit
fi
#
rm pres_temp_4D_wr.o
#
mv a.out pres_temp_4D_wr
./pres_temp_4D_wr > pres_temp_4D_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pres_temp_4D_wr."
  exit
fi
rm pres_temp_4D_wr
#
echo "Program output written to pres_temp_4D_wr_output.txt"
