#!/bin/bash
#
cp ~/include/netcdf.h .
gcc -c -g pres_temp_4D_rd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pres_temp_4D_rd.c."
  exit
fi
rm compiler.txt
rm netcdf.h
#
gcc pres_temp_4D_rd.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pres_temp_4D_rd.o."
  exit
fi
#
rm pres_temp_4D_rd.o
#
mv a.out pres_temp_4D_rd
./pres_temp_4D_rd > pres_temp_4D_rd_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pres_temp_4D_rd."
  exit
fi
rm pres_temp_4D_rd
#
echo "Program output written to pres_temp_4D_rd_output.txt"
