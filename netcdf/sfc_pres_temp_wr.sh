#!/bin/bash
#
cp ~/include/netcdf.h .
gcc -c -g sfc_pres_temp_wr.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sfc_pres_temp_wr.c."
  exit
fi
rm compiler.txt
rm netcdf.h
#
gcc sfc_pres_temp_wr.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sfc_pres_temp_wr.o."
  exit
fi
#
rm sfc_pres_temp_wr.o
#
mv a.out sfc_pres_temp_wr
./sfc_pres_temp_wr > sfc_pres_temp_wr_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sfc_pres_temp_wr."
  exit
fi
rm sfc_pres_temp_wr
#
echo "Program output written to sfc_pres_temp_wr_output.txt"
