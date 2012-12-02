#!/bin/bash
#
cp ~/include/netcdf.h .
gcc -c -g sfc_pres_temp_rd.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sfc_pres_temp_rd.c."
  exit
fi
rm compiler.txt
rm netcdf.h
#
gcc sfc_pres_temp_rd.o -L$HOME/lib/$ARCH -lnetcdf -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sfc_pres_temp_rd.o."
  exit
fi
#
rm sfc_pres_temp_rd.o
#
mv a.out sfc_pres_temp_rd
./sfc_pres_temp_rd > sfc_pres_temp_rd_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sfc_pres_temp_rd."
  exit
fi
rm sfc_pres_temp_rd
#
echo "Program output written to sfc_pres_temp_rd_output.txt"
