#!/bin/bash
#
gcc -c rand_defect.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rand_defect.c."
  exit
fi
rm compiler.txt
#
gcc rand_defect.o
if [ $? -ne 0 ]; then
  echo "Errors linking rand_defect.o."
  exit
fi
#
rm rand_defect.o
#
mv a.out rand_defect
./rand_defect >  rand_defect_output.txt
./rand_defect >> rand_defect_output.txt
./rand_defect >> rand_defect_output.txt
./rand_defect >> rand_defect_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rand_defect."
  exit
fi
rm rand_defect
#
echo "Program output written to rand_defect_output.txt"
