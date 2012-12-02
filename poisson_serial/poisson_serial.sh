#!/bin/bash
#
gcc -c poisson_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling poisson_serial.c"
  exit
fi
rm compiler.txt
#
gcc poisson_serial.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading poisson_serial.o"
  exit
fi
rm poisson_serial.o
#
mv a.out poisson_serial
./poisson_serial > poisson_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running poisson_serial"
  exit
fi
rm poisson_serial
#
echo "Program output written to poisson_serial_output.txt"
