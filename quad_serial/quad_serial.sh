#!/bin/bash
#
gcc -c quad_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_serial.c."
  exit
fi
rm compiler.txt
#
gcc quad_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking quad_serial.o."
  exit
fi
#
rm quad_serial.o
#
mv a.out quad_serial
./quad_serial > quad_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quad_serial."
  exit
fi
rm quad_serial
#
echo "Program output written to quad_serial_output.txt"
