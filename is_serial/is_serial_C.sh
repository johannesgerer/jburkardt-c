#!/bin/bash
#
cp npbparams_C.h npbparams.h
#
gcc -c is_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling is_serial.c."
  exit
fi
rm compiler.txt
#
gcc -c c_print_results.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c_print_results.c."
  exit
fi
rm compiler.txt
#
gcc -c c_timers.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c_timers.c."
  exit
fi
rm compiler.txt
#
gcc -c wtime.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime.c."
  exit
fi
rm compiler.txt
#
gcc is_serial.o c_print_results.o c_timers.o wtime.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading is_serial.o + c_print_results.o + c_timers.o + wtime.o"
  exit
fi
#
rm npbparams.h
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/is_serial_C
#
echo "Executable installed as ~/binc/$ARCH/is_serial_C"
