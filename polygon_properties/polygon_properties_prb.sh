#!/bin/bash
#
gcc -c -I/$HOME/include polygon_properties_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_properties_prb.c"
  exit
fi
rm compiler.txt
#
gcc polygon_properties_prb.o /$HOME/libc/$ARCH/polygon_properties.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polygon_properties_prb.o."
  exit
fi
#
rm polygon_properties_prb.o
#
mv a.out polygon_properties_prb
./polygon_properties_prb > polygon_properties_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polygon_properties_prb."
  exit
fi
rm polygon_properties_prb
#
echo "Program output written to polygon_properties_prb_output.txt"
