#!/bin/bash
#
cp simplex_coordinates.H /$HOME/include
#
g++ -c -g -I /$HOME/include simplex_coordinates.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_coordinates.C."
  exit
fi
rm compiler.txt
#
mv simplex_coordinates.o ~/libcpp/$ARCH/simplex_coordinates.o
#
echo "Library installed as ~/libcpp/$ARCH/simplex_coordinates.o"
