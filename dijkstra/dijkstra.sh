#!/bin/bash
#
gcc -c dijkstra.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dijkstra.c."
  exit
fi
rm compiler.txt
#
gcc dijkstra.o
if [ $? -ne 0 ]; then
  echo "Errors linking dijkstra.o."
  exit
fi
#
rm dijkstra.o
#
mv a.out dijkstra
./dijkstra > dijkstra_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dijkstra."
  exit
fi
rm dijkstra
#
echo "Program output written to dijkstra_output.txt"
