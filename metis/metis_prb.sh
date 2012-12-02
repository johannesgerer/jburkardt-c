#!/bin/bash
#
#  Get the graph file.
#
cp ../../data/metis_graph/4elt.graph .
#
gcc -c -g -I. metis_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling metis_prb.c."
  exit
fi
rm compiler.txt
#
gcc -c -g -I. io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling io.c."
  exit
fi
rm compiler.txt
#
gcc metis_prb.o io.o -L/$HOME/libc/$ARCH -lmetis -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading metis_prb.o."
  exit
fi
#
rm metis_prb.o
rm io.o
#
mv a.out metis_prb
./metis_prb 4elt.graph > metis_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running metis_prb."
  exit
fi
rm metis_prb
rm 4elt.graph
#
echo "Program output written to metis_prb_output.txt"
