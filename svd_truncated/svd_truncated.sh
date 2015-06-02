#!/bin/bash
#
gcc -c -g -I$HOME/include svd_truncated.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_truncated.c"
  exit
fi
rm compiler.txt
#
gcc svd_truncated.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_truncated.o"
  exit
fi
#
rm svd_truncated.o
#
chmod ugo+x a.out
mv a.out svd_truncated
./svd_truncated > svd_truncated_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running svd_truncated"
  exit
fi
rm svd_truncated
#
echo "Program output written to svd_truncated_output.txt"
