#!/bin/bash
#
cp test_opt_con.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_opt_con.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_opt_con.c"
  exit
fi
rm compiler.txt
#
mv test_opt_con.o ~/libc/$ARCH/test_opt_con.o
#
echo "Library installed as ~/libc/$ARCH/test_opt_con.o"
