#!/bin/bash
#
gcc -c -g -I/$HOME/include fem2d_pack_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_pack_prb.c"
  exit
fi
rm compiler.txt
#
gcc fem2d_pack_prb.o /$HOME/libc/$ARCH/fem2d_pack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_pack_prb.o."
  exit
fi
#
rm fem2d_pack_prb.o
#
mv a.out fem2d_pack_prb
./fem2d_pack_prb > fem2d_pack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_pack_prb."
  exit
fi
rm fem2d_pack_prb
#
echo "Program output written to fem2d_pack_prb_output.txt"
