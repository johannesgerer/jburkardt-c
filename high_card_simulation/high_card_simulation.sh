#!/bin/bash
#
cp high_card_simulation.h /$HOME/include
#
gcc -c -g -I/$HOME/include high_card_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling high_card_simulation.c"
  exit
fi
rm compiler.txt
#
mv high_card_simulation.o ~/libc/$ARCH/high_card_simulation.o
#
echo "Library installed as ~/libc/$ARCH/high_card_simulation.o"
