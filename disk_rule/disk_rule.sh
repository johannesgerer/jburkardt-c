#!/bin/bash
#
cp disk_rule.h /$HOME/include
#
gcc -c -g -I/$HOME/include disk_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_rule.c"
  exit
fi
rm compiler.txt
#
mv disk_rule.o ~/libc/$ARCH/disk_rule.o
#
echo "Library installed as ~/libc/$ARCH/disk_rule.o"
