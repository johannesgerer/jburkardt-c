#!/bin/bash
#
cp line_nco_rule.h /$HOME/include
#
gcc -c -I/$HOME/include line_nco_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_nco_rule.c"
  exit
fi
#
mv line_nco_rule.o ~/libc/$ARCH/line_nco_rule.o
#
echo "Library installed as ~/libc/$ARCH/line_nco_rule.o"
