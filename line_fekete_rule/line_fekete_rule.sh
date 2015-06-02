#!/bin/bash
#
cp line_fekete_rule.h /$HOME/include
#
gcc -c -I/$HOME/include line_fekete_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_fekete_rule.c"
  exit
fi
#
mv line_fekete_rule.o ~/libc/$ARCH/line_fekete_rule.o
#
echo "Library installed as ~/libc/$ARCH/line_fekete_rule.o"
