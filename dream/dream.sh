#!/bin/bash
#
cp dream.h /$HOME/include
cp dream_user.h /$HOME/include
#
gcc -c -I/$HOME/include dream.c
if [ $? -ne 0 ]; then
  echo "Errors compiling dream.c"
  exit
fi
#
mv dream.o ~/libc/$ARCH/dream.o
#
echo "Library installed as ~/libc/$ARCH/dream.o"
