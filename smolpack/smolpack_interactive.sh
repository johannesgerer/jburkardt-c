#!/bin/bash
#
gcc -c -g smolpack_interactive.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling smolpack_interactive.c."
  exit
fi
rm compiler.txt
#
gcc smolpack_interactive.o -L/$HOME/libc/$ARCH -lsmolpack -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading smolpack_interactive.o."
  exit
fi
#
rm smolpack_interactive.o
#
mv a.out /$HOME/binc/$ARCH/smolpack_interactive
#
echo "Executable installed as ~/binc/$ARCH/smolpack_interactive."
