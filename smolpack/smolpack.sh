#!/bin/bash
#
cp smolpack.h /$HOME/include
#
gcc -c -g ccsmolyak.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ccsmolyak.c."
  exit
fi
rm compiler.txt
#
gcc -c -g smolyak.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling smolyak.c."
  exit
fi
rm compiler.txt
#
ar qc libsmolpack.a *.o
rm *.o
mv libsmolpack.a ~/libc/$ARCH
#
echo "Library installed as ~/libc/ARCH/libsmolpack.a"
