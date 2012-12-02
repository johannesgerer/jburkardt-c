#!/bin/bash
#
gcc -c -g -I/$HOME/include csparse_demo1.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csparse_demo1.c."
  exit
fi
rm compiler.txt
#
gcc csparse_demo1.o $HOME/libc/$ARCH/csparse.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading csparse_demo1.o + csparse.o."
  exit
fi
#
rm csparse_demo1.o
#
chmod ugo+x a.out
mv a.out csparse_demo1
#
./csparse_demo1 < kershaw.st > csparse_demo1_output.txt
rm csparse_demo1
#
echo "Program output written to csparse_demo1_output.txt."
