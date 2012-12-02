#!/bin/bash
#
gcc -c -g -I/$HOME/include csparse_demo2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csparse_demo2.c."
  exit
fi
rm compiler.txt
#
gcc -c -g -I/$HOME/include csparse_demo.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csparse_demo.c."
  exit
fi
rm compiler.txt
#
gcc csparse_demo2.o csparse_demo.o $HOME/libc/$ARCH/csparse.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading csparse_demo2.o + csparse.o."
  exit
fi
#
rm csparse_demo2.o
rm csparse_demo.o
#
chmod ugo+x a.out
mv a.out csparse_demo2
#
./csparse_demo2 < kershaw.st > csparse_demo2_output.txt
rm csparse_demo2
#
echo "Program output written to csparse_demo2_output.txt."
