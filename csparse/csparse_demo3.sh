#!/bin/bash
#
gcc -c -g -I/$HOME/include csparse_demo3.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling csparse_demo3.c."
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
gcc csparse_demo3.o csparse_demo.o $HOME/libc/$ARCH/csparse.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading csparse_demo3.o + csparse.o."
  exit
fi
#
rm csparse_demo3.o
rm csparse_demo.o
#
chmod ugo+x a.out
mv a.out csparse_demo3
#
./csparse_demo3 < kershaw.st > csparse_demo3_output.txt
rm csparse_demo3
#
echo "Program output written to csparse_demo3_output.txt."
