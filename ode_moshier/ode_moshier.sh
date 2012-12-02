#!/bin/bash
#
cp ode_moshier.h /$HOME/include
#
gcc -c -g adams3.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling adams3.c"
  exit
fi
rm compiler.txt
#
gcc -c -g rungek.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rungek.c"
  exit
fi
rm compiler.txt
#
ar qc libode_moshier.a *.o
rm *.o
mv libode_moshier.a ~/libc/$ARCH/libode_moshier.a
#
echo "Library installed as ~/libc/$ARCH/libode_moshier.a"
