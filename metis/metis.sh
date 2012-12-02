#!/bin/bash
#
cp metis.h /$HOME/include
#
for FILE in `ls -1 *.c`;
do
  gcc -c -g -I. $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
#
ar qc libmetis.a *.o
rm *.o
#
mv libmetis.a ~/libc/$ARCH
#
echo "Library installed as ~/libc/$ARCH/libmetis.a"

