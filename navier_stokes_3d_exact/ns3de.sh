#!/bin/bash
#
cp ns3de.h /$HOME/include
#
gcc -c -I/$HOME/include ns3de.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ns3de.c"
  exit
fi
#
mv ns3de.o ~/libc/$ARCH/ns3de.o
#
echo "Library installed as ~/libc/$ARCH/ns3de.o"
