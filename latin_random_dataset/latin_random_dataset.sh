#!/bin/bash
#
gcc -c -I$HOME/include latin_random_dataset.c
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random_dataset.c"
  exit
fi
#
gcc latin_random_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_random_dataset.o."
  exit
fi
#
rm latin_random_dataset.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/latin_random_dataset
#
echo "Executable installed as ~/binc/$ARCH/latin_random_dataset"
