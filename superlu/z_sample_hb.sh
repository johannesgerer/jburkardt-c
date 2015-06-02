#!/bin/bash
#
#  Compile
#
gcc -c -I/$HOME/include z_sample_hb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling z_sample_hb.c"
  exit
fi
#
#  Link and load
#
#gcc z_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran z_sample_hb.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
#gcc z_sample_hb.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading z_sample_hb.o"
  exit
fi
rm z_sample_hb.o
mv a.out z_sample_hb
#
#  Run
#
./z_sample_hb < sample_cua.txt > z_sample_hb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running z_sample_hb"
  exit
fi
rm z_sample_hb
#
#  Terminate.
#
echo "Program output written to z_sample_hb_output.txt"
