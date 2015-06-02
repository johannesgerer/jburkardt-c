#!/bin/bash
#
#  Compile
#
gcc -c -I/$HOME/include z_sample_st.c
if [ $? -ne 0 ]; then
  echo "Errors compiling z_sample_st.c"
  exit
fi
#
#  Link and load
#
#gcc z_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
gfortran z_sample_st.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -lblas
#gcc z_sample_st.o -L$HOME/lib/$ARCH \
#  -L/$HOME/libc/$ARCH -lsuperlu_4.3 -lm -framework veclib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading z_sample_st.o"
  exit
fi
rm z_sample_st.o
mv a.out z_sample_st
#
#  Run
#
./z_sample_st < sample_cst.txt > z_sample_st_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running z_sample_st"
  exit
fi
rm z_sample_st
#
#  Terminate.
#
echo "Program output written to z_sample_st_output.txt"
