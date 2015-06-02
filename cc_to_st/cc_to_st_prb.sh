#!/bin/bash
#
gcc -c -I/$HOME/include cc_to_st_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_to_st_prb.c"
  exit
fi
#
gcc -o cc_to_st_prb cc_to_st_prb.o /$HOME/libc/$ARCH/cc_to_st.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cc_to_st_prb.o."
  exit
fi
#
rm cc_to_st_prb.o
#
./cc_to_st_prb > cc_to_st_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cc_to_st_prb."
  exit
fi
rm cc_to_st_prb
#
echo "Program output written to cc_to_st_prb_output.txt"
