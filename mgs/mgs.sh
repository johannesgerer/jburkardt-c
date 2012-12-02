#!/bin/bash
#
gcc -c mgs.c
if [ $? -ne 0 ]; then
  echo "Errors compiling mgs.c"
  exit
fi
#
echo "The mgs.c file was compiled."
