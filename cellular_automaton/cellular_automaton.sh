#!/bin/bash
#
gcc -c -g cellular_automaton.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cellular_automaton.c"
  exit
fi
rm compiler.txt
#
gcc cellular_automaton.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cellular_automaton.o"
  exit
fi
rm cellular_automaton.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/cellular_automaton
#
echo "Program installed as ~/binc/$ARCH/cellular_automaton"
