#!/bin/bash
#
gcc -c -g triangulation_node_to_element.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_node_to_element.c."
  exit
fi
rm compiler.txt
#
gcc triangulation_node_to_element.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_node_to_element.o."
  exit
fi
#
rm triangulation_node_to_element.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangulation_node_to_element
#
echo "Executable installed as ~/binc/$ARCH/triangulation_node_to_element"
