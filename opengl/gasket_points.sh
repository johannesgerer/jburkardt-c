#!/bin/bash
#
gcc -c gasket_points.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gasket_points.c"
  exit
fi
rm compiler.txt
#
gcc gasket_points.o -lm -framework GLUT -framework OpenGL
# gcc gasket_points.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gasket_points.o"
  exit
fi
#
rm gasket_points.o
mv a.out ~/binc/$ARCH/gasket_points
#
echo "Executable installed as ~/binc/$ARCH/gasket_points"
