#!/bin/bash
#
gcc -c lights_out_opengl.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lights_out_opengl.c"
  exit
fi
rm compiler.txt
#
gcc lights_out_opengl.o -framework OpenGL -framework GLUT
if [ $? -ne 0 ]; then
  echo "Errors linking lights_out_opengl.o"
  exit
fi
#
rm lights_out_opengl.o
mv a.out ~/binc/$ARCH/lights_out_opengl
#
echo "Executable installed as ~/binc/$ARCH/lights_out_opengl"
