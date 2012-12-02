#!/bin/bash
#
gcc -c teapot.c
if [ $? -ne 0 ]; then
  echo "Errors compiling the program."
  exit
fi
#
#  This is how to access the OpenGL and GLUT libraries under Mac OSX.
#
gcc teapot.o -framework GLUT -framework OpenGL
#
#  This is how to access the OpenGL and GLUT libraries under a standard installation.
#
#gcc teapot.o -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the program."
  exit
fi
#
rm teapot.o
mv a.out teapot
echo "Executable installed as teapot"
