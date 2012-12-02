#!/bin/bash
#
gcc -c planet.c
if [ $? -ne 0 ]; then
  echo "Errors compiling the program."
  exit
fi
#
#  This is how to access the OpenGL and GLUT libraries under Mac OSX.
#
gcc planet.o -framework GLUT -framework OpenGL
#
#  This is how to access the OpenGL and GLUT libraries under a standard installation.
#
#gcc planet.o -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the program."
  exit
fi
#
rm planet.o
mv a.out planet
echo "Executable installed as planet"
