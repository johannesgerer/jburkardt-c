#!/bin/bash
#
gcc -c hello.c
if [ $? -ne 0 ]; then
  echo "Errors compiling the program."
  exit
fi
#
#  This is how to access the OpenGL and GLUT libraries under Mac OSX.
#
gcc hello.o -framework GLUT -framework OpenGL
#
#  This is how to access the OpenGL and GLUT libraries under a standard installation.
#
#gcc hello.o -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the program."
  exit
fi
#
rm hello.o
mv a.out hello
echo "Executable installed as hello"
