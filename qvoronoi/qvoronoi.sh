#!/bin/bash
#
gcc -c qvoronoi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qvoronoi.c."
  exit
fi
rm compiler.txt
#
gcc -c geom.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geom.c."
  exit
fi
rm compiler.txt
#
gcc -c geom2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geom2.c."
  exit
fi
rm compiler.txt
#
gcc -c global.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling global.c."
  exit
fi
rm compiler.txt
#
gcc -c io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling io.c."
  exit
fi
rm compiler.txt
#
gcc -c mem.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mem.c."
  exit
fi
rm compiler.txt
#
gcc -c merge.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling merge.c."
  exit
fi
rm compiler.txt
#
gcc -c poly.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poly.c."
  exit
fi
rm compiler.txt
#
gcc -c poly2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poly2.c."
  exit
fi
rm compiler.txt
#
gcc -c qhull.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qhull.c."
  exit
fi
rm compiler.txt
#
gcc -c qset.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qset.c."
  exit
fi
rm compiler.txt
#
gcc -c stat.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stat.c."
  exit
fi
rm compiler.txt
#
gcc -c user.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling user.c."
  exit
fi
rm compiler.txt
#
gcc qvoronoi.o geom.o geom2.o global.o io.o mem.o merge.o poly.o poly2.o qhull.o qset.o stat.o user.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading the object files."
  exit
fi
#
rm *.o
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/qvoronoi
#
echo "Executable installed as ~/binc/$ARCH/qvoronoi"
