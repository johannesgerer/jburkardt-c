#!/bin/bash
#
gcc peri1d.c pdblas.c -lblas
mv a.out peri1d
echo "The executable program has been created as peri1d."

