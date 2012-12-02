#!/bin/bash
#
#  Get the mesh file.
#
cp ../../data/metis_mesh/metis.mesh .
#
#  Run the program, requesting a partition into 8 parts.
#
partnmesh metis.mesh 8 > partnmesh_prb_output.txt
#
#  Get rid of the mesh file.
#
rm metis.mesh
#
echo "Program output written to partnmesh_prb_output.txt"
