#!/bin/bash
#
#  Get the mesh file.
#
cp ../../data/metis_mesh/metis.mesh .
#
#  Run the program, requesting a partition into 8 parts.
#
partdmesh metis.mesh 8 > partdmesh_prb_output.txt
#
#  Get rid of the mesh file.
#
rm metis.mesh
#
echo "Program output written to partdmesh_prb_output.txt"
