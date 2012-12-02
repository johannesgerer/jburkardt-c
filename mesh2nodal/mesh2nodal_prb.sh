#!/bin/bash
#
#  Get the mesh file.
#
cp ../../data/metis_mesh/metis.mesh .
#
#  Run the program.
#
mesh2nodal metis.mesh > mesh2nodal_prb_output.txt
#
#  Get rid of the mesh file.
#
rm metis.mesh
#
echo "Program output written to mesh2nodal_prb_output.txt"
