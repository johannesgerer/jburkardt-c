#!/bin/bash
#
#  Get the mesh file.
#
cp ../../data/metis_mesh/metis.mesh .
#
#  Run the program.
#
./mesh2dual metis.mesh > mesh2dual_prb_output.txt
#
#  Get rid of the mesh file.
#
rm metis.mesh
#
echo "Program output written to mesh2dual_prb_output.txt"
