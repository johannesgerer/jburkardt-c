#!/bin/bash
#
#  Get the graph file.
#
cp ../../data/metis_graph/4elt.graph .
#
#  Have the program check the file.
#
graphchk 4elt.graph > graphchk_prb_output.txt
#
#  Get rid of the graph file.
#
rm 4elt.graph
#
echo "Program output written to graphchk_prb_output.txt"
