#!/bin/bash
#
#  Get the graph file.
#
cp ../../data/metis_graph/4elt.graph .
#
#  Run the program, requesting a reordering.
#
onmetis 4elt.graph > onmetis_prb_output.txt
#
#  Get rid of the graph file.
#
rm 4elt.graph
#
echo "Program output written to onmetis_prb_output.txt"
