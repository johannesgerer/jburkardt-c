#!/bin/bash
#
#  Get the graph file.
#
cp ../../data/metis_graph/4elt.graph .
#
#  Run the program, requesting a partition into 8 parts.
#
kmetis 4elt.graph 8 > kmetis_prb_output.txt
#
#  Get rid of the graph file.
#
rm 4elt.graph
#
echo "Program output written to kmetis_prb_output.txt"
