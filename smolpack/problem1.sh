#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem1_output.txt ]; then
  rm problem1_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 1 1 5 123456789 2 >  problem1_output.txt
smolpack_interactive 1 2 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 3 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 4 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 5 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 6 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 7 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 8 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 9 5 123456789 2 >> problem1_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 1 10 1 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 2 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 3 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 4 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 5 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 6 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 7 123456789 1 >> problem1_output.txt
smolpack_interactive 1 10 8 123456789 1 >> problem1_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 1 10 1 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 2 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 3 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 4 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 5 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 6 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 7 123456789 2 >> problem1_output.txt
smolpack_interactive 1 10 8 123456789 2 >> problem1_output.txt
#
echo "smolpack_interactive run on problem1."
