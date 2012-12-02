#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem4_output.txt ]; then
  rm problem4_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 4 1 5 123456789 2 >  problem4_output.txt
smolpack_interactive 4 2 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 3 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 4 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 5 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 6 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 7 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 8 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 9 5 123456789 2 >> problem4_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 4 10 1 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 2 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 3 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 4 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 5 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 6 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 7 123456789 1 >> problem4_output.txt
smolpack_interactive 4 10 8 123456789 1 >> problem4_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 4 10 1 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 2 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 3 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 4 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 5 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 6 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 7 123456789 2 >> problem4_output.txt
smolpack_interactive 4 10 8 123456789 2 >> problem4_output.txt
#
echo "smolpack_interactive run on problem4."
