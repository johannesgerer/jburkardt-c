#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem7_output.txt ]; then
  rm problem7_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 7 1 5 123456789 2 >  problem7_output.txt
smolpack_interactive 7 2 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 3 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 4 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 5 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 6 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 7 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 8 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 9 5 123456789 2 >> problem7_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 7 10 1 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 2 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 3 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 4 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 5 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 6 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 7 123456789 1 >> problem7_output.txt
smolpack_interactive 7 10 8 123456789 1 >> problem7_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 7 10 1 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 2 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 3 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 4 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 5 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 6 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 7 123456789 2 >> problem7_output.txt
smolpack_interactive 7 10 8 123456789 2 >> problem7_output.txt
#
echo "smolpack_interactive run on problem7."
