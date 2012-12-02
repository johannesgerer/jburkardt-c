#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem2_output.txt ]; then
  rm problem2_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 2 1 5 123456789 2 >  problem2_output.txt
smolpack_interactive 2 2 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 3 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 4 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 5 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 6 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 7 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 8 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 9 5 123456789 2 >> problem2_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 2 10 1 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 2 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 3 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 4 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 5 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 6 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 7 123456789 1 >> problem2_output.txt
smolpack_interactive 2 10 8 123456789 1 >> problem2_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 2 10 1 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 2 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 3 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 4 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 5 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 6 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 7 123456789 2 >> problem2_output.txt
smolpack_interactive 2 10 8 123456789 2 >> problem2_output.txt
#
echo "smolpack_interactive run on problem2."
