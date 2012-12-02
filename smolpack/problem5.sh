#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem5_output.txt ]; then
  rm problem5_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 5 1 5 123456789 2 >  problem5_output.txt
smolpack_interactive 5 2 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 3 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 4 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 5 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 6 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 7 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 8 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 9 5 123456789 2 >> problem5_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 5 10 1 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 2 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 3 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 4 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 5 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 6 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 7 123456789 1 >> problem5_output.txt
smolpack_interactive 5 10 8 123456789 1 >> problem5_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 5 10 1 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 2 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 3 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 4 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 5 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 6 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 7 123456789 2 >> problem5_output.txt
smolpack_interactive 5 10 8 123456789 2 >> problem5_output.txt
#
echo "smolpack_interactive run on problem5."
