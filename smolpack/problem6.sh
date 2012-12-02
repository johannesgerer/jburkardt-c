#!/bin/bash
#
#  Delete old output file.
#
if [ -f problem6_output.txt ]; then
  rm problem6_output.txt
fi
#
#  Run in dimensions 1:9 using standard Clenshaw Curtis rule.
#
smolpack_interactive 6 1 5 123456789 2 >  problem6_output.txt
smolpack_interactive 6 2 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 3 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 4 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 5 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 6 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 7 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 8 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 9 5 123456789 2 >> problem6_output.txt
#
#  Run in dimension 10 using delayed Clenshaw Curtis rule.
#
smolpack_interactive 6 10 1 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 2 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 3 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 4 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 5 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 6 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 7 123456789 1 >> problem6_output.txt
smolpack_interactive 6 10 8 123456789 1 >> problem6_output.txt
#
#  Run in dimension 10 using standard Clenshaw Curtis rule.
#
smolpack_interactive 6 10 1 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 2 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 3 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 4 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 5 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 6 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 7 123456789 2 >> problem6_output.txt
smolpack_interactive 6 10 8 123456789 2 >> problem6_output.txt
#
echo "smolpack_interactive run on problem6."
