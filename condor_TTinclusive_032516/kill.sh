#!/bin/tcsh

foreach i (`seq 746 766`)
    condor_rm 546$i
end
#EOF
