#!/bin/tcsh -f
#
#  Strip "unknown branch" error lines from RA2bZinvAnalysis output
#

if ( $#argv < 1 ) then
  echo "Usage: tidyOut <output file name>"
  exit 1
endif

set OUTPUT = $1
set TIDYOUTPUT = ${OUTPUT}b
grep -v "unknown branch" ${OUTPUT} > ${TIDYOUTPUT}
touch -r ${OUTPUT} ${TIDYOUTPUT}
mv -f ${TIDYOUTPUT} ${OUTPUT}
