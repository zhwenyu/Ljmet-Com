#!/bin/bash
#
# Template of the c-shell script for submitting a CONDOR job
#
# Need to be submitted as arguments from condor .jdl file:
#    - CMSSWBASE       - local release base ($CMSSW_BASE)
#    - OUTPUT_DIR      - condor work and output dir
#    - PREFIX          - some generic name for the set of jobs (like ttbar, wjets)
#    - CONDOR_PROCESS  - condor job number (ranges from 0 to N-1 jobs)
#    - CONDOR_CLUSTER  - condor cluster (for recordkeeping)
# Arguments = INPUTTAR OUTPUT_DIR PREFIX $(process)
INPUTTAR=$1
OUTPUT_DIR=$2
PREFIX=$3
JOBID=$4

tar -xf ${INPUTTAR}.tar
cd ${INPUTTAR}
scram b ProjectRename

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -

sed -i "s|CONDOR_RELBASE|$PWD/$INPUTTAR|" ${PREFIX}_${JOBID}.py
ljmet ${PREFIX}_${JOBID}.py

# copy output to eos
echo "xrdcp .root output for condor"
for FILE in *.root
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTPUT_DIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done
