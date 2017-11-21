#!/bin/bash


# testfile=testMC.py

testfile=testMC_singlelep.py
# testfile=testData_singlelep.py

START=$(date +%s.%N)
echo ""

echo "--------------------"
echo "Running $testfile"
echo "--------------------"
echo ""

ljmet $testfile

END=$(date +%s.%N)
DIFF=$(echo "($END - $START)" | bc)
DIFFmin=$(echo "($END - $START)/60." | bc)
DIFFhr100=$(echo "($END - $START)/6./6. " | bc)
msg=" seconds has elapsed"
msgmin=" minutes has elapsed"
msghr=" hr has elapsed for 100x more events"
out=$DIFF$msg
outmin=$DIFFmin$msgmin
outhr100=$DIFFhr100$msghr
echo $out
echo $outmin
echo $outhr100

