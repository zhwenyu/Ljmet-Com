#!/bin/bash

START=$(date +%s.%N)
echo ""
echo "running testData.py"
ljmet testData.py

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
msg=" seconds has elapsed"
out=$DIFF$msg
echo $out
