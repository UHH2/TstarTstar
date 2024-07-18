#!/bin/bash

for channel in mu ele
do
    for region in SR VR
    do
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\"\) > log_bgest_${channel}_${region}.txt
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\",\"btagging_totalUp\"\) >> log_bgest_${channel}_${region}.txt
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\",\"btagging_totalDown\"\) >> log_bgest_${channel}_${region}.txt
    done
done