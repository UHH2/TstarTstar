#!/bin/bash

for year in UL18 UL17 UL16postVFP UL16preVFP
do
    for channel in mu ele
    do
        for region in SignalRegion ValidationRegion
        do
            root -l -b -q combinePDFRMSforST.C\(\"$year\",\"$channel\",\"$region\"\)
        done
    done
done