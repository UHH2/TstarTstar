#!/bin/bash

root -l -q cutflowPlot.C
root -l -q cutflowPlot.C\(\"_mu\"\)
root -l -q cutflowPlot.C\(\"_mu\_lowpt\"\)
root -l -q cutflowPlot.C\(\"_mu\_highpt\"\)
root -l -q cutflowPlot.C\(\"_ele\"\)
root -l -q cutflowPlot.C\(\"_ele\_lowpt\"\)
root -l -q cutflowPlot.C\(\"_ele\_highpt\"\)
