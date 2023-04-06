#!/bin/bash
#loop over samples
denomstep=After2D_mu
for sample in TstarTstar\_M-700 TstarTstar\_M-800 TstarTstar\_M-900 TstarTstar\_M-1000 TstarTstar\_M-1100 TstarTstar\_M-1200 TstarTstar\_M-1300 TstarTstar\_M-1400 TstarTstar\_M-1500 TstarTstar\_M-1600
do
    #loop over triggers
#    for step in SemiLepTTBarMatchGENRECO_triggerSingleLeptonEle SemiLepTTBarMatchGENRECO_triggerSingleJet_ele SemiLepTTBarMatchGENRECO_triggerHT_ele SemiLepTTBarMatchGENRECO_triggerPFHT_ele
    for step in triggerSingleLeptonMu
    do
        #loop over hists
  for hist in pt_mu eta_mu pt_ele eta_ele pt_jet1 pt_ak8jet1
	do
	    root -l -b -q TriggerEffPlots.C\(\"uhh2.AnalysisModuleRunner.MC.${sample}.root\",\"${sample}\",\"${step}\",\"${hist}\",\"${denomstep}\"\)
	done
    done
done
