#!/usr/bin/env python3

# this very nice script was adapted from Christopher's version

import os
import sys
import csv
from collections import OrderedDict
import argparse
from itertools import permutations
from termcolor import colored
import shutil

# to include CrossSectionHelper:
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets'))
from CrossSectionHelper import MCSampleValuesHelper
helper = MCSampleValuesHelper()

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/TstarTstar/config'))
from constants import _YEARS


class configContainer:

   '''Container class for configurating XML files'''

   userName = str()
   uhh2Dir = str()
   userMail = str()
   yearVars = dict()
   used_samples = OrderedDict()

   def __init__(self):

      self.userName = os.environ.get('USER')
      self.uhh2Dir = os.environ.get('CMSSW_BASE')+'/src/UHH2/'
      self.userMail = '@'.join([self.userName, 'mail.desy.de']) # Avoid spam due to public code on GitHub

      self.yearVars['yearVersions'] = { # e.g. 'v3' in case  of 2016v3
         'UL16preVFP': '',
         'UL16postVFP': '',
         'UL17': '',
         'UL18': '',
      }

      # Set these values such that there are no more than 2,500 jobs per Preselection. This way, you can submit two Preselections in parallel to avoid going over 5,000 jobs (current user limit for NAF)
      self.yearVars['preselFileSplit'] = {
         'UL16preVFP': '50',
         'UL16postVFP': '50',
         'UL17': '50',
         'UL18': '50',
      }

      self.yearVars['targetLumis'] = {
         'UL16preVFP': _YEARS['UL16preVFP'].get('lumi_pb'),
         'UL16postVFP': _YEARS['UL16postVFP'].get('lumi_pb'),
         'UL17': _YEARS['UL17'].get('lumi_pb'),
         'UL18': _YEARS['UL18'].get('lumi_pb'),
      }

      self.yearVars['lumiFiles'] = {
         'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP_normtag.root',
         'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16postVFP_normtag.root',
         'UL17': self.uhh2Dir+'common/UHH2-data/UL17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root',
         'UL18': self.uhh2Dir+'common/UHH2-data/UL18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root',
      }

      self.yearVars['pileupFiles'] = {
         'mc': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyMCPileupHistogram_UL16preVFP.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyMCPileupHistogram_UL16postVFP.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyMCPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyMCPileupHistogram_UL18.root',
         },
         'data': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyDataPileupHistogram_UL16preVFP.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyDataPileupHistogram_UL16postVFP.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18.root',
         },
         'dataUp': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyDataPileupHistogram_UL16preVFP_72383.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyDataPileupHistogram_UL16postVFP_72383.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17_72383.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18_72383.root',
         },
         'dataDown': {
            'UL16preVFP': self.uhh2Dir+'common/UHH2-data/UL16preVFP/MyDataPileupHistogram_UL16preVFP_66017.root',
            'UL16postVFP': self.uhh2Dir+'common/UHH2-data/UL16postVFP/MyDataPileupHistogram_UL16postVFP_66017.root',
            'UL17': self.uhh2Dir+'common/UHH2-data/UL17/MyDataPileupHistogram_UL17_66017.root',
            'UL18': self.uhh2Dir+'common/UHH2-data/UL18/MyDataPileupHistogram_UL18_66017.root',
         },
      }

      self.yearVars['HOTVRSFs'] = {
         'UL16preVFP': '',
         'UL16postVFP': '',
         'UL17': self.uhh2Dir+'HOTVR/data/TopTaggingScaleFactors_RunIISummer19UL17_PUPPIv15.root',
         'UL18': self.uhh2Dir+'HOTVR/data/TopTaggingScaleFactors_RunIISummer19UL18_PUPPIv15.root',
      }

      # everything will be passed through until the end, so just put the new ones!
      self.additionalBranches = {
         'Preselection': "",
         'Selection': "trigger_decision trigger_decision_ele MC_isfake2017B MC_isfake2016B is_muevt evt_weight is_highpt weight_pu weight_pu_up weight_pu_down prefiringWeight prefiringWeightDown prefiringWeightUp",
         'Analysis': "neutrino is_btagevent weight_sfmu_id weight_sfmu_id_down weight_sfmu_id_up weight_sfmu_iso weight_sfmu_iso_down weight_sfmu_iso_up weight_sfmu_trigger weight_sfmu_trigger_up weight_sfmu_trigger_down weight_sfelec_id weight_sfelec_id_down weight_sfelec_id_up weight_sfelec_reco weight_sfelec_reco_down weight_sfelec_reco_up weight_sfelec_trigger weight_sfelec_trigger_up weight_sfelec_trigger_down weight_btagdisc_central weight_btagdisc_lf_up weight_btagdisc_lf_down weight_btagdisc_hf_up weight_btagdisc_hf_down weight_btagdisc_hfstats1_up weight_btagdisc_hfstats1_down weight_btagdisc_hfstats2_up weight_btagdisc_hfstats2_down weight_btagdisc_lfstats1_up weight_btagdisc_lfstats1_down weight_btagdisc_lfstats2_up weight_btagdisc_lfstats2_down weight_btagdisc_cferr1_up weight_btagdisc_cferr1_down weight_btagdisc_cferr2_up weight_btagdisc_cferr2_down  ST_AK4 ST_HOTVR weight_murmuf_upup weight_murmuf_upnone weight_murmuf_noneup weight_murmuf_nonedown weight_murmuf_downnone weight_murmuf_downdown",
         'DNN': "ST_weight DNN_Inputs TstarTstar_Hyp_gHOTVR TstarTstar_Hyp_gAK4",
         'DNN_datadriven': "",
      }

      self.additionalBranches["Selection"] = self.additionalBranches["Selection"]+" "+self.additionalBranches["Preselection"]
      self.additionalBranches["Analysis"] = self.additionalBranches["Analysis"]+" "+self.additionalBranches["Selection"]
      self.additionalBranches["DNN"] = self.additionalBranches["DNN"]+" "+self.additionalBranches["Analysis"]
      self.additionalBranches["DNN_datadriven"] = self.additionalBranches["DNN_datadriven"]+" "+self.additionalBranches["DNN"]
      self.additionalBranches["DNN_datadriven_variation"] = self.additionalBranches["DNN_datadriven"]+" "+self.additionalBranches["DNN"]

      self.additionalBranches["TriggerEff"] = self.additionalBranches["Selection"]

      self.systematics = list()

   @staticmethod
   def read_database_106X_v2(years: list, groups: list, triggerEff = False):

      if(triggerEff): from sample_db_triggerEff import samplesDict
      else: from sample_db import samplesDict
      used_samples = OrderedDict()
      for year in years:
          temp = OrderedDict()
          for group in groups:
             temp[group] = list()
             for k, v in samplesDict.items():
                use_me_year = v.get('years') == None or year in v.get('years', [])

                # block to group config files
                use_me_group = False
                if k.startswith("DATA"):
                    if (group == "DATA"): use_me_group = True
                elif k.startswith("TstarTstarToTgammaTgamma"):
                    if group == "Signal": use_me_group = True
                elif k.startswith("TstarTstarToTgluonTgluon_Spin32"):
                    if group == "Signal_Spin32": use_me_group = True
                else:
                    if group == "SM": use_me_group = True
                # end grouping block

                if (use_me_year and use_me_group):

                   sample_entity = sampleEntity((k, v, year,))
                   if not os.path.isfile(sample_entity.xmlPath):
                      # mabe this was split manually?
                      print(colored('XML for sample  '+sample_entity.nickName+' ('+year+')  does not exist. Skipping this sample', 'red'))
                      continue
                   temp[group].append(sample_entity)
          used_samples[year] = temp
      configContainer.used_samples = used_samples

class sampleEntity:

   '''Container to hold information about a data or MC sample, as read from CSV database'''

   def __init__(self, csvRow):

      k, v, year = csvRow
      self.year = year
      self.is_data = k.startswith('DATA_')
      self.nickName = k
      if(self.is_data): self.nickName = k[5:] # removing the DATA_ prefix
      self.n_das = None
      self.n_pnfs = None
      self.pnfs_sum_of_weights = helper.get_nevt(v['db_name'], '13TeV', self.year)
      self.xs_gen = None
      self.xs_theo = None
      self.lumi = 1. if self.is_data else helper.get_lumi(v['db_name'], '13TeV', self.year, kFactor=v.get('kfac', False), Corrections=v.get('corr', False))
      self.xsection = self.pnfs_sum_of_weights / self.lumi
      self.xmlPath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets', helper.get_xml(v['db_name'], '13TeV', self.year))


class xmlCreator:

   '''Creates XML files for SFrame'''

   def __init__(self, step: str, year: str, group: str, variationDirections, isTriggerEff):

      confCon = configContainer()
      self.uhh2Dir = confCon.uhh2Dir
      self.userName = confCon.userName
      self.userMail = confCon.userMail
      self.yearVars = confCon.yearVars
      self.sample_list = confCon.used_samples[year][group]
      self.systematics = confCon.systematics
      self.additionalBranches = confCon.additionalBranches
      self.jersmear_direction = variationDirections["jersmear"]
      self.jecsmear_direction = variationDirections["jecsmear"]
      self.isTriggerEff = isTriggerEff

      if step not in ['Preselection', 'Selection', 'Analysis', 'DNN', "DNN_datadriven", "DNN_datadriven_variation", "TriggerEff"]:
         sys.exit('Given value of argument "selection" not valid. Abort.')
      self.step = step
      if (step == "Preselection"):
          self.previousFolder = "none"
          self.analysisModule = "TstarTstarPreselectionModule"
      elif (step == "Selection"):
          self.previousFolder = "Preselection"
          self.analysisModule = "TstarTstarSelectionModule"
      elif (step == "Analysis"):
          self.previousFolder = "Selection"
          self.analysisModule = "TstarTstarAnalysisModule"
      elif (step == "DNN" or step == "DNN_datadriven" or step == "DNN_datadriven_variation"):
          self.previousFolder = "Analysis"
          self.analysisModule = "TstarTstarDNNModule"
      elif (step == "TriggerEff"):
          self.previousFolder = "Selection_TriggerEff"
          self.analysisModule = "TstarTstarTriggerSFModule"
      if(self.isTriggerEff and not step == "TriggerEff"): self.previousFolder += "_TriggerEff"

      if year not in ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']:
         sys.exit('Given value of argument "year" not valid. Abort.')
      self.year = year
      self.yearVersion = self.yearVars['yearVersions'][year]

      self.outputDirBase = '/nfs/dust/cms/user/flabe/TstarTstar/data/'
      if not os.path.isdir(self.outputDirBase):
         sys.exit('Warning: Make sure to create output directory via "ln -s". Abort.')
      self.outputDirJERJER = ''
      if (not self.jersmear_direction == "nominal"): self.outputDirJERJER = "/JER_"+self.jersmear_direction+"/"
      elif (not self.jecsmear_direction == "nominal"): self.outputDirJERJER = "/JEC_"+self.jecsmear_direction+"/"
      self.inputDirJERJER = ''
      if not (step == "DNN_datadriven" or step == "DNN_datadriven_variation"):
          if (not self.jersmear_direction == "nominal"): self.inputDirJERJER = "/JER_"+self.jersmear_direction+"/"
          elif (not self.jecsmear_direction == "nominal"): self.inputDirJERJER = "/JEC_"+self.jecsmear_direction+"/"

      self.xmlFileName = '_'.join(['parsedConfigFile', self.step, self.year, group])+'.xml'
      if(self.isTriggerEff): self.xmlFilePathBase = self.uhh2Dir+'TstarTstar/config/'+'_'.join(['config', self.step, "TriggerEff", self.year])+'/'
      else: self.xmlFilePathBase = self.uhh2Dir+'TstarTstar/config/'+'_'.join(['config', self.step, self.year])+'/'
      os.makedirs(self.xmlFilePathBase, exist_ok=True)
      self.xmlFilePath = self.xmlFilePathBase+self.xmlFileName
      if(self.isTriggerEff): self.workdirName = '_'.join(['workdir', self.step, "TriggerEff", self.year, group])
      else: self.workdirName = '_'.join(['workdir', self.step, self.year, group])

      shutil.copyfile(self.uhh2Dir+"/TstarTstar/config/JobConfig.dtd", self.xmlFilePathBase+"/JobConfig.dtd")

      self.write_xml_successful = False


   def write_xml(self):

      with open(self.xmlFilePath, 'w') as file:
         file.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
         file.write('''\n''')
         file.write('''<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd"[\n''')
         file.write('''\n''')
         file.write('''<!ENTITY TargetLumi "'''+str(self.yearVars['targetLumis'][self.year])+'''">\n''')
         if not (self.step == "Preselection"):
            if not (self.step == "Selection"): file.write('''<!ENTITY INPUTdir "'''+(self.outputDirBase+self.previousFolder+'/'+self.year+self.inputDirJERJER)+'''">\n''')
            else: file.write('''<!ENTITY INPUTdir "'''+(self.outputDirBase+self.previousFolder+'/'+self.year)+'''">\n''')
            file.write('''<!ENTITY INPUTfilename "uhh2.AnalysisModuleRunner">\n''')
         if(self.isTriggerEff):file.write('''<!ENTITY OUTPUTdir "'''+(self.outputDirBase+self.step+'_TriggerEff/'+self.year)+'''">\n''')
         else: file.write('''<!ENTITY OUTPUTdir "'''+(self.outputDirBase+self.step+'/'+self.year+self.outputDirJERJER)+'''">\n''')
         file.write('''<!ENTITY b_Cacheable "False">\n''')
         file.write('''<!ENTITY NEVT "-1">\n''')
         file.write('''<!ENTITY YEARsuffix "_'''+self.year+self.yearVersion+'''">\n''')
         file.write('''<!ENTITY PROOFdir "/nfs/dust/cms/user/'''+self.userName+'''/.proof2">\n''')
         file.write('''\n''')
         for s in self.sample_list:
            if not (self.step == "Preselection"):
               file.write('''<!ENTITY '''+s.nickName+''' "&INPUTdir;/&INPUTfilename;'''+('.DATA.' if s.is_data else '.MC.')+s.nickName+'''&YEARsuffix;.root">\n''')
            else:
               file.write('''<!ENTITY '''+s.nickName+''' SYSTEM "'''+s.xmlPath+'''">\n''')
         file.write('''\n''')
         file.write(''']>\n''')
         file.write('''\n''')
         file.write('''<!--\n''')
         file.write('''<ConfigParse NEventsBreak="'''+('500000' if not (self.step == "Preselection") else '0')+'''" FileSplit="'''+('0' if not (self.step == "Preselection") else self.yearVars['preselFileSplit'][self.year])+'''" AutoResubmit="5"/>\n''')
         file.write('''<ConfigSGE RAM="4" DISK="3" Mail="'''+self.userMail+'''" Notification="as" Workdir="'''+self.workdirName+'''"/>\n''')
         file.write('''-->\n''')
         file.write('''\n''')
         file.write('''<!-- OutputLevel controls which messages are printed; set to VERBOSE or DEBUG for more verbosity, to WARNING or ERROR for less -->\n''')
         file.write('''<JobConfiguration JobName="ExampleCycleJob" OutputLevel="INFO">\n''')
         file.write('''<Library Name="libSUHH2TstarTstar"/>\n''')
         file.write('''<Package Name="SUHH2TstarTstar.par"/>\n''')
         file.write('''<Cycle Name="uhh2::AnalysisModuleRunner" OutputDirectory="&OUTPUTdir;/" PostFix="" TargetLumi="&TargetLumi;">\n''')
         file.write('''\n''')
         for s in self.sample_list:
             if(self.step == "Preselection"):
                 file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> &'''+s.nickName+'''; <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
             else:
                 file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> <In FileName="&'''+s.nickName+''';" Lumi="0.0"/> <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
         file.write('''\n''')
         file.write('''<UserConfig>\n''')
         file.write('''\n''')
         file.write('''<Item Name="PrimaryVertexCollection" Value="offlineSlimmedPrimaryVertices"/>\n''')
         file.write('''<Item Name="METName" Value="slimmedMETsPuppi"/>\n''')
         file.write('''<Item Name="ElectronCollection" Value="slimmedElectronsUSER"/>\n''')
         file.write('''<Item Name="PhotonCollection" Value="slimmedPhotonsUSER"/>\n''')
         file.write('''<Item Name="MuonCollection" Value="slimmedMuonsUSER"/>\n''')
         file.write('''<Item Name="JetCollection" Value="jetsAk4Puppi"/>\n''')
         file.write('''<Item Name="GenJetCollection" Value="slimmedGenJets"/>\n''')
         file.write('''<Item Name="TopJetCollection" Value="hotvrPuppi"/>\n''')
         file.write('''<Item Name="GenTopJetCollection" Value="hotvrGen"/>\n''')
         file.write('''<Item Name="GenParticleCollection" Value="GenParticles"/>\n''')
         file.write('''<Item Name="GenInfoName" Value="genInfo"/>\n''')
         file.write('''\n''')
         if not (self.additionalBranches[self.step] == ""):
             file.write('''<Item Name="additionalBranches" Value="'''+self.additionalBranches[self.step]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="lumi_file" Value="'''+self.yearVars['lumiFiles'][self.year]+'''"/>\n''')
         file.write('''<Item Name="lumihists_lumi_per_bin" Value="500."/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="pileup_directory" Value="'''+self.yearVars['pileupFiles']['mc'][self.year]+'''"/>\n''')
         file.write('''<Item Name="pileup_directory_data" Value="'''+self.yearVars['pileupFiles']['data'][self.year]+'''"/>\n''')
         file.write('''<Item Name="pileup_directory_data_up" Value="'''+self.yearVars['pileupFiles']['dataUp'][self.year]+'''"/>\n''')
         file.write('''<Item Name="pileup_directory_data_down" Value="'''+self.yearVars['pileupFiles']['dataDown'][self.year]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Tell AnalysisModuleRunner NOT to use the MC event weight from SFrame; rather let MCLumiWeight (called via CommonModules) calculate the MC event weight. The MC event weight assigned by MCLumiWeight is InputData.Lumi / Cycle.TargetLumi. -->\n''')
         file.write('''<Item Name="use_sframe_weight" Value="false"/>\n''')
         file.write('''<Item Name="AnalysisModule" Value="'''+self.analysisModule+'''"/>\n''')
         file.write('''<Item Name="uhh2Dir" Value="'''+self.uhh2Dir+'''"/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="jersmear_direction" Value="'''+ self.jersmear_direction+'''"/>\n''')
         file.write('''<Item Name="jecsmear_direction" Value="'''+ self.jecsmear_direction+'''"/>\n''')
         if (self.step == "Preselection"):
             file.write('''<Item Name="ScaleVariationMuR" Value = "central" />''')
             file.write('''<Item Name="ScaleVariationMuF" Value = "central" />''')
         elif(self.step == "Selection"):
             file.write('''<!-- scale factor configuration -->\n''')
             file.write('''<Item Name="HOTVRTopTagSFs" Value="'''+self.yearVars['HOTVRSFs'][self.year]+'''"/>\n''')
             file.write('''<Item Name="SF_path" Value="/nfs/dust/cms/user/flabe/TstarTstar/CMSSW_10_2_17/src/UHH2/TstarTstar/factors/" />\n''')
             file.write('''<Item Name="NLOCorrections" Value = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/sfs/nnlo" />''')
         elif(self.step == "DNN_datadriven"):
             file.write('''<Item Name="use_data_for" Value="background_extrapolation"/>\n''')
         elif(self.step == "DNN_datadriven_variation"):
             file.write('''<Item Name="use_data_for" Value="background_extrapolation_variation"/>\n''')
         if(self.isTriggerEff): file.write('''<Item Name="IsTriggerSFMeasurement" Value="True"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Switch for debugging of the central AnalysisModule -->\n''')
         file.write('''<Item Name="debug" Value="false"/>\n''')
         file.write('''\n''')
         file.write('''</UserConfig>\n''')
         file.write('''\n''')
         file.write('''</Cycle>\n''')
         file.write('''</JobConfiguration>\n''')

      self.write_xml_successful = True

      print('Created '+self.xmlFilePath)

      return self.xmlFilePath

   def write_bash_scripts(self):

       scriptFilePath_sframe_batch = os.path.join(self.xmlFilePathBase, 'run_all_sframe_batch.sh')
       with open(scriptFilePath_sframe_batch, 'w') as outfile:
          outfile.write('#!/bin/bash\n')
          newline_base = 'sframe_batch.py $1 '
          outfile.write(newline_base+self.xmlFilePath+'\n')
          for systXmlFilePath in self.systXmlFilePaths:
             outfile.write(newline_base+systXmlFilePath+'\n')
       print('Created '+scriptFilePath_sframe_batch)

       scriptFilePath_run_local = os.path.join(self.xmlFilePathBase, 'run_all_local.sh')
       with open(scriptFilePath_run_local, 'w') as outfile:
          outfile.write('#!/bin/bash\n')
          newline_base = 'python run_local.py $1 '
          outfile.write(newline_base+self.workdirName+'\n')
          for systWorkdirName in self.systWorkdirNames:
             outfile.write(newline_base+systWorkdirName+'\n')
       print('Created '+scriptFilePath_run_local)
       # copy_run_local_command = 'cp -n '+os.path.join(self.xmlFilePathBase, '..', 'run_local.py')+' '+os.path.join(self.xmlFilePathBase, '.')
       link_run_local_command = 'ln -s '+os.path.abspath(os.path.join(self.xmlFilePathBase, '..', 'run_local.py'))+' '+os.path.join(self.xmlFilePathBase, 'run_local.py')
       os.system(link_run_local_command)



if __name__=='__main__':

   selections = ['Preselection', 'Selection', 'Analysis', 'DNN', 'DNN_datadriven', 'DNN_datadriven_variation', 'TriggerEff']
   years = ['UL16preVFP','UL16postVFP','UL17', 'UL18']
   jersmear_directions = ['nominal', 'up', 'down']
   jecsmear_directions = ['nominal', 'up', 'down']

   if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
   parser = argparse.ArgumentParser()
   parser.add_argument('--all', action='store_true', help='Create XML files for all selections and years.')
   parser.add_argument('-s', '--selections', choices=selections, nargs='*', default=[])
   parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
   parser.add_argument('-a', '--auto-complete', action='store_true', help='Auto-complete arguments if not all arguments for selections and years are given.')
   parser.add_argument('--jer', choices=jersmear_directions, nargs=1, default = ["nominal"])
   parser.add_argument('--jec', choices=jecsmear_directions, nargs=1, default = ["nominal"])
   parser.add_argument('--triggerEff', action='store_true')

   args = parser.parse_args(sys.argv[1:])

   if(args.all == True):
      if(len(args.selections) + len(args.years) != 0):
         sys.exit('Not allowed to use "--all" option jointly with manually given selection or year argument. Exit.')
      args.selections = selections
      args.years = years
      if(args.auto_complete):
         print('Warning: You already specified "--all". Therefore, "--auto-complete" will not have any effect.')
   else:
      if args.auto_complete:
         print('Auto-completing arguments.')
         if not args.selections: args.selections = selections
         if not args.years: args.years = years
      else:
         for p in permutations([args.selections, args.years]):
            if p[0] and not p[1]:
               sys.exit('You specified arguments for at least one of the two options: "--selections" or "--years", but not for both of them. Also, you did not specify "--auto-complete" to compensate for this. Exit.')

   print('Going to create XML files for:')
   print('  Selections: '+', '.join(str(x) for x in args.selections))
   print('  Years: '+', '.join(str(x) for x in args.years))

   groups = ['DATA', 'Signal', 'Signal_Spin32', 'SM']

   configContainer.read_database_106X_v2(args.years, groups, args.triggerEff)

   # setting jec, jer directions
   directions = {}
   directions["jecsmear"] = args.jec[0]
   directions["jersmear"] = args.jer[0]

   for selection in args.selections:
      for year in args.years:
         if(selection == "DNN_datadriven" or selection == "DNN_datadriven_variation"):
            x = xmlCreator(selection, year, "DATA", directions, False)
            x.write_xml()
         else:
            for group in groups:
               x = xmlCreator(selection, year, group, directions, args.triggerEff)
               x.write_xml()
