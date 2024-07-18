#!/usr/bin/env python3

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

   '''Container class config stuff'''

   used_samples = OrderedDict()

   def __init__(self):
       self.uhh2Dir = os.environ.get('CMSSW_BASE')+'/src/UHH2/'


   @staticmethod
   def read_database_106X_v2(years: list, groups: list):

      from sample_db import samplesDict
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
                elif k.startswith("TstarTstar"):
                    if group == "Signal": use_me_group = True
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
      self.nickName = k
      self.pnfs_sum_of_weights = helper.get_nevt(v['db_name'], '13TeV', self.year)
      self.xs = helper.get_xs(v['db_name'], '13TeV', self.year)
      self.xs_source = helper.get_xs(v['db_name'], '13TeV', self.year, "Source")
      self.br = helper.get_br(v['db_name'], '13TeV', self.year)
      self.br_source = helper.get_br(v['db_name'], '13TeV', self.year, "Source")
      self.lumi = helper.get_lumi(v['db_name'], '13TeV', self.year, kFactor=v.get('kfac', False), Corrections=v.get('corr', False))
      self.xsection = self.pnfs_sum_of_weights / self.lumi
      self.xsprint = self.xs * self.br
      self.xmlPath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets', helper.get_xml(v['db_name'], '13TeV', self.year))
      self.xmlSource = helper.get_xml(v['db_name'], '13TeV', self.year, "Source")

class tableCreator:

   def __init__(self, year: str, group: str):

      confCon = configContainer()
      self.uhh2Dir = confCon.uhh2Dir
      self.sample_list = confCon.used_samples[year][group]

      if year not in ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']:
         sys.exit('Given value of argument "year" not valid. Abort.')
      self.year = year

      self.xmlFileName = '_'.join(['table', self.year, group])+'.txt'
      self.xmlFilePathBase = self.uhh2Dir+'TstarTstar/config/LaTeX/'
      os.makedirs(self.xmlFilePathBase, exist_ok=True)
      self.xmlFilePath = self.xmlFilePathBase+self.xmlFileName

      self.write_table_successful = False


   def write_table(self):

      with open(self.xmlFilePath, 'w') as file:

          file.write('''\\begin{center}\n''')
          file.write('''\\begin{tabular}{ | l | c | c | c |}\n''')
          file.write('''\\hline\n''')
          file.write('''sample & $N_{\\mathrm{evt}}$ & $\\sigma [\\mathrm{pb}]$ & references  \\\\ \n''')
          file.write('''\\hline\n''')
          campaign = ""
          for s in self.sample_list:
              file.write('''\\texttt{''')
              file.write(s.xmlSource.split("/")[1].replace("_", "\\_"))
              campaign_new = s.xmlSource.split("/")[2].replace("_", "\\_")
              if(campaign_new[-2] == "v"): campaign_new = campaign_new[:-1]+"*"
              if(campaign == ""): campaign = campaign_new
              elif(not campaign == campaign_new):
                  print("ERROR: non-constant campaign found!")
                  print("Had: "+campaign)
                  print("Found: "+s.xmlSource.split("/")[2].replace("_", "\\_"))
              file.write('''} & $'''+str(s.pnfs_sum_of_weights)+'''$ & $'''+str(s.xsprint)+'''$ & ''')
              if( not s.xs_source == ""):
                  if(s.xs_source[:4] == "http"): file.write('''$\\sigma$ taken from \\href{'''+s.xs_source+'''}{here}.''')
                  else: file.write('''$\\sigma$ obtained using '''+s.xs_source+'''.''')
              else: file.write(''' - ''')
              file.write(''' \\\\ \n''')
          file.write('''\\hline\n''')
          file.write('''\\end{tabular}\n''')
          file.write('''\\captionof{table}{Overview of MC background samples for the '''+self.year+''' era. ''')
          file.write('''All samples are part of the campaign \\texttt{'''+campaign+'''} and are taken from MINIAOD.}\n''')
          file.write('''\\label{tab:MCbkg_'''+self.year+'''}\n''')
          file.write('''\\end{center}\n''')

      self.write_table_successful = True

      print('Created '+self.xmlFilePath)

      return self.xmlFilePath

if __name__=='__main__':

   years = ['UL16preVFP','UL16postVFP','UL17', 'UL18']

   parser = argparse.ArgumentParser()
   parser.add_argument('-y', '--years', choices=years, nargs='*', default=years)

   args = parser.parse_args(sys.argv[1:])

   print('Going to create LaTeX tables for:')
   print('  Years: '+', '.join(str(x) for x in args.years))

   configContainer.read_database_106X_v2(args.years, ["SM"])

   for year in args.years:
       x = tableCreator(year, "SM")
       x.write_table()
