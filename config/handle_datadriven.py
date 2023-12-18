import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--option', required=True)
args = parser.parse_args()

years = ["UL16preVFP", "UL16postVFP", "UL17", "UL18"]

for year in years:
    os.system("sframe_batch.py -" + args.option + "  config_DNN_datadriven_" + year + "/parsedConfigFile_DNN_datadriven_" + year + "_DATA.xml")