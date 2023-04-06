#!/bin/bash

sframe_batch.py -s Selection_TstarTstar_2016.xml
sframe_batch.py -s Selection_TstarTstar_2017.xml
sframe_batch.py -s Selection_TstarTstar_2018.xml

sframe_batch.py -s Selection_SM_2016.xml
sframe_batch.py -s Selection_SM_2017.xml
sframe_batch.py -s Selection_SM_2018.xml

sframe_batch.py -s Selection_addSM_2016.xml
sframe_batch.py -s Selection_addSM_2017.xml
sframe_batch.py -s Selection_addSM_2018.xml

sframe_batch.py -s Selection_Data_2016.xml
sframe_batch.py -s Selection_Data_2017.xml
sframe_batch.py -s Selection_Data_2018.xml
