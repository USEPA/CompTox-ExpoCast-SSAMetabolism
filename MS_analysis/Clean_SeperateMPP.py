# -*- coding: utf-8 -*-
"""
Created on Wed May 19 13:01:26 2021

@author: MBOYCE
"""

import os
import pandas as pd
import numpy as np
import csv
import glob
    
sep_Value = 6
conditions = ['Pellet', 'Super', 'Gluc']
hours = [0, 1, 4]
sampleTypes = ['Blank', 'Sample']
replicates = 3

def CleanAndSeperate(directory):
    save_dir = directory+'\\WebApp_Input'
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    
    for file in glob.glob(directory + r'\MPP\*.txt'):
        print(file)
        fileName = file.rsplit(r'\\', 1)[-1]
        if 'pos' in fileName:
            mode = 'pos'
        else:
            mode = 'neg'

        df = pd.read_csv(file, skiprows = 4, sep = '\t')
        df['Annotations'] = df['Annotations'].apply(modify_annotation)
        df['Compound'] = [x.replace(' ','') for x in df['Compound']]
        df['Compound'] = [x+'_mfg' if 'mfg=' in y else x for x,y in zip(df['Compound'], df['Annotations'])]
        df = df[(df['Formula'] != '> limit') & (df['Formula'] != '<none>')]
            
        dataDF, infoDF = seperateDF(df)
        dataDF = renameCols(dataDF)
        combineSave(infoDF, dataDF, mode, save_dir)

def modify_annotation(x):
    try:
        split_x = x.split(", tgt=", 2)[:2]
        new_x = split_x[0] + ", tgt=" + split_x[1] + " ]"
    except:
        new_x = np.nan
    return (new_x)

def seperateDF(df):
    dataColumns = [col for col in df.columns if 'raw' in col]
    infoColumns = [col for col in df.columns if 'raw' not in col]
    return df[dataColumns], df[infoColumns]

def renameCols(df):
    col = []
    for condition in conditions:
        for hour in hours:
            for sampleType in sampleTypes:
                for repNum in range(replicates):
                    if(sampleType == 'Blank'):
                        col += [sampleType+condition[0]+str(hour)+str(repNum)]
                    else:
                        col += [condition+str(hour)+str(repNum)]
    df.columns = col
    return df

def combineSave(infoDF, dataDF, mode, save_directory):
    name_dict = {}
    idx = 0
    for condition in conditions:
        for hour in hours:
            name_dict[idx] = condition+"_"+str(hour)+"_"+mode
            idx = idx + 1
            
    sampleCount = len(dataDF.columns)
    for block in range(int(sampleCount/sep_Value)):
        lowerIdx = block * sep_Value
        upperIdx = lowerIdx + sep_Value
        cols = dataDF.columns[lowerIdx:upperIdx]
        tmpDF = pd.concat([infoDF, dataDF[cols]], axis = 1)
        columnOrder = ['Compound']+dataDF[cols].columns.tolist()+['Alignment Value','Annotations','Compound Name','CompoundAlgo','Formula','Frequency','Ionization mode','Mass','MS1 Composite Spectrum', 'Retention Time','Score (DB)']
        tmpDF[columnOrder].to_csv(save_directory + '\\' + name_dict[block]+'.csv', index = False)
        
########################################################################