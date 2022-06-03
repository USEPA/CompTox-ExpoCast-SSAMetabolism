# -*- coding: utf-8 -*-
"""
Created on Tue May 01 12:34:34 2018

@author: HALGhoul
"""

import time
import os
import glob
import pandas as pd
import mgf_parser_v24_AC_v2 as mg

def FeatureCompare(CFMID_Results_DF, Merged_Feature_DF):
    matched_features = []
    for mass in CFMID_Results_DF['MASS']:
        diffSum = Merged_Feature_DF[['Negative','Positive']].apply(lambda x: abs(x - mass)).sum(axis = 1)
        matched_features.append(Merged_Feature_DF.loc[diffSum == diffSum.min(), 'MPP'].values[0])
    return matched_features
    

def MergeCFMIDResults(rootdir, featureDF):
    mergedDF = pd.DataFrame()
    dataPaths = [rootdir + r'\Negative_Results', rootdir + r'\Positive_Results']
    for path in dataPaths:
        for filePath in glob.glob(os.path.join(path, '*.xlsx')):
            df = pd.read_excel(filePath)[['MASS','DTXCID','energy_sum','MASS_in_MGF']]
            df['Feature_MASS'] = FeatureCompare(df, featureDF)
            try:
                mergedDF = mergedDF.merge(df, on = ['DTXCID','MASS', 'Feature_MASS'], how = 'outer')
            except:
                mergedDF = df
    return mergedDF

path = r'/home/mboyce/Documents/MSMS_Haloperidol/MSData'
os.chdir(f'{path}/MSMS_query')

# This loops through all mgf files and generates csv files
for infile in glob.glob( os.path.join(f'{path}/MGF/Negative', '*.mgf') ): # Only reads mgf files in a directory

    if not os.path.exists(infile.rsplit('.',1)[0]+r'.csv'):   
        print("MGF to CSV: current file is: " + infile)
        mg.parseMGF(infile) # parse the file

#time_dict = {'SQL':0, 'scoring':0}

fpcdl = pd.read_csv(os.getcwd()+'/FeatureList_Negative.csv') # This is the list of masses to input that it will look for in the mgf/csv file
# CFMID search each csv file generated from above
for infile in glob.glob( os.path.join(path, '*.csv') ): # Reads all CSV files in directory
    print("CSV search/score: current file is: " + infile)
    t0=time.monotonic()
    #dfcfmid = mg.compare_mgf_df(infile,infile,10,0.02,POSMODE=True,filtering=False)
    mg.compare_mgf_df(infile,infile,10,0.02,POSMODE=False,filtering=False, mass_feature_list = fpcdl)
    t1=time.monotonic()
    print("time to Process is " + str(t1-t0))
    #print("Total time processing is " + str(round((t1-t_start)/60)) + " minutes")
    
#import featureList
positiveFeatures = pd.read_csv(f'{path}/MGF/FeatureList_Positive.csv')
negativeFeatures = pd.read_csv(f'{path}/MGF/FeatureList_Negative.csv')
MPPFeatures = pd.read_csv(f'{path}/MGF/FeatureList_MPP.csv')

All_Features_DF = pd.DataFrame(zip(positiveFeatures['MASS'], negativeFeatures['MASS'], MPPFeatures['MASS']),
                        columns = ['Positive', 'Negative', 'MPP'])




