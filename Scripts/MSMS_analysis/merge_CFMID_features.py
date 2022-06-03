# -*- coding: utf-8 -*-
"""
Merge CFMID Data



Created on Wed Jul 14 23:14:54 2021
@author: MBOYCE
"""
import pandas as pd
import os
import glob


#1) Set root location
#2) Open files and store them, merge on x
#3) Store merged data with mergedDF from NTA webapp

rootDir = r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Haloperidol_CaseStudy'



CFMID_path = rootDir + r'\CFMID_Results'

saveDir = rootDir + '\CFMID_Merged'
if not os.path.exists(saveDir):
    os.makedirs(saveDir)

filterDir = rootDir + r'\Filtered_features'
#import featureList

def FeatureCompare(CFMID_Results_DF, FeatureList):
    matched_features = []
    for mass in CFMID_Results_DF['MASS']:
        diffSum = FeatureList['Mass'].apply(lambda x: abs(x - mass))
        matched_features.append(FeatureList.loc[diffSum == diffSum.min(), 'Mass'].values[0])
    return matched_features
    

def MergeCFMIDResults(rootdir, featureDF):
    mergedDF = pd.DataFrame()
    dataPaths = [rootdir + r'\Negative', rootdir + r'\Positive']
    for path in dataPaths:
        for filePath in glob.glob(os.path.join(path, '*.xlsx')):
            df = pd.read_excel(filePath)[['MASS','DTXCID','FORMULA','SMILES','energy_sum','MASS_in_MGF']]
            df['Feature_MASS'] = FeatureCompare(df, featureDF)
            try:
                mergedDF = mergedDF.merge(df, on = ['DTXCID','MASS','FORMULA','SMILES', 'Feature_MASS'], how = 'outer')
            except:
                mergedDF = df
    return mergedDF

MPPFeatures = pd.read_csv(filterDir + r'\FeatureList_subset.csv')
merged_CFMID_DF = MergeCFMIDResults(CFMID_path, MPPFeatures)

merged_CFMID_DF['energy_sum_MEDIAN'] = merged_CFMID_DF[[x for x in merged_CFMID_DF.columns if 'energy_sum' in x]].median(axis = 1)
groups = merged_CFMID_DF.groupby(by = 'Feature_MASS').max()
merged_CFMID_DF['quotient_SCORE'] = [x['energy_sum_MEDIAN']/groups.loc[x['Feature_MASS'],'energy_sum_MEDIAN'] for idx,x in merged_CFMID_DF.iterrows()]
merged_CFMID_DF = merged_CFMID_DF[['MASS', 'DTXCID','FORMULA','SMILES','Feature_MASS','energy_sum_MEDIAN', 'quotient_SCORE']].sort_values(by = 'Feature_MASS')
merged_CFMID_DF.to_csv(saveDir + '\MergedCFMIDResults.csv', index = False)
merged_CFMID_DF[merged_CFMID_DF['quotient_SCORE'] > 0.8].sort_values(by = ['Feature_MASS','quotient_SCORE'], ascending = [True, False]).to_csv(saveDir + '\MergedCFMIDResults_quotientCutOff.csv', index = False)

