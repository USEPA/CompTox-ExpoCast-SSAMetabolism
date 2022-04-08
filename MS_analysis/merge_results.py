# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:28:52 2021

@author: MBOYCE
"""


import os
import pandas as pd
import numpy as np

def MergeResults(directory):
    fileDir = directory + '\\WebApp_Results'
    files = os.listdir(fileDir)
    mergedDF = None
    dataFiles = [file for file in files if r'full.csv' in file]
    for file in dataFiles:
        df = pd.read_csv(fileDir + '\\'+ file)
        df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        notation = df.columns[df.columns.str.contains('BlankSub_Median_')==True][0].split('_', -1)[-1]        
        df.rename(columns = {df.columns[df.columns.str.contains('BlankSub_Max_Median_')==True][0] : 
                             df.columns[df.columns.str.contains('BlankSub_Max_Median_')==True][0]+'_'+notation}, 
                  inplace = True)
        if mergedDF is None:
            mergedDF = df
            continue
        #Combine results for each adduct format then drop column from merging DF
        mergedDF = pd.merge(mergedDF, df, how = "outer", on = ['Feature_ID', 'Formula', 'Score', 'Mass', 'Retention_Time'])

    hasAdductCols = mergedDF.columns[mergedDF.columns.str.contains(pat ='Has_Adduct')==True]
    mergedDF['Has_Adduct_or_Loss'] = np.where(mergedDF[hasAdductCols].sum(axis = 1) > 0, 1, 0)
    mergedDF.drop(labels = ['Has_Adduct_or_Loss_x', 'Has_Adduct_or_Loss_y'], inplace = True, axis = 1)
    
    isAdductCols = mergedDF.columns[mergedDF.columns.str.contains(pat ='Is_Adduct')==True]
    mergedDF['Is_Adduct_or_Loss'] = np.where(mergedDF[isAdductCols].sum(axis = 1) > 0, 1, 0)
    mergedDF.drop(labels = ['Is_Adduct_or_Loss_x', 'Is_Adduct_or_Loss_y'], inplace = True, axis = 1)
    
    infoAdductCols = mergedDF.columns[mergedDF.columns.str.contains(pat ='Adduct_or_Loss_Info')==True]
    mergedDF[infoAdductCols].fillna('').sum(axis = 1)
    mergedDF['Adduct_or_Loss_Info'] = mergedDF[infoAdductCols].fillna('').sum(axis = 1)
    mergedDF.drop(labels = ['Adduct_or_Loss_Info_x', 'Adduct_or_Loss_Info_y'], inplace = True, axis = 1)
    
    saveDir = directory + "\WebApp_Combined"
    
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    
    mergedDF.rename(columns = lambda x: x.replace('_x', '').replace('_y',''), inplace = True)
    mergedDF.to_csv(saveDir + '\\WebApp_Results_combined.csv', index = False)
            