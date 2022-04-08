# -*- coding: utf-8 -*-
"""
Created on Mon May 24 15:45:31 2021

@author: MBOYCE
"""

import pandas as pd
import numpy as np
import os
import csv
import time
import logging
import traceback
import shutil
import json
#import gridfs
#import pymongo
from datetime import datetime
#from dask.distributed import Client, LocalCluster, fire_and_forget
import re
from operator import itemgetter
from difflib import SequenceMatcher
from itertools import groupby


BLANKS = ['MB_', 'blank', 'blanks', 'BLANK', 'Blank']

minimum_rt = 0.0
rt_accuracy = 0.05
mass_accuracy_units = 'ppm'
mass_accuracy_tr = 5
rt_accuracy_tr = 0.1

sample_to_blank = 3
min_replicate_hits = 3
max_replicate_cv = 1.3

#REP_NUM = 3
HBR = 3.0 # High_Blank_Ratio condition
HMR = 1.5 # High_Mid_Ratio condition
SCORE = 90 # formula match is 90

FILENAMES = {'stats': ['stats_pos', 'stats_neg'],
             'tracers': ['tracers_pos', 'tracers_neg'],
             'tracer_plots': ['tracer_plot_pos', 'tracer_plot_neg'],
             #'cleaned': ['cleaned_pos', 'cleaned_neg'],
             #'flags': ['flags_pos', 'flags_neg'],
             #'combined': 'combined',
             'mpp_ready': ['for_stats_full', 'for_stats_reduced'],
             #'dashboard': 'dashboard_search',
             'toxpi': ['final_output_full', 'final_output_reduced']
             }

def RunNTA(directory):

    mainDir = directory + "\WebApp_Input"
    saveDir = directory + "\WebApp_Results"
    
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    fileType = ['Pellet', 'Super', 'Gluc']
    Hour = [0,1,4]
    fileNum = 0
    
    for filetype in fileType:
        for hour in Hour:
            df0 = fn_read_data(mainDir + r'\\'+filetype+"_"+str(hour)+'_neg.csv', 0)
            df1 = fn_read_data(mainDir + r'\\'+filetype+"_"+str(hour)+'_pos.csv', 1)
            
            filePrefix = filetype + '_' +str(hour)
            
            dfs = [df0, df1]
            tracer_df = []
            
            void_filter_dfs = filter_void_volume(dfs, minimum_rt)
            dup_filter_dfs = filter_duplicates(void_filter_dfs)
            stat_dfs = calc_statistics(dup_filter_dfs)
            #check_tracers(stat_dfs, tracer_df)
            clean_dfs = clean_features(stat_dfs)
            create_flags(clean_dfs)
            combined_df, MPP_df = combine_modes(clean_dfs)
            MPP_df['Has_Adduct_or_Loss'] = combined_df['Has_Adduct_or_Loss']
            MPP_df['Is_Adduct_or_Loss'] = combined_df['Is_Adduct_or_Loss']
            MPP_df['Adduct_or_Loss_Info'] = combined_df['Adduct_or_Loss_Info']
            MPP_df.to_csv(saveDir + r'\\'+filePrefix+ r'_for_stats_full.csv')
            reduced_file(MPP_df).to_csv(saveDir + r'\\' + filetype + "_" + str(hour) + r'_for_stats_reduced.csv')
            
            fileNum += 1


        
def fn_fix_names(df,index): # parse the Dataframe into a numpy array
        #df.columns = df.columns.str.replace(': Log2','') #log specific code
        df.columns = df.columns.str.replace(' ','_')
        df.columns = df.columns.str.replace('\([^)]*\)','')
        df['Compound'] = df['Compound'].str.replace("\ Esi.*$","")
        if 'Ionization_mode' in df.columns:
            df.rename(columns = {'Ionization_mode':'Ionization_Mode'},inplace=True)
        #df.drop(['CompositeSpectrum','Compound_Name'],axis=1)
        df.drop(['Compound_Name'],axis=1)
        Headers = fn_parse_headers(df,index)
        Abundance = [item for sublist in Headers for item in sublist if len(sublist)>1]    
        Samples= [x for x in Abundance]
        NewSamples = common_substrings(Samples)
        df.drop([col for col in df.columns if 'Spectrum' in col], axis=1,inplace=True)
        for i in range(len(Samples)):
            df.rename(columns = {Samples[i]:NewSamples[i]},inplace=True)
        #df = df
        return df
    
def common_substrings(ls=None):
    match  = SequenceMatcher(None,ls[0],ls[len(ls)-1]).find_longest_match(0,len(ls[0]),0,len(ls[len(ls)-1]))
    common = ls[0][match.a: match.a + match.size]
    #print((" ********* " + common))
    lsnew = list()
    for i in range(len(ls)):
        if len(common) > 3:
            lsnew.append(ls[i].replace(common,''))
        else:
            lsnew.append(ls[i])
            #print ls
    return lsnew

def score(df):  # Get score from annotations.
    regex = "^.*=(.*) \].*$"  # a regex to find the score looking for a pattern of "=something_to_find ]"
    if "Annotations" in df:
        if df.Annotations.isnull().all():  # make sure there isn't a totally blank Annotations column
            df['Score'] = None
            return df
        if df.Annotations.str.contains('overall=').any():
            df['Score'] = df.Annotations.str.extract(regex, expand=True).astype('float64')
    elif "Score" in df:
        pass
    else:
        df['Score'] = None
    return df

def fn_read_data(file,index):  # read a csv file into a DataFrame
        ext = os.path.splitext(file)[1]
        #print(ext)
        if ext == '.tsv':
            df = pd.read_csv(file,sep='\t',comment='#',na_values= 1 | 0)
        if ext == '.csv':
            df = pd.read_csv(file,comment='#',na_values= 1 | 0)
        df = fn_fix_names(df,index)
        return df

def fn_parse_headers(df,index): #group headers into a group of samples
        #global df
        headers = [[],[]]
        headers[index] = df.columns.values.tolist()
        countS=0
        countD=0
        new_headers = [[],[]]
        New_Headers = [None,None]
        Headers = [None,None]
        groups = [None,None]
        for s in range(0,len(headers[index])-1):
            #print headers[s],headers[s+1],list(set(str(headers[s])) - set(str(headers[s+1])))
            if 'blank' or 'Blank' or 'MB' in headers[index][s]:
                if differences(str(headers[index][s]),str(headers[index][s+1])) < 2: #3 is more common
                            countS += 1
                if differences(str(headers[index][s]),str(headers[index][s+1])) >= 2:
                            countD += 1
                            countS = countS + 1
            else:
                if differences(str(headers[index][s]),str(headers[index][s+1])) < 2: #2 is more common
                            countS += 1
                if differences(str(headers[index][s]),str(headers[index][s+1])) >= 2:

                            countD += 1
                            countS = countS + 1
                    #print "These are different "
            if "_Flags" in headers[index][s]:
                break
            new_headers[index].append([headers[index][countS],countD])
            new_headers[index].sort(key = itemgetter(1))
            groups[index] = groupby(new_headers[index], itemgetter(1))
            New_Headers[index] = [[item[0] for item in data] for (key, data) in groups[index]] 
        Headers[index] = New_Headers[index]
        #print((Headers[1]))
        return Headers[index]
    
    
def filter_void_volume(dfs, min_rt):
    dfs = [df.loc[df['Retention_Time'] > min_rt].copy() for index, df in enumerate(dfs)]
    return dfs

def filter_duplicates(dfs):
    dfs = [task_fun_duplicates(df) for df in dfs]
    return dfs

def task_fun_duplicates(df, mass_cutoff=0.005, rt_cutoff=0.05):
    """
    Drop features that are deemed to be duplicates. Duplicates are defined as two features whose differences in
    both mass and retention time are less than the defined cutoffs.
    :param df: A dataframe of MS feature data
    :param mass_cutoff: Mass differences below this number are possible duplicates. Units of Da.
    :param rt_cutoff: Retention time differences below this number are possible duplicates. Units of mins.
    :return: A dataframe with duplicates removed
    """
    df_new = df.copy()
    samples_df = df.filter(like='Sample', axis=1)
    df_new['all_sample_mean'] = samples_df.mean(axis=1)  # mean intensity across all samples
    df_new.sort_values(by=['all_sample_mean'], inplace=True, ascending=False)
    df_new.reset_index(drop=True, inplace=True)
    mass = df_new['Mass'].to_numpy()
    rts = df_new['Retention_Time'].to_numpy()
    masses_matrix = np.reshape(mass, (len(mass), 1))
    rts_matrix = np.reshape(rts, (len(rts), 1))
    diff_matrix_mass = masses_matrix - masses_matrix.transpose()
    diff_matrix_rt = rts_matrix - rts_matrix.transpose()
    duplicates_matrix = np.where((abs(diff_matrix_mass) <= mass_cutoff) & (abs(diff_matrix_rt) <= rt_cutoff),1,0)
    np.fill_diagonal(duplicates_matrix, 0)
    row_sums = np.sum(duplicates_matrix, axis=1)  # gives number of duplicates for each df row
    duplicates_matrix_lower = np.tril(duplicates_matrix)  # lower triangle of matrix
    lower_row_sums = np.sum(duplicates_matrix_lower, axis=1)
    to_keep = df_new[(row_sums == 0) | (lower_row_sums == 0)].copy()
    to_keep.sort_values(by=['Mass'], inplace=True)
    to_keep.reset_index(drop=True, inplace=True)
    to_keep = to_keep.drop(['all_sample_mean'], axis=1).copy()
    return to_keep

def calc_statistics(dfs, mass_accuracy = 10, mass_accuracy_units = 'ppm', rt_accuracy = 0.05):
    ppm = mass_accuracy_units == 'ppm'
    dfs = [taskfun_statistics(df) for df in dfs]
    dfs[0] = taskfun_assign_feature_id(dfs[0])
    dfs[1] = taskfun_assign_feature_id(dfs[1], start=len(dfs[0].index)+1)
    dfs[0] = taskfun_adduct_identifier(dfs[0], mass_accuracy, rt_accuracy, ppm,
                                             ionization='positive', id_start=1)
    dfs[1] = taskfun_adduct_identifier(dfs[1], mass_accuracy, rt_accuracy, ppm,
                                             ionization='negative', id_start=len(dfs[0].index)+1)
    #mongo_save(dfs[0], FILENAMES['stats'][0])
    #mongo_save(dfs[1], FILENAMES['stats'][1])
    return dfs

def taskfun_statistics(df_in):
    """
    # Calculate Mean,Median,STD,CV for every feature in a sample of multiple replicates
    :param df_in: the dataframe to calculate stats for
    :return: a new dataframe including the stats
    """
    df = df_in.copy()
    all_headers = taskfun_parse_headers(df)
    abundance = [item for sublist in all_headers for item in sublist if len(sublist) > 1]
    df = score(df)
    filter_headers= ['Compound','Ionization_Mode','Score','Mass','Retention_Time','Frequency'] + abundance
    df = df[filter_headers].copy()
    for the_list in all_headers:
        REP_NUM = len(the_list)
        if REP_NUM > 1:
            for i in range(0, REP_NUM):
                # match finds the indices of the largest common substring between two strings
                match = SequenceMatcher(None, the_list[i], the_list[i+1]).find_longest_match(0, len(the_list[i]),0, len(the_list[i+1]))
                df['Mean_'+ str(the_list[i])[match.a:match.a +  match.size]] = df[the_list[i:i + REP_NUM]].mean(axis=1).round(0)
                df['Median_'+ str(the_list[i])[match.a:match.a +  match.size]] = df[the_list[i:i + REP_NUM]].median(axis=1,skipna=True).round(0)
                df['STD_'+ str(the_list[i])[match.a:match.a +  match.size]] = df[the_list[i:i + REP_NUM]].std(axis=1,skipna=True).round(0)
                df['CV_'+ str(the_list[i])[match.a:match.a +  match.size]] = (df['STD_'+ str(the_list[i])[match.a:match.a +  match.size]]/df['Mean_' + str(the_list[i])[match.a:match.a + match.size]]).round(4)
                df['N_Abun_'+ str(the_list[i])[match.a:match.a +  match.size]] = df[the_list[i:i + REP_NUM]].count(axis=1).round(0)
                break
    df.sort_values(['Mass', 'Retention_Time'], ascending=[True, True], inplace=True)
    df['Rounded_Mass'] = df['Mass'].round(0)
    return df

def taskfun_parse_headers(df_in):
    '''
    A function to group the dataframe's column headers into sets of similar names which represent replicates
    :param df_in: the dataframe of features
    :return: a list of groups of column labels
    '''
    df = df_in.copy()
    headers = df.columns.values.tolist()
    countS = 0
    countD = 0
    new_headers = []
    for s in range(0,len(headers)-1):
        if 'blank' or 'Blank' or 'MB' in headers[s]:
            if differences(str(headers[s]),str(headers[s+1])) < 2: #3 is more common
                countS += 1
            if differences(str(headers[s]),str(headers[s+1])) >= 2:
                countD += 1
                countS = countS + 1
        else:
            if differences(str(headers[s]),str(headers[s+1])) < 2: #2 is more common
                countS += 1
            if differences(str(headers[s]),str(headers[s+1])) >= 2:

                countD += 1
                countS = countS + 1
            #print "These are different "
        if "_Flags" in headers[s]:
            break
        new_headers.append([headers[countS],countD])
        new_headers.sort(key = itemgetter(1))
    groups = groupby(new_headers, itemgetter(1))
    new_headers_list = [[item[0] for item in data] for (key, data) in groups]
    return new_headers_list

def taskfun_assign_feature_id(df_in, start = 1):
    """
    A function to assign unique feature ids to a nta dataset
    :param df_in: the dataframe to assign ids to
    :param start: assign ids starting at this integer
    :return: returns the new df with unique feature ids added
    """
    df = df_in.copy()
    row_nums = list(range(0, len(df.index)))
    to_assign = [x + start for x in row_nums]
    df.insert(0, 'Feature_ID', to_assign.copy())
    return df

def taskfun_adduct_identifier(df_in, Mass_Difference, Retention_Difference, ppm, ionization, id_start = 1):  # TODO optimize memory usage
    """
    Label features which could have adduct or loss products in the feature list, or may be adduct of loss products of
    other features. This information is added to the dataframe in three columns, Has_Adduct_or_Loss - true/false,
    Is_Adduct_or_Loss - true/false, and Adduct_or_Less_Info which points to the feature number of the associated
    adduct/loss or parent feature and gives the type of adduct/loss.
    :param df_in: A dataframe of MS features
    :param Mass_Difference: Mass differences below this number are possible duplicates. Units of ppm or Da based on the
    ppm parameter.
    :param Retention_Difference: Retention time differences below this number are possible duplicates. Units of mins.
    :param ppm: True if mass differences are given in ppm, otherwise False and units are Da.
    :param ionization: 'positive' if data are from positive ionization mode, otherwise 'negative'.
    :param id_start: The first feature id in the dataset (defaults to 1)
    :return: A Dataframe where adduct info is given in three new columns.
    """
    df = df_in.copy()
    mass = df['Mass'].to_numpy()
    rts = df['Retention_Time'].to_numpy()
    masses_matrix = np.reshape(mass, (len(mass), 1))
    rts_matrix = np.reshape(rts, (len(rts),1))
    diff_matrix_mass = masses_matrix - masses_matrix.transpose()
    diff_matrix_rt = rts_matrix - rts_matrix.transpose()
    pos_adduct_deltas = {'Na': 22.989218, 'K': 38.963158, 'NH4': 18.033823}
    neg_adduct_deltas = {'Cl': 34.969402, 'Br': 78.918885, 'HCO2': 44.998201, 'CH3CO2': 59.013851, 'CF3CO2': 112.985586}
    neutral_loss_deltas= {'H2O': -18.010565, 'CO2': -43.989829}
    proton_mass = 1.007276
    if ionization == "positive":
        # we observe Mass+(H+) and Mass+(Adduct)
        possible_adduct_deltas = {k: v - proton_mass for (k,v) in pos_adduct_deltas.items()}
    else:
        # we observe Mass-(H+) and Mass+(Adduct)
        possible_adduct_deltas = {k: v + proton_mass for (k,v) in neg_adduct_deltas.items()}
    possible_adduct_deltas.update(neutral_loss_deltas)  # add our neutral losses
    df['Has_Adduct_or_Loss'] = 0
    df['Is_Adduct_or_Loss'] = 0
    df['Adduct_or_Loss_Info'] = ""
    unique_adduct_number = np.zeros(len(df.index))
    for a_name, delta in sorted(possible_adduct_deltas.items()):
        is_adduct_diff = abs(diff_matrix_mass - delta)
        has_adduct_diff = abs(diff_matrix_mass + delta)
        if ppm:
            is_adduct_diff = (is_adduct_diff/masses_matrix)*10**6
            has_adduct_diff = (has_adduct_diff/masses_matrix)*10**6
        is_adduct_matrix = np.where((is_adduct_diff < Mass_Difference) & (abs(diff_matrix_rt) < Retention_Difference), 1, 0)
        has_adduct_matrix = np.where((has_adduct_diff < Mass_Difference) & (abs(diff_matrix_rt) < Retention_Difference), 1, 0)
        np.fill_diagonal(is_adduct_matrix, 0)  # remove self matches
        np.fill_diagonal(has_adduct_matrix, 0)  # remove self matches
        row_num = len(mass)
        is_id_matrix = np.tile(np.arange(row_num),row_num).reshape((row_num,row_num)) + id_start
        has_id_matrix = is_id_matrix.transpose()
        is_adduct_number = is_adduct_matrix * is_id_matrix
        is_adduct_number_flat = np.max(is_adduct_number, axis=1) # if is adduct of multiple, keep highest # row
        is_adduct_number_flat_index = np.where(is_adduct_number_flat > 0, is_adduct_number_flat -1, 0)
        is_adduct_of_adduct = np.where((is_adduct_number_flat > 0) &
                                       (df['Is_Adduct_or_Loss'][pd.Series(is_adduct_number_flat_index-id_start).clip(lower=0)] > 0), 1, 0)
        is_adduct_number_flat[is_adduct_of_adduct == 1] = 0
        has_adduct_number = has_adduct_matrix * is_id_matrix
        has_adduct_number_flat = np.max(has_adduct_number, axis=1)  # these will all be the same down columns
        unique_adduct_number = np.where(has_adduct_number_flat != 0, has_adduct_number_flat, is_adduct_number_flat).astype(int)
        #unique_adduct_number = np.where(unique_adduct_number == 0, unique_adduct_number_new, unique_adduct_number)
        df['Has_Adduct_or_Loss'] = np.where((has_adduct_number_flat > 0) & (df['Is_Adduct_or_Loss'] == 0),
                                            df['Has_Adduct_or_Loss']+1, df['Has_Adduct_or_Loss'])
        df['Is_Adduct_or_Loss'] = np.where((is_adduct_number_flat > 0) & (df['Has_Adduct_or_Loss'] == 0), 1, df['Is_Adduct_or_Loss'])
        # new_cols = ['unique_{}_number'.format(a_name), 'has_{}_adduct'.format(a_name), 'is_{}_adduct'.format(a_name)]
        #unique_adduct_number_str =
        df['Adduct_or_Loss_Info'] = np.where((has_adduct_number_flat > 0) & (df['Is_Adduct_or_Loss'] == 0),
                                             df['Adduct_or_Loss_Info'] + unique_adduct_number.astype(str) + "({});".format(a_name), df['Adduct_or_Loss_Info'])
        df['Adduct_or_Loss_Info'] = np.where((is_adduct_number_flat > 0) & (df['Has_Adduct_or_Loss'] == 0),
                                             df['Adduct_or_Loss_Info'] + unique_adduct_number.astype(str) + "({});".format(a_name), df['Adduct_or_Loss_Info'])
    return df

def check_tracers(dfs, tracer_df):

    ppm = mass_accuracy_units_tr == 'ppm'
    tracer_dfs_out = [fn.check_feature_tracers(df, tracer_df, mass_accuracy_tr, rt_accuracy_tr, ppm) for index, df in enumerate(dfs)]
    tracer_dfs_out = [format_tracer_file(df) for df in tracer_dfs_out]
    tracer_plots_out = [create_tracer_plot(df) for df in tracer_dfs_out]
    mongo_save(tracer_dfs_out[0], FILENAMES['tracers'][0])
    mongo_save(tracer_dfs_out[1], FILENAMES['tracers'][1])
    mongo_save(tracer_plots_out[0], FILENAMES['tracer_plots'][0], df=Falselse)
    mongo_save(tracer_plots_out[1], FILENAMES['tracer_plots'][1], df=False)
    return

def clean_features(dfs, sample_to_blank = 3, min_replicate_hits = 3, max_replicate_cv = 1.3):
    controls = [sample_to_blank, min_replicate_hits, max_replicate_cv]
    dfs = [task_fun_clean_features(df, controls) for index, df in enumerate(dfs)]
    dfs = [fn_Blank_Subtract(df, index) for index, df in enumerate(dfs)]  # subtract blanks from medians
    #self.mongo_save(self.dfs[0], FILENAMES['cleaned'][0])
    #self.mongo_save(self.dfs[1], FILENAMES['cleaned'][1])
    return dfs

def task_fun_clean_features(df, controls):  # a method that drops rows based on conditions
    Abundance=  df.columns[df.columns.str.contains(pat ='N_Abun_')].tolist()
    blanks = ['MB','mb','mB','Mb','blank','Blank','BLANK']
    Median = df.columns[df.columns.str.contains(pat ='Median_')].tolist()
    Median_Samples = [md for md in Median if not any(x in md for x in blanks)]
    Median_High = [md for md in Median if 'C' in md]
    Median_Mid = [md for md in Median if 'B' in md]
    Median_Low = [md for md in Median if 'A' in md]
    Median_MB = [md for md in Median if any(x in md for x in blanks)]
    N_Abun_High = [N for N in Abundance if 'C' in N]
    N_Abun_MB = [N for N in Abundance if any(x in N for x in blanks)]
    N_Abun_Samples = [N for N in Abundance if not any(x in N for x in blanks)]
    #N_Abun_MB= [N for N in Abundanceif 'MB' in N]
    CV =  df.columns[df.columns.str.contains(pat ='CV_')].tolist()
    CV_Samples= [C for C in CV if not any(x in C for x in blanks)]
    #set medians where feature abundance is less than some cutoff to nan
    df['AnySamplesDropped'] = np.nan
    for median,N in zip(Median_Samples,N_Abun_Samples):
        #print((str(median) + " , " +str(N)))
        df.loc[df[N] < controls[1], median] = np.nan
        df.loc[df[N] < controls[1], 'AnySamplesDropped'] = 1
    # remove all features where the abundance is less than some cutoff in all samples
    df.drop(df[(df[N_Abun_Samples] < controls[1]).all(axis=1)].index, inplace=True)
    df.drop(df[(df[CV_Samples] > controls[2]).all(axis=1)].index, inplace=True)
    # blank out samples that do not meet the CV cutoff
    cv_not_met = df[CV_Samples] > controls[2]
    m = df[Median_Samples].copy()
    cv_not_met.columns = m.columns
    df[Median_Samples] = m.mask(cv_not_met)
    #find the median of all samples and select features where median_samples/ median_blanks >= cutoff
    df['Max_Median_ALLSamples'] = df[Median_Samples].max(axis=1,skipna=True).round(0)
    df['SampletoBlanks_ratio'] = df['Max_Median_ALLSamples'].astype('float')/df[Median_MB[0]].astype('float')
    df = df[(df[N_Abun_MB[0]] == 0) | (df['SampletoBlanks_ratio'] >= controls[0])].copy()
    return df

def fn_Blank_Subtract(df,index):
    """
    Calculate the median blank intensity for each feature and subtract that value from each sample's median value for
    that feature
    """
    Abundance = [[],[]]
    Headers = [0,0]
    Blanks = [[],[]]
    Median = [[],[]]
    Headers[index] = fn_parse_headers(df,index)
    blanks_str = BLANKS
    Abundance[index] = [item for sublist in Headers[index] for item in sublist if (len(sublist)>1) & (not any(x in item for x in blanks_str))]
    # On with the agony of subtracting the MB median from Samples
    Blanks[index] = df.columns[(df.columns.str.contains(pat ='MB_|blank|blanks|BLANK|Blank')) &
                               (df.columns.str.contains(pat='Mean|Median|CV|STD|N_Abun|ratio') == False)].tolist()
    Median[index] = df.columns[(df.columns.str.contains(pat ='Median_')==True) & (df.columns.str.contains(pat ='MB|blank|blanks|BLANK|Blank')==False)].tolist()
    df['Median_ALLMB'] = df[Blanks[index]].median(axis=1,skipna=True).round(0).fillna(0)  # instead using median calc in statistics
    #df[Abundance[index]] = df[Abundance[index]].sub(df['Median_ALLMB'],axis=0)
    for median in Median[index]:
        df["BlankSub_"+str(median)] = df[median].sub(df['Median_ALLMB'],axis=0)
        df["BlankSub_"+str(median)] = df["BlankSub_"+str(median)].clip(lower=0).replace({0:np.nan})
    #df[Abundance[index]] = df[Abundance[index]].clip(lower=0).replace({0:np.nan})
    return df


def create_flags(dfs):
    dfs = [fn_flags(df) for df in dfs]
    #self.mongo_save(self.dfs[0], FILENAMES['flags'][0])
    #self.mongo_save(self.dfs[1], FILENAMES['flags'][1])
    
def fn_flags(df): # a method to develop required flags
    df['Neg_Mass_Defect'] = np.where((df.Mass - df.Mass.round(0)) < 0 , '1','0')
    df['Halogen'] = np.where(df.Compound.str.contains('F|l|r|I'),'1','0')
    df['Formula_Match'] = np.where(df.Score != df.Score,'0','1') #check if it does not have a score
    df['Formula_Match_Above90'] = np.where(df.Score >= SCORE,'1','0')
    df['X_NegMassDef_Below90'] = np.where(((df.Score < SCORE) & (df.Neg_Mass_Defect == '1') & (df.Halogen == '1')),'1','0')
    df['For_Dashboard_Search'] = np.where(((df.Formula_Match_Above90 == '1') | (df.X_NegMassDef_Below90 == '1')) , '1', '0')
    df.sort_values(['Formula_Match','For_Dashboard_Search','Formula_Match_Above90','X_NegMassDef_Below90'],ascending=[False,False,False,False],inplace=True)
    #df.to_csv('input-afterflag.csv', index=False)
    #print df1
    df.sort_values('Compound',ascending=True,inplace=True)
    return df

def combine_modes(dfs):
    df_combined = fn_combine(dfs[0], dfs[1])
    #self.mongo_save(self.df_combined, FILENAMES['combined'])
    mpp_ready = fn_MPP_Ready(df_combined)
    #mongo_save(mpp_ready, FILENAMES['mpp_ready'][0])
    #mongo_save(reduced_file(mpp_ready), FILENAMES['mpp_ready'][1])  # save the reduced version
    return df_combined, mpp_ready
    
def fn_combine(df1,df2):
    #Headers = [[],[]]
        #Headers[0] = parse_headers(df1,0)
        #Headers[1] = parse_headers(df2,1)
    #print("##############")
    Abundance=[[],[]]
    Abundance[0] = df1.columns.values.tolist()
    Abundance[1] = df2.columns.values.tolist()
    #diff = append_headers(Abundance[0],Abundance[1])
    #print len(df1.columns.values.tolist())
    #for i in range(len(Abundance[0])):
    #    #print (Abundance[0][i],Abundance[1][i])
    #    df1.rename(columns = {Abundance[0][i]:new_headers[i]},inplace=True)
    #    df2.rename(columns = {Abundance[1][i]:new_headers[i]},inplace=True)
    #print df1.columns.values.tolist()
    #print(" ||||___|||| - - - - - - ")
    #print df2.columns.values.tolist()
    #df1[list(set(Abundance[0])-set(diff))]    = np.nan
    #df2[list(set(Abundance[1])-set(diff))]    = np.nan
    dfc = pd.concat([df1,df2], sort=True) #fixing pandas FutureWarning
    dfc = dfc.reindex(columns = df1.columns)
    columns = dfc.columns.values.tolist()
    #print((str(len(columns)) + " ##### " + str(len(df1.columns.values.tolist())) + " #### " + str(len(df2.columns.values.tolist()))))
    dfc = pd.merge(dfc,df2,suffixes=['','_x'],on='Compound',how='left')
    dfc = pd.merge(dfc,df1,suffixes=['','_y'],on='Compound',how='left')

    # create new flags
    dfc = dfc.drop_duplicates(subset=['Compound','Mass','Retention_Time','Score'])
    dfc['Both_Modes'] = np.where(((abs(dfc.Mass_x-dfc.Mass_y)<=0.005) & (abs(dfc.Retention_Time_x-dfc.Retention_Time_y)<=1)),'1','0')
    dfc['N_Compound_Hits'] = dfc.groupby('Compound')['Compound'].transform('size')
    Median_list =  dfc.columns[(dfc.columns.str.contains(pat ='Median_')==True)\
                 & (dfc.columns.str.contains(pat ='MB|blank|blanks|BlankSub|_x|_y')==False)].tolist()
    #print(Median_list)
    dfc['N_Abun_Samples'] = dfc[Median_list].count(axis=1,numeric_only=True)
    dfc['Median_Abun_Samples'] = dfc[Median_list].median(axis=1,skipna=True).round(0)
    dfc['One_Mode_No_Isomers'] = np.where(((dfc.Both_Modes == '0') & (dfc.N_Compound_Hits == 1)),'1','0')
    dfc['One_Mode_Isomers'] = np.where(((dfc.Both_Modes == '0') & (dfc.N_Compound_Hits > 1)),'1','0')
    dfc['Two_Modes_No_Isomers'] = np.where(((dfc.Both_Modes == '1') & (dfc.N_Compound_Hits == 2)),'1','0')
    dfc['Two_Modes_Isomers'] = np.where(((dfc.Both_Modes == '1') & (dfc.N_Compound_Hits > 2)),'1','0')
    dfc['Est_Chem_Count'] = None #Default to non-type
    dfc.loc[dfc['One_Mode_No_Isomers'] == '1','Est_Chem_Count'] = 1
    dfc.loc[dfc['One_Mode_Isomers'] == '1','Est_Chem_Count'] = dfc['N_Compound_Hits']
    dfc.loc[(dfc['Two_Modes_No_Isomers'] == '1') | (dfc['Two_Modes_Isomers'] == '1'),'Est_Chem_Count'] = dfc['N_Compound_Hits']/2
    columns.extend(('Both_Modes','N_Compound_Hits','N_Abun_Samples','Median_Abun_Samples','One_Mode_No_Isomers','One_Mode_Isomers','Two_Modes_No_Isomers',
            'Two_Modes_Isomers','Est_Chem_Count'))
    dfc = dfc[columns].sort_values(['Compound'],ascending=[True])

    #dft.reset_index() 
    #dft.dropna(inplace=True)
    return dfc

def fn_MPP_Ready(dft, directory='',file=''):
    #dft = dft.rename(columns = {'Compound':'Formula','Retention_Time':'RT'})
    #dft['Compound Name'] = dft['Formula']
    dft = dft.rename(columns = {'Compound':'Formula'})
    Headers = fn_parse_headers(dft,0)
    raw_samples= [item for sublist in Headers for item in sublist if (len(sublist) > 2) & ('BlankSub' not in item)]
    blank_subtracted_medians = dft.columns[dft.columns.str.contains(pat='BlankSub')].tolist()
    #Blanks = dft.columns[dft.columns.str.contains(pat ='MB_')].tolist()
    #Samples = [x for x in Abundance if x not in Blanks]
    #NewSamples = common_substrings(Samples)
    #for i in range(len(Samples)):
    #    dft.rename(columns = {Samples[i]:NewSamples[i]},inplace=True)
    #columns = dft.columns.values.tolist()
    #dft = dft.reindex(columns=Columns)
    #print dft
    #dft.to_csv(directory+'/'+file+'_MPP_Ready.csv', index=False)
    dft = dft[['Feature_ID','Formula','Score', 'Mass','Retention_Time'] + raw_samples + blank_subtracted_medians]
    #dft.to_csv(directory+'/'+'Data_Both_Modes_MPP_Ready.csv', index=False)
    return dft

def connect_to_mongo_gridfs(address):
    db = pymongo.MongoClient(host=address).nta_storage
    print("Connecting to mongodb at {}".format(address))
    fs = gridfs.GridFS(db)
    return fs

def mongo_save(file, step="", df=True):
    if df:
        to_save = file.to_json(orient='split')
    else:
        to_save = file
    id = "00000000" + "_" + step
    gridfs = connect_to_mongo_gridfs(None)
    gridfs.put(to_save, _id=id, encoding='utf-8', project_name = 'test_project')
    
def reduced_file(df_in):
    df = df_in.copy()
    headers = fn_parse_headers(df, 0)
    keeps_str = ['MB_', 'blank', 'blanks', 'BLANK', 'Blank', 'Median', 'Sub']
    to_drop = [item for sublist in headers for item in sublist if
                        (len(sublist) > 1) & (not any(x in item for x in keeps_str))]
    to_drop.extend(df.columns[(df.columns.str.contains(pat ='CV_|N_Abun_|Mean_|STD_')==True)].tolist())
    to_drop.extend(df.columns[(df.columns.str.contains(pat ='Median_') == True) &
                              (df.columns.str.contains(pat ='MB|blank|blanks|BLANK|Blank|Sub')==False)].tolist())
    if 'Median_ALLMB' in df.columns.values.tolist():
        to_drop.extend(['Median_ALLMB'])
    df.drop(to_drop, axis=1, inplace=True)
    return df

def differences(s1,s2): #find the number of different characters between two strings (headers)
        s1 = re.sub(re.compile(r'\([^)]*\)'),'',s1)
        s2 = re.sub(re.compile(r'\([^)]*\)'),'',s2)
        count = sum(1 for a, b in zip(s1, s2) if a != b) + abs(len(s1) - len(s2))
        return count