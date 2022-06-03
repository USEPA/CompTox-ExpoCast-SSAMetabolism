# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 17:40:09 2021

@author: MBOYCE
"""
import os
import sys

script_dir = r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Python Scripts'
sys.path.append(script_dir)
import Clean_SeperateMPP as CleanSeperate_script
import NTA_summary_compression as NTA_script
import merge_results as merge_script

roots = [
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Haloperidol_CaseStudy',
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\CP-122721_CaseStudy',
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Celecoxib_CaseStudy',
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Dapsone_CaseStudy',
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Curcumin_CaseStudy',
    r'L:\Lab\NCCT_ExpoCast\ExpoCast2020\SSA-Metabolism\CaseStudy\Sulindac_CaseStudy'
    ]

for root_dir in roots:
    CleanSeperate_script.CleanAndSeperate(root_dir)
    NTA_script.RunNTA(root_dir)
    merge_script.MergeResults(root_dir)
