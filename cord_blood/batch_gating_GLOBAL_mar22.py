#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:17:31 2020

@author: antton
"""

import aligater as ag
import pandas as pd
import numpy as np
from math import inf
from styleframe import StyleFrame, utils 
from random import randint
from datetime import date
import os
import sys


#Get file list
#Get repeats, store paths to them in 'repeats_filepaths' variable

path_to_excel = '/home/antton/TFM/data/repeats_info.xlsx'
path_to_files = "/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs" 

# # Extracting info from manually modified file


def only_cells_with_blue_background(cell):
    ## If cell has right color (non-white, non-red), write "repeat" in it and return. Else, leave alone.
    
    if cell.style.bg_color not in {utils.colors.white, '00000000'} and cell.style.bg_color not in {utils.colors.red, 'FFD9D9D9'}:
        cell.value = 'repeat'
    return cell

def get_filepaths(maxDate=None, maxSampleNum=None, includeRepeats=False):
    ## Create a list with the exact filepath to each sample file of interest.
    ## This includes genotyped samples and samples that have been measured twice.
    if maxDate is None:
        maxDate = date(2019,9,9) #last date we are interested in
    
    # Repeated measurements. Specific file-names obtained from excel file from Natsumi
    repeats_filepaths = []  # Empty array that will contain paths to every repeat sample .fcs  
    if includeRepeats:
        sf = StyleFrame.read_excel(path_to_excel, read_style=True, use_openpyxl_styles=False)
        sf_marked = StyleFrame(sf.applymap(only_cells_with_blue_background))  # All blue cells have been noted
        repeats_df = sf_marked.data_df  # DataFrame identical to the excel file with NaN for all values except for blue cells
        for index, row in repeats_df.iterrows():  # Iterate on the DataFrame row by row
            if row.notes !='avoid' and row.notes != 'poor viability':  # Filter out crossed-out cells
                rep_locs = np.where(row.values == 'repeat')[0]  # Array with column number of blue cells
                for column_index in rep_locs:
                    col_name = repeats_df.columns[column_index]  # Get actual column names (part of folder name) from index
                    for folder in os.listdir(path_to_files):
                        if folder.startswith(str(col_name)):  # Find the appropriate folder using the column name
                            path_into_folder = path_to_files+'/'+folder
                            for filename in os.listdir(path_into_folder):
                                if filename[3:].startswith(str(row.ID)):  # Find the specific file in folder using row ID
                                    repeats_filepaths.append(path_into_folder+'/'+filename)  # Save complete filepath
    
    #                                 
    relevant_folder_paths = []
    folderpaths = "/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/"
    all_folders = sorted(os.listdir(folderpaths))  # Take ALL the folder names in 'fcs'
    if '.DS_Store' in all_folders: all_folders.remove('.DS_Store') #damn mac hidden folders...

    for folder in all_folders:  # For each subfolder
        folder_date = date(int('20'+folder[:2]),int(folder[2:4]),int(folder[4:6]))  # Take date from folder name
        if folder_date > maxDate: # If we exceed the max date, stop
            break
        relevant_folder_paths.append(path_to_files+"/"+folder)  # Add the complete path to the subfolders
    
    genotyped_filepaths = []  # List of all filepath to input into the experiment object.
    for subfolder in relevant_folder_paths:
        #ls_path = subfolder.replace(' ','/ ')
        for file in os.listdir(subfolder):
            if file.endswith('.fcs'):  # For each fcs file in the folder...
                if maxSampleNum:
                    sampleNum = file.split('.')[0].split(' ')[1]
                    if '-' in sampleNum:
                        sampleNum = sampleNum.split('-')[0]
                    try:
                        sampleNum = int(sampleNum)
                    except ValueError:
                        continue
                    if sampleNum <= maxSampleNum:
                        genotyped_filepaths.append(subfolder+"/"+file)
                else:
                    genotyped_filepaths.append(subfolder+"/"+file)  # If there is no maxSampleNum, just add everything 
                
    final_filepaths = genotyped_filepaths + list(set(repeats_filepaths) - set(genotyped_filepaths)) # Fuse litsts. Unique paths only.  
    print('All repeat measurements: ', len(repeats_filepaths))
    print('All genotyped files: ', len(genotyped_filepaths))
    print('All relevant files: ', len(final_filepaths))
    return final_filepaths

def out_folder_list(): # List of all output folders for the images
    folder_list = ["/home/antton/TFM/output/plots/01-Viable",
        "/home/antton/TFM/output/plots/02-Singlet",
        "/home/antton/TFM/output/plots/03-PBMCs",
        "/home/antton/TFM/output/plots/04-CD45high",        
        "/home/antton/TFM/output/plots/05-CD45posCD34pos",
        "/home/antton/TFM/output/plots/06-CD3negBACKGATE",
        "/home/antton/TFM/output/plots/07-CD4posBACKGATE",
        "/home/antton/TFM/output/plots/08-CD19negBACKGATE",       
        "/home/antton/TFM/output/plots/09-CD14negBACKGATE",
        
        "/home/antton/TFM/output/plots/10-CD16CD56/10A-CD16neg_56negBACKGATE",
        "/home/antton/TFM/output/plots/10-CD16CD56/10B-CD16pospos_56posposBACKGATE",
        "/home/antton/TFM/output/plots/10-CD16CD56/10C-CD16neg_56posposBACKGATE",
     
        "/home/antton/TFM/output/plots/11-LinnegCD34pos",   
        "/home/antton/TFM/output/plots/12-CD38neg",    
        "/home/antton/TFM/output/plots/13-HSC_MLP_MPP",
        "/home/antton/TFM/output/plots/14-B_NK",
        "/home/antton/TFM/output/plots/15-CD10pos(MLP)",
        "/home/antton/TFM/output/plots/16-CMP_GMP_MEP"]
    
    return folder_list

def get_blacklist():
    blacklist = ['/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/181203 CB/C7 98.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190401 CB/D3 333.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190610 CB/C7 529.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190701 CB/C6 616.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190715 CB/D6 651.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190909 CB/E8 848.fcs']

    return blacklist

def gateFullDataset(my_sample, save_images=True, *args, **kwargs):
    #Setup
    ag.agconf.ag_verbose=False
    sample_filepath = my_sample.filePath
    
    if sample_filepath in get_blacklist():
        return my_sample
    
    date_plate = ag.getFileName(ag.getParent(sample_filepath))
    sampleName = ag.getFileName(sample_filepath)
    
    event_cap = 2500000
    if len(my_sample()) > event_cap:  ## Cap samples with too many events. Take only the first 3M
        my_sample.fcsDF = my_sample.fcsDF.iloc[0:event_cap]
        sys.stderr.write("Event cap reached. Downsapmled to "+str(event_cap)+" events.\n")
    #Put inconsistent marker names into variables
    markers = my_sample().columns.tolist()
    for marker in markers[6:]:
        if marker.startswith('CD45'):
            CD45 = marker
        if marker.startswith('7AAD') or marker.startswith('CD235'):
            marker_7AAD = marker
    
    if save_images:
        fileName = 'notNone'  # Dirty coding
    else:
        fileName = None
    
    ## Gate 1:  Viable/non-viable separation through 7AAD
    livedead_step1 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=214000, parentGate=None, orientation='vertical', population='lower')
    livedead_step2 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=-800, parentGate=livedead_step1, orientation='horisontal', population='upper')   
    halfcut = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=150000, parentGate=livedead_step2, orientation='vertical', population='upper')    
    ylim_back = ag.densityDelimitation(my_sample, xCol= marker_7AAD, parentGate=halfcut, interval=[300,1100], limit_threshold=0.2, direction='right',scale='bilog',T=200)    
    if ylim_back == inf:  # Failsafe in case of no limit
        ylim_back = 900    
    if fileName:
        fileName="/home/antton/TFM/output/plots/01-Viable/"+date_plate+"-"+sampleName+"-Viable.jpeg"    
    livedead_final = ag.gateThreshold(my_sample, name="Viable", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='logish', thresh=ylim_back, parentGate=livedead_step2, orientation='horisontal', population='lower', filePlot=fileName)
    

    ## Gate 2:  Singlets
    singlets_step1 = ag.gateThreshold(my_sample, name="remove_high_clutter", xCol='FSC 488/10-A', yCol='FSC 488/10-H', scale='linear', thresh=210000, parentGate=livedead_final,  orientation='vertical', population='lower')
    if fileName:
        fileName="/home/antton/TFM/output/plots/02-Singlet/"+date_plate+"-"+sampleName+"-Singlet.jpeg"
    singlets = ag.gatePC(my_sample, name="Singlets", xCol='FSC 488/10-A', yCol='FSC 488/10-H', center='centroid', adjustAngle=3,widthScale=2.5,heightScale=3.5, parentGate=singlets_step1, filePlot=fileName)


    ## Gate 3: PBMC cells  
    
    halfcut_middle = ag.gateThreshold(my_sample, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=120000, parentGate=singlets, orientation='vertical', population='upper')

    ylim_bot = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_middle, interval=[10000,20000], limit_threshold=0.05, direction='left',scale='linear')
    if ylim_bot == inf:
        ylim_bot = 10000
    PBMC_step1 = ag.gateCorner(my_sample, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=85000, yThresh=ylim_bot, xOrientation='lower', yOrientation='lower', Outer=True, parentGate=singlets)
    
    halfcut_tail = ag.gateThreshold(my_sample, name="right_tail",  xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=180000, parentGate=singlets, orientation='vertical', population='upper')
    ylim_top = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000], limit_threshold=0.2, direction='right',scale='linear') + 25000
    if ylim_top == inf:
        ylim_top = 140000

    PBMC = ag.horisontalPath(my_sample, name="PBMC",
                        xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                        startY=ylim_bot+30000, endY=ylim_top, xboundaries=[55000,140000],
                        yboundaries=[ylim_bot+25000,ylim_top+5000],
                        leftRight=True , direction='both',
                        maxStep=2, phi=0.1, bins=100, sigma=1,
                        scale='linear', parentGate=PBMC_step1)
    
    if fileName:
        fileName="/home/antton/TFM/output/plots/03-PBMCs/"+date_plate+"-"+sampleName+"-PBMCs.jpeg"
    
    ag.backGate(my_sample, population=PBMC, background_population=singlets, xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', T=200, markersize=0.1,  filePlot=fileName)
        
    my_sample.update(ag.AGgate(PBMC, None,'FSC 488/10-A', 'SSC 488/10-A', "PBMC"), QC=False)

    
    ## Gate 4: CD45+ AND CD34+ out of PBMC cells, "Declutter gate" to get rid of CD45- debree 
    xlim_arbitrary = 11000
    right_half = ag.gateThreshold(my_sample, name="right_half", xCol='CD34 PE-Cy7-A' , yCol=CD45, scale='bilog', T=200, thresh=xlim_arbitrary, parentGate=PBMC, orientation='vertical', population='upper')    
        # Get highest density point of CD34 blob to use as reference
    mean, median, sigma, maxVal = ag.axisStats(my_sample(), xCol=CD45, vI=right_half())
    upper_half = ag.gateThreshold(my_sample, name="upper_half", xCol='CD34 PE-Cy7-A', yCol=CD45, scale='bilog', T=200, thresh=maxVal, parentGate=PBMC, orientation='horisontal', population='upper')    
    lowr_halfof_upper_half = ag.gateThreshold(my_sample, name="lower_half_of_upper_half", xCol='CD34 PE-Cy7-A', yCol=CD45, scale='bilog', T=200, thresh=maxVal+3000, parentGate=upper_half, orientation='horisontal', population='lower')   
    xlim_middle = ag.valleySeek(my_sample, xCol='CD34 PE-Cy7-A', interval=[800, 11000], require_local_min=True, scale='bilog', T=200, parentGate=lowr_halfof_upper_half) -1000        
    if fileName:
        fileName="/home/antton/TFM/output/plots/04-CD45high/"+date_plate+"-"+sampleName+"-CD45high.jpeg"
    CD45_high = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal-2900), theta=-40, scale='bilog', T=200, population='upper', parentGate=PBMC, filePlot=fileName)
    
    
    ## Gate 5: CD34+ ratio relative to CD45+
    CD34pos = ag.gateTiltedLine(my_sample, name="CD34+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='lower', parentGate=CD45_high)
    if fileName:
        fileName="/home/antton/TFM/output/plots/05-CD45posCD34pos/"+date_plate+"-"+sampleName+"-CD45posCD34pos.jpeg"
    CD45pos = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='upper', parentGate=CD45_high, filePlot=fileName)

    my_sample.update(ag.AGgate(CD34pos, CD45pos,'CD34 PE-Cy7-A', CD45, "CD34+"), QC=False)


    ## Gate 6: CD3 out of CD45+
        ## Gate 6A: CD3-        
    ylim_middle = ag.valleySeek(my_sample, xCol='CD3 APC-H7-A', interval=[-1000, 7000], require_local_min=True, scale='bilog', T=200, parentGate=CD45_high)    
    CD3neg = ag.gateThreshold(my_sample, name="CD3-", xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD45_high, orientation='horisontal', population='lower')    
    if fileName:
        fileName="/home/antton/TFM/output/plots/06-CD3negBACKGATE/"+date_plate+"-"+sampleName+"-CD3negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD3neg, background_population=PBMC, xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)  
        ## Gate 6B: CD3+ 
    CD3pos = ag.gateThreshold(my_sample, name="CD3+", xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD45_high, orientation='horisontal', population='upper')
    
    my_sample.update(ag.AGgate(CD3pos, CD45pos, CD45, 'CD3 APC-H7-A', "CD3+"), QC=False)
    
    ##Gate 7: CD4, CD8 out of CD3+
        ## Gate 7A: CD4+
    xlim_middle = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True, scale='bilog', T=200, parentGate=CD3pos)
    cd4cd8_step1 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=xlim_middle, parentGate=CD3pos, orientation='vertical', population='upper')    
    CD4pos_final = ag.gatePC(my_sample, name="CD4+", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', center='centroid', adjustAngle=3,widthScale=2.5, scale='bilog', T=200, heightScale=3.5, parentGate=cd4cd8_step1)    
    if fileName:
        fileName="/home/antton/TFM/output/plots/07-CD4posBACKGATE/"+date_plate+"-"+sampleName+"-CD4posBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD4pos_final, background_population=CD3pos, xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName) 
   
    my_sample.update(ag.AGgate(CD4pos_final, CD3pos,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD4+_1"), QC=False)
        ##Gate 7B: CD8+
    ylim_middle = ag.valleySeek(my_sample, xCol='CD8 PerCP-Cy5.5-A', interval=[200, 2000], require_local_min=True, parentGate=CD3pos)
    cd4cd8_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD3pos, orientation='horisontal', population='upper')    
    xlim_middle_upper = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True, scale='bilog', T=200, parentGate=cd4cd8_step2)   
    cd4cd8_step3 = ag.gateThreshold(my_sample, name="separate_middle_x_axis", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=xlim_middle_upper, parentGate=cd4cd8_step2, orientation='vertical', population='lower')   
    CD8pos_final = ag.gatePC(my_sample, name="CD8+", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', center='centroid', adjustAngle=0,widthScale=2, scale='bilog', T=200, heightScale=3.5, parentGate=cd4cd8_step3)
    
    my_sample.update(ag.AGgate(CD8pos_final, CD3pos,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD8+"), QC=False)
    my_sample.update(ag.AGgate(CD4pos_final, CD8pos_final,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD4+_2"), QC=False)

    ## Gate 8: CD19 out of CD3-
        ## Gate 8A: CD19-
    ylim_middle_cd19 = ag.valleySeek(my_sample, xCol='CD19 PE-Texas Red-A', interval=[200, 10000], require_local_min=True, scale='bilog', T=200, parentGate=CD3neg)
    CD19neg = ag.gateThreshold(my_sample, name="CD19-", xCol='CD3 APC-H7-A' , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horisontal', population='lower')
    if fileName:
        fileName="/home/antton/TFM/output/plots/08-CD19negBACKGATE/"+date_plate+"-"+sampleName+"-CD19negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD19neg, background_population=CD3neg, xCol='CD3 APC-H7-A'  , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName) 
        ## Gate 8B: CD19+
    CD19pos = ag.gateThreshold(my_sample, name="CD19+", xCol='CD3 APC-H7-A' , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horisontal', population='upper', filePlot=fileName)
            
    my_sample.update(ag.AGgate(CD19pos, CD45pos, 'CD3 APC-H7-A', 'CD19 PE-Texas Red-A', "CD19+"), QC=False)

    ## Gate 9: CD14 out of CD3-/CD19-    
        ## Gate 7A: CD14- out of CD19-
    ylim_middle_cd14 = ag.valleySeek(my_sample, xCol='CD14 (BV) 605-A', interval=[0, 3000], require_local_min=True, parentGate=CD19neg, scale='bilog')
    CD14neg = ag.gateThreshold(my_sample, name="CD14-", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A', scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg, orientation='horisontal', population='lower')
    if fileName:
        fileName="/home/antton/TFM/output/plots/09-CD14negBACKGATE/"+date_plate+"-"+sampleName+"-CD14negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD14neg, background_population=CD19neg, xCol='CD16 (BV) 786-A' , yCol='CD14 (BV) 605-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)     
        ## Gate 7B: CD14+ out of CD19-        
    CD14pos = ag.gateThreshold(my_sample, name="CD14+", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A', scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg, orientation='horisontal', population='upper')
        
    my_sample.update(ag.AGgate(CD14pos, CD45pos, 'CD19 PE-Texas Red-A', 'CD14 (BV) 605-A', "CD14+"), QC=False)

    ## Gate 10: CD16 and CD56 values out of CD14-

    aprox_xlim_middle = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True, scale='bilog', T=200, parentGate=CD14neg)    
    left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle, parentGate=CD14neg, orientation='vertical', population='lower')    
    aprox_ylim_middle = ag.valleySeek(my_sample, xCol='CD56 (BV) 650-A', interval=[700, 1300], require_local_min=True, scale='bilog', T=200, parentGate=CD14neg)    
    lower_left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_ylim_middle, parentGate=left_half, orientation='horisontal', population='lower')    
    main_ylim_middle = ag.densityDelimitation(my_sample, xCol='CD56 (BV) 650-A', parentGate=lower_left_half, interval=[0,2000], limit_threshold=0.07, direction='right',scale='linear')    
    if main_ylim_middle == inf:
        main_ylim_middle = aprox_ylim_middle    
    lower_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg, orientation='horisontal', population='lower')    
    upper_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg, orientation='horisontal', population='upper')    
    aprox_xlim_middle2 = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True, scale='bilog', T=200, parentGate=lower_half)    
    CD16neg_56neg_step1 = ag.gateThreshold(my_sample, name="lower_left_side", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2, parentGate=lower_half, orientation='vertical', population='lower')    
    main_xlim_middle = ag.densityDelimitation(my_sample, xCol='CD16 (BV) 786-A', parentGate=CD16neg_56neg_step1, interval=[0,2000], limit_threshold=0.07, direction='right',scale='linear')
    if main_xlim_middle == inf:
        main_xlim_middle = aprox_xlim_middle    
    CD16neg_56neg_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle, parentGate=CD16neg_56neg_step1, orientation='vertical', population='lower')
        ## CD16-CD56-, a.k.a. LIN NEGATIVE
    CD16neg_56neg_final = ag.gatePC(my_sample, name="CD16-CD56-", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', center='centroid', adjustAngle=0,widthScale=2, scale='bilog', T=200, heightScale=3.5, parentGate=CD16neg_56neg_step2)
    if fileName:
        fileName="/home/antton/TFM/output/plots/10-CD16CD56/10A-CD16neg_56negBACKGATE/"+date_plate+"-"+sampleName+"-CD16neg_56negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD16neg_56neg_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
    central_section = ag.horisontalPath(my_sample, name="hor_path", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', population='lower', startY=3000, endY=12000,  xboundaries=[0,500000], yboundaries=[2000,15000], leftRight=True , direction='both', maxStep=2, phi=0.1, bins=100, sigma=1, scale='bilog', T=200, parentGate=upper_half)   
         ## CD16++CD56++
    CD16pospos_56pospos_final = ag.gateThreshold(my_sample, name="CD16+CD56+", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle, parentGate=central_section, orientation='vertical', population='upper')
    if fileName:
        fileName="/home/antton/TFM/output/plots/10-CD16CD56/10B-CD16pospos_56posposBACKGATE/"+date_plate+"-"+sampleName+"-CD16pospos_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16pospos_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
    
    my_sample.update(ag.AGgate(CD16pospos_56pospos_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16++CD56++"), QC=False)
        ## CD16-CD56++    
    CD16neg_56pospos_final = ag.horisontalPath(my_sample, name="CD16-CD56++", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', population='upper', startY=3000, endY=12000, xboundaries=[0,500000], yboundaries=[2000,15000], leftRight=True , direction='both', maxStep=2, phi=0.1, bins=100, sigma=1, scale='bilog', T=200, parentGate=upper_half)
    if fileName:
        fileName="/home/antton/TFM/output/plots/10-CD16CD56/10C-CD16neg_56posposBACKGATE/"+date_plate+"-"+sampleName+"-CD16neg_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16neg_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
        
    my_sample.update(ag.AGgate(CD16neg_56pospos_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16-CD56++"), QC=False)
        ## CD16++CD56-
    CD16pospos_56neg_final = ag.gateThreshold(my_sample, name="CD16++CD56-", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2, parentGate=lower_half, orientation='vertical', population='upper')
 
    my_sample.update(ag.AGgate(CD16pospos_56neg_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16++CD56-"), QC=False)
    
    ## Gate 11: Lin-CD34+ out of CD16-CD56-/CD14-/CD19-/CD3-/CD45+
    
    if fileName:
        fileName="/home/antton/TFM/output/plots/11-LinnegCD34pos/"+date_plate+"-"+sampleName+"-LinnegCD34pos.jpeg"
    CD34pos_step1 = ag.gateTiltedLine(my_sample, name="tilted_gate_cd34", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='lower', parentGate=CD16neg_56neg_final, filePlot=fileName)
    linnegCD34pos = ag.gatePC(my_sample, name="Lin-CD34+", xCol='CD34 PE-Cy7-A' , yCol=CD45, center='centroid', adjustAngle=3,widthScale=2.5, heightScale=3.5, parentGate=CD34pos_step1)
    
    my_sample.update(ag.AGgate(linnegCD34pos, CD45pos,'CD34 PE-Cy7-A' , CD45,"linneg_cd34pos"), QC=True, xlim=[-300,200000], ylim=[-500,100000])

    ## Gate 12: CD38 out of Lin-CD34+ (PLACEHOLDER until new CD38 gate is developed)
    
    gated_df = my_sample.fcsDF.loc[CD34pos()].copy()
    cd38_values_list = sorted(list(gated_df['CD38 (BV) 421-A']))
    if cd38_values_list:  #If the list is not empty
        index = int(len(cd38_values_list)*0.3)
        ylim_middle = cd38_values_list[index]
    else:
        return my_sample

        ## Gate 11A: CD38-
    if fileName:
        fileName="/home/antton/TFM/output/plots/12-CD38neg/"+date_plate+"-"+sampleName+"-CD38neg.jpeg"
    CD38neg= ag.gateThreshold(my_sample, name="CD38-", xCol='CD34 PE-Cy7-A' , yCol='CD38 (BV) 421-A',  scale='bilog', T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horisontal', population='lower', filePlot=fileName)
        ## Gate 11B: CD38+
    CD38pos = ag.gateThreshold(my_sample, name="CD38+", xCol='CD34 PE-Cy7-A' , yCol='CD38 (BV) 421-A',  scale='bilog', T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horisontal', population='upper')
    
    ## Gate 13: HSCs, MPPs and MLPs out of CD38-   
    ylim = ag.densityDelimitation(my_sample, xCol='CD90 PE (R-phycoerythrin)-A', interval=[200,5000], limit_threshold=0.5, direction='right',scale='bilog',T=200, parentGate=CD38neg)    
    xlim = ag.densityDelimitation(my_sample, xCol='CD45RA FITC-A', interval=[-100,300], limit_threshold=0.05, direction='right',scale='bilog',T=200, parentGate=CD38neg)   
    if ylim == inf:  # Plan B...
        right_hand= ag.gateThreshold(my_sample, name="remove_clutter_3", xCol='CD45RA FITC-A', yCol='CD90 PE (R-phycoerythrin)-A', scale='bilog', T=200, thresh=xlim, parentGate=CD38neg, orientation='vertical', population='upper')        
        ylim = ag.densityDelimitation(my_sample, xCol='CD90 PE (R-phycoerythrin)-A', interval=[200,5000], limit_threshold=0.2, direction='right',scale='bilog',T=200, parentGate=right_hand)
        if ylim == inf:
            ylim = 5000  # Plan C...
    if fileName:
        fileName="/home/antton/TFM/output/plots/13-HSC_MLP_MPP/"+date_plate+"-"+sampleName+"-HSC_MLP_MPP.jpeg"
    HSC, irrelevant, MLP, MPP = ag.quadGate(my_sample, names=['HSC', 'NA', 'MLP', 'MPP'], xCol='CD45RA FITC-A', yCol='CD90 PE (R-phycoerythrin)-A', xThresh=xlim, yThresh=ylim, scale='bilog', T=200, parentGate=CD38neg, filePlot=fileName)

    my_sample.update(ag.AGgate(HSC, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"HSC_1"), QC=False)
    my_sample.update(ag.AGgate(HSC, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"HSC_2"), QC=False)
    my_sample.update(ag.AGgate(MLP, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MLP_1"), QC=False)
    my_sample.update(ag.AGgate(MLP, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MLP_2"), QC=False)
    my_sample.update(ag.AGgate(MPP, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MPP_1"), QC=False)
    my_sample.update(ag.AGgate(MPP, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MPP_2"), QC=False)
    
    ## Gate 14: non-BNK (CD10neg) and B-NK out of CD38+
    ylim_cd10 = ag.densityDelimitation(my_sample, xCol='CD10 APC (Allophycocyanin)-A', parentGate=halfcut, interval=[-100,800], limit_threshold=0.2, direction='right',scale='bilog',T=200)
    
        ## B-NK progenitors, CD10+?
    if fileName:
        fileName="/home/antton/TFM/output/plots/14-B_NK/"+date_plate+"-"+sampleName+"-B_NK.jpeg"        
    BNK_prog = ag.gateCorner(my_sample, name="BNK", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', xThresh=200, yThresh=ylim_cd10, xOrientation='upper', yOrientation='upper', Outer=False, scale='bilog', T=200, parentGate=CD38pos, filePlot=fileName)
    
    my_sample.update(ag.AGgate(BNK_prog, CD38pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_1"), QC=False) 
    my_sample.update(ag.AGgate(BNK_prog, CD34pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_2"), QC=False)
    my_sample.update(ag.AGgate(BNK_prog, linnegCD34pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_3"), QC=False)

        ## CD10-
    nonBNK= ag.gateThreshold(my_sample, name="CD10-", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=CD38pos, orientation='horisontal', population='lower')

    ## Gate 15: CD10 out of MLP
    if fileName:
        fileName="/home/antton/TFM/output/plots/15-CD10pos(MLP)/"+date_plate+"-"+sampleName+"-CD10pos(MLP).jpeg"     
    CD10pos_MLP = ag.gateThreshold(my_sample, name="cd10+_MLP", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=MLP, orientation='horisontal', population='upper',filePlot=fileName)
    CD10neg_MLP = ag.gateThreshold(my_sample, name="cd10-_MLP", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=MLP, orientation='horisontal', population='lower')

    my_sample.update(ag.AGgate(CD10pos_MLP, MLP,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_1"), QC=False)
    my_sample.update(ag.AGgate(CD10pos_MLP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_2"), QC=False)
    my_sample.update(ag.AGgate(CD10pos_MLP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_3"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, MLP,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_1"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_2"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_3"), QC=False)    
    
    ## Gate 16: CNPs, GMPs and MEPs out of non-BNK CD38+ cells
    ylim = ag.densityDelimitation(my_sample, xCol='CD135 (BV) 711-A', interval=[-100,800], limit_threshold=0.5, direction='left',scale='bilog',T=200, parentGate=nonBNK)
    
    xlim_cut_tail = ag.densityDelimitation(my_sample, xCol='CD45RA FITC-A', interval=[-100,800], limit_threshold=0.2, direction='right',scale='bilog',T=200, parentGate=nonBNK)
    if fileName:
        fileName="/home/antton/TFM/output/plots/16-CMP_GMP_MEP/"+date_plate+"-"+sampleName+"-CMP_GMP_MEP.jpeg"
    CMP, GMP, irrelevant, MEP = ag.quadGate(my_sample, names=['CMP', 'GMP', 'NA', 'MEP'], xCol='CD45RA FITC-A', yCol='CD135 (BV) 711-A', xThresh=xlim_cut_tail, yThresh=ylim, scale='bilog', T=200,  parentGate=nonBNK, filePlot=fileName)
    
    my_sample.update(ag.AGgate(CMP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CMP_1"), QC=False)
    my_sample.update(ag.AGgate(CMP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CMP_2"), QC=False)
    my_sample.update(ag.AGgate(MEP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"MEP_1"), QC=False)
    my_sample.update(ag.AGgate(MEP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"MEP_2"), QC=False)
    my_sample.update(ag.AGgate(GMP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"GMP_1"), QC=False)
    my_sample.update(ag.AGgate(GMP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"GMP_2"), QC=False)
    
    return my_sample

if __name__ == '__main__':
    

    filepaths = get_filepaths(maxDate=date(2021,3,12), maxSampleNum=2020)  # Include only samples with number below maxSampleNum

    out_folders = out_folder_list()
    if not os.path.exists("/home/antton/TFM/output/plots/10-CD16CD56/"):
            os.mkdir("/home/antton/TFM/output/plots/10-CD16CD56/")
    for folder in out_folders:
        if not os.path.exists(folder):
            os.mkdir(folder)
    #target_file =  ['/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/181029 CB 60min/B8 27-2.fcs']
    #Define experiment object
    #print(filepaths.index('/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190610 CB/C7 529.fcs'))
    #filepaths = ['/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190211 CB/A5 162.fcs', '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/200210 CB/E1 1435.fcs']
    CB_exp=ag.AGExperiment(filepaths, filters=['fcs'], mask=['30min','45min','Neg','test'], experiment_name="cord_blood_experiment", flourochrome_area_filter=True, QC=True, QCbins=128)
    CB_exp.apply(gateFullDataset)
    CB_exp.printExperiment("/home/antton/TFM/output/aligater_output_test/cblood_test_"+str(date.today())+".txt")
