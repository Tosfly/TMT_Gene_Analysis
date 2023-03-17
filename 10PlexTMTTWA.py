#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:03:13 2022

@author: yi-zhiwang
"""
import os
os.chdir("/home/ywd617/data_analysis")
import numpy
import math
import pandas as pd
import seaborn as sns
import statistics
#######################################################################################################################################################################
## Input/Output codes
def excel_reader(file_name,sheet_Num):
    import xlrd
    data = xlrd.open_workbook(file_name)
    table = data.sheets()[sheet_Num]
    return table

def GN_single(name):
    if name.find('GN=') > 0:
        s = name.split('GN=')[1]
        e = s.find(' ')
        return s[:e]
    else:
        return name.split(' OS=')[0]    
    
def ComCsv(filename,data):
##generate a new csv file 
    import csv
    with open((filename + '.csv'),'w') as f:
        f_csv = csv.writer(f,lineterminator='\n')
        for i in data:
            f_csv.writerow(i)
    return print('The data is generated as csv file')

def Ex2Dat(filename,number):
    test = excel_reader(filename,number)
    temp = []
    for i in range(len(test.col_values(0))):
        temp.append(test.row_values(i))
    return(temp)

def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split('>')[1:]

MouseData_raw = read_FASTA_strings('MouseDatabase.fasta')
MouseData = [x for x in MouseData_raw if not 'Reverse' in x]

def GeneItem(entry):
##find a gene in a entry
##entry = read_FASTA_strings('MouseDatabase.fasta') [0]
    if 'GN=' in entry:
        Fullname = entry.split('GN=')[1]
        temp = Fullname.split(' ')[0]
        return temp.split('\n')[0]
    else:
        return 'Uncharacterized protein X'

def ID2GN(Locus):
    Temp = [x for x in MouseData if Locus in x][0]
    if Temp.find('GN=') > 0:
        s = Temp.split('GN=')[1]
        if ' ' in s:
            e = s.find(' ')
            return s[:e]
        elif '\n' in s:
            return s.split('\n')[0]
        else:
            return Locus
    else:
        return Locus

#############################################################################################################
# #Load Excel files
TMT = Ex2Dat('AlixcKO.xlsx',0)
#################################################################################################################################
#Extract all isoforms under the same genename
def IDPeptideExtr(Locus,rawdata):
##find a ID of an isoform, and extract all peptides
    Isoform = [x for x in rawdata if x[1] == Locus][0]
    SpecCount = Isoform[2]
    IDIndex = rawdata.index(Isoform)
    while rawdata[IDIndex][0] == 'P':
        IDIndex += 1
    End = int(IDIndex + SpecCount)
    return rawdata[IDIndex:End]
        
def PeptidesFor(Locus,rawdata,group):
    Peptides = IDPeptideExtr(Locus,rawdata)
##Extract peptide information:sequece, Normal_TMT value, scan number
    if group == 'SPNs':
        Head = [['SEQUENCE','dSPN','dSPN','dSPN','dSPN','dSPN','iSPN','iSPN','iSPN','iSPN','iSPN','ScanNum']]
    elif group == 'cKO':
        Head = [['SEQUENCE','KO','KO','KO','KO','KO','WT','WT','WT','WT','WT','ScanNum']]
    for i in Peptides:
        Head.append([i[2],i[4],i[6],i[8],i[10],i[12],i[14],i[16],i[18],i[20],i[22],i[28]])
    return Head

def GeneFinder(genename,rawdata,group):
##e.g.,rawdata = DvA_Raw
    AllIDs = [x[1] for x in rawdata if GN_single(x[61]) == genename]
    Temp = []
    for i in AllIDs:
        Temp += PeptidesFor(i,rawdata,group)
    if group == 'SPNs':
        Head = ['SEQUENCE','dSPN','dSPN','dSPN','dSPN','dSPN','iSPN','iSPN','iSPN','iSPN','iSPN','ScanNum']
    elif group == 'cKO':
        Head = ['SEQUENCE','KO','KO','KO','KO','KO','WT','WT','WT','WT','WT','ScanNum']
    return [Head] + [x for x in Temp if not x == Head]
        
def PeptideCleanUp(genename,rawdata,group):
## remove multiple counted peptide by ScanNum
    Allpeptides_raw = GeneFinder(genename,rawdata,group)
    Allpeptides = [x for x in Allpeptides_raw if type(x[0]) == str] 
    ScanNum = list(set([x[11] for x in Allpeptides[1:]]))
    Temp = [Allpeptides[0]]
    for i in ScanNum:
        Peptide = [x for x in Allpeptides if x[11] == i][0]
        Temp.append(Peptide)
    return Temp

def ANOVAFormater(genename,rawdata,group):
## format data into format for two-way ANOVA analysis
    PeptideClean = PeptideCleanUp(genename,rawdata,group)
    Subset = PeptideClean[1:]
    Head = [['Genotype'] + PeptideClean[0][1:11]]
    for i in range(len(Subset)):
        Index = ['P' + str(i + 1)]
        Head.append(Index + Subset[i][1:11])
    return numpy.array(Head).transpose()
    ##################################################################################################################################
##Two-way ANOVA test on peptide levels
def PepFold(genename,rawdata,group):
##calculate the all isoform changes of one protein
    data = ANOVAFormater(genename,rawdata,group)
    Genotype_A = [statistics.median(x[1:].astype(float)) for x in data[1:6]]
    Genotype_B = [statistics.median(x[1:].astype(float)) for x in data[6:]]
    return sum(Genotype_A)/sum(Genotype_B)

def GeneExFold(rawdata,group):
##Calculate the fold changes of all gene products
    GeneListRaw = [x for x in rawdata if x[0] == 'P']
    GeneListCon = [x for x in GeneListRaw if not 'contaminant' in x[1]]
    GeneList = [x for x in GeneListCon if not 'Reverse' in x[1]]
    GeneName = list(set([GN_single(x[61]) for x in GeneList if (x[2] > 1) and (x[3] > 1)]))
    Temp = [['GeneName','Fold_HET/WT']]
    for i in GeneName:
        Temp.append([i] + [PepFold(i,rawdata,group)])
    return Temp
    
    
def TwoWayANOVA(genename,rawdata,group):
##TwoWayANOVA for each gene
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    import scipy.stats as stats
    from bioinfokit.analys import stat 
    data = ANOVAFormater(genename,rawdata,group)
    d = pd.DataFrame(data[1:],columns = data[0])
    d[data[0][1:]] = d[data[0][1:]].apply(pd.to_numeric)
    d_melt = pd.melt(d, id_vars=['Genotype'], value_vars=data[0][1:])
    d_melt.columns = ['Genotype', 'P', 'value']
    model = ols('value ~ C(Genotype) + C(P) + C(Genotype):C(P)', data=d_melt).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    res = stat()
    res.anova_stat(df=d_melt, res_var='value', anova_model='value~C(Genotype)+C(P)+C(Genotype):C(P)')
    w, pvalue = stats.shapiro(res.anova_model_out.resid)
    return [anova_table["PR(>F)"][0],anova_table["PR(>F)"][1],anova_table["PR(>F)"][2],pvalue]
    
def TMTANOVA(rawdata,group):
##Calculate the ANOVA p-value for each gene:
    GeneListRaw = [x for x in rawdata if x[0] == 'P']
    GeneListCon = [x for x in GeneListRaw if not 'contaminant' in x[1]]
    GeneList = [x for x in GeneListCon if not 'Reverse' in x[1]]
    GeneName = list(set([GN_single(x[61]) for x in GeneList if (x[2] > 2) and (x[3] > 2)]))
    Temp = [['GeneName','C(Genotype)','C(P)','C(Genotype):C(P)','Shapiro-Wilk','Fold_HET/WT']]
    for i in GeneName:
        Temp.append([i] + TwoWayANOVA(i,rawdata,group) + [PepFold(i,rawdata,group)])
    return Temp

TMT_ANOVA = TMTANOVA(TMT,'cKO')





























# -*- coding: utf-8 -*-

