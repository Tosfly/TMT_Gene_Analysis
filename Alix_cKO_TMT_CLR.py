#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 16:03:13 2022

@author: yi-zhiwang
"""
import os
os.chdir("/Users/yi-zhiwang/Projects/cLTP/Raw_Data/Alix_cKO")
import numpy
import math
import pandas as pd
import seaborn as sns
import statistics
#from numba import jit, cuda
#######################################################################################################################################################################
## Input/Output codes
def excel_reader(file_name,sheet_Num):
    import xlrd
    data = xlrd.open_workbook(file_name)
    table = data.sheets()[sheet_Num]
    return table

def GN_select(names):
    new = []
    for i in names:
        if i.find('GN=') > 0:
           s = i.split('GN=')[1]
           e = s.find(' ')
           new.append(s[:e])
        else:
            new.append(i)
    return new


def GN_single(name):
    if name.find('GN=') > 0:
        s = name.split('GN=')[1]
        e = s.find(' ')
        return s[:e]
    else:
        return name.split(' OS=')[0]


def report(name,data):
##name = 'name.txt'
    name = name + '.csv'
    new = open(name,'w')
    for i in data:
        new.write(str(i)+'\n')
    return new.close()

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
        
def PeptidesFor(Locus,rawdata):
    Peptides = IDPeptideExtr(Locus,rawdata)
##Extract peptide information:sequece, Normal_TMT value, scan number
    Head = [['SEQUENCE','KO','KO','KO','KO','KO','WT','WT','WT','WT','WT','ScanNum','TMT_purity','Signal-noise']]
    for i in Peptides:
        Head.append([i[2],i[3],i[5],i[7],i[9],i[11],i[13],i[15],i[17],i[19],i[21],i[28],i[24],i[25]])
    return Head

def GeneFinder(genename,rawdata):
##e.g.,rawdata = DvA_Raw
    AllIDs_raw = [x[1] for x in rawdata if GN_single(x[61]) == genename]
    AllIDs_com = [x for x in AllIDs_raw if not 'Reverse_' in x]
    AllIDs = [x for x in AllIDs_com if not 'contaminant_' in x]  
    Temp = []
    for i in AllIDs:
        Temp += PeptidesFor(i,rawdata)
    Head = ['SEQUENCE','KO','KO','KO','KO','KO','WT','WT','WT','WT','WT','ScanNum','TMT_purity','Signal-noise']
    return [Head] + [x for x in Temp if not x == Head]
        
def ReplaceZero(Rowdata):
##replace 0 with 1
    Temp = []
    for i in Rowdata:
        if i == 0:
            Temp += [1]
        else:
            Temp += [i]
    return Temp

def Close(Rowdata):
##Close the data
    ReporterIon = Rowdata[1:11]
    Total = sum(ReporterIon)
    if Total != 0:
        return [Rowdata[0]] + [x/Total for x in ReporterIon] + Rowdata[11:]
    else:
        return 'NA'

def CLRrow(Rowdata):
##CLR on rowdata
    CloseData = Close(Rowdata)
    if CloseData != 'NA':
        ReporterIon = CloseData[1:11]
        GeoMean = sum([math.log2(x) for x in ReporterIon])/10
        return [Rowdata[0]] + [(math.log2(x) - GeoMean) for x in ReporterIon] + Rowdata[11:]
    else:
        return 'NA'
        
def PeptideCleanUp(genename,rawdata):
## remove multiple counted peptide by ScanNum
    Allpeptides_raw = GeneFinder(genename,rawdata)
    Allpeptides = [x for x in Allpeptides_raw if type(x[0]) == str] 
    ScanNum = list(set([x[11] for x in Allpeptides[1:]]))
    Temp = []
    for i in ScanNum:
        Peptide = [x for x in Allpeptides if x[11] == i][0]
        Temp.append(Peptide)
    NoZero = [ReplaceZero(x) for x in Temp]
    return [Allpeptides[0]] + [CLRrow(x) for x in NoZero]

def GeneExpress(genename,rawdata):
    GeneData = PeptideCleanUp(genename,rawdata)[1:]
    cKO1 = sum([x[1] for x in GeneData])
    cKO2 = sum([x[2] for x in GeneData])
    cKO3 = sum([x[3] for x in GeneData])
    cKO4 = sum([x[4] for x in GeneData])
    cKO5 = sum([x[5] for x in GeneData])
    WT1 = sum([x[6] for x in GeneData])
    WT2 = sum([x[7] for x in GeneData])
    WT3 = sum([x[8] for x in GeneData])
    WT4 = sum([x[9] for x in GeneData])
    WT5 = sum([x[10] for x in GeneData])
    return [genename,cKO1,cKO2,cKO3,cKO4,cKO5,WT1,WT2,WT3,WT4,WT5]

#@jit(target_backend='cuda')

def ValidGene(dataset):
    AllGeneNames = list(set([GN_single(x[61]) for x in dataset if not x[61] == '']))
    AllIDs_ori = [x[1] for x in dataset if GN_single(x[61]) in AllGeneNames]
    AllIDs_com = [x for x in AllIDs_ori if not 'Reverse_' in x]
    AllIDs = [x for x in AllIDs_com if not 'contaminant_' in x]
    return list(set([GN_single(x[61]) for x in dataset if x[1] in AllIDs]))
                
def TMTRefor(rawdata,genelist):
    Temp = [['GeneName','cKO1','cKO2','cKO3','cKO4','cKO5','WT1','WT2','WT3','WT4','WT5']]
    for i in genelist:    
        Info = GeneExpress(i, rawdata)
        if sum(Info[1:]) != 0:
            Temp.append(Info)
    return Temp

AllGeneNames = ValidGene(TMT)

report('AllGeneNames',AllGeneNames)

DataReform_1 = TMTRefor(TMT,AllGeneNames[:1000])
ComCsv('CLR_AlixcKOvsW_1',DataReform_1)

DataReform_2 = TMTRefor(TMT,AllGeneNames[1000:2000])
ComCsv('CLR_AlixcKOvsW_2',DataReform_2)

DataReform_3 = TMTRefor(TMT,AllGeneNames[2000:3000])
ComCsv('CLR_AlixcKOvsW_3',DataReform_3)

DataReform_4 = TMTRefor(TMT,AllGeneNames[3000:4000])
ComCsv('CLR_AlixcKOvsW_4',DataReform_4)

DataReform_5 = TMTRefor(TMT,AllGeneNames[4000:5000])
ComCsv('CLR_AlixcKOvsW_5',DataReform_5)

DataReform_6 = TMTRefor(TMT,AllGeneNames[5000:6000])
ComCsv('CLR_AlixcKOvsW_6',DataReform_6)

DataReform_7 = TMTRefor(TMT,AllGeneNames[6000:7000])
ComCsv('CLR_AlixcKOvsW_7',DataReform_7)

DataReform_8 = TMTRefor(TMT,AllGeneNames[7000:])
ComCsv('CLR_AlixcKOvsW_8',DataReform_8)
        
            





























