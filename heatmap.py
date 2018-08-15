# -*- coding: utf-8 -*-
"""
Spyder Editor
Analysis of THP sampels from brian_wilhem
This is a temporary script file.
"""



#from __future__ import print_function
import sys
from commands import *
import os
import re
import numpy as np
import pandas as pd
from pandas import Series, DataFrame 
import csv as csv
from itertools import groupby
import fnmatch
import subprocess
import glob2
from scipy.stats import mstats
from math import exp, expm1, log
import subprocess
import string
import datetime
#import sqlite3 as lite
import MySQLdb as my


#def f_drop_mysql_table(f_table_id):    
#    sql_cmd="""DROP TABLE *"""
#    return (sql_cmd)    

def f_create_mysql_table(f_table_id):    
    if f_table_id in ['circexp','htseq','featurecount']:
        sql_cmd="""CREATE TABLE `"""+f_table_id+"""` (
          `Project_ID` varchar(100) DEFAULT NULL,
          `Project_Type` varchar(100) DEFAULT NULL,
          `Sample` varchar(100) DEFAULT NULL,
          `File_name` varchar(100) NOT NULL,
          `Location` varchar(512) NOT NULL,
           PRIMARY KEY (`file_name`,`Location`)
        )"""
        print "Performing query: "+sql_cmd
        
    elif f_table_id in ['chipseq','deseq']:
        sql_cmd="""CREATE TABLE `"""+f_table_id+"""` (
          `Project_ID` varchar(100) DEFAULT NULL,
          `Project_Type` varchar(100) DEFAULT NULL,
          `Sample` varchar(100) NOT NULL,
          `Location` varchar(512) NOT NULL,
           PRIMARY KEY (`Sample`,`Location`)
        )"""
        print "Performing query: "+sql_cmd
     
    elif f_table_id in ['raw_fastq']:
        sql_cmd="""CREATE TABLE `"""+f_table_id+"""` (
          `Project_ID` varchar(100) DEFAULT NULL,
          `Project_Type` varchar(100) DEFAULT NULL,
          `File_name` varchar(100) NOT NULL,
          `Location` varchar(512) NOT NULL,
           PRIMARY KEY (`file_name`,`Location`)
        )"""
        print "Performing query: "+sql_cmd
        
    else: #bam and bwa alignment
        sql_cmd="""CREATE TABLE `"""+f_table_id+"""` (
          `Project_ID` varchar(100) DEFAULT NULL,
          `Project_Type` varchar(100) DEFAULT NULL,
          `Sample` varchar(100) NOT NULL,
          `File_name` varchar(100) NOT NULL,
          `Protocol` varchar(100) DEFAULT NULL,
          `Algorithm` varchar(10) DEFAULT NULL,
          `Location` varchar(512) NOT NULL,
           PRIMARY KEY (`file_name`,`Location`)
        )"""
        print "Performing query: "+sql_cmd
        
    return (sql_cmd)    

    
def f_load_data_mysql(f_data):   
    #sql_cmd="""LOAD DATA INFILE '~/estorage/data/chipseq.txt' 
    table_name=f_data.split(".")[0]
    sql_cmd="""LOAD DATA LOCAL INFILE '/home/dbeck/estorage/data/"""+f_data+"""' 
    INTO TABLE """+table_name+""" 
    FIELDS TERMINATED BY '\t'
    LINES TERMINATED BY '\n'
    IGNORE 1 LINES;"""
    print "Performing query: "+sql_cmd
    return (sql_cmd)    


def f_project_id(path):
    analysis_type=['star','circexplorer2','htseq','featurecount','bwa']
    project_id=[]
    for i in range(len(path)):
        row=path['Project_Type'][i]
        if row=="public" or path['location'][i].split('/')[2] in analysis_type:
            project_id.append(path['location'][i].split('/')[1])    
#        elif row=="labdata":
#            project_id.append(path['location'][i].split('/')[1]+"/"+path['location'][i].split('/')[2]+"/"+path['location'][i].split('/')[3])
        else: 
            #project_id.append(path['location'][i].split('/')[1]+"/"+path['location'][i].split('/')[2])
            project_id.append(path['location'][i].split('/')[2]) # for labdata
     
    return project_id       


list_files = pd.read_csv("file_structure_Diego_database.txt", sep="\t", header = None , low_memory=False)
list_files.rename(columns={0: 'location'}, inplace=True)

#1. Raw data
raw_fastq=list_files[list_files['location'].str.contains("fastq")].reset_index(drop=True)
raw_fastq['file_name']=raw_fastq['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
raw_fastq['temp']=raw_fastq['file_name'].apply(lambda x: x.split(".")[len(x.split("."))-1])

raw_fastq=raw_fastq[((raw_fastq['temp']=='fastq') | (raw_fastq['temp']=='fq'))] 
raw_fastq=raw_fastq[raw_fastq['file_name'].str.contains("\.")]


raw_fastq['Project_Type']=raw_fastq['location'].apply(lambda x: x.split("/")[0])
raw_fastq['Project_ID']=raw_fastq['location'].apply(lambda x: x.split("/")[1])
raw_fastq=raw_fastq[['Project_ID','Project_Type','file_name','location']]
raw_fastq.to_csv('raw_fastq.txt',index=False,header=True,sep='\t')   

#Experiments
#2. chipseq
chipseq=list_files[list_files['location'].str.contains("PeakCallers")].reset_index(drop=True)
chipseq['location']=chipseq['location'].apply(lambda x: x.split("PeakCallers/")[0])+"PeakCallers/"
chipseq=chipseq[chipseq['location'].str.contains("/PeakCallers/")].reset_index(drop=True)
chipseq['Sample']=chipseq['location'].apply(lambda x: x.split("/")[len(x.split("/"))-3])

chipseq['Project_Type']=chipseq['location'].apply(lambda x: x.split("/")[0])
chipseq['Project_ID']=f_project_id(chipseq)
#project_id
chipseq=chipseq[['Project_ID','Project_Type','Sample','location']].drop_duplicates().reset_index(drop=True)
chipseq.to_csv('chipseq.txt',index=False,header=True,sep='\t')   

#alignment
#3. star 
star_bam=list_files[list_files['location'].str.contains("star")].reset_index(drop=True)
star_bam['file_name']=star_bam['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
star_bam['Sample']=star_bam['location'].apply(lambda x: x.split("/")[len(x.split("/"))-2])
star_bam['temp']=star_bam['file_name'].apply(lambda x: x.split(".")[len(x.split("."))-1])

star_bam=star_bam[(star_bam['temp']=='bam') | (star_bam['temp']=='sam')].reset_index(drop=True) #
star_bam=star_bam[star_bam['file_name'].str.contains("\.")].reset_index(drop=True)

star_bam['Algorithm']="STAR"
star_bam['Protocol']="RNA-seq"
star_bam['Project_Type']=star_bam['location'].apply(lambda x: x.split("/")[0])
star_bam['Project_ID']=f_project_id(star_bam)

star_bam=star_bam[['Project_ID','Project_Type','Sample','file_name','Protocol','Algorithm','location']]
star_bam.to_csv('star_bam.txt',index=False,header=True,sep='\t')   

#4. bwa
bwa_bam=list_files[(list_files['location'].str.contains("bwa")) | (list_files['location'].str.contains("PeakCallers"))].reset_index(drop=True)
bwa_bam['Protocol']="ChIP-seq_Chrom-Access"
bwa_bam['file_name']=bwa_bam['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
bwa_bam['Sample']=bwa_bam['location'].apply(lambda x: x.split("/")[len(x.split("/"))-2])
bwa_bam['temp']=bwa_bam['file_name'].apply(lambda x: x.split(".")[len(x.split("."))-1])

bwa_bam=bwa_bam[(bwa_bam['temp']=='bam') | (bwa_bam['temp']=='sam')].reset_index(drop=True) #
bwa_bam=bwa_bam[bwa_bam['file_name'].str.contains("\.")].reset_index(drop=True)

bwa_bam['Algorithm']="BWA"

bwa_bam['Project_Type']=bwa_bam['location'].apply(lambda x: x.split("/")[0])
bwa_bam['Project_ID']=f_project_id(bwa_bam)

bwa_bam=bwa_bam[['Project_ID','Project_Type','Sample','file_name','Protocol','Algorithm','location']]
bwa_bam.to_csv('test',index=False,header=True,sep='\t')   
bwa_bam.to_csv('bwa_bam.txt',index=False,header=True,sep='\t')   

#5. deseq
deseq=list_files[list_files['location'].str.contains("DiffExpr")].reset_index(drop=True)
deseq['location']=deseq['location'].apply(lambda x: x.split("DiffExpr",)[0])
deseq['Sample']=deseq['location'].apply(lambda x: x.split("/")[len(x.split("/"))-2])

#deseq['temp']=deseq['location'].apply(lambda x: x.split("DiffExpr",)[0])

#deseq=deseq[-deseq['temp'].str.contains('deseq')].reset_index(drop=True)

#deseq['sample']=deseq['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])

#star_bam=star_bam[star_bam['temp']=='bam'].reset_index(drop=True)
#star_bam=star_bam[star_bam['file_name'].str.contains("\.")]
deseq['Project_Type']=deseq['location'].apply(lambda x: x.split("/")[0])
deseq['Project_ID']=f_project_id(deseq)
deseq=deseq[['Project_ID','Project_Type','Sample','location']].drop_duplicates().reset_index(drop=True)
deseq.to_csv('deseq.txt',index=False,header=True,sep='\t')   

#6. circexplores2
#circexp=list_files[list_files['location'].str.contains("circexplorer2")].reset_index(drop=True)
circexp=list_files[list_files['location'].str.contains("circ_fusion.txt")].reset_index(drop=True)

circexp['file_name']=circexp['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
circexp['Sample']=circexp['location'].apply(lambda x: x.split("/")[len(x.split("/"))-3])

circexp['Project_Type']=circexp['location'].apply(lambda x: x.split("/")[0])
circexp['Project_ID']=f_project_id(circexp)
circexp=circexp[['Project_ID','Project_Type','Sample','file_name','location']]
circexp.to_csv('circexp.txt',index=False,header=True,sep='\t')   

#7. htseq
htseq=list_files[list_files['location'].str.contains(".htseq.count")].reset_index(drop=True)
htseq=htseq[-htseq['location'].str.contains("summary")].reset_index(drop=True)

htseq['file_name']=htseq['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
htseq['Sample']=htseq['file_name'].apply(lambda x: x.split(".")[0])

htseq['Project_Type']=htseq['location'].apply(lambda x: x.split("/")[0])
htseq['Project_ID']=f_project_id(htseq)
htseq=htseq[['Project_ID','Project_Type','Sample','file_name','location']].drop_duplicates().reset_index(drop=True)
htseq.to_csv('htseq.txt',index=False,header=True,sep='\t')   

#8. featurecount
featurecount=list_files[(list_files['location'].str.contains("featurecount")) & (list_files['location'].str.contains(".count"))].reset_index(drop=True)
featurecount=featurecount[-featurecount['location'].str.contains("summary")].reset_index(drop=True)

featurecount['file_name']=featurecount['location'].apply(lambda x: x.split("/")[len(x.split("/"))-1])
featurecount['Sample']=featurecount['file_name'].apply(lambda x: x.split(".")[0])

featurecount['Project_Type']=featurecount['location'].apply(lambda x: x.split("/")[0])
featurecount['Project_ID']=f_project_id(featurecount)
featurecount=htseq[['Project_ID','Project_Type','Sample','file_name','location']].drop_duplicates().reset_index(drop=True)
featurecount.to_csv('featurecount.txt',index=False,header=True,sep='\t')   



print "Opening mysql connection to update tables..." 

db = my.connect(host="138.25.210.11",
user="root",
passwd="bloodBioinf",
db="BloodGenomics_DB",
local_infile = 1
)

#print(db)
cursor = db.cursor()


list_tables=['raw_fastq','chipseq','star_bam','bwa_bam','deseq','circexp','htseq','featurecount']
for i in range(len(list_tables)):
    table_id=list_tables[i]
    #cursor.execute("""DROP TABLE `"""+table_id+"""`""")
    sql_cmd=f_create_mysql_table("""DROP TABLE IF EXISTS `"""+table_id+"""`;""")
    cursor.execute(sql_cmd)
    
    print "Updating table "+table_id+" ..." 
    sql_cmd=f_create_mysql_table(table_id)
    cursor.execute(sql_cmd)
    sql_cmd=f_load_data_mysql(table_id+".txt")
    cursor.execute(sql_cmd)
    print "Complete..." 

#
print "Saving changes..." 
db.commit()
print "Closing connection..." 
db.close()



 
