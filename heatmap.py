# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 14:36:56 2018
Before: Heatmap generated from normalised data of 199 cells
After: Heatmap generated from normalised data of 207 samples. HThen, 8 samples wer eremoved manually.


@author: diego
"""
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
from collections import Counter





#out_dir=sys.argv[3]


print "for Seurat2"

in_folder="./"
samples2exlude=['MCD4899_A37','MCD4899_A50','MCD4899_A73','MCD4899_A72','MCD4899_A14','MCD4899_A206','MCD4899_A26','MCD4899_A54']
cluster_info="./exported_clusters_without_lowseq_outliers_199.csv" # seurat 2
out_dir="./"

#
#print "-option 1"
#preffix="Seurat2_"
#normalised_count="../../C_vs_VT_207/DiffExpr/DeSeq/deSeq_normalized.txt" # seurat 2
#
#print "-option 2"
#preffix="Seurat2_CPM_"
#normalised_count="./207/logcountsCPM_normalisation.txt" # seurat 2

print "-option 3"
preffix="Seurat2_Scran_"
normalised_count="./207/Scran_normalisation_size30.txt" # seurat 2

print "Renaming columns from normalised data..."

clusters_ali= pd.read_csv(cluster_info, sep=",", header = 0 , low_memory=False)    
normal_genes= pd.read_csv(normalised_count, sep="\t", header = 0 , low_memory=False)    

print "Filtering out 8 samples..."
normal_genes_transp=normal_genes.transpose()
normal_genes_transp=normal_genes_transp[normal_genes_transp.index.isin(clusters_ali['Unnamed: 0'])]

print "Updating matrix of normalised data..."
normal_genes=normal_genes_transp.transpose()

print "Renaming columns (samples) of matrix..."
clusters_ali['Unnamed: 0']=clusters_ali['Unnamed: 0'].apply(lambda x: x.split(".")[3])
clusters_ali['x']=map(str,list(clusters_ali['x']))
clusters_ali['x']="c_"+clusters_ali['x']

normal_genes.columns=list(clusters_ali['Unnamed: 0'])
   

    
seurat_genes=['0_1','0_2','0_3','1_2','1_3','2_3']

for i in range(len(seurat_genes)):
    print str(i)
    seurat_cluster= pd.read_csv(in_folder+"seurat_clusters_"+seurat_genes[i]+".csv", sep=",", header = 0 , low_memory=False)
    #subprocess.call("cp ./Seurat_Diff/seurat_clusters_"+seurat_genes[i]+".csv ./new_heatmap/", shell=True)     
    
    seurat_cluster=seurat_cluster[['Unnamed: 0','p_val','avg_logFC','p_val_adj']]
    seurat_cluster=seurat_cluster[seurat_cluster['p_val']<0.05]
    seurat_cluster.rename(columns={'p_val': 'p_val_'+seurat_genes[i], 'avg_logFC': 'avg_logF_'+seurat_genes[i], 'p_val_adj': 'p_val_adj_'+seurat_genes[i]}, inplace=True)
    
    if i==0:
        merge_table=seurat_cluster.copy()
    else:
        merge_table=pd.merge(merge_table,seurat_cluster,how="outer",on=['Unnamed: 0'])


    

#freq_genes=pd.DataFrame()
print "Defining the genes differentially expressed in at least one of the clusters..."

#for FC_thr in [0,0.5,1,1.5,2,2.5,3]:
for FC_thr in [1.5]:
    print "Generating heatmap for fold-change > "+str(FC_thr)    
    merge_table_eval=merge_table[(merge_table['avg_logF_0_1'].abs()>FC_thr)|(merge_table['avg_logF_0_2'].abs()>FC_thr)|(merge_table['avg_logF_0_3'].abs()>FC_thr)|(merge_table['avg_logF_1_2'].abs()>FC_thr)|(merge_table['avg_logF_1_3'].abs()>FC_thr)|(merge_table['avg_logF_2_3'].abs()>FC_thr)].reset_index(drop=True)
    
    #merge_table_eval.fillna(-100, inplace=True)
    diffExp_genes=merge_table_eval[['Unnamed: 0','avg_logF_0_1','p_val_0_1','avg_logF_0_2','p_val_0_2','avg_logF_0_3','p_val_0_3','avg_logF_1_2','p_val_1_2','avg_logF_1_3','p_val_1_3','avg_logF_2_3','p_val_2_3']]
    diffExp_genes.to_csv(out_dir+preffix+"Lisa_DiffExp_Genes_clustering_"+str(FC_thr),index=False,header=True,sep='\t')              
    
    print "Getting normalised count matrix for differential expressed genes..."
    expression_heatmap=normal_genes[normal_genes.index.isin(diffExp_genes['Unnamed: 0'])] 
    print "Renaming columns using cluster assignation..."
    expression_heatmap.columns=list(clusters_ali['x'])
    print "Generating output..."
    #normal_genes=DataFrame.transpose(normal_genes)
    expression_heatmap.to_csv(out_dir+preffix+"Lisa_heatmap_"+str(FC_thr),index=True,header=True,sep='\t')              
    
#    
#
#print "2 PART: FILTERING NORMALISED DATA FROM GENES OF INTEREST..."
#
#clusters_ali= pd.read_csv(cluster_info, sep=",", header = 0 , low_memory=False)    
#clusters_ali['Unnamed: 0']=clusters_ali['Unnamed: 0'].apply(lambda x: x.split(".")[3])
#clusters_ali['x']=map(str,list(clusters_ali['x']))
#clusters_ali['x']="c_"+clusters_ali['x']
#
#normal_genes= pd.read_csv(normalised_count, sep="\t", header = 0 , low_memory=False)    
#normal_genes.columns=list(clusters_ali['Unnamed: 0'])
#genes2filter= pd.read_csv("GENES_INTEREST.txt", sep="\t", header = 0 , low_memory=False)   
#
#
#expression_heatmap=normal_genes[(normal_genes.index.isin(genes2filter['LIST1'])) | (normal_genes.index.isin(genes2filter['LIST2']))]
#expression_heatmap.columns=list(clusters_ali['x'])
#expression_heatmap.to_csv(out_dir+"Lisa_Expression_heatmap_FULL_GENES_INTEREST",index=True,header=True,sep='\t')              
#
#
#print "3 PART: WHICH GENES FROM THE LIST OF INTERESET ARE DIFF EXPRESSED..."
#total_diffexp_genes= pd.read_csv(out_dir+"Lisa_DiffExp_Genes_from_cluster_0", sep="\t", header = 0 , low_memory=False)    
#
#total_diffexp_genes=total_diffexp_genes[(total_diffexp_genes['Unnamed: 0'].isin(genes2filter['LIST1'])) |  (total_diffexp_genes['Unnamed: 0'].isin(genes2filter['LIST2']))]
#total_diffexp_genes.to_csv(out_dir+"Lisa_DiffExp_GENES_INTEREST_0",index=False,header=True,sep='\t')              
#
#
#
#expression_heatmap=normal_genes[normal_genes.index.isin(total_diffexp_genes['Unnamed: 0'])] 
#expression_heatmap.columns=list(clusters_ali['x'])
#expression_heatmap.to_csv(out_dir+"Lisa_heatmap_only_diffExp_GENES_INTEREST",index=True,header=True,sep='\t')              





