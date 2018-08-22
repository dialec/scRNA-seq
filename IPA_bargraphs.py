# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 14:36:56 2018

@author: diego
"""
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
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
import math

#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt


input_genes="0_Input.txt"
canon_paths="1_CanonicalPathways.txt"
diseas="2_Diseases.txt"
mol_cellFunc="2_MolecularCellularFunctions.txt"
syst_devel="2_SystemDevelopment.txt"
ups_reg="3_UpstreamRegulators.txt"
gene_net="4_GeneNetworks.txt"
overl_bloodGenNet="5_OverlapBloodGenomicsNetworks.txt"

input_table_summary=sys.argv[1]
print input_table_summary
variable_x=sys.argv[2]
output_figure_svg=input_table_summary.split(".")[0]+".svg"
output_figure_png=input_table_summary.split(".")[0]+".png"


#input_table_summary="3_UpstreamRegulators.txt"
#output_figure=input_table_summary.split(".")[0]+".svg"
#variable_x="p-value of overlap"

DOWN_scale=['FF0000','FF1C1C','FF3838','FF5555','FF7171','FF8D8D','FFAAAA','FFC6C6','FFE2E2','FFFFFF']
UP_scale=['0000FF','1C1CFF','3838FF','5555FF','7171FF','8D8DFF','AAAAFF','C6C6FF','E2E2FF','FFFFFF']

#color_ratio=[-1,-0.95,-0.90,-0.85,-0.80,-0.75,-0.70,-0.65,-0.60,-0.65,0.5,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1]
#color_scale=['FF0000','FF1600','FF2C01','FF4201','FF5802','FF6E02','FF8403','FF9A03','FFB004','FFC604','FFDC05','E5D004','CCC404','B2B803','99AC03','7FA002','669402','4C8801','337C01','197000','006400']

color_ratio=[-1.00,-0.90,-0.80,-0.70,-0.60,0.5,0.60,0.70,0.80,0.90,1.00]
color_scale=['#FF0000','#FF2C01','#FF5802','#FF8403','#FFB004','#FFDC05','#CCC404','#99AC03','#669402','#337C01','#006400']
             
color_ratio=[-1.00,-0.90,-0.80,-0.70,-0.60,0.5,0.60,0.70,0.80,0.90,1.00]
color_scale=['#810E26','#9A3D44','#B36D62','#CC9D80','#E5CD9E','#FEFDBD','#CCE0A2','#9AC388','#68A66E','#368954','#056C3A']             
             
#color_ratio=[-1.00,-0.90,-0.80,-0.70,-0.60,0.5,0.60,0.70,0.80,0.90,1.00]
#color_scale=['#941229','#A84046','#BD6F63','#D29D80','#E7CC9D','#FCFBBA','#D2DFA4','#A9C48E','#7FA879','#568D63','#2D724E']

#color_ratio=[-1.00,-0.90,-0.80,-0.70,-0.60,0.5,0.60,0.70,0.80,0.90,1.00]
#color_scale=['#FF0000','#F02424','#E14848','#D26C6C','#C39090','#B5B5B5','#9090C3','#6C6CD2','#4848E1','#2424F0','#0000FF']


             
color_df=pd.DataFrame({'color_ratio':color_ratio,'RGB':color_scale})


print "Executing tables...."
suffix="_summary.txt" 
#input_tables=[canon_paths.split(".")[0], diseas.split(".")[0], mol_cellFunc.split(".")[0], syst_devel.split(".")[0],  ups_reg.split(".")[0], gene_net.split(".")[0], overl_bloodGenNet.split(".")[0]]
summary_tables=[canon_paths.split(".")[0]+suffix,  diseas.split(".")[0]+suffix, mol_cellFunc.split(".")[0]+suffix, syst_devel.split(".")[0]+suffix, ups_reg.split(".")[0]+suffix, gene_net.split(".")[0]+suffix, overl_bloodGenNet.split(".")[0]+suffix]




#print "Summary "+summary_tables[0]
if input_table_summary!="3_UpstreamRegulators.txt":
    table_summary=pd.read_csv(input_table_summary, sep="\t", header = 0 , low_memory=False)
else:
    table_summary=pd.read_csv(input_table_summary, sep="\t", header = 1 , low_memory=False)


print "Sorting by fold change and retrieveing top10..."
if "-log" in variable_x:
    table_summary=table_summary.sort_values([variable_x],ascending=[False]).reset_index(drop=True)
else:
    table_summary=table_summary.sort_values([variable_x],ascending=[True]).reset_index(drop=True)

table_summary=table_summary[:10][:]
table_summary['flag'] = np.where(table_summary['UpGenes']>=table_summary['DownGenes'], 'UP', 'DOWN')
table_summary['color_ratio'] = np.where(table_summary['flag']=="UP", table_summary['UpGenes']/table_summary['Total_Genes'], -table_summary['DownGenes']/table_summary['Total_Genes'])

table_summary['color_ratio']= table_summary['color_ratio'].apply(lambda x:np.round(x,decimals=1))
table_summary=pd.merge(table_summary,color_df,on=['color_ratio']).sort_values([table_summary.columns[1]],ascending=[False]).reset_index(drop=True)

print "Splitting name in two lines..."
table_summary['len_name']=table_summary[table_summary.columns[0]].apply(lambda x: len(x.split(" "))/2)
table_summary['name_split']=table_summary[table_summary.columns[0]].apply(lambda x: x.split(" "))
new_name_col=[]
for i in  range(len(table_summary)):
    new_name=table_summary['name_split'][i][0]
    for j in range(len(table_summary['name_split'][i])):
        print "row:"+ str(i)+" space:"+str(j)
        
        if j>0:
        
            if len(table_summary['name_split'][i])<7:
            
                if j==table_summary['len_name'][i]:
                    new_name=new_name+"\n"+table_summary['name_split'][i][j]
                else:
                    new_name=new_name+" "+table_summary['name_split'][i][j]
    
            
            else:
            
                if j==4 :
                    new_name=new_name+"\n"+table_summary['name_split'][i][j]
                else:
                    new_name=new_name+" "+table_summary['name_split'][i][j]
    
    new_name_col.append(new_name)    
            
        

table_summary['mult_line']=new_name_col
    


print "Plotting graph ..."
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt; plt.rcdefaults()
from pylab import *
print "working"
##
#objects = tuple(table_summary[table_summary.columns[0]])
objects = tuple(table_summary['mult_line'])
y_pos = np.arange(len(objects))
p_value = table_summary[variable_x]

#fig, axes = plt.subplots(figsize=(20,8))
fig, axes = plt.subplots(figsize=(20,6))
#fig, ax = plt.subplots()
#axes.set_ylabel( multialignment='center')
axes.set_xlabel(variable_x, fontsize=16)
#axes.set_xticklabels(fontsize=16)
axes.set_yticks(y_pos)
axes.set_yticklabels(objects, fontsize=14)

#table_summary['x_val']=table_summary['Total_Genes']*table_summary['color_ratio'].apply(lambda x: abs(x))
for y, x in enumerate(table_summary['Total_Genes'].apply(lambda x: int(x))):
    plt.annotate(str(x), xy=(0.5, y), va='center', size=15)

#axes.barh(y_pos, p_value,  align='center', color='green', ecolor='black')
axes.barh(y_pos, p_value,  align='center', color=table_summary['RGB'], edgecolor='black' )
axes.set_title(input_table_summary.split(".")[0], fontsize=18)
axes.tick_params(labelsize=12)
plt.gca().invert_yaxis()


ax1 = fig.add_axes([0.92, 0.11, 0.025, 0.77])  #[left, bottom, width, height] 

cmap = mpl.cm.RdYlGn
norm = mpl.colors.Normalize(vmin=-1, vmax=1) 

#bounds = [0.5, -0.75, 0.50, 0.75, 1]
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,  orientation='vertical', ticks=[0,0.25,0.5,0.75,1])
ax1.set_yticklabels(['100%', '75%','50%','75%','100%'])
#ax1.label_params(labelsize=14)
ax1.tick_params(labelsize=12)


#cb1.set_label('<-------------  Downregulated  ------------------->|<-----------------   Upregulated   ----------------->\n\nGene Regulation in Pathways')#, fontsize=12
cb1.set_label('<--------  Downregulated  -------->|<---------   Upregulated   --------->\n\nGene Regulation in Pathways')#, fontsize=12



#fig.colorbar(color_ratio)
#ax = fig.add_subplot()

#x = np.linspace(0, 2*np.pi, 400)
#y = np.sin(x**2)
#
#axes[1].scatter(x, y)
#ax.annotate('local max', xy=(2, 1), xytext=(3, 1.5),  arrowprops=dict(facecolor='black', shrink=0.05),)
#ax.set_ylim(-2,2)

#axes.set_xlim((min(p_value), max(p_value)))



fig.savefig(output_figure_svg, dpi=300, format='svg')
fig.savefig(output_figure_png, dpi=300, format='png')
print "working"
###
#np.random.seed(19680801)
#
#color_range=np.array([-1.00,-0.90,-0.80,-0.70,-0.60,0.5,0.60,0.70,0.80,0.90,1.00])
#plt.subplot(211)
#
##plt.imshow(np.random.random((100,100)), cmap=plt.cm.RdYlBu)
##plt.imshow(color_range, cmap=plt.cm.RdYlBu)
##plt.subplot(212)
##plt.imshow(np.random.random((100, 100)), cmap=plt.cm.RdYlBu)
##plt.imshow(color_range, cmap=plt.cm.RdYlBu)
#
#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
##cax = plt.axes([-1.00,-0.60,0.5,0.8])
#plt.colorbar(cax=cax, cmap=plt.cm.RdYlBu)
#    
#plt.savefig("test")
