# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 09:54:00 2022

@author: 18811
"""


import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import pickle
from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri
from matplotlib_venn import venn2,venn3,venn3_circles


numpy2ri.activate()
pandas2ri.activate()
robjects.r('''source('E:/Code/python/WuQiong/limma.R')''') 

do_limma = robjects.globalenv['do_limma']


def foldchange(con_value,drug_value,log2=True):
    con_mean = con_value.mean(axis=1).values
    drug_mean = drug_value.mean(axis=1).values
    fc = drug_mean/con_mean
    if log2: fc = np.log2(fc)
    return fc


def t_test(con_value,drug_value):
    p_value = []
    for i in range(con_value.shape[0]):
        p = stats.ttest_ind(drug_value.iloc[i,:].values,con_value.iloc[i,:].values)[1]
        p_value.append(p)
    return p_value


def limma(con_value,drug_value):
    limma_input = pd.concat([con_value,drug_value],axis=1)
    group_list = np.array(['vehicle']*len(con_value.columns)+['drug']*len(drug_value.columns))
    pvalue = do_limma(limma_input, group_list).reshape(len(limma_input),)
    return pvalue
    

def cal_adj_pvalue(pvalue,method='fdr_bh',adjust=True,log10=True):
    if adjust: adj_pvalue = sm.stats.multipletests(pvalue, alpha = 0.05, method=method)[1]
    else: adj_pvalue = pvalue
    if log10: log_adj_pvalue = -np.log10(adj_pvalue)
    return adj_pvalue,log_adj_pvalue


def volcano(table,vechicle,drug,name,method,fc_thres = 1, pv_thres = 0.05):
    con_value,drug_value = table.loc[:,vechicle],table.loc[:,drug]
    table_ = pd.concat([table.iloc[:,:2],table.loc[:,vechicle+drug]],axis=1)
    fc = foldchange(con_value,drug_value)
    if method == 't_test': p = t_test(con_value,drug_value)
    if method == 'limma': p = limma(np.log2(con_value),np.log2(drug_value))
    adj_p,log_p = cal_adj_pvalue(p)
    
    table_['log2(FC)'] = fc
    table_['pvalue'] = p
    table_['adj_pvalue'] = adj_p
    table_['-log10(pvalue)'] = -np.log10(p)
    table_['-log10(adj.pvalue)'] = log_p
    
    sig = np.where((table_['adj_pvalue']<np.array(pv_thres))&(abs(table_['log2(FC)'])>np.array(fc_thres)))[0]
    sig_up = np.where((table_['adj_pvalue']<np.array(pv_thres))&((table_['log2(FC)'])>np.array(fc_thres)))[0]
    sig_down = np.where((table_['adj_pvalue']<np.array(pv_thres))&((table_['log2(FC)'])<-np.array(fc_thres)))[0]
    sig_gene = table_.iloc[sig,:].reset_index(drop=True)
    
    plt.figure(dpi = 900,figsize=(5,5)) 
    plt.scatter(fc, log_p, color = '#ABABA6',  marker = '.', s=150)
    plt.axvline(x = fc_thres,ls = '--', color = 'black')
    plt.axvline(x = -fc_thres,ls = '--', color = 'black')
    plt.axhline(y = -np.log10(pv_thres), ls = '--', color = 'black')
    plt.xlabel('log2(FC) ({:}/Vechicle)'.format(name),fontsize=14,labelpad=10)
    plt.ylabel('-log10 P-value',fontsize=14,labelpad=10)
    plt.scatter(fc[sig_up], log_p[sig_up], color = '#FA8260', edgecolors='#C85A40', marker = '.', s=150, linewidths=0.5) 
    plt.scatter(fc[sig_down], log_p[sig_down], color = '#4D8FD1', edgecolors='#315E94', marker = '.', s=150, linewidths=0.5) 
   
    """ 
    for k in sig:
        plt.text(fc[k]+0.1, log_p[k]+0.1, table_.loc[k,'PG.Genes'].split(';')[0], color = 'red', fontsize=7)
      """ 
    plt.title('{:}-{:}'.format(name,method),fontsize=16,pad=10)
    plt.xlim(-3,3)
    plt.savefig('E:/Code/R_input_files/{:}-{:}.pdf'.format(name,method),format='pdf', dpi=900, bbox_inches='tight')
    plt.show()
    return sig_gene,table_


table = pd.read_csv('E:/Code/R_input_files/test.csv')

name = 'Durg vs DMSO'

vechicle = list(table.iloc[:,2:5].columns) 
drug = list(table.iloc[:,5:8].columns) 

sig_prot,table_all = volcano(table,vechicle,drug,name,'t_test')
sig_prot = sig_prot.sort_values(by='-log10(pvalue)',ascending=False).reset_index(drop=True)
sig_prot.to_excel('E:/Code/R_input_files/t_test_{:}_sig_prot.xlsx'.format(name),index=False)


sig_prot,table_all = volcano(table,vechicle,drug,name,'limma')
sig_prot = sig_prot.sort_values(by='-log10(pvalue)',ascending=False).reset_index(drop=True)
sig_prot.to_excel('E:/Code/R_input_files/limma_{:}_sig_prot.xlsx'.format(name),index=False)

    
    


