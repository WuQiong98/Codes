# -*- coding: utf-8 -*-
"""
Created on Tue May  6 20:25:55 2025

@author: 18811
"""

import numpy as np
from os import listdir
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import math
from tqdm import tqdm
from sklearn.metrics import r2_score


def data_clean(table):
    table = table[table['# PSMs']>=1].reset_index(drop=True)
    for i in range(len(table)):
        if 'GN=' not in str(table.loc[i,'Description']):
            table.loc[i,'Description'] = 'nan'
        else:
            table.loc[i,'Description'] = table.loc[i,'Description'].split('GN=')[1].split(' ')[0]
    krt = ['KRT{:}'.format(i) for i in range(1,201)]
    table = table[~table['Description'].isin(krt)].reset_index(drop=True)
    for i in range(16):
        table = table[~table.iloc[:,i+3].isna()]
        table = table[~(table.iloc[:,i+3]==0)]
    table = table.reset_index(drop=True) 
    return table

def boxplot(table, title, log2=False, div=True):
    table_clean = table.iloc[:,-16:]
    
    if log2: 
        table_clean = np.log2(table_clean)
        ylabel = 'log2 Intensity'
    elif div:
        for i in range(0, 16, 8):
            table_clean.iloc[:,i:i+8] = table_clean.iloc[:,i:i+8].div(table_clean.iloc[:,i], axis=0)
        ylabel = 'Soluble Fraction'
        
    fig,ax = plt.subplots(figsize=(8,4),dpi=300)
    plt.rcParams['font.sans-serif'] = 'Arial'
    sns.boxplot(data=table_clean,
                width=0.5,
                showfliers=False,orient='v')
    plt.xticks(rotation=45,fontsize=12,ha='right')
    plt.yticks(fontsize=12)
    plt.ylabel(ylabel,fontsize=16,labelpad=5)
    plt.title(title,fontsize=16,pad=5)
    fig.autofmt_xdate()
   
    plt.savefig(f'E:/Code/R_input_files/boxplot_{title}.pdf', format='pdf', bbox_inches='tight', dpi=900)
    plt.show()

def pca(drug, table, title):
    data = table.iloc[:,2:].reset_index(drop=True).T
    pca = PCA(n_components=0.99, whiten=True)
    reduced_data = pd.DataFrame(pca.fit_transform(data)).iloc[:,:2]
    reduced_data.columns = ['PCA_1', 'PCA_2']
    rep = data.index.str.split('-', expand=True)
    temps = data.index.str.split('_', expand=True)
    reduced_data['replicate'] = [r[1] for r in rep]
    reduced_data['tempreture'] = [r[1] for r in temps]
    
    x,y = pca.explained_variance_ratio_[:2]
    
    plt.figure(dpi=300,figsize=(4,4))
    plt.rcParams['font.sans-serif'] = 'Arial'
    sns.scatterplot(data=reduced_data,x='PCA_1', y='PCA_2', 
                    hue='tempreture', style='replicate')
    plt.legend(bbox_to_anchor=[1,1],fontsize=10)
    plt.xlabel('PCA_1: {:.2f}%'.format(x*100),fontsize=14,labelpad=5)
    plt.ylabel('PCA_2: {:.2f}%'.format(y*100),fontsize=14,labelpad=5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(drug,fontsize=16,pad=5)
    plt.ylim(-2,4)
    
    plt.savefig(f'E:/Code/R_input_files/PCA_{title}.pdf', format='pdf', bbox_inches='tight', dpi=900)
    plt.show()

def curve_template0(t, a, b, p ):
    return (1-p)/ (1+ np.exp(-a/(t+273.15) + b)) + p

def cal_tm(y,popt):
    t = (popt[0]/(popt[1]-math.log(((1-popt[2])/(y-popt[2]))-1)))-273.15
    return round(t,2)


col_list = [3, 4, 10]+list(range(18, 34))
group = ['DMSO', 'Drug']
temp = [37, 44, 47, 50, 53, 55, 59, 66]
temp_arr = np.array(temp)
tem = np.linspace(37,66,1000)

path = r"E:/Code/R_input_files/TPP-WQ"
summary = {}
for name in listdir(path):
    if 'xlsx' not in name: continue
    table = pd.read_excel(os.path.join(path, f'{name}'), usecols=col_list)
    table = data_clean(table)
    table.columns = ['Accession','GeneSymbol','PSMs']+[f'{g}_{t}' for g in group for t in temp]
    summary.setdefault(name.split('.')[0], table)


total = {}
for name, table in summary.items():
   
    boxplot(table, f"prenormalized_{name.split('.')[0]}")
    table_ = table.copy()
    table_.drop(columns=['PSMs'],inplace=True)
    for i in range(2, 18, 8):
        table_.iloc[:,i:i+8] = table_.iloc[:,i:i+8].div(table_.iloc[:,i], axis=0)
    table_.columns = ['Accession','GeneSymbol']+[f'{g}_{t}_{name}' for g in group for t in temp]
    if name.split('-')[0] not in total.keys():
        total[name.split('-')[0]] = table_.copy()
    else:
        t = total[name.split('-')[0]]
        t = pd.merge(t, table_, on=['Accession','GeneSymbol'])
        total.update({name.split('-')[0]:t})

for drug, table in total.items():
    
    pca(drug, table, title=f'prenormalized_{drug}')
    

drugs = list(total.keys())
norm_tables = {}
for drug in drugs:
    table = pd.DataFrame()
    for i in range(1,4):
        t1 = summary[f'{drug}-{i}'].copy()
        t1.columns = ['Accession','GeneSymbol',f'{drug}-{i}_PSMs']+[f'{drug}-{i}_{g}_{t}' for g in group for t in temp]
        if table.empty: table = t1.copy()
        else: table = pd.merge(table, t1, on=['Accession','GeneSymbol'], how='outer')
            
    norm_table = table[['Accession','GeneSymbol']]
    psm = table[[f'{drug}-1_PSMs',f'{drug}-2_PSMs',f'{drug}-3_PSMs']].sum(axis=1)
    table = table.replace(np.nan,0)
    for g in group:
        for t in temp:
            norm_table[f'{g}_{t}'] = (table[f'{drug}-1_PSMs']*table[f'{drug}-1_{g}_{t}']+
                                      table[f'{drug}-2_PSMs']*table[f'{drug}-2_{g}_{t}']+
                                      table[f'{drug}-3_PSMs']*table[f'{drug}-3_{g}_{t}'])/psm
    norm_tables.setdefault(drug, norm_table)

final_table = {}
for name, table in norm_tables.items():
    table_ = table.copy()
    for i in range(2, 18, 8):
        table_.iloc[:,i:i+8] = table_.iloc[:,i:i+8].div(table_.iloc[:,i], axis=0)
    
    for i in range(2, 10):
        m = table_.iloc[:,[i, i+8]].median(axis=0)
        value = m.median()/m
        table_.iloc[:,i] = table_.iloc[:,i]*value[0]
        table_.iloc[:,i+8] = table_.iloc[:,i+8]*value[1]
    
    median = table_.median().values[:8]
    popt0, pcov0 = curve_fit(curve_template0, temp, median, maxfev=10000)
    yvals = curve_template0(temp_arr,*popt0)
    m = list(yvals/median)
    for i in range(2, 18, 8):
        table_.iloc[:,i:i+8] = table_.iloc[:,i:i+8]*m
    final_table.setdefault(name, table_)
    
  
    title = f'normalized_{name}'
    
    fig,ax = plt.subplots(figsize=(8,4),dpi=300)
    plt.rcParams['font.sans-serif'] = 'Arial'
    sns.boxplot(data=table_,
                width=0.5,
                showfliers=False,orient='v')
    plt.xticks(rotation=45,fontsize=12,ha='right')
    plt.yticks(fontsize=12)
    plt.ylabel('Soluble Fraction',fontsize=16,labelpad=5)
    plt.title(name,fontsize=16,pad=5)
    fig.autofmt_xdate()
    plt.savefig(f'E:/Code/R_input_files/boxplot_{title}.pdf', format='pdf', bbox_inches='tight', dpi=900)
    plt.show()
    
    
    pca(name, table_, title=f'normalized_{name}')


for drug, table in final_table.items():
    table_ = table.copy()
    feat = []
    for i in tqdm(table_.index):
        values = table_.iloc[i,2:].values
        a, b, p, r, tm = [], [], [], [], []
        for v in range(0, 16, 8):
            try:
                popt0, pcov0 = curve_fit(curve_template0, temp, values[v:v+8], maxfev=10000)
                tm.append(cal_tm(0.5, popt0))
                y_pre = curve_template0(temp_arr, *popt0)
                a.append(popt0[0])
                b.append(popt0[1])
                p.append(popt0[2])
                r.append(round(r2_score(values[v:v+8], y_pre),4))
            except:
                tm.append(np.nan)
                a.append(np.nan)
                b.append(np.nan)
                p.append(np.nan)
                r.append(np.nan)
        feat.append(tuple((a[0], b[0], p[0], r[0], a[1], b[1], p[1], r[1], round(tm[1]-tm[0],2))))
    feat = pd.DataFrame(feat, columns=[f'{a}_{b}' for a in ['DMSO','Drug'] for b in ['a','b','p','R2']]+['ΔTm'])
    table_ = pd.concat([table_, feat], axis=1)
    table_.to_csv(os.path.join(path, f'files/{drug}.csv'), index=False, encoding="gbk")


drug2tar = {'MTX':['DHFR'], 'Pan':['HDAC1','HDAC2'], 'Ral':['TYMS']}
colors = ['steelblue', 'darkorange']
for drug, tar in drug2tar.items():
    table = final_table[drug]
    for t in tar:
        values = table[table['GeneSymbol']==t].iloc[0,2:].values
        plt.figure(figsize=(4,4),dpi=300)
        plt.rcParams['font.sans-serif'] = 'Arial'
        tm = []
        for i in range(0, 16, 8):
            popt0, pcov0 = curve_fit(curve_template0, temp, values[i:i+8], maxfev=10000)
            yvals = curve_template0(tem, *popt0)
            plt.plot(temp, values[i:i+8], 'o',markersize=5,color=colors[i//8])
            plt.plot(tem, yvals,color=colors[i//8])
            tm.append(cal_tm(0.5, popt0))
        plt.xticks(temp, temp, fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('Temperature', fontsize=16)
        plt.ylabel('Soluble Fraction', fontsize=16)
        plt.title(f'{drug}-{t}', fontsize=16)
        plt.text(54, 0.95, f'Tm={round(tm[1]-tm[0],2)} °C', fontsize=14)
        
        plt.savefig(f'E:/Code/R_input_files/{drug}-{t}.pdf', format='pdf', dpi=900, bbox_inches='tight')
        plt.show()
        




