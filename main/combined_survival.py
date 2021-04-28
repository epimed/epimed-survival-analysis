from service import survivalService, figureService

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

def getGroupName(groupNumberGenes, n):
    for groupName in groupNumberGenes:
        nmin = groupNumberGenes[groupName][0]
        nmax = groupNumberGenes[groupName][1]
        if (n>=nmin and n<=nmax):
            return groupName

def getGroupColors(groupNumberGenes):
    nbGroups = len(groupNumberGenes)
    groupColors = {'P1' : 'royalblue', 'P2' : 'black', 'Not expressed' : 'royalblue', 'Expressed' : 'black'}
    if (nbGroups==3):
        groupColors = {'P1' : 'royalblue', 'P2' : 'crimson', 'P3' : 'black'}
    if (nbGroups==4):
        groupColors = {'P1' : 'royalblue', 'P2' : 'orange', 'P3' : 'crimson', 'P4' : 'black'}
    if (nbGroups==5):
        groupColors = {'P1' : 'royalblue', 'P2' : 'orange', 'P3' : 'crimson', 'P4' : 'darkviolet', 'P5' : 'black'}
    if (nbGroups==6):
        groupColors = {'P1' : 'cyan', 'P2' : 'royalblue', 'P3' : 'orange', 'P4' : 'crimson', 'P5' : 'darkviolet', 'P6' : 'black'}
    if (nbGroups==7):
        groupColors = {'P1' : 'cyan', 'P2' : 'royalblue', 'P3' : 'orange', 'P4' : 'red', 'P5' :'crimson', 'P6' : 'darkviolet', 'P7' : 'black'}
    return  groupColors


def generateGroupLabel(groupNumberGenes, group):
    nmin = groupNumberGenes[group][0]
    nmax = groupNumberGenes[group][1]
    if (nmin<nmax):
        return str(nmin) + '-' + str(nmax)
    else:
        return str(nmin)


arialNarrowFont = figureService.createArialNarrowFont()
regular, medium, small, tiny = figureService.createFontSizes(regular=22)
todayFormat = figureService.createTodayFormat()

max_survival_time = 120
group_colors = {'Low': 'royalblue', 'High': 'crimson'}

# Data file

dataset = 'TCGA-COAD'
threshold_col = 'threshold' # Attention! The column name is different depending on dataset. For GSE39582, it is test_threshold.

'''
dataset = 'GSE39582'
threshold_col = 'test_threshold'
'''

combined_features = ['ERFE@151176', 'CCDC154@645811', 'LINC01356@100996702', 'ODF3L2@284451', 'SLC6A1@6529']

# groupNumberGenes = {'P1': [0, 0], 'P2': [1, len(combined_features)]}  
groupNumberGenes = {'P1': [0, 0], 'P2': [1, 1], 'P3': [2, len(combined_features)]} 
# groupNumberGenes = {'P1': [0, 0], 'P2' : [1, 1], 'P3' : [2, 2], 'P4': [3, len(combined_features)]}


thresholds = pd.read_csv('../data/thresholds.csv', sep=';')
thresholds.index = thresholds['gene_symbol'] + '@' + thresholds['gene'].apply('{:.0f}'.format)
thresholds = thresholds.loc[combined_features, threshold_col]
print(thresholds.head())

data = pd.read_csv('../data/expression_data_' + dataset + '.csv', sep=';', index_col='id_sample')

features = list(data.columns)
features.remove('time')
features.remove('event')


for id_sample in data.index:
    
    # Ajust time and event to the follow-up stopped at max_survival_time
    shifted_time, shifted_event = survivalService.shiftToMaxSurvivalTime(data.loc[id_sample, 'time'], data.loc[id_sample, 'event'], max_survival_time)
    data.loc[id_sample, 'time'] = shifted_time
    data.loc[id_sample, 'event'] = shifted_event
    
    # Calculate number of activated gene above threshold and the corresponding group
    n = 0
    for feature in combined_features:
        if (data.loc[id_sample, feature]) > thresholds.loc[feature]:
            n = n + 1
    data.loc[id_sample, 'n'] = n
    data.loc[id_sample, 'group'] = getGroupName(groupNumberGenes, n)
    
print(data.head())

# Logrank test
logrank = multivariate_logrank_test(data['time'], data['group'], data['event'])  

# Cox model    
cph = CoxPHFitter()
cox = data[['n', 'time', 'event']]
cph.fit(cox, duration_col='time', event_col='event', show_progress=False)
print(cph.summary)
cox_pvalue = cph.summary.p['n']
cox_hr = cph.summary['exp(coef)']['n']


# FIGURE

nMax = int(data['n'].max()) 
sortedGroups = sorted(groupNumberGenes.keys())
lastGroup = sortedGroups[len(sortedGroups)-1]
groupNumberGenes[lastGroup][1] = nMax

gecText = str(len(combined_features)) + '-GEC'

if (len(combined_features)<10):
    gene_symbols = [f.split('@')[0] for f in combined_features]
    gecText = gecText + ': ' + ', '.join(sorted(gene_symbols))

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
fig.text(0.5, 1.18, dataset + ' (n=' + str(data.shape[0]) +')', ha='center', va='center', color='grey',  fontsize=regular, fontweight='normal', **arialNarrowFont)
fig.text(0.5, 1.1, gecText, ha='center', va='center', color='grey', fontsize=medium, fontweight='bold', **arialNarrowFont)

groupColors = getGroupColors(groupNumberGenes)

kmf = KaplanMeierFitter()
for group in sorted(groupNumberGenes.keys()):
    selection = data[data['group']==group]
    nSamples = selection.shape[0]
    labelText = generateGroupLabel(groupNumberGenes, group) + ' (n=' + str(nSamples) +')'
    kmf.fit(selection['time'], selection['event'], label=labelText)
    kmf.plot(ax=ax, ci_show=False, show_censors=True, color=groupColors[group], linewidth=3)
    print(group, 'Samples', nSamples, list(selection.index))

cox_text = 'cox p-value = ' +  '{:.1e}'.format(cox_pvalue) + ' ' + figureService.getSignificanceSymbol(cox_pvalue)
logrank_text = 'logrank p-value = ' + '{:.1e}'.format(logrank.p_value) + ' ' + figureService.getSignificanceSymbol(logrank.p_value)
hr_text = 'hazard ratio = ' + '{:.1f}'.format(cox_hr)

ax.set_title(cox_text + '\n' + logrank_text + '\n' + hr_text, fontsize=medium, **arialNarrowFont) 
    
# L = ax.legend(fontsize=small, loc='center', bbox_to_anchor=(0.5, -0.25), ncol=1)
L = ax.legend(fontsize=small, loc='lower left')
plt.setp(L.texts, **arialNarrowFont)

ax.set_ylim([-0.05, 1.05])
ax.set_xlim([0.0 - 0.05*max_survival_time, max_survival_time + 0.05*max_survival_time])

step = 10.0 * np.ceil(max_survival_time/100.0)
xticks = np.arange(0, max_survival_time + step, step)
ax.set_xticks(xticks)
ax.set_xticklabels(xticks, **arialNarrowFont)
ax.set_xlabel('Time in months', fontsize=regular, **arialNarrowFont)

yticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ax.set_yticks(yticks)
ax.set_yticklabels(yticks, **arialNarrowFont)
ax.set_ylabel('Overall survival', fontsize=regular, **arialNarrowFont)
ax.tick_params(axis='both', labelsize=medium)

plt.close(fig)

file_prefix = 'kaplan_meier_' +  dataset + '_' + str(len(combined_features)) + '_features_' + str(len(groupNumberGenes.keys())) + '_groups'
figureService.savefigWithResolution(fig, '../figures/combined_survival/', file_prefix, dpi=100, ext='png')  
    
