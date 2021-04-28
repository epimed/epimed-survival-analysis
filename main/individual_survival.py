from service import survivalService, figureService

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

arialNarrowFont = figureService.createArialNarrowFont()
regular, medium, small, tiny = figureService.createFontSizes(regular=22)
todayFormat = figureService.createTodayFormat()

max_survival_time = 120
group_colors = {'Low': 'royalblue', 'High': 'crimson'}

# Data file

'''
dataset = 'TCGA-COAD'
threshold_col = 'threshold' # Attention! The column name is different depending on dataset. For GSE39582, it is test_threshold.
'''

dataset = 'GSE39582'
threshold_col = 'test_threshold'


thresholds = pd.read_csv('../data/thresholds.csv', sep=';')
thresholds.index = thresholds['gene_symbol'] + '@' + thresholds['gene'].apply('{:.0f}'.format)
print(thresholds.head())

data = pd.read_csv('../data/expression_data_' + dataset + '.csv', sep=';', index_col='id_sample')
print(data.head())

features = list(data.columns)
features.remove('time')
features.remove('event')

# Ajust time and event to the follow-up stopped at max_survival_time
for id_sample in data.index:
    shifted_time, shifted_event = survivalService.shiftToMaxSurvivalTime(data.loc[id_sample, 'time'], data.loc[id_sample, 'event'], max_survival_time)
    data.loc[id_sample, 'time'] = shifted_time
    data.loc[id_sample, 'event'] = shifted_event

print(data.head())

# features = ['ERFE@151176']

cph = CoxPHFitter()
kmf = KaplanMeierFitter()

for feature in features:
    threshold = thresholds.loc[feature, threshold_col]
    print(dataset, feature, 'threshold', threshold)
    
    # Prepare DataFrame for Cox model: covariate = expression of the gene
    cox_expression = data[[feature, 'time', 'event']]
    
    # Cox model for expression
    cph.fit(cox_expression, duration_col='time', event_col='event', show_progress=False)
    cox_pvalue_expression = cph.summary.p[feature]
    cox_hr_expression = cph.summary['exp(coef)'][feature]
    
    
    # Prepare DataFrame for Cox model: covariate = binary group: "low" - under the threshold, "high" - over the threshold 
    cox_group =  cox_expression.copy()
    cox_group['group'] = 0.0
    cox_group.loc[cox_group[feature]>threshold, 'group'] = 1.0
    cox_group = cox_group[['group', 'time', 'event']]
    
    # Cox model for groups
    cph.fit(cox_group, duration_col='time', event_col='event', show_progress=False)
    cox_pvalue_group = cph.summary.p['group']
    cox_hr_group = cph.summary['exp(coef)']['group']
    
    # Logrank test for groups
    logrank = multivariate_logrank_test(cox_group['time'], cox_group['group'], cox_group['event'])
    
    
    # FIGURE
        
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    
    groups = { 'Low': cox_group[cox_group['group']==0.0], 'High': cox_group[cox_group['group']==1.0]}
    
    for group_name in groups.keys():
        group = groups[group_name]
        labelText = group_name + ' (n=' + str(group.shape[0]) +')'
        kmf.fit(group['time'], group['event'], label=labelText)
        kmf.plot(ax=ax, ci_show=False, show_censors=True, color=group_colors[group_name], linewidth=3)
        # kmf.plot(ax=ax)
    
    gene_symbol = feature.split('@')[0]
    cox_text = 'cox p-value = ' +  '{:.1e}'.format(cox_pvalue_expression) + ' ' + figureService.getSignificanceSymbol(cox_pvalue_expression)
    logrank_text = 'logrank p-value = ' + '{:.1e}'.format(logrank.p_value) + ' ' + figureService.getSignificanceSymbol(logrank.p_value)
    hr_text = 'HR between groups = ' + '{:.1f}'.format(cox_hr_group)
    
    fig.text(0.5, 1.1, gene_symbol + ' - ' + dataset, ha='center', va='center', fontweight='bold', color='grey', fontsize=regular, **arialNarrowFont)
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
    
    file_prefix = 'kaplan_meier_' + feature + '_' + dataset
    figureService.savefigWithResolution(fig, '../figures/individual_survival/', file_prefix, dpi=100, ext='png')  
    