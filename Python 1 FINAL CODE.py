# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:08:26 2024
final revision
@author: ryanj
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Load the CSV file
file_path = r"C:\Users\ryanj\Project_1\source\38,102,104, 84.csv"
data = pd.read_csv(file_path)

# Calculate FFF/FFI
data['FFF/FFI'] = abs(1 - (data['fill_factor_initial'] / data['fill_factor_final']))

# Function to calculate Cohen's d for effect size
def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)
    return (np.mean(x) - np.mean(y)) / pooled_std

# Define a function to perform the analysis for a given metric
def analyze_metric(data, metric, reaction1, reaction2):
    # Extract data for the specified reactions from the 'reaction' column
    group1 = data[data['reaction'] == reaction1][metric]
    group2 = data[data['reaction'] == reaction2][metric]

    # Convert data to numeric, handling non-numeric values
    group1 = pd.to_numeric(group1, errors='coerce').dropna()
    group2 = pd.to_numeric(group2, errors='coerce').dropna()

    # Summary statistics
    group1_summary = group1.describe()
    group2_summary = group2.describe()

    # Boxplots
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    sns.boxplot(y=group1)
    plt.title(f'{reaction1} {metric} Boxplot')
    plt.subplot(1, 2, 2)
    sns.boxplot(y=group2)
    plt.title(f'{reaction2} {metric} Boxplot')
    plt.show()

    # Log transformation
    group1_log = np.log1p(group1)
    group2_log = np.log1p(group2)
    group1_log_skewness = group1_log.skew()
    group2_log_skewness = group2_log.skew()

    # T-test and Mann-Whitney U test
    t_stat_log, p_value_log = stats.ttest_ind(group1_log, group2_log)
    u_stat, p_value = stats.mannwhitneyu(group1, group2)
    effect_size = cohen_d(group1, group2)

    # Violin plots
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    sns.violinplot(y=group1)
    plt.title(f'{reaction1} {metric} Violin Plot')
    plt.subplot(1, 2, 2)
    sns.violinplot(y=group2)
    plt.title(f'{reaction2} {metric} Violin Plot')
    plt.show()

    return {
        'group1_summary': group1_summary,
        'group2_summary': group2_summary,
        'group1_log_skewness': group1_log_skewness,
        'group2_log_skewness': group2_log_skewness,
        't_stat_log': t_stat_log,
        'p_value_log': p_value_log,
        'u_stat': u_stat,
        'p_value': p_value,
        'effect_size': effect_size
    }

# Perform the analysis for each metric for NK-R102 vs 38 and NK-R104 vs 38
metrics = ['doublets', 'FFF/FFI', 'total_droplets']
comparisons = [
    ('NK-R102', 'NK-R38'),
    ('NK-R104', 'NK-R38'),
    ('NK-R102', 'NK-R104'),
    ('NK-R102', 'NK-R84'),
    ('NK-R104', 'NK-R84')
]

results = {f"{reaction1}_vs_{reaction2}": {metric: analyze_metric(data, metric, reaction1, reaction2) for metric in metrics} for reaction1, reaction2 in comparisons}

# Compile results into a DataFrame with individual columns for summary statistics
def compile_results(results):
    compiled_data = []
    for comparison, metrics in results.items():
        for metric, result in metrics.items():
            row = {
                'Comparison': comparison,
                'Metric': metric,
                'Log Skewness Group1': result['group1_log_skewness'],
                'Log Skewness Group2': result['group2_log_skewness'],
                'T-test p-value (log)': result['p_value_log'],
                'Mann-Whitney U test p-value': result['p_value'],
                'Effect size (Cohen\'s d)': result['effect_size']
            }
            row.update({f'Group1 {stat}': value for stat, value in result['group1_summary'].items()})
            row.update({f'Group2 {stat}': value for stat, value in result['group2_summary'].items()})
            compiled_data.append(row)
    return pd.DataFrame(compiled_data)

compiled_results = compile_results(results)

# Save the results to an Excel file
output_file = r"C:\Users\ryanj\Project_1\output\final_compiled_results.xlsx"
with pd.ExcelWriter(output_file) as writer:
    compiled_results.to_excel(writer, index=False, sheet_name='Results')

# Display the DataFrame
print(compiled_results)