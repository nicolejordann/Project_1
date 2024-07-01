# -*- coding: utf-8 -*-
"""
Created on Jun 27 15:10:11 2024

@author: ryanj
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import numpy as np

# Load the CSV file
file_path = r"C:\Users\ryanj\Project_1\source\38,102,104, 84.csv"

data = pd.read_csv(file_path)

# Calculate FFF/FFI
data['FFF/FFI'] = abs(1 - (data['fill_factor_initial'] / data['fill_factor_final']))

# Define a function to perform the analysis for a given metric
def analyze_metric(data, metric, reaction1, reaction2):
    # Extract data for the specified reactions from the 'reaction' column
    group1 = data[data['reaction'] == reaction1][metric]
    group2 = data[data['reaction'] == reaction2][metric]

    # Convert data to numeric, handling non-numeric values
    group1 = pd.to_numeric(group1, errors='coerce')
    group2 = pd.to_numeric(group2, errors='coerce')

    # Drop missing values
    group1 = group1.dropna()
    group2 = group2.dropna()

    # Step 1: EDA
    # Summary statistics
    group1_summary = group1.describe()
    group2_summary = group2.describe()

    # Plot boxplots
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    sns.boxplot(y=group1)
    plt.title(f'{reaction1} {metric} Boxplot')

    plt.subplot(1, 2, 2)
    sns.boxplot(y=group2)
    plt.title(f'{reaction2} {metric} Boxplot')

    plt.show()

    # Step 2: Data Transformation
    # Apply log transformation to handle skewness
    group1_log = np.log1p(group1)
    group2_log = np.log1p(group2)

    # Check skewness after transformation
    group1_log_skewness = group1_log.skew()
    group2_log_skewness = group2_log.skew()

    # Step 3: Hypothesis Testing
    # Conduct t-test after log transformation
    t_stat_log, p_value_log = stats.ttest_ind(group1_log, group2_log)

    # Conduct Mann-Whitney U test on original data
    u_stat, p_value = stats.mannwhitneyu(group1, group2)

    # Step 4: Effect Size Calculation
    # Calculate Cohen's d for effect size
    def cohen_d(x, y):
        nx = len(x)
        ny = len(y)
        dof = nx + ny - 2
        pooled_std = np.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)
        return (np.mean(x) - np.mean(y)) / pooled_std

    effect_size = cohen_d(group1, group2)

    # Step 5: Advanced Visualization
    # Plot violin plots
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    sns.violinplot(y=group1)
    plt.title(f'{reaction1} {metric} Violin Plot')

    plt.subplot(1, 2, 2)
    sns.violinplot(y=group2)
    plt.title(f'{reaction2} {metric} Violin Plot')

    plt.show()

    return group1_summary, group2_summary, group1_log_skewness, group2_log_skewness, t_stat_log, p_value_log, u_stat, p_value, effect_size

# Perform the analysis for each metric for NK-R102 vs 38 and NK-R104 vs 38
metrics = ['doublets', 'FFF/FFI', 'total_droplets']
results_102_vs_38 = {}
results_104_vs_38 = {}
results_102_vs_104 = {}
results_102_vs_84 = {}
results_104_vs_84 = {}

for metric in metrics:
    results_102_vs_38[metric] = analyze_metric(data, metric, 'NK-R102', 'NK-R38')
    results_104_vs_38[metric] = analyze_metric(data, metric, 'NK-R104', 'NK-R38')
    results_102_vs_104[metric] = analyze_metric(data, metric, 'NK-R102', 'NK-R104')
    results_102_vs_84[metric] = analyze_metric(data, metric, 'NK-R102', 'NK-R84')
    results_104_vs_84[metric] = analyze_metric(data, metric, 'NK-R104', 'NK-R84')

results_102_vs_38, results_104_vs_38, results_102_vs_104, results_102_vs_84, results_104_vs_84

# Compile results into a DataFrame
def compile_results(results):
    compiled_data = []
    for metric, result in results.items():
        compiled_data.append({
            'Metric': metric,
            'Group1 Summary': result[0].to_dict(),
            'Group2 Summary': result[1].to_dict(),
            'Log Skewness Group1': result[2],
            'Log Skewness Group2': result[3],
            'T-test p-value (log)': result[5],
            'Mann-Whitney U test p-value': result[7],
            'Effect size (Cohen\'s d)': result[8]
        })
    return pd.DataFrame(compiled_data)

df_102_vs_38 = compile_results(results_102_vs_38)
df_104_vs_38 = compile_results(results_104_vs_38)
df_102_vs_104 = compile_results(results_102_vs_104)
df_102_vs_84 = compile_results(results_102_vs_84)
df_104_vs_84 = compile_results(results_104_vs_84)

# Display the DataFrame
print("NK-R102 vs NK-R38:")
print(df_102_vs_38)
print("\nNK-R104 vs NK-R38:")
print(df_104_vs_38)
print("\nNK-R102 vs NK-R104:")
print(df_102_vs_104)
print("\nNK-R102 vs NK-R84:")
print(df_102_vs_84)
print("\nNK-R104 vs NK-R84:")
print(df_104_vs_84)
