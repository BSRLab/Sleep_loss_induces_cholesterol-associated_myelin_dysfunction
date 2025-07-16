#!/usr/bin/env python
# coding: utf-8

# In[1]:


import mne
import numpy as np
import matplotlib.pyplot as plt
from mne_connectivity import spectral_connectivity_epochs
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy import signal



raw= mne.io.read_raw_edf("C:/.../EDF_example.edf", preload=True)

freq_low = 0.15  # Lower frequency limit in Hz
freq_high = 40.0  # Upper frequency limit in Hz


raw=raw.filter(freq_low, freq_high, fir_design='firwin')
raw

data1, times = raw.copy().pick_channels(['PARIRTAL_RIGHT'])[0]

data2, _ = raw.copy().pick_channels(['PARIETAL_LEFT'])[0]


times

min_length = min(len(data1), len(data2))
data1 = data1[:min_length]
data2 = data2[:min_length]


print(data1.shape)
print(data2.shape)


data1


data2



sfreq = raw.info['sfreq']
sfreq


window_size = 4
window_samples = int(window_size * sfreq)


num_windows = (len(data1[0]) - window_samples) // window_samples
num_windows=num_windows+1


print("num_windows", num_windows)
print("window_samples",window_samples)

df_stage_info = pd.read_excel("C:/.../STAGING_example.xlsx")


window_indices = np.arange(num_windows)


num_windows_data = df_stage_info[df_stage_info['EpochNo'].isin(window_indices)]

# Separate the data into three groups based on the 'Stage' (NR, R, W)
stage_nr_data = num_windows_data[num_windows_data['Stage'] == 'NR']
stage_r_data = num_windows_data[num_windows_data['Stage'] == 'R']
stage_w_data = num_windows_data[num_windows_data['Stage'] == 'W']

num_windows_data['StartIdx'] = num_windows_data['EpochNo'] * window_samples

num_windows_data


threshold=-1
cross_correlation_nr = []
cross_correlation_r = []
cross_correlation_w = []
time_values_nr = []
time_values_r = []
time_values_w = []
time_shift_nr = []
time_shift_r = []
time_shift_w = []
time_delay_nr = []
time_delay_r = []
time_delay_w = []


for i, row in num_windows_data[num_windows_data['Stage'] == 'NR'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1=np.mean(data1_window)
        mean2=np.mean(data2_window)
        cross_corr_matrix = np.correlate(data1_window, data2_window, mode='full')

        time_shift = np.argmax(cross_corr_matrix) - len(data1_window) + 1  # Calculate time shift
        time_delay = time_shift / sfreq  # Calculate time delay in seconds
        cross_corr = cross_corr_matrix[np.argmax(cross_corr_matrix)]

        if cross_corr >= threshold:
            cross_correlation_nr.append(cross_corr)
        

            time_values_nr.append(row['EpochNo'])  # Use 'EpochNo' directly as x-axis values
            time_shift_nr.append(time_shift)
            time_delay_nr.append(time_delay)

        

for i, row in num_windows_data[num_windows_data['Stage'] == 'R'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1=np.mean(data1_window)
        mean2=np.mean(data2_window)
        cross_corr_matrix = np.correlate(data1_window, data2_window, mode='full')
        time_shift = np.argmax(cross_corr_matrix) - len(data1_window) + 1  # Calculate time shift
        time_delay = time_shift / sfreq  # Calculate time delay in seconds
        cross_corr = cross_corr_matrix[np.argmax(cross_corr_matrix)]
                # Check if the correlation value is above the threshold
        if cross_corr >= threshold:
            cross_correlation_r.append(cross_corr)
            time_values_r.append(row['EpochNo'])  # Use 'EpochNo' directly as x-axis values
            time_shift_r.append(time_shift)
            time_delay_r.append(time_delay)

for i, row in num_windows_data[num_windows_data['Stage'] == 'W'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1=np.mean(data1_window)
        mean2=np.mean(data2_window)
        cross_corr_matrix = np.correlate(data1_window, data2_window, mode='full')
        time_shift = np.argmax(cross_corr_matrix) - len(data1_window) + 1  # Calculate time shift
        time_delay = time_shift / sfreq  # Calculate time delay in seconds
        cross_corr =  cross_corr_matrix[np.argmax(cross_corr_matrix)]
        if cross_corr >= threshold:
            cross_correlation_w.append(cross_corr)
       
            time_values_w.append(row['EpochNo'])  # Use 'EpochNo' directly as x-axis values
            time_shift_w.append(time_shift)
            time_delay_w.append(time_delay)

# Calculate the mean cross-correlation for each stage
mean_cross_correlation_nr = np.mean(cross_correlation_nr)
mean_cross_correlation_r = np.mean(cross_correlation_r)
mean_cross_correlation_w = np.mean(cross_correlation_w)

mean_time_shift_nr = np.mean(time_shift_nr)
mean_time_shift_r = np.mean(time_shift_r)
mean_time_shift_w = np.mean(time_shift_w)

mean_time_delay_nr = np.mean(time_delay_nr)
mean_time_delay_r = np.mean(time_delay_r)
mean_time_delay_w = np.mean(time_delay_w)


print(data2_window/mean2)

plt.hist(cross_correlation_nr)
len(cross_correlation_nr)

print(mean_cross_correlation_nr)

plt.figure(figsize=(10, 6))

# Plot for Stage NR
plt.subplot(311)
plt.vlines(time_values_nr, ymin=0, ymax=cross_correlation_nr, color='blue', alpha=0.7)
plt.axhline(mean_cross_correlation_nr, color='red', linestyle='dashed', label='Mean')
plt.xlabel('Epoch Number')
plt.ylabel('Cross-Correlation')
plt.title('Stage NR')
#plt.ylim(0, 0.002) eventually for the same x axis
plt.grid(False)
plt.legend()

# Plot for Stage R
plt.subplot(312)
plt.vlines(time_values_r, ymin=0, ymax=cross_correlation_r,color='purple' , alpha=0.7)
plt.axhline(mean_cross_correlation_r, color='red', linestyle='dashed', label='Mean')
plt.xlabel('Epoch Number')
plt.ylabel('Cross-Correlation')
plt.title('Stage R')
plt.grid(False)
plt.legend()

# Plot for Stage W
plt.subplot(313)
plt.vlines(time_values_w, ymin=0, ymax=cross_correlation_w,color='green' , alpha=0.7)
plt.axhline(mean_cross_correlation_w, color='red', linestyle='dashed', label='Mean')
plt.xlabel('Epoch Number')
plt.ylabel('Cross-Correlation')
plt.title('Stage W')
plt.grid(False)
plt.legend()

plt.tight_layout()
plt.show()


std_cross_correlation_nr = np.std(cross_correlation_nr)
std_cross_correlation_r = np.std(cross_correlation_r)
std_cross_correlation_w = np.std(cross_correlation_w)

data = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Cross-Correlation': [mean_cross_correlation_nr, mean_cross_correlation_r, mean_cross_correlation_w],
    'Standard Deviation': [std_cross_correlation_nr, std_cross_correlation_r, std_cross_correlation_w]
}

result_table = pd.DataFrame(data)


print(result_table)



num_data_points_nr = len(cross_correlation_nr)
num_data_points_r = len(cross_correlation_r)
num_data_points_w = len(cross_correlation_w)


print("Number of data points for Stage NR:", num_data_points_nr)
print("Number of data points for Stage R:", num_data_points_r)
print("Number of data points for Stage W:", num_data_points_w)
total_data_points = num_data_points_nr + num_data_points_r + num_data_points_w


print("Total number of data points for all stages:", total_data_points)

###pearsonr
threshold=0.50
correlation_nr = []
correlation_r = []
correlation_w = []
time_values_nr = []
time_values_r = []
time_values_w = []

for i, row in num_windows_data[num_windows_data['Stage'] == 'NR'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1 = np.mean(data1_window)
        mean2 = np.mean(data2_window)
        cross_corr_matrix = np.corrcoef(data1_window, data2_window)
        cross_corr = cross_corr_matrix[0, 1]

        if cross_corr >= threshold:
            correlation_nr.append(cross_corr)
            time_values_nr.append(row['EpochNo'])



for i, row in num_windows_data[num_windows_data['Stage'] == 'R'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1 = np.mean(data1_window)
        mean2 = np.mean(data2_window)
        cross_corr_matrix = np.corrcoef(data1_window, data2_window)
        cross_corr = cross_corr_matrix[0, 1]

        if cross_corr >= threshold:
            correlation_r.append(cross_corr)
            time_values_r.append(row['EpochNo'])


for i, row in num_windows_data[num_windows_data['Stage'] == 'W'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]
        mean1 = np.mean(data1_window)
        mean2 = np.mean(data2_window)
        cross_corr_matrix = np.corrcoef(data1_window, data2_window)
        cross_corr = cross_corr_matrix[0, 1]

        if cross_corr >= threshold:
            correlation_w.append(cross_corr)
            time_values_w.append(row['EpochNo'])



mean_correlation_nr = np.mean(correlation_nr)
mean_correlation_r = np.mean(correlation_r)
mean_correlation_w = np.mean(correlation_w)

print(mean_correlation_nr)
print(mean_correlation_r)
print(mean_correlation_w)
plt.hist(correlation_nr)

mean_cross_correlation_nr = np.mean(cross_correlation_nr)
mean_cross_correlation_r = np.mean(cross_correlation_r)
mean_cross_correlation_w = np.mean(cross_correlation_w)

std_cross_correlation_nr = np.std(cross_correlation_nr)
std_cross_correlation_r = np.std(cross_correlation_r)
std_cross_correlation_w = np.std(cross_correlation_w)

mean_correlation_nr = np.mean(correlation_nr)
mean_correlation_r = np.mean(correlation_r)
mean_correlation_w = np.mean(correlation_w)

std_correlation_nr = np.std(correlation_nr)
std_correlation_r = np.std(correlation_r)
std_correlation_w = np.std(correlation_w)


data = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Cross-Correlation': [mean_cross_correlation_nr, mean_cross_correlation_r, mean_cross_correlation_w],
    'Standard Deviation (Cross-Correlation)': [std_cross_correlation_nr, std_cross_correlation_r, std_cross_correlation_w],
        'Mean Correlation': [mean_correlation_nr, mean_correlation_r, mean_correlation_w],
    'Standard Deviation (Correlation)': [std_correlation_nr, std_correlation_r, std_correlation_w]
}

result_table = pd.DataFrame(data)

# Display the table
print(result_table)
result_table

