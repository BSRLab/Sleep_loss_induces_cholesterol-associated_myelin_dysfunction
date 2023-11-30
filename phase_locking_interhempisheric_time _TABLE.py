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
from scipy.signal import hilbert


# Data import
raw= mne.io.read_raw_edf("C:/Users/Utente/Desktop/FABRIZIO/file_name", preload=True)# Import data file into workspace


# Define the frequency range for the band-pass filter
freq_low = 0.5  # Lower frequency limit in Hz
freq_high = 9.0  # Upper frequency limit in Hz

raw=raw.filter(freq_low, freq_high, fir_design='firwin')


data1, times = raw.copy().pick_channels(['PARIETAL_RIGHT'])[0]
#data1, times = raw.copy().pick_channels(['FRONTAL_RIGHT'])[0]
#data2, _ = raw.copy().pick_channels(['FRONTAL_LEFT'])[0]
data2, _ = raw.copy().pick_channels(['PARIETAL_LEFT'])[0]


times


min_length = min(len(data1), len(data2))
data1 = data1[:min_length]
data2 = data2[:min_length]


print(data1.shape)
print(data2.shape)


data1


data2

# Sampling frequency
sfreq = raw.info['sfreq']
sfreq



window_size = 4
window_samples = int(window_size * sfreq)

num_windows = (len(data1[0]) - window_samples) // window_samples
num_windows=num_windows+1


print("num_windows", num_windows)
print("window_samples",window_samples)


df_stage_info = pd.read_excel("C:/Users/Utente/Desktop/Ffile_name_staging.xlsx")

window_indices = np.arange(num_windows)


num_windows_data = df_stage_info[df_stage_info['EpochNo'].isin(window_indices)]

# Separate the data into three groups based on the 'Stage' (NR, R, W)
stage_nr_data = num_windows_data[num_windows_data['Stage'] == 'NR']
stage_r_data = num_windows_data[num_windows_data['Stage'] == 'R']
stage_w_data = num_windows_data[num_windows_data['Stage'] == 'W']


num_windows_data['StartIdx'] = num_windows_data['EpochNo'] * window_samples
num_windows_data



#####repeat for stage R and W

len(stage_nr_data)

plv_nr = []
plv_r = []
plv_w = []
time_values_nr = []
time_values_r = []
time_values_w = []

for i, row in num_windows_data[num_windows_data['Stage'] == 'NR'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]

        # Apply Hilbert transform to extract the analytic signal
        analytic_signal1 = hilbert(data1_window)
        analytic_signal2 = hilbert(data2_window)

        # Compute the phase difference between the two analytic signals
        phase_diff = np.angle(analytic_signal1) - np.angle(analytic_signal2)

        # Calculate the Phase Locking Value (PLV)
        plv = np.abs(np.mean(np.exp(1j * phase_diff)))

        plv_nr.append(plv)
        time_values_nr.append(row['EpochNo'])

mean_plv_nr = np.mean(plv_nr)
std_plv_nr = np.std(plv_nr)
#mean_plv_r = np.mean(plv_r)
#mean_plv_w = np.mean(plv_w)

print("mean_plv_nr",mean_plv_nr)
print("std_plv_nr",std_plv_nr)
len(plv_nr)


# Data import FOR REC 2
raw= mne.io.read_raw_edf("C:/Users/Utente/Desktop/file_name_rec2.edf", preload=True)# Import data file into workspace

freq_low = 0.5  # Lower frequency limit in Hz
freq_high = 9.0  # Upper frequency limit in Hz


raw=raw.filter(freq_low, freq_high, fir_design='firwin')

data1, times = raw.copy().pick_channels(['PARIETAL_RIGHT'])[0]
data2, _ = raw.copy().pick_channels(['PARIETAL_LEFT'])[0]


min_length = min(len(data1), len(data2))
data1 = data1[:min_length]
data2 = data2[:min_length]

print(data1.shape)
print(data2.shape)


data1

data2

# Sampling frequency
sfreq = raw.info['sfreq']
sfreq

window_size = 4
window_samples = int(window_size * sfreq)


num_windows = (len(data1[0]) - window_samples) // window_samples
num_windows=num_windows+1

print("num_windows", num_windows)
print("window_samples",window_samples)


df_stage_info = pd.read_excel("C:/Users/Utente/Desktop/file_name_staging_rec2.xlsx")

window_indices = np.arange(num_windows)

num_windows_data = df_stage_info[df_stage_info['EpochNo'].isin(window_indices)]

stage_nr_data = num_windows_data[num_windows_data['Stage'] == 'NR']
stage_r_data = num_windows_data[num_windows_data['Stage'] == 'R']
stage_w_data = num_windows_data[num_windows_data['Stage'] == 'W']

num_windows_data['StartIdx'] = num_windows_data['EpochNo'] * window_samples

num_windows_data

len(stage_nr_data)

plv_nr_rec2 = []
plv_r = []
plv_w = []
time_values_nr = []
time_values_r = []
time_values_w = []

for i, row in num_windows_data[num_windows_data['Stage'] == 'NR'].iterrows():
    start_idx = row['StartIdx']
    end_idx = start_idx + window_samples
    if end_idx <= len(data1[0]):
        data1_window = data1[0, start_idx:end_idx]
        data2_window = data2[0, start_idx:end_idx]

        # Apply Hilbert transform to extract the analytic signal
        analytic_signal1 = hilbert(data1_window)
        analytic_signal2 = hilbert(data2_window)

        # Compute the phase difference between the two analytic signals
        phase_diff = np.angle(analytic_signal1) - np.angle(analytic_signal2)

        # Calculate the Phase Locking Value (PLV)
        plv = np.abs(np.mean(np.exp(1j * phase_diff)))

        plv_nr_rec2.append(plv)
        time_values_nr.append(row['EpochNo'])



# Calculate the mean PLV for each stage
mean_plv_nr_rec2 = np.mean(plv_nr_rec2)
std_plv_nr_rec2 = np.std(plv_nr_rec2)
#mean_plv_r = np.mean(plv_r)
#mean_plv_w = np.mean(plv_w)


print("mean_plv_nr_rec2",mean_plv_nr_rec2)
print("std_plv_nr_rec2",std_plv_nr_rec2)
len(plv_nr_rec2)

