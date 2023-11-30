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


# In[2]:


# Data import
raw= mne.io.read_raw_edf("C:/Users/Utente/Desktop/file_name", preload=True)# Import data file into workspace


# In[3]:


# Define the frequency range for the band-pass filter
freq_low = 0.15  # Lower frequency limit in Hz
freq_high = 40.0  # Upper frequency limit in Hz

# Apply band-pass filter to the data
raw.filter(freq_low, freq_high, fir_design='firwin')


# In[4]:


#raw


# In[5]:


data1, times = raw.copy().pick_channels(['PARIETAL_RIGHT'])[0]
data2, _ = raw.copy().pick_channels(['PARIETAL_LEFT'])[0]


# In[6]:


times


# In[7]:


min_length = min(len(data1), len(data2))
data1 = data1[:min_length]
data2 = data2[:min_length]


# In[8]:


print(data1.shape)
print(data2.shape)



# Sampling frequency
sfreq = raw.info['sfreq']
sfreq


window_size = 4  # Time interval for cross-correlation in seconds
window_samples = int(window_size * sfreq)  # Convert window size to samples

# Calculate the number of windows (number of epochs)
num_windows = (len(data1[0]) - window_samples) // window_samples
num_windows=num_windows+1

print("num_windows", num_windows)
print("window_samples",window_samples)

df_stage_info = pd.read_excel("C:/Users/Utente/Desktop/file_name_staging.xlsx")

window_indices = np.arange(num_windows)


num_windows_data = df_stage_info[df_stage_info['EpochNo'].isin(window_indices)]

# Separate the data into three groups based on the 'Stage' (NR, R, W)
stage_nr_data = num_windows_data[num_windows_data['Stage'] == 'NR']
stage_r_data = num_windows_data[num_windows_data['Stage'] == 'R']
stage_w_data = num_windows_data[num_windows_data['Stage'] == 'W']

num_windows_data['StartIdx'] = num_windows_data['EpochNo'] * window_samples



num_windows_data

time_values_nr = []
time_values_r = []
time_values_w = []


# In[20]:


#########coeherence frequency domain##########


# In[21]:


window_size = 4  # Time interval for coherence in seconds
window_samples = int(window_size * sfreq)  # Convert window size to samples

# Calculate the number of windows
num_windows = (len(data1[0]) - window_samples) // window_samples
num_windows=num_windows+1

# Initialize arrays to store magnitude squared coherence and epoch numbers
coherence_values = np.zeros(num_windows)
epoch_numbers = np.arange(num_windows)

# Calculate magnitude squared coherence for each window
for i in range(num_windows):
    start_idx = i * window_samples
    end_idx = start_idx + window_samples
    
    # Extract the data within the window
    data1_window = data1[0, start_idx:end_idx]
    data2_window = data2[0, start_idx:end_idx]
    
    # Calculate Power Spectral Density (PSD) for both channels
    Pxx, f = plt.mlab.psd(data1_window, NFFT=512, Fs=sfreq, scale_by_freq=True)
    Pyy, f = plt.mlab.psd(data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Cross Spectral Density (CSD) between the two channels
    Pxy, f = plt.mlab.csd(data1_window, data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Coherence
    coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    coherence_values[i] = np.mean(coherence)
    



# In[22]:


coherence_values_nr = coherence_values[num_windows_data['Stage'] == 'NR']
coherence_values_r = coherence_values[num_windows_data['Stage'] == 'R']
coherence_values_w = coherence_values[num_windows_data['Stage'] == 'W']


mean_coherence_nr = np.mean(coherence_values_nr)
mean_coherence_r = np.mean(coherence_values_r)
mean_coherence_w = np.mean(coherence_values_w)


std_coherence_nr = np.std(coherence_values_nr)
std_coherence_r = np.std(coherence_values_r)
std_coherence_w = np.std(coherence_values_w)


data = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr, mean_coherence_r, mean_coherence_w],
    'Standard Deviation': [std_coherence_nr, std_coherence_r, std_coherence_w]
}

coherence_table = pd.DataFrame(data)
print(coherence_table)
print(len(coherence_values_nr))

indices_nr = num_windows_data[num_windows_data['Stage'] == 'NR'].index


coherence_values_nr = np.zeros(len(indices_nr))
epoch_numbers_nr = np.arange(len(indices_nr))


for i, idx in enumerate(indices_nr):
    start_idx = idx * window_samples
    end_idx = start_idx + window_samples


    data1_window = data1[0, start_idx:end_idx]
    data2_window = data2[0, start_idx:end_idx]

    # Calculate Power Spectral Density (PSD) for both channels
    Pxx, f = plt.mlab.psd(data1_window, NFFT=512, Fs=sfreq, scale_by_freq=True)
    Pyy, _ = plt.mlab.psd(data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Cross Spectral Density (CSD) between the two channels
    Pxy, _ = plt.mlab.csd(data1_window, data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)


    coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    coherence_values_nr[i] = np.mean(coherence)

# Calculate the mean coherence for the 'NR' stage
mean_coherence_nr = np.mean(coherence_values_nr)

# Plot the coherence spectrum for the 'NR' stage
plt.plot(f, coherence)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Coherence')
plt.title('Coherence Spectrum for Stage NR')
plt.xlim(0, 50)  # Set the x-axis limit between 0 and 50


plt.axvline(x=0.5, color='blue', linestyle='--', label='f = 0.5')
plt.axvline(x=4, color='blue', linestyle='--', label='f = 4')
plt.axvline(x=9, color='blue', linestyle='--', label='f = 9')


plt.axhline(mean_coherence_nr, color='r', linestyle='--', label='Mean Coherence (NR)')

plt.legend()
plt.show()

len(coherence)


indices_r = num_windows_data[num_windows_data['Stage'] == 'R'].index


coherence_values_r = np.zeros(len(indices_r))
epoch_numbers_r = np.arange(len(indices_r))


for i, idx in enumerate(indices_r):
    start_idx = idx * window_samples
    end_idx = start_idx + window_samples


    data1_window = data1[0, start_idx:end_idx]
    data2_window = data2[0, start_idx:end_idx]

    Pxx, f = plt.mlab.psd(data1_window, NFFT=512, Fs=sfreq, scale_by_freq=True)
    Pyy, _ = plt.mlab.psd(data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Cross Spectral Density (CSD) between the two channels
    Pxy, _ = plt.mlab.csd(data1_window, data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Coherence
    coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    coherence_values_r[i] = np.mean(coherence)


mean_coherence_r = np.mean(coherence_values_r)


plt.plot(f, coherence,color='purple')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Coherence')
plt.title('Coherence Spectrum for Stage R')
plt.xlim(0, 50)

plt.axvline(x=0.5, color='blue', linestyle='--', label='f = 0.5')
plt.axvline(x=4, color='blue', linestyle='--', label='f = 4')
plt.axvline(x=9, color='blue', linestyle='--', label='f = 9')


plt.axhline(mean_coherence_r, color='r', linestyle='--', label='Mean Coherence (R)')

plt.legend()
plt.show()



indices_w = num_windows_data[num_windows_data['Stage'] == 'W'].index


coherence_values_w = np.zeros(len(indices_w))
epoch_numbers_w = np.arange(len(indices_w))


for i, idx in enumerate(indices_w):
    start_idx = idx * window_samples
    end_idx = start_idx + window_samples


    data1_window = data1[0, start_idx:end_idx]
    data2_window = data2[0, start_idx:end_idx]


    Pxx, f = plt.mlab.psd(data1_window, NFFT=512, Fs=sfreq, scale_by_freq=True)
    Pyy, _ = plt.mlab.psd(data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Cross Spectral Density (CSD) between the two channels
    Pxy, _ = plt.mlab.csd(data1_window, data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)

    # Calculate Coherence
    coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy)
    coherence_values_w[i] = np.mean(coherence)


mean_coherence_w = np.mean(coherence_values_w)

plt.plot(f, coherence,color='green')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Coherence')
plt.title('Coherence Spectrum for Stage W')
plt.xlim(0, 50)

plt.axvline(x=0.5, color='blue', linestyle='--', label='f = 0.5')
plt.axvline(x=4, color='blue', linestyle='--', label='f = 4')
plt.axvline(x=9, color='blue', linestyle='--', label='f = 9')


plt.axhline(mean_coherence_w, color='r', linestyle='--', label='Mean Coherence (W)')

plt.legend()
plt.show()

delta_band = (0.5, 4)
theta_band = (4, 9)
alpha_band = (9, 12)
sigma_band = (12, 15)
beta_band = (15, 25)
gamma_band = (25, 50)
mean_band=(0,50)


mean_coherence_delta = np.zeros(num_windows)
mean_coherence_theta = np.zeros(num_windows)
mean_coherence_alpha = np.zeros(num_windows)
mean_coherence_sigma = np.zeros(num_windows)
mean_coherence_beta = np.zeros(num_windows)
mean_coherence_gamma = np.zeros(num_windows)
mean_coherence_50 = np.zeros(num_windows)


for i in range(num_windows):
    start_idx = i * window_samples
    end_idx = start_idx + window_samples


    data1_window = data1[0, start_idx:end_idx]
    data2_window = data2[0, start_idx:end_idx]


    Pxx, f = plt.mlab.psd(data1_window, NFFT=512, Fs=sfreq, scale_by_freq=True)
    Pyy, _ = plt.mlab.psd(data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)


    Pxy, _ = plt.mlab.csd(data1_window, data2_window, NFFT=512, Fs=sfreq, scale_by_freq=True)


    coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy)


    mean_coherence_delta[i] = np.mean(coherence[(f >= delta_band[0]) & (f <= delta_band[1])])
    mean_coherence_theta[i] = np.mean(coherence[(f >= theta_band[0]) & (f <= theta_band[1])])
    mean_coherence_alpha[i] = np.mean(coherence[(f >= alpha_band[0]) & (f <= alpha_band[1])])
    mean_coherence_sigma[i] = np.mean(coherence[(f >= sigma_band[0]) & (f <= sigma_band[1])])
    mean_coherence_beta[i] = np.mean(coherence[(f >= beta_band[0]) & (f <= beta_band[1])])
    mean_coherence_gamma[i] = np.mean(coherence[(f >= gamma_band[0]) & (f <= gamma_band[1])])
    mean_coherence_50[i] = np.mean(coherence[(f >= mean_band[0]) & (f <= mean_band[1])])


    coherence_values[i] = np.mean(coherence)

len(coherence_values)


# In[29]:


coherence_values_nr_50 = mean_coherence_50[num_windows_data['Stage'] == 'NR']
coherence_values_r_50 = mean_coherence_50[num_windows_data['Stage'] == 'R']
coherence_values_w_50 = mean_coherence_50[num_windows_data['Stage'] == 'W']

mean_coherence_nr_50 = np.mean(coherence_values_nr_50)
mean_coherence_r_50 = np.mean(coherence_values_r_50)
mean_coherence_w_50 = np.mean(coherence_values_w_50)

std_coherence_nr_50 = np.std(coherence_values_nr_50)
std_coherence_r_50 = np.std(coherence_values_r_50)
std_coherence_w_50 = np.std(coherence_values_w_50)


data_50 = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_50, mean_coherence_r_50, mean_coherence_w_50],
    'Standard Deviation': [std_coherence_nr_50, std_coherence_r_50, std_coherence_w_50]
}

coherence_table_50 = pd.DataFrame(data_50)
print(coherence_table_50)
print(len(coherence_values_nr_50))

coherence_values_nr_delta = mean_coherence_delta[num_windows_data['Stage'] == 'NR']
coherence_values_r_delta = mean_coherence_delta[num_windows_data['Stage'] == 'R']
coherence_values_w_delta = mean_coherence_delta[num_windows_data['Stage'] == 'W']


mean_coherence_nr_delta = np.mean(coherence_values_nr_delta)
mean_coherence_r_delta = np.mean(coherence_values_r_delta)
mean_coherence_w_delta = np.mean(coherence_values_w_delta)


std_coherence_nr_delta = np.std(coherence_values_nr_delta)
std_coherence_r_delta = np.std(coherence_values_r_delta)
std_coherence_w_delta = np.std(coherence_values_w_delta)


data_delta = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_delta, mean_coherence_r_delta, mean_coherence_w_delta],
    'Standard Deviation': [std_coherence_nr_delta, std_coherence_r_delta, std_coherence_w_delta]
}

coherence_table_delta = pd.DataFrame(data_delta)
print(coherence_table_delta)
print(len(coherence_values_nr_delta))

coherence_values_nr_theta = mean_coherence_theta[num_windows_data['Stage'] == 'NR']
coherence_values_r_theta = mean_coherence_theta[num_windows_data['Stage'] == 'R']
coherence_values_w_theta = mean_coherence_theta[num_windows_data['Stage'] == 'W']


mean_coherence_nr_theta = np.mean(coherence_values_nr_theta)
mean_coherence_r_theta = np.mean(coherence_values_r_theta)
mean_coherence_w_theta = np.mean(coherence_values_w_theta)


std_coherence_nr_theta = np.std(coherence_values_nr_theta)
std_coherence_r_theta = np.std(coherence_values_r_theta)
std_coherence_w_theta = np.std(coherence_values_w_theta)


data_theta = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_theta, mean_coherence_r_theta, mean_coherence_w_theta],
    'Standard Deviation': [std_coherence_nr_theta, std_coherence_r_theta, std_coherence_w_theta]
}

coherence_table_theta = pd.DataFrame(data_theta)
print(coherence_table_theta)
print(len(coherence_values_nr_theta))

coherence_values_nr_sigma = mean_coherence_sigma[num_windows_data['Stage'] == 'NR']
coherence_values_r_sigma = mean_coherence_sigma[num_windows_data['Stage'] == 'R']
coherence_values_w_sigma = mean_coherence_sigma[num_windows_data['Stage'] == 'W']


mean_coherence_nr_sigma = np.mean(coherence_values_nr_sigma)
mean_coherence_r_sigma = np.mean(coherence_values_r_sigma)
mean_coherence_w_sigma = np.mean(coherence_values_w_sigma)

std_coherence_nr_sigma = np.std(coherence_values_nr_sigma)
std_coherence_r_sigma = np.std(coherence_values_r_sigma)
std_coherence_w_sigma = np.std(coherence_values_w_sigma)

data_sigma = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_sigma, mean_coherence_r_sigma, mean_coherence_w_sigma],
    'Standard Deviation': [std_coherence_nr_sigma, std_coherence_r_sigma, std_coherence_w_sigma]
}

coherence_table_sigma = pd.DataFrame(data_sigma)
print(coherence_table_sigma)
print(len(coherence_values_nr_sigma))

coherence_values_nr_alpha = mean_coherence_alpha[num_windows_data['Stage'] == 'NR']
coherence_values_r_alpha = mean_coherence_alpha[num_windows_data['Stage'] == 'R']
coherence_values_w_alpha = mean_coherence_alpha[num_windows_data['Stage'] == 'W']


mean_coherence_nr_alpha = np.mean(coherence_values_nr_alpha)
mean_coherence_r_alpha = np.mean(coherence_values_r_alpha)
mean_coherence_w_alpha = np.mean(coherence_values_w_alpha)

std_coherence_nr_alpha = np.std(coherence_values_nr_alpha)
std_coherence_r_alpha = np.std(coherence_values_r_alpha)
std_coherence_w_alpha = np.std(coherence_values_w_alpha)

data_alpha = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_alpha, mean_coherence_r_alpha, mean_coherence_w_alpha],
    'Standard Deviation': [std_coherence_nr_alpha, std_coherence_r_alpha, std_coherence_w_alpha]
}

coherence_table_alpha = pd.DataFrame(data_alpha)
print(coherence_table_alpha)
print(len(coherence_values_nr_alpha))

coherence_values_nr_beta = mean_coherence_beta[num_windows_data['Stage'] == 'NR']
coherence_values_r_beta = mean_coherence_beta[num_windows_data['Stage'] == 'R']
coherence_values_w_beta = mean_coherence_beta[num_windows_data['Stage'] == 'W']


mean_coherence_nr_beta = np.mean(coherence_values_nr_beta)
mean_coherence_r_beta = np.mean(coherence_values_r_beta)
mean_coherence_w_beta = np.mean(coherence_values_w_beta)

std_coherence_nr_beta = np.std(coherence_values_nr_beta)
std_coherence_r_beta = np.std(coherence_values_r_beta)
std_coherence_w_beta = np.std(coherence_values_w_beta)

data_beta = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_beta, mean_coherence_r_beta, mean_coherence_w_beta],
    'Standard Deviation': [std_coherence_nr_beta, std_coherence_r_beta, std_coherence_w_beta]
}

coherence_table_beta = pd.DataFrame(data_beta)
print(coherence_table_beta)
print(len(coherence_values_nr_beta))


coherence_values_nr_gamma = mean_coherence_gamma[num_windows_data['Stage'] == 'NR']
coherence_values_r_gamma = mean_coherence_gamma[num_windows_data['Stage'] == 'R']
coherence_values_w_gamma = mean_coherence_gamma[num_windows_data['Stage'] == 'W']


mean_coherence_nr_gamma = np.mean(coherence_values_nr_gamma)
mean_coherence_r_gamma = np.mean(coherence_values_r_gamma)
mean_coherence_w_gamma = np.mean(coherence_values_w_gamma)

std_coherence_nr_gamma = np.std(coherence_values_nr_gamma)
std_coherence_r_gamma = np.std(coherence_values_r_gamma)
std_coherence_w_gamma = np.std(coherence_values_w_gamma)

data_gamma = {
    'Stage': ['NR', 'R', 'W'],
    'Mean Coherence': [mean_coherence_nr_gamma, mean_coherence_r_gamma, mean_coherence_w_gamma],
    'Standard Deviation': [std_coherence_nr_gamma, std_coherence_r_gamma, std_coherence_w_gamma]
}

coherence_table_gamma = pd.DataFrame(data_gamma)
print(coherence_table_gamma)
print(len(coherence_values_nr_gamma))

len(mean_coherence_alpha)

len(coherence_values_r)
plt.figure(figsize=(10, 8))

plt.subplot(311)
plt.axhline(mean_coherence_nr, color='red', linestyle='dashed', label='Mean Coherence')
plt.xlabel('Epoch Number')
plt.ylabel('Coherence')
plt.title('Coherence for Stage NR')
plt.legend()
plt.ylim(0, 1)
plt.vlines(epoch_numbers[num_windows_data['Stage'] == 'NR'], 0, coherence_values_nr, colors='blue', linewidth=1)

plt.subplot(312)
plt.axhline(mean_coherence_r, color='red', linestyle='dashed', label='Mean Coherence')
plt.xlabel('Epoch Number')
plt.ylabel('Coherence')
plt.title('Coherence for Stage R')
plt.legend()
plt.ylim(0, 1)
plt.vlines(epoch_numbers[num_windows_data['Stage'] == 'R'], 0, coherence_values_r, colors='purple', linewidth=1)

plt.subplot(313)
plt.axhline(mean_coherence_w, color='red', linestyle='dashed', label='Mean Coherence')
plt.xlabel('Epoch Number')
plt.ylabel('Coherence')
plt.title('Coherence for Stage W')
plt.legend()
plt.ylim(0, 1)
plt.vlines(epoch_numbers[num_windows_data['Stage'] == 'W'], 0, coherence_values_w, colors='green', linewidth=1)

plt.tight_layout()
plt.show()

