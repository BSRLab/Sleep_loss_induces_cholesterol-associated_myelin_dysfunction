
"Sleep loss impairs myelin function by altering cholesterol metabolism in oligodendroglia"
Reyila Simayi, Eleonora Ficiarà, Oluwatomisin Faniyan, Antonio Cerdán Cerdá, Amina Aboufares El Alaoui, Rosamaria Fiorini, Adele Cutignano, Fabiana Piscitelli, Pablo Andújar, Alexandra Santos, Federico Del Gallo, Luisa de Vivo, Silvia De Santis, Michele Bellesi
bioRxiv 2023.11.27.568716; doi: https://doi.org/10.1101/2023.11.27.568716

Codes:
- Python Codes used to calculate the interhemispheric synchronization in the manuscript
- "LME_statistics_intensity.R" R code used for analysis based on linear mixed models 
- "analysis_outline_DE.R" R code used for differential expression analysis 
- "sleep_architecture_eeg.m" MATLAB code to calculate sleep architecture in EEG analysis 

SOFTWARE REQUIREMENTS
Python ≥ 3.10

Core libraries to install: mne 1.8.0, mneconnectivity 0.7.0, numpy 1.26.4, pandas 2.3.1,
matplotlib 3.8.0, scipy 1.11.4

R = 4.3.2
Required packages:
readxl_1.4.5 lmerTest_3.1-3 lme4_1.1-37 Matrix_1.6-5

DEMO
Input for the Python and R scripts
EDF file that contains the channels RIGHT and LEFT
Excel file with the columns EpochNo and Stage (NR, R, W)
Coherence_FREQ_DOMAIN_BANDS.py
o Expected output: tables reporting mean ± SD coherence for each sleep stage and for
every frequency band.
Correlation_interhemispheric_time_domain_TABLE.py
o Expected output: for each stage, a table with the mean ± SD of two metrics:
Crosscorrelation and Pearson correlation.
 Phase_locking_interhemispheric_time_TABLE.py
o Expected output: for stage NR, a table with the mean ± SD PhaseLocking Value
(PLV).

LME_statistics_intensity.R
o Input: the Excel file containing the columns Intensity (numeric), Condition
(categorical) and Rat (categorical).
o Expected output:
- likelihoodratio ANOVA table comparing two randomintercept mixed
models (null vs. + Condition)
- diagnostic plots: residual QQ plot, residualvsfit scatter, residual
histogram.
