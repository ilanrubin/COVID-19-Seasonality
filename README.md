# Seasonal forcing and waning immunity drive the sub-annual periodicity of the COVID-19 epidemic

This repository includes the code used in the following manuscript https://www.medrxiv.org/content/10.1101/2025.03.05.25323464v2.

### Wavelet Analysis
Wavelet analysis is in the directory wavelets.

*county_wavelet_analysis.R* - What the user will likely want to run. Runs the wavelet analysis over all FIPS (county and county equivalents) in the United States and plots the results.

*wavelet_function.R* - Contains the functions that run the wavelet analysis

*peaks_function.R* - Contains the peak finding algorithm.

*plot_function.R* - A simple function to plot cases.

The rest of the files load the necessary data. The data directory contains a few files that are not easily downloaded via an api (e.g., they are downloaded as a .zip and have to be uncompressed first), but all are still publically available.

### SIRS Model with Waning Immunity
An individual-based implementation of the SIRS model with waning immunity is in the directory SIRS.

*run_model_LHS.R* - What the user will likely want to run. Is a script to do Latine Hyopercube Sampling and run the model over defined parameter ranges.

*SIRS_with_waning.R* - Contains the model.

*peaks_function.R* - Contains the peak finding algorithm.

*plot_LHS.R* - Plots the LHS results.

*parameters.dat* - The default parameterization that gets read into the model.
