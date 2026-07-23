# FastICA-Ultrasound
Blind source separation of ultrasound signals using FastICA.

## Overview

This repository contains an R implementation for blind source separation of biomedical ultrasound signals using Independent Component Analysis (FastICA).

The workflow includes signal preprocessing, FastICA-based source separation, linear combination of the estimated independent components, quantitative validation, and graphical visualization using both experimental and simulated ultrasound signals.

## Workflow

1. Load ultrasound signals.
2. Signal preprocessing:
   - Butterworth low-pass filter
   - Moving average filter
   - Kalman filter
   - Savitzky–Golay filter
3. Signal normalization and centering.
4. Blind source separation using FastICA.
5. Linear combination of the estimated independent components.
6. Performance evaluation:
   - Pearson correlation
   - RMSE
   - Windowed stability analysis
   - Local correlation analysis
7. Time-domain and frequency-domain visualization.

## Requirements

R

Required packages:

- fastICA
- ggplot2
- gridExtra
- signal

## Input files

- mezcla1.txt
- mezcla2.txt
- originalUS.txt
- originalRF.txt

## Output

The script generates:

- Independent components
- Linearly combined signals
- Quantitative performance metrics
- Diagnostic plots
- FFT analysis



Department of Statistics and Operations Research

University of Granada
