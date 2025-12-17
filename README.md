# ChronoPhot-visualization-module
Visualization module for long-term fiber photometry data, generating SVG plots of corrected signals, ΔF/F, and Z-score traces.

# ChronoPhot-View

ChronoPhot-View is a Python-based visualization module for long-term fiber photometry data.
This script is designed to generate publication-ready figures and to facilitate quality control
of extended photometry recordings.

## Overview

The code processes preprocessed fiber photometry data files and generates SVG plots of
raw signals, tdTomato-corrected signals, ΔF/F, and Z-score–normalized traces across long
recording periods. It is optimized for long-term recordings (e.g., several hours) and
automatically excludes initial zero or baseline segments when specified.

ChronoPhot-View is intended to be used as a companion visualization tool to batch analysis
pipelines, enabling inspection and presentation of photometry signals after preprocessing.

## Features

- Visualization of long-term fiber photometry recordings
- tdTomato-based regression correction for motion and bleaching
- Signal detrending and normalization
- ΔF/F and Z-score trace plotting
- Automatic generation of SVG figures suitable for publication
- Batch processing of multiple CSV files

## Input

- CSV files containing fiber photometry signals
- Signal channels for calcium-dependent fluorescence and control (tdTomato)

## Output

- SVG plots of:
  - Raw and corrected fluorescence signals
  - ΔF/F traces
  - Z-score–normalized signals
- Output files are generated automatically for each input dataset

## Requirements

- Python 3.x
- numpy
- pandas
- matplotlib
- scipy

## Usage

Edit the input and output directory paths in the script as needed, then run:

```bash
python chronophot_view.py
