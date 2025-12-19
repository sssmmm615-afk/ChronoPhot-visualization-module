# Fiber Photometry Control vs Test Batch Analysis

## Introduction
This repository provides a Python-based batch analysis pipeline for comparing Control and Test groups in long-term fiber photometry experiments. The script is designed for lock-in demodulated CSV files and enables fully automated preprocessing, normalization, and group-level visualization suitable for publication and data sharing.

The pipeline processes raw photometry data on a per-animal basis, applies standardized correction and normalization steps, and generates both individual outputs and group mean ± SEM plots for direct comparison between experimental conditions.

---

## Data Organization
The recommended folder structure is as follows:


Each CSV file is assumed to correspond to a single animal and to contain:
- A time column (in seconds)
- A calcium-dependent signal channel (e.g., 465 nm, GFP)
- A reference or isosbestic channel (e.g., 405 nm)

The script automatically detects the header row and standardizes column names based on common keywords (e.g., `time`, `gfp`, `tomato`).  
If your CSV files use different labels (such as `405`, `isosbestic`, or `470`), the keyword mapping in the script may need to be adjusted accordingly.

---

## Data Processing Workflow
For each animal in the Control and Test folders, the pipeline performs the following steps in a batch-oriented manner.

### 1. Photobleaching Correction
Photobleaching correction is applied independently to both fluorescence channels.  
A linear trend is estimated using two predefined windows: an early baseline window near the start of the recording and a late window anchored to the end of the recording. The trend defined by the mean signals in these windows is subtracted from the raw trace to correct slow signal drift over long durations.

### 2. Motion Correction
Motion correction is performed using linear regression.  
The reference channel (e.g., 405 nm) is fitted to the signal channel (e.g., 465 nm), and the fitted component is removed from the signal channel. This procedure reduces motion-related and shared noise while preserving calcium-dependent dynamics.

### 3. Z-score Normalization
After photobleaching and motion correction, the corrected signal is converted to a Z-score–normalized trace using a fixed global baseline interval (default: 25–35 minutes after recording onset).  
The same baseline definition is applied to all animals and to both experimental groups.

Photobleaching correction and Z-score normalization are treated as conceptually independent steps.

---

## Baseline Strategy
Z-score normalization relies on a single, predefined global baseline interval. This conservative approach is particularly well suited for long-term recordings and avoids window-specific or adaptive baseline recalculation.

As a result, differences observed between Control and Test groups reflect genuine signal dynamics rather than differences arising from normalization criteria.

---

## Group-level Analysis and Visualization
After per-animal processing, Z-score traces within each group are binned to 1-second resolution and merged by time.  
For each time point, the pipeline computes the group mean and the standard error of the mean (SEM) across animals. NaN values are ignored when calculating the mean; SEM is computed using the total number of animals in each group.

The final output is an SVG figure showing mean ± SEM Z-score traces for Control and Test groups over the specified long-term time window. This figure is suitable for direct inclusion in manuscripts or supplementary materials.

---

## Output Files
The pipeline automatically generates the following outputs:

### Per-animal processed CSV files
Stored separately for Control and Test groups, these files contain:
- Photobleaching-corrected fluorescence signals  
- Motion-corrected signal traces  
- Z-score–normalized signals  

### Group-level SVG figure
An SVG file displaying mean ± SEM Z-score traces comparing Control and Test groups.

All results are written to the `output/` directory and organized by experimental group.

---

## Reproducibility and Configuration
All key parameters—including baseline intervals, photobleaching windows, and plotting ranges—are explicitly defined and can be adjusted via command-line arguments.

The pipeline requires no interactive steps or manual intervention, ensuring that analyses can be rerun under identical conditions and supporting transparent, reproducible research workflows.

---

## Intended Use
This pipeline is particularly suited for:
- Long-duration fiber photometry recordings  
- Group comparisons involving delayed or sustained neural responses  
- Hormonal, metabolic, or pharmacological intervention studies  
- Reproducible batch processing for publication-quality analysis  

This repository is intended to serve as a transparent and reproducible record of the analysis pipeline used for Control vs Test comparisons in long-term fiber photometry experiments and to support data sharing alongside manuscripts or collaborative projects.

