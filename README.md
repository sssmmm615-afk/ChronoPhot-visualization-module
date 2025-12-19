# Fiber Photometry Control vs Test Batch Analysis

## Introduction

This repository provides a Python-based batch analysis pipeline for comparing **Control** and **Test** groups in long-term fiber photometry experiments. The script is designed for **lock-in demodulated CSV files** and enables fully automated preprocessing, normalization, and group-level visualization suitable for publication and data sharing.

The pipeline processes raw photometry data on a per-animal basis, applies standardized correction and normalization steps, and generates both individual outputs and group mean ± SEM plots for direct comparison between experimental conditions.

---

## Data Organization

The recommended folder structure is as follows:


Each CSV file is assumed to correspond to a single animal and to contain a time column (in seconds) together with two fluorescence channels: a calcium-dependent signal channel (e.g., 465 nm / GFP) and a reference or isosbestic channel (e.g., 405 nm / Tomato).  
The script automatically detects the header row and standardizes column names based on common keywords.

---

## Data Processing Workflow

For each animal in the Control and Test folders, the pipeline performs the following steps in a batch-oriented manner.

First, **photobleaching correction** is applied independently to both fluorescence channels. A linear trend is estimated using an early recording window and a late window anchored to the end of the recording, and this trend is subtracted from the raw signal to correct slow signal drift over long durations.

Next, **motion correction** is performed using linear regression. The reference channel (405 nm) is fitted to the signal channel (465 nm), and the fitted component is removed from the signal channel to reduce motion-related and shared noise while preserving calcium-dependent dynamics.

After these corrections, the motion-corrected signal is converted to a **Z-score–normalized trace** using a **fixed global baseline interval** (default: 25–35 minutes after recording onset). The same baseline definition is applied to all animals and both experimental groups.

---

## Baseline Strategy

Z-score normalization relies on a single, predefined global baseline interval. This conservative approach is well suited for long-term recordings and avoids window-specific or adaptive baseline recalculation. As a result, differences observed between Control and Test groups reflect genuine signal dynamics rather than changes in normalization criteria.

Photobleaching correction and Z-score normalization are treated as conceptually independent steps.

---

## Group-level Analysis and Visualization

After per-animal processing, Z-score traces from each group are aligned to a common time axis and binned to one-second resolution. For each time point, the pipeline computes the **group mean and standard error of the mean (SEM)** across animals.

The final output is an SVG figure showing **mean ± SEM Z-score traces** for Control and Test groups over the specified long-term time window. This figure is suitable for direct inclusion in manuscripts or supplementary materials.

---

## Output Files

The pipeline generates the following outputs automatically:

- **Per-animal processed CSV files**  
  Containing corrected fluorescence signals and Z-score–normalized traces.

- **Group-level SVG figure**  
  Mean ± SEM Z-score traces comparing Control and Test groups.

All results are written to the `output/` directory, organized by group.

---

## Reproducibility and Configuration

All key parameters—including baseline intervals, photobleaching windows, and plotting ranges—are explicitly defined and can be adjusted via command-line arguments. The pipeline requires no interactive steps or manual intervention, ensuring that analyses can be rerun under identical conditions.

This repository is intended to serve as a transparent and reproducible record of the analysis pipeline used for Control vs Test comparisons in long-term fiber photometry experiments, and to support data sharing alongside manuscripts or collaborative projects.

---

## Intended Use

This pipeline is particularly suited for:

- Long-duration fiber photometry recordings
- Group comparisons with delayed or sustained neural responses
- Hormonal, metabolic, or pharmacological intervention studies
- Reproducible batch processing for publication-quality analysis
