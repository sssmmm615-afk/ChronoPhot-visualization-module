ChronoPhot-View
Overview

ChronoPhot-View is a Python-based visualization module for long-term fiber photometry data. It generates publication-ready SVG figures and facilitates quality control of extended recordings, typically spanning several hours. This tool is designed to visualize signals obtained from long-duration photometry experiments and to support inspection of signal stability and long-term dynamics after preprocessing.

What This Tool Does

ChronoPhot-View loads preprocessed fiber photometry CSV files and produces continuous visualizations of fluorescence signals across extended recording periods. It supports visualization of corrected fluorescence signals, ΔF/F traces, and Z-score–normalized signals, enabling users to evaluate baseline stability, slow signal drift, and long-term fluctuations in neural activity. By presenting multiple signal representations within a unified visualization workflow, the module facilitates systematic quality control and efficient generation of figures for long-term recordings.

Design Philosophy

This module is optimized for long-duration recordings rather than trial-aligned or event-triggered analyses. It is intended to be used as a companion visualization tool for batch photometry analysis pipelines, focusing on inspection and presentation of already preprocessed data. ChronoPhot-View supports batch processing of multiple datasets and automatically generates scalable vector graphics (SVG) figures suitable for direct inclusion in publications. Optional exclusion of initial zero or baseline segments is supported to improve the interpretability of long-term traces.

Intended Use

ChronoPhot-View focuses exclusively on visualization and post-processing inspection rather than raw signal correction or preprocessing. As such, it complements existing photometry analysis pipelines by providing a lightweight and reproducible solution for inspecting, validating, and presenting fiber photometry data following preprocessing.

License

This software is distributed under the MIT License.

```bash
python chronophot_view.py
