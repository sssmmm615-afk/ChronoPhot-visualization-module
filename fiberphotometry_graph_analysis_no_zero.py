#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Photometry Control vs Test (portable GitHub-ready; paper-matched)

- Batch process raw CSVs in Control/Test folders
- Photobleaching correction
- Motion correction (405 -> 465 regression)
- Z-score normalization using a global baseline interval
- Saves per-animal processed CSVs and group mean±SEM SVG

Paper-matched note (important for reproducibility):
- This script intentionally passes pandas Series into mean()/std() where appropriate,
  so that std() uses pandas default ddof=1 (matching typical "paper-used" scripts).
- Avoids .values for Z-score and photobleaching steps to prevent ddof=0 differences.

Default folder layout (recommended):
repo/
  photometry_control_vs_test.py
  control/   (put raw CSVs here)
  test/      (put raw CSVs here)
  output/    (auto-generated)

Run:
  python photometry_control_vs_test.py

Or specify folders:
  python photometry_control_vs_test.py --control ./my_control --test ./my_test --out ./my_output
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Defaults (can be overridden by CLI)
# -----------------------------
DEFAULT_BASELINE_INTERVAL_GLOBAL = (1500, 2100)  # 25–35 min baseline for Z-score
DEFAULT_BASELINE_PRE = (100, 600)                # photobleaching: pre interval
DEFAULT_BASELINE_POST_START = 500                # photobleaching: post interval start offset from end
DEFAULT_BASELINE_POST_END = 0                    # photobleaching: post interval end offset from end

DEFAULT_PLOT_LOWER_SEC = 2700    # 45 min
DEFAULT_PLOT_UPPER_SEC = 24300   # 6 h 45 min


# -----------------------------
# Correction / normalization functions
# -----------------------------
def correct_photobleaching(ts, ys, pre_interval, post_interval):
    """
    Linear photobleaching correction using means from two windows.

    IMPORTANT:
    - Designed to behave like typical paper-used implementations using pandas Series.
    - No fallback paths that change outputs silently.
    """
    pre_mask = (ts >= pre_interval[0]) & (ts <= pre_interval[1])
    post_mask = (ts >= post_interval[0]) & (ts <= post_interval[1])
    if pre_mask.sum() < 5 or post_mask.sum() < 5:
        raise ValueError("Baseline intervals too short for photobleaching correction.")

    pre_mean = ys[pre_mask].mean()
    post_mean = ys[post_mask].mean()
    slope = (post_mean - pre_mean) / (post_interval[1] - pre_interval[0])
    intercept = pre_mean - slope * pre_interval[0]
    return ys - (slope * ts + intercept)


def correct_motion(fluo465, fluo405):
    """
    Motion correction by linear regression (405 fitted to 465).

    NOTE:
    - Uses numpy lstsq; passing Series or ndarray yields equivalent results in practice.
    """
    if len(fluo405) < 10:
        raise ValueError("Not enough data points for motion correction.")
    A = np.vstack([np.asarray(fluo405), np.ones_like(np.asarray(fluo405))]).T
    coeffs, _, _, _ = np.linalg.lstsq(A, np.asarray(fluo465), rcond=None)
    fitted = A @ coeffs
    corrected = np.asarray(fluo465) - (fitted - np.mean(fitted))
    return np.asarray(fluo405), corrected


def transform_to_zscore(ts, ys, baseline_interval=(0, 60)):
    """
    Z-score normalization using a global baseline interval.

    CRITICAL (paper-matched):
    - If ys is a pandas Series, .std() uses ddof=1 by default (matching many paper scripts).
    - To preserve that behavior, call this with pandas Series (not .values).
    """
    mask = (ts >= baseline_interval[0]) & (ts <= baseline_interval[1])
    if mask.sum() < 5:
        raise ValueError("Baseline interval too short for Z-score.")

    mu = ys[mask].mean()
    sd = ys[mask].std()  # pandas default ddof=1 when ys is Series
    if sd == 0:
        raise ValueError("Cannot Z-score (std=0).")
    return (ys - mu) / sd


# -----------------------------
# IO helpers
# -----------------------------
def detect_header_line(filepath, max_lines=20):
    """
    Attempts to find the header line containing 'Time' and either 'gfp' or 'tomato' keywords.
    """
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f.readlines()[:max_lines]):
            low = line.lower()
            if ("time" in low) and (("gfp" in low) or ("tomato" in low)):
                return i
    return None


def load_and_standardize_csv(filepath):
    """
    Loads CSV and maps columns to: time, F-465 (GFP), AF-405 (tdTomato)
    Returns standardized df with columns ['time','F-465','AF-405'].
    """
    header_line = detect_header_line(filepath)
    if header_line is None:
        raise ValueError("Header line not detected (missing Time/GFP/Tomato keywords).")

    df = pd.read_csv(filepath, skiprows=header_line, low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    col_map = {}
    for col in df.columns:
        low = col.lower()
        if "time" in low:
            col_map[col] = "time"
        elif "gfp" in low:
            col_map[col] = "F-465"
        elif "tomato" in low:
            col_map[col] = "AF-405"

    df = df.rename(columns=col_map)

    required = ["time", "F-465", "AF-405"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Required columns missing: {missing}")

    df = df[required].dropna()
    return df


# -----------------------------
# ① Per-animal processing
# -----------------------------
def process_folder(input_folder, output_folder,
                   baseline_interval_global,
                   baseline_pre_interval,
                   baseline_post_start_offset,
                   baseline_post_end_offset,
                   nan_ratio_threshold=0.2):
    files = glob.glob(os.path.join(input_folder, "*.csv"))
    if not files:
        print(f"⚠ No CSV files found: {input_folder}")
        return 0

    os.makedirs(output_folder, exist_ok=True)

    n_saved = 0
    for fpath in files:
        fname = os.path.basename(fpath)
        print(f"\n📌 Processing: {fname}")

        try:
            df = load_and_standardize_csv(fpath)
        except Exception as e:
            print(f"⚠ Load/header error -> skip: {e}")
            continue

        if df.empty or len(df) < 10:
            print("⚠ Too few rows -> skip")
            continue

        max_time = df["time"].max()
        pre_interval = baseline_pre_interval
        post_interval = (max_time - baseline_post_start_offset, max_time - baseline_post_end_offset)

        # Photobleaching correction (paper-matched: pass Series, not .values)
        try:
            df["fluo465-pbc"] = correct_photobleaching(df["time"], df["F-465"], pre_interval, post_interval)
            df["fluo405-pbc"] = correct_photobleaching(df["time"], df["AF-405"], pre_interval, post_interval)
        except Exception as e:
            print(f"⚠ Photobleaching correction failed -> skip: {e}")
            continue

        nan_ratio = df[["fluo465-pbc", "fluo405-pbc"]].isna().mean().mean()
        if nan_ratio > nan_ratio_threshold:
            print(f"⚠ Too many NaNs after photobleaching correction ({nan_ratio:.2%}) -> skip")
            continue

        # Motion correction
        try:
            fluo405_maf, fluo465_mac = correct_motion(df["fluo465-pbc"], df["fluo405-pbc"])
            # Store as Series aligned to df index
            df["fluo405-maf"] = pd.Series(fluo405_maf, index=df.index)
            df["fluo465-mac"] = pd.Series(fluo465_mac, index=df.index)
        except Exception as e:
            print(f"⚠ Motion correction failed -> skip: {e}")
            continue

        # Z-score (paper-matched: pass Series, not .values; preserves ddof=1 behavior)
        try:
            df["fluo465-zsc"] = transform_to_zscore(df["time"], df["fluo465-mac"],
                                                   baseline_interval=baseline_interval_global)
        except Exception as e:
            print(f"⚠ Z-score failed -> skip: {e}")
            continue

        animal_name = os.path.splitext(fname)[0]
        out_csv = os.path.join(output_folder, f"{animal_name}-phmtry.csv")
        df.to_csv(out_csv, index=False)
        print(f"✅ Saved: {out_csv}")
        n_saved += 1

    return n_saved


# -----------------------------
# ② Group mean±SEM plot
# -----------------------------
def load_group_data(folder, plot_lower_sec, plot_upper_sec):
    files = glob.glob(os.path.join(folder, "*-phmtry.csv"))
    if not files:
        raise ValueError(f"No processed '*-phmtry.csv' files found: {folder}")

    all_df = []
    for fpath in files:
        try:
            df = pd.read_csv(fpath, usecols=["time", "fluo465-zsc"])
        except Exception as e:
            print(f"⚠ Read processed CSV failed -> skip: {fpath} ({e})")
            continue

        upper_limit = min(plot_upper_sec, df["time"].max())
        df = df[(df["time"] >= plot_lower_sec) & (df["time"] <= upper_limit)]
        if df.empty or len(df) < 10:
            print(f"⚠ Not enough data in plot window -> skip: {os.path.basename(fpath)}")
            continue

        # bin to 1-second and average
        df["time"] = df["time"].round().astype(int)
        df = df.groupby("time", as_index=False).mean(numeric_only=True)

        animal_name = os.path.basename(fpath).split("-")[0]
        df = df.rename(columns={"fluo465-zsc": animal_name})
        all_df.append(df)

    if not all_df:
        raise ValueError(f"No valid processed files in folder: {folder}")

    merged = all_df[0]
    for d in all_df[1:]:
        merged = pd.merge(merged, d, on="time", how="outer")

    merged = merged.drop_duplicates(subset=["time"]).sort_values("time").reset_index(drop=True)
    return merged


def plot_control_vs_test(control_data, test_data, out_svg):
    control_times = control_data["time"].values
    control_vals = control_data.drop(columns=["time"]).values
    control_mean = np.nanmean(control_vals, axis=1)
    control_sem = np.nanstd(control_vals, axis=1) / np.sqrt(control_vals.shape[1])

    test_times = test_data["time"].values
    test_vals = test_data.drop(columns=["time"]).values
    test_mean = np.nanmean(test_vals, axis=1)
    test_sem = np.nanstd(test_vals, axis=1) / np.sqrt(test_vals.shape[1])

    os.makedirs(os.path.dirname(out_svg), exist_ok=True)

    plt.figure(figsize=(12, 6))
    plt.title("Z-score Mean ± SEM (Control vs Test)", fontsize=14)
    plt.xlabel("Time (s)", fontsize=12)
    plt.ylabel("Z-score", fontsize=12)

    plt.plot(control_times, control_mean, label="Control", lw=2)
    plt.fill_between(control_times, control_mean - control_sem, control_mean + control_sem, alpha=0.3)

    plt.plot(test_times, test_mean, label="Test", lw=2)
    plt.fill_between(test_times, test_mean - test_sem, test_mean + test_sem, alpha=0.3)

    plt.legend()
    plt.tight_layout()
    plt.savefig(out_svg, format="svg")
    plt.close()


# -----------------------------
# Main
# -----------------------------
def main():
    repo_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Fiber photometry Control vs Test pipeline (GitHub-ready; paper-matched).")
    parser.add_argument("--control", default=os.path.join(repo_dir, "control"),
                        help="Folder containing raw Control CSV files (default: ./control)")
    parser.add_argument("--test", default=os.path.join(repo_dir, "test"),
                        help="Folder containing raw Test CSV files (default: ./test)")
    parser.add_argument("--out", default=os.path.join(repo_dir, "output"),
                        help="Output base folder (default: ./output)")

    parser.add_argument("--baseline_start", type=float, default=DEFAULT_BASELINE_INTERVAL_GLOBAL[0])
    parser.add_argument("--baseline_end", type=float, default=DEFAULT_BASELINE_INTERVAL_GLOBAL[1])

    parser.add_argument("--pre_start", type=float, default=DEFAULT_BASELINE_PRE[0])
    parser.add_argument("--pre_end", type=float, default=DEFAULT_BASELINE_PRE[1])
    parser.add_argument("--post_start_offset", type=float, default=DEFAULT_BASELINE_POST_START)
    parser.add_argument("--post_end_offset", type=float, default=DEFAULT_BASELINE_POST_END)

    parser.add_argument("--plot_lower", type=float, default=DEFAULT_PLOT_LOWER_SEC)
    parser.add_argument("--plot_upper", type=float, default=DEFAULT_PLOT_UPPER_SEC)

    args = parser.parse_args()

    control_out_dir = os.path.join(args.out, "Control")
    test_out_dir = os.path.join(args.out, "Test")
    result_dir = os.path.join(args.out, "Control_vs_Test")
    out_svg = os.path.join(result_dir, "Control_vs_Test.svg")

    os.makedirs(args.out, exist_ok=True)

    baseline_interval_global = (args.baseline_start, args.baseline_end)
    baseline_pre_interval = (args.pre_start, args.pre_end)

    print("\n=== Control group processing ===")
    n_control = process_folder(
        input_folder=args.control,
        output_folder=control_out_dir,
        baseline_interval_global=baseline_interval_global,
        baseline_pre_interval=baseline_pre_interval,
        baseline_post_start_offset=args.post_start_offset,
        baseline_post_end_offset=args.post_end_offset,
    )

    print("\n=== Test group processing ===")
    n_test = process_folder(
        input_folder=args.test,
        output_folder=test_out_dir,
        baseline_interval_global=baseline_interval_global,
        baseline_pre_interval=baseline_pre_interval,
        baseline_post_start_offset=args.post_start_offset,
        baseline_post_end_offset=args.post_end_offset,
    )

    if n_control == 0 or n_test == 0:
        print("\n⚠ Not enough processed files to plot (need at least 1 in each group).")
        return

    print("\n=== Group mean±SEM & plotting ===")
    try:
        control_data = load_group_data(control_out_dir, args.plot_lower, args.plot_upper)
        test_data = load_group_data(test_out_dir, args.plot_lower, args.plot_upper)
        plot_control_vs_test(control_data, test_data, out_svg)
        print(f"✅ Saved SVG: {out_svg}")
    except Exception as e:
        print(f"⚠ Plotting failed: {e}")


if __name__ == "__main__":
    main()
