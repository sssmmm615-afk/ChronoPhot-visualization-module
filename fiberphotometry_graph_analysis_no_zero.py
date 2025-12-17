import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# === 設定 ===
sampling_rate = 2  # Hz
baseline_duration_sec = 30 * 60
total_duration_sec = 6 * 60 * 60
baseline_points = baseline_duration_sec * sampling_rate
total_points = total_duration_sec * sampling_rate
window_size = 301  # 退色補正用の移動中央値ウィンドウ

# === フォルダ設定 ===
desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
input_folder = os.path.join(desktop_path, 'Dric_CSV')
main_output_folder = os.path.join(desktop_path, 'Dric_6h_graph')
os.makedirs(main_output_folder, exist_ok=True)

# === グラフ描画関数（色付き）===
def plot_with_colors_and_save(x, y_data_dict, title, xlabel, ylabel, filename):
    plt.figure()
    for label, (y, color) in y_data_dict.items():
        plt.plot(x, y, label=label, color=color)
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename, format='svg')
    plt.close()

# === CSVファイル処理 ===
csv_files = glob.glob(os.path.join(input_folder, '*.csv'))

for file_path in csv_files:
    filename = os.path.basename(file_path)
    base_name = os.path.splitext(filename)[0]
    sample_folder = os.path.join(main_output_folder, base_name)
    os.makedirs(sample_folder, exist_ok=True)

    try:
        df = pd.read_csv(file_path, header=1)
        gfp_col = [col for col in df.columns if 'gfp' in col.lower()]
        tdt_col = [col for col in df.columns if 'tdtomato' in col.lower() or 'red' in col.lower() or 'tomato' in col.lower()]
        if not gfp_col or not tdt_col:
            print(f"⚠ カラムが見つかりません（{filename}）")
            continue

        gfp_all = df[gfp_col[0]].astype(float).values[:total_points]
        tdt_all = df[tdt_col[0]].astype(float).values[:total_points]

        # === 0分値（index 0）を除外 ===
        gfp_all = gfp_all[1:]
        tdt_all = tdt_all[1:]
        time_axis_min = np.arange(1, len(gfp_all) + 1) / (sampling_rate * 60)

        # 形状整形
        gfp = gfp_all.reshape(-1, 1)
        tdt = tdt_all.reshape(-1, 1)

        # 移動中央値による退色補正
        gfp_smooth = pd.Series(gfp.flatten()).rolling(window=window_size, center=True, min_periods=1).median().values
        gfp_detrended = gfp.flatten() - gfp_smooth

        # tdTomatoで回帰補正（ベースライン範囲のみで学習）
        model = LinearRegression()
        model.fit(tdt[:baseline_points], gfp_detrended[:baseline_points].reshape(-1, 1))
        gfp_fitted = model.predict(tdt)
        gfp_corrected = gfp_detrended.reshape(-1, 1) - gfp_fitted

        # Z-score
        baseline = gfp_corrected[:baseline_points]
        z_score = (gfp_corrected - np.mean(baseline)) / np.std(baseline)
        z_score = z_score.flatten()

        # ΔF/Fと正規化ΔF/F
        baseline_fluo = np.mean(gfp[:baseline_points])
        delta_f_over_f = (gfp.flatten() - baseline_fluo) / baseline_fluo
        delta_f_over_f_z = (delta_f_over_f - np.mean(delta_f_over_f[:baseline_points])) / np.std(delta_f_over_f[:baseline_points])

        # === グラフ出力 ===
        plot_with_colors_and_save(
            time_axis_min,
            {'GFP (raw)': (gfp.flatten(), 'green'), 'tdTomato': (tdt.flatten(), 'red')},
            f"{base_name} - Signal vs Control",
            'Time (min)', 'Fluorescence',
            os.path.join(sample_folder, "signal_vs_control.svg")
        )

        plot_with_colors_and_save(
            time_axis_min,
            {'GFP Detrended': (gfp_detrended, 'green'), 'Fitted GFP from tdTomato': (gfp_fitted.flatten(), 'red')},
            f"{base_name} - Fitted Control",
            'Time (min)', 'Signal',
            os.path.join(sample_folder, "fitted_control.svg")
        )

        plot_with_colors_and_save(
            time_axis_min,
            {'Delta F/F': (delta_f_over_f, 'green'), 'Normalized Delta F/F': (delta_f_over_f_z, 'blue')},
            f"{base_name} - Delta F/F",
            'Time (min)', 'Value',
            os.path.join(sample_folder, "deltaF.svg")
        )

        plot_with_colors_and_save(
            time_axis_min,
            {'Z-score': (z_score, 'black')},
            f"{base_name} - Z-score over 6h",
            'Time (min)', 'Z',
            os.path.join(sample_folder, "zscore.svg")
        )

        print(f"✅ グラフ出力完了（0分除外）: {base_name}")

    except Exception as e:
        print(f"❌ エラー（{filename}）: {e}")
