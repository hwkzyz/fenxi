from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
FIG_DIR = OUT / "figures_step07j_bundle_ablation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

DIRECT_TREND = OUT / "Step07J_NestedStaticWarp_VPFullWave_Trend_20250527_B1_S136_main_gaptilt.csv"
NODIRECT_TREND = OUT / "Step07J_NestedStaticWarp_VPFullWave_Trend_20250527_B1_S136_main_gaptilt_nodirect.csv"


def load_case(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    keep = pd.DataFrame(
        {
            "window": df["window_id"].astype(int),
            "case": label,
            "points": df["phase_safe_point_count"].astype(float),
            "rmse_mV": df["tilt_rmse_mV"].astype(float),
            "amplitude_mm": df["tilt_amplitude_mm"].astype(float),
            "EO": df["tilt_EO"].astype(int),
        }
    )
    return keep


def main() -> None:
    direct = load_case(DIRECT_TREND, "Direct bundle")
    nodirect = load_case(NODIRECT_TREND, "No direct bundle")
    data = pd.concat([direct, nodirect], ignore_index=True)

    wide = direct.merge(
        nodirect,
        on="window",
        suffixes=("_direct", "_nodirect"),
    )
    wide["point_ratio_direct_to_nodirect"] = wide["points_direct"] / wide["points_nodirect"]
    wide["rmse_reduction_mV"] = wide["rmse_mV_nodirect"] - wide["rmse_mV_direct"]
    wide["rmse_reduction_percent"] = 100 * wide["rmse_reduction_mV"] / wide["rmse_mV_nodirect"]
    compare_csv = OUT / "Step07J_BundleAblation_WindowCompare_20250527_B1_S136.csv"
    wide.to_csv(compare_csv, index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    colors = {
        "Direct bundle": "#0072B2",
        "No direct bundle": "#D55E00",
    }
    markers = {
        "Direct bundle": "o",
        "No direct bundle": "s",
    }

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(7.0, 5.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.12},
        constrained_layout=False,
    )

    for label, group in data.groupby("case", sort=False):
        axes[0].plot(
            group["window"],
            group["points"] / 1000,
            marker=markers[label],
            markersize=4.5,
            linewidth=1.6,
            color=colors[label],
            label=label,
        )
        axes[1].plot(
            group["window"],
            group["rmse_mV"],
            marker=markers[label],
            markersize=4.5,
            linewidth=1.6,
            color=colors[label],
            label=label,
        )

    axes[0].set_ylabel("Valid points per window (x10^3)")
    axes[1].set_ylabel("Weighted RMSE (mV)")
    axes[1].set_xlabel("Sliding window index")

    axes[0].set_ylim(0, max(data["points"]) / 1000 * 1.12)
    axes[1].set_ylim(0, max(data["rmse_mV"]) * 1.18)
    axes[1].set_xticks(np.arange(1, int(data["window"].max()) + 1, 1))

    for ax in axes:
        ax.grid(axis="y", color="#D0D0D0", linewidth=0.55, alpha=0.8)
        ax.set_xlim(0.5, data["window"].max() + 0.5)

    axes[0].legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.18))

    mean_direct_pts = direct["points"].mean()
    mean_nodirect_pts = nodirect["points"].mean()
    mean_direct_rmse = direct["rmse_mV"].mean()
    mean_nodirect_rmse = nodirect["rmse_mV"].mean()
    ratio = mean_direct_pts / mean_nodirect_pts
    rmse_drop = mean_nodirect_rmse - mean_direct_rmse

    text = (
        f"Mean valid points: {mean_direct_pts:,.0f} vs {mean_nodirect_pts:,.0f} "
        f"({ratio:.2f}x)\n"
        f"Mean RMSE: {mean_direct_rmse:.2f} vs {mean_nodirect_rmse:.2f} mV "
        f"({rmse_drop:.2f} mV lower)"
    )
    axes[0].text(
        0.02,
        0.66,
        text,
        transform=axes[0].transAxes,
        va="top",
        ha="left",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#B0B0B0"},
    )

    fig.suptitle(
        "Step07J bundle-source ablation, 20250527 B1 S136",
        y=0.985,
        fontsize=10,
        fontweight="bold",
    )
    fig.align_ylabels(axes)
    fig.subplots_adjust(top=0.89, left=0.12, right=0.98, bottom=0.10)

    base = FIG_DIR / "Step07J_BundleAblation_ValidPoints_RMSE_20250527_B1_S136"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(base.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)

    print(f"Saved figure base: {base}")
    print(f"Saved compare CSV: {compare_csv}")
    print(
        "Mean points direct/no-direct: "
        f"{mean_direct_pts:.1f}/{mean_nodirect_pts:.1f}; "
        f"mean RMSE {mean_direct_rmse:.3f}/{mean_nodirect_rmse:.3f} mV"
    )


if __name__ == "__main__":
    main()
