"""
Reproduce Catoni et al. 2018 (NAR, DMRcaller paper) Figure 3 style on Muna's
SMA WGBS data. Two-panel figure per null model (loose + strict modes).

Inputs: parameter_benchmark_{label_swap,archie_scramble,stratified_scramble}.csv
Outputs: PDF + PNG into results/dmr_benchmark/plots_v2/
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

WORKSPACE = Path("results/dmr_benchmark")
OUT = WORKSPACE / "plots_v2"
OUT.mkdir(exist_ok=True, parents=True)

files = {
    "Label-swap":             "parameter_benchmark_label_swap.csv",
    "Read count permutation": "parameter_benchmark_archie_scramble.csv",
    "Stratified scramble":    "parameter_benchmark_stratified_scramble.csv",
}

frames = []
for label, fname in files.items():
    df = pd.read_csv(WORKSPACE / fname)
    df["null_model"] = label
    frames.append(df)
df = pd.concat(frames, ignore_index=True)

df["window_size"] = df["window_size"].astype(int)
df = df[~((df["method"] == "noise_filter") & (df["kernel"].isin(["uniform", "epanechnicov"])))]

method_label = {
    "bins": "DMRcaller-B",
    "neighbourhood": "DMRcaller-NB",
    "noise_filter": "DMRcaller-NF",
}
df["method_label"] = df["method"].map(method_label)

df["cov_real_mb"]  = df["n_real"]      * df["window_size"] / 1e6
df["cov_scram_mb"] = df["n_scrambled"] * df["window_size"] / 1e6
df["cov_diff_mb"]  = df["cov_real_mb"] - df["cov_scram_mb"]

COL = {
    "DMRcaller-B":  "#0F8C7E",
    "DMRcaller-NB": "#2E5DA8",
    "DMRcaller-NF": "#C2410C",
}
MARKER = {"DMRcaller-B": "s", "DMRcaller-NB": "^", "DMRcaller-NF": "D"}
WINDOWS = [100, 200, 300, 500, 1000, 2000]
COVERAGE_METHODS = ["DMRcaller-B", "DMRcaller-NF"]


def panel_A(ax, sub, title):
    methods_present = [m for m in COVERAGE_METHODS if m in sub["method_label"].unique()]
    for m in methods_present:
        d = sub[sub["method_label"] == m].sort_values("window_size")
        ax.plot(d["window_size"], d["cov_real_mb"],
                color=COL[m], marker=MARKER[m], markersize=7, markerfacecolor="white",
                markeredgewidth=1.6, linewidth=1.7, linestyle="-", zorder=3)
    for m in methods_present:
        d = sub[sub["method_label"] == m].sort_values("window_size")
        ax.plot(d["window_size"], d["cov_scram_mb"],
                color=COL[m], marker=MARKER[m], markersize=6, markerfacecolor="white",
                markeredgewidth=1.4, linewidth=1.4, linestyle="--", zorder=2, alpha=0.85)
    _format_axis(ax, "DMR genome coverage (Mb)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=10)


def panel_B(ax, sub, title):
    methods_present = [m for m in COVERAGE_METHODS if m in sub["method_label"].unique()]
    ax.axhline(0, linestyle=":", color="#888", linewidth=1, zorder=1)
    for m in methods_present:
        d = sub[sub["method_label"] == m].sort_values("window_size")
        ax.plot(d["window_size"], d["cov_diff_mb"],
                color=COL[m], marker=MARKER[m], markersize=7, markerfacecolor="white",
                markeredgewidth=1.6, linewidth=1.7, linestyle="-", zorder=3)
    _format_axis(ax, "Delta coverage (real minus scrambled, Mb)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=10)


def _format_axis(ax, ylabel):
    ax.set_xscale("log")
    ax.set_xticks(WINDOWS)
    ax.set_xticklabels([str(w) for w in WINDOWS])
    ax.set_xlabel("bin/window size (bp)", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(True, which="major", color="#dadada", linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#666")
    ax.spines["bottom"].set_color("#666")


def make_legend(ax_a, ax_b):
    handles_A = []
    for m in COVERAGE_METHODS:
        handles_A.append(Line2D([0], [0], color=COL[m], marker=MARKER[m],
                                markersize=7, markerfacecolor="white",
                                markeredgewidth=1.4, linewidth=1.6,
                                label=f"{m} (real)"))
    for m in COVERAGE_METHODS:
        handles_A.append(Line2D([0], [0], color=COL[m], marker=MARKER[m],
                                markersize=6, markerfacecolor="white",
                                markeredgewidth=1.2, linewidth=1.4, linestyle="--",
                                label=f"{m} (scrambled)"))
    ax_a.legend(handles=handles_A, loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False, fontsize=8.5, handlelength=2.4)

    handles_B = [Line2D([0], [0], color=COL[m], marker=MARKER[m],
                        markersize=7, markerfacecolor="white",
                        markeredgewidth=1.4, linewidth=1.6,
                        label=f"{m} (delta)") for m in COVERAGE_METHODS]
    ax_b.legend(handles=handles_B, loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                frameon=False, fontsize=8.5, handlelength=2.4)


def render_panels(null_model_label, mode, fname_tag):
    sub = df[(df["null_model"] == null_model_label) & (df["mode"] == mode)].copy()
    if sub.empty:
        return None

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.7), constrained_layout=False)
    plt.subplots_adjust(left=0.07, right=0.80, top=0.84, bottom=0.20, wspace=0.65)

    titleA = (f"A   CpG DMRs - ASO_VPA vs ASO_CTRL  ({mode})\n"
              f"null model: {null_model_label}")
    titleB = (f"B   CpG DMRs difference (real minus scrambled)  ({mode})\n"
              f"null model: {null_model_label}")
    panel_A(axA, sub, titleA)
    panel_B(axB, sub, titleB)
    make_legend(axA, axB)

    fig.text(0.07, 0.02,
             "DMRcaller-NB excluded from coverage panels: count is constant w.r.t. window size. Shown in ratio plot.",
             fontsize=8, color="#555", ha="left")

    pdf_path = OUT / f"radu_panels_{fname_tag}_{mode}.pdf"
    png_path = OUT / f"radu_panels_{fname_tag}_{mode}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return pdf_path


def render_ratio(mode):
    sub = df[df["mode"] == mode].copy()
    collapsed = (sub["n_real"] == 0) & (sub["n_scrambled"] == 0)
    sub.loc[collapsed, "ratio"] = 0.0
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["ratio"])

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), constrained_layout=False, sharey=False)
    plt.subplots_adjust(left=0.07, right=0.93, top=0.82, bottom=0.18, wspace=0.38)
    null_order = ["Label-swap", "Read count permutation", "Stratified scramble"]
    for i, null in enumerate(null_order):
        ax = axes[i]
        s = sub[sub["null_model"] == null]
        for m in list(COL):
            d = s[s["method_label"] == m].sort_values("window_size")
            if d.empty:
                continue
            ax.plot(d["window_size"], d["ratio"], color=COL[m],
                    marker=MARKER[m], markersize=7, markerfacecolor="white",
                    markeredgewidth=1.5, linewidth=1.6, label=m)
        ax.axhline(1.0, linestyle=":", color="#888", linewidth=1)
        ax.set_xscale("log")
        ax.set_xticks(WINDOWS)
        ax.set_xticklabels([str(w) for w in WINDOWS], fontsize=9)
        ax.set_xlabel("bin/window size (bp)", fontsize=10)
        ax.set_ylabel("Signal/noise ratio (real / scrambled)" if i == 0 else "")
        ax.set_title(null, loc="left", fontsize=11, fontweight="bold", pad=8)
        ax.grid(True, which="major", color="#dadada", linewidth=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if i == 2:
            ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    fig.suptitle(f"Signal/noise across null models   -   {mode} threshold",
                 fontsize=12, fontweight="bold", x=0.07, ha="left", y=0.96)
    pdf_path = OUT / f"ratio_all_{mode}.pdf"
    png_path = OUT / f"ratio_all_{mode}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return pdf_path


def render_threshold_comparison(null_model_label, fname_tag):
    sub = df[df["null_model"] == null_model_label].copy()
    collapsed = (sub["n_real"] == 0) & (sub["n_scrambled"] == 0)
    sub.loc[collapsed, "ratio"] = 0.0
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=["ratio"])
    if sub.empty:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=False)
    plt.subplots_adjust(left=0.08, right=0.94, top=0.82, bottom=0.18, wspace=0.32)
    for i, mode in enumerate(("loose", "strict")):
        ax = axes[i]
        s = sub[sub["mode"] == mode]
        for m in list(COL):
            d = s[s["method_label"] == m].sort_values("window_size")
            if d.empty:
                continue
            ax.plot(d["window_size"], d["ratio"], color=COL[m],
                    marker=MARKER[m], markersize=7, markerfacecolor="white",
                    markeredgewidth=1.5, linewidth=1.6, label=m)
        ax.axhline(1.0, linestyle=":", color="#888", linewidth=1)
        ax.set_xscale("log")
        ax.set_xticks(WINDOWS)
        ax.set_xticklabels([str(w) for w in WINDOWS], fontsize=9)
        ax.set_xlabel("bin/window size (bp)", fontsize=10)
        ax.set_ylabel("Signal/noise ratio (real / scrambled)" if i == 0 else "")
        ax.set_title(f"{mode} threshold", loc="left", fontsize=11, fontweight="bold", pad=8)
        ax.grid(True, which="major", color="#dadada", linewidth=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if i == 1:
            ax.legend(loc="upper left", fontsize=8.5, frameon=False)
    fig.suptitle(f"Threshold sensitivity   -   {null_model_label}",
                 fontsize=12, fontweight="bold", x=0.08, ha="left", y=0.96)
    pdf_path = OUT / f"threshold_comparison_{fname_tag}.pdf"
    png_path = OUT / f"threshold_comparison_{fname_tag}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    return pdf_path


def main():
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["axes.unicode_minus"] = False

    out = []
    for null in files:
        tag = null.lower().replace(" ", "_").replace("-", "_")
        for mode in ("loose", "strict"):
            r = render_panels(null, mode, tag)
            if r:
                out.append(r)
        r = render_threshold_comparison(null, tag)
        if r:
            out.append(r)
    for mode in ("loose", "strict"):
        r = render_ratio(mode)
        if r:
            out.append(r)
    print(f"\nGenerated {len(out)} files in {OUT}")
    for o in out:
        print("  -", o.name)


if __name__ == "__main__":
    main()
