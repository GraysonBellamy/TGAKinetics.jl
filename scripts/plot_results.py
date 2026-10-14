"""Plot TGAKinetics output: measured and model mass and mass-loss rate
against temperature, one colour per dataset.

    python plot_results.py results/              # fit_*.csv  -> results/fit.png
    python plot_results.py results/ --data       # data_*.csv -> results/data.png (from `check`)
    python plot_results.py results/simulation/   # sim_*.csv  -> sim.png

Needs numpy and matplotlib.
"""

import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np


def read(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder")
    ap.add_argument("--data", action="store_true", help="plot the processed data written by `check`")
    args = ap.parse_args()

    prefix = "data" if args.data else None
    if prefix is None:
        prefix = "fit" if glob.glob(os.path.join(args.folder, "fit_*.csv")) else "sim"
    natural = lambda f: [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", f)]
    files = sorted(glob.glob(os.path.join(args.folder, f"{prefix}_*.csv")), key=natural)
    if not files:
        raise SystemExit(f"no {prefix}_*.csv files in {args.folder}")

    fig, (ax_m, ax_r) = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True)
    for k, path in enumerate(files):
        name = os.path.basename(path)[len(prefix) + 1 : -4]
        d = read(path)
        c = f"C{k}"
        if prefix == "data":
            T = d["temperature_smoothed_K"]
            ax_m.plot(T, d["mass"], ".", ms=1.5, color=c, alpha=0.4)
            ax_m.plot(T, d["mass_smoothed"], "-", color=c, lw=1.2, label=name)
            ax_r.plot(T, d["mlr_per_s"], "-", color=c, lw=1.2, label=name)
            continue
        T = d["temperature_K"]
        if np.isfinite(d["mass_exp"]).any():
            fit = d["in_fit_range"] > 0
            ax_m.plot(T[fit], d["mass_exp"][fit], "o", ms=2.5, mfc="none", color=c, alpha=0.6)
            ax_r.plot(T[fit], d["mlr_exp_per_s"][fit], "o", ms=2.5, mfc="none", color=c, alpha=0.6)
        ax_m.plot(T, d["mass_model"], "-", color=c, lw=1.5, label=name)
        ax_r.plot(T, d["mlr_model_per_s"], "-", color=c, lw=1.5, label=name)

    ax_m.set_ylabel("Normalized mass (-)")
    ax_r.set_ylabel("Mass-loss rate (1/s)")
    ax_r.set_xlabel("Temperature (K)")
    ax_m.legend(frameon=False)
    if prefix != "data":
        ax_r.text(0.02, 0.95, "circles: measured, lines: model", transform=ax_r.transAxes, va="top", fontsize=8)
    fig.tight_layout()
    out = os.path.join(args.folder, f"{prefix}.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main()
