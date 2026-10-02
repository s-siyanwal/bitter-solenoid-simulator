"""Write results/emulation_compare.json and figures/fig6_emulation.png.

Does not touch results/results.json or the published optimum.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "package"))

from bittersim import BitterDesign
from bittersim.emulation import emulate


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = emulate(BitterDesign(), realizations=50, seed=1, homogeneity=True)
    results = os.path.join(root, "results")
    figures = os.path.join(root, "figures")
    if not os.path.isdir(results):
        os.makedirs(results)
    if not os.path.isdir(figures):
        os.makedirs(figures)
    path = os.path.join(results, "emulation_compare.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = ["B0_T", "V_V", "P_W", "T_hot_C", "ppm"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    xs = list(range(len(labels)))
    for i, key in enumerate(labels):
        p = out["percentiles"][key]
        nom = out["nominal_summary"][key]
        ax.errorbar([i], [p["p50"] / nom], yerr=[[(p["p50"] - p["p05"]) / nom], [(p["p95"] - p["p50"]) / nom]],
                    fmt="o", color="#1d3557", capsize=4)
    ax.axhline(1.0, color="#888", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylabel("value / continuum nominal")
    ax.set_title("Mesoscopic emulation, not molecular dynamics (seed 1, n=50)")
    fig.tight_layout()
    fig.savefig(os.path.join(figures, "fig6_emulation.png"), dpi=120)
    print(path)
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "update_report.py")
    if os.path.isfile(script):
        os.system("%s %s" % (sys.executable, script))


if __name__ == "__main__":
    main()
