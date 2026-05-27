import argparse
import copy
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SIM_DIR = ROOT / "sim"
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

from config import SIMULATION_PARAMS
from main import simulate_batch
from simulated_annealing import SEROTYPES, chamfer_distance


def parse_args():
    parser = argparse.ArgumentParser(description="Generate visual comparisons for new simulator output.")
    parser.add_argument("--replicas", type=int, default=48, help="Number of simulated subjects.")
    parser.add_argument("--base-seed", type=int, default=12345, help="Base seed for reproducibility.")
    parser.add_argument("--duration-days", type=int, default=None, help="Override duration_days.")
    parser.add_argument(
        "--dataset",
        default="sim/data/train_VCN7-Tf_fit.csv",
        help="Clinical dataset CSV path.",
    )
    parser.add_argument(
        "--out-dir",
        default="images/new_model",
        help="Output directory for generated figures.",
    )
    parser.add_argument(
        "--params-json",
        default=None,
        help="Optional JSON file with simulation params to override SIMULATION_PARAMS.",
    )
    return parser.parse_args()


def _series_to_points(series):
    return np.array([(row[0], row[1]) for row in series], dtype=float)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = copy.deepcopy(SIMULATION_PARAMS)
    if args.params_json is not None:
        with open(args.params_json, "r", encoding="utf-8") as f:
            loaded_params = json.load(f)
        params.update(loaded_params)
    if args.duration_days is not None:
        params["duration_days"] = args.duration_days

    clinical = pd.read_csv(args.dataset)
    sim_summary = simulate_batch(params, n_replicas=args.replicas, base_seed=args.base_seed)
    runs = sim_summary["runs"]

    metrics = {}
    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    axes = axes.flatten()

    for idx, (serotype, cols) in enumerate(SEROTYPES.items()):
        ax = axes[idx]
        real_points = np.column_stack((clinical[cols[0]].values, clinical[cols[1]].values))
        sim_points = np.array([(run[serotype][0], run[serotype][1]) for run in runs], dtype=float)

        real_centroid = np.mean(real_points, axis=0)
        sim_centroid = np.mean(sim_points, axis=0)

        chamfer = chamfer_distance(real_points, sim_points)
        centroid_l2 = float(np.linalg.norm(real_centroid - sim_centroid))
        metrics[serotype] = {
            "chamfer": float(chamfer),
            "centroid_l2": centroid_l2,
            "real_centroid": real_centroid.tolist(),
            "sim_centroid": sim_centroid.tolist(),
        }

        ax.scatter(real_points[:, 0], real_points[:, 1], color="orange", alpha=0.75, label="real")
        ax.scatter(sim_points[:, 0], sim_points[:, 1], color="royalblue", alpha=0.65, label="sim")
        ax.scatter(real_centroid[0], real_centroid[1], color="crimson", marker="X", s=90, label="real centroid")
        ax.scatter(sim_centroid[0], sim_centroid[1], color="navy", marker="X", s=90, label="sim centroid")
        ax.set_title(f"sp{serotype} | Chamfer={chamfer:.3f} | C-L2={centroid_l2:.3f}")
        ax.set_xlabel("pre")
        ax.set_ylabel("post")
        ax.grid(alpha=0.25)

    axes[7].axis("off")
    axes[8].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Clinical vs New Simulator Output by Serotype", fontsize=16)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    grid_path = out_dir / "serotypes_grid.png"
    fig.savefig(grid_path, dpi=180)
    plt.close(fig)

    # one image per serotype
    for serotype, cols in SEROTYPES.items():
        fig_s, ax_s = plt.subplots(figsize=(6, 5))
        real_points = np.column_stack((clinical[cols[0]].values, clinical[cols[1]].values))
        sim_points = np.array([(run[serotype][0], run[serotype][1]) for run in runs], dtype=float)
        ax_s.scatter(real_points[:, 0], real_points[:, 1], color="orange", alpha=0.75, label="real")
        ax_s.scatter(sim_points[:, 0], sim_points[:, 1], color="royalblue", alpha=0.65, label="sim")
        ax_s.set_title(
            f"sp{serotype} | Chamfer={metrics[serotype]['chamfer']:.3f} | "
            f"C-L2={metrics[serotype]['centroid_l2']:.3f}"
        )
        ax_s.set_xlabel("pre")
        ax_s.set_ylabel("post")
        ax_s.grid(alpha=0.25)
        ax_s.legend()
        fig_s.tight_layout()
        fig_s.savefig(out_dir / f"sp{serotype}.png", dpi=180)
        plt.close(fig_s)

    metrics_path = out_dir / "metrics.json"
    with metrics_path.open("w", encoding="utf-8") as f:
        json.dump(
            {
                "replicas": args.replicas,
                "base_seed": args.base_seed,
                "duration_days": params["duration_days"],
                "metrics": metrics,
            },
            f,
            indent=2,
        )

    print(f"Saved grid: {grid_path}")
    print(f"Saved metrics: {metrics_path}")
    for st, m in metrics.items():
        print(f"sp{st}: chamfer={m['chamfer']:.4f}, centroid_l2={m['centroid_l2']:.4f}")


if __name__ == "__main__":
    main()
