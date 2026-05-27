import argparse
import json
from pathlib import Path


SEROTYPES = ["1", "5", "6B", "14", "18C", "19F", "23F"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Show per-serotype plasma/memory factors from annealing log."
    )
    parser.add_argument(
        "--log",
        default="misc/annealing_log.jsonl",
        help="Path to annealing jsonl log file.",
    )
    parser.add_argument(
        "--iter",
        type=int,
        default=None,
        dest="target_iter",
        help="Specific iteration to read. Default: last iteration.",
    )
    return parser.parse_args()


def _read_target_row(log_path: Path, target_iter):
    selected = None
    with log_path.open("r", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            if target_iter is None:
                selected = row
            elif row.get("iter") == target_iter:
                selected = row
                break
    return selected


def main():
    args = parse_args()
    log_path = Path(args.log)
    if not log_path.exists():
        print(f"log not found: {log_path}")
        return

    row = _read_target_row(log_path, args.target_iter)
    if row is None:
        if args.target_iter is None:
            print("log is empty")
        else:
            print(f"iteration {args.target_iter} not found")
        return

    params = row.get("params_candidate", {})
    plasma = params.get("plasma_production_factor", {})
    memory = params.get("memory_production_factor", {})

    print(f"iter: {row.get('iter')}")
    print(f"best_energy: {row.get('best_energy')}")
    print("serotype factors:")
    print("serotype, plasma_factor, memory_factor")
    for st in SEROTYPES:
        p = plasma.get(st, "NA") if isinstance(plasma, dict) else "NA"
        m = memory.get(st, "NA") if isinstance(memory, dict) else "NA"
        print(f"{st}, {p}, {m}")


if __name__ == "__main__":
    main()
