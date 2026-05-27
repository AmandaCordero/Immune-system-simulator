import argparse
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Show serotype objective scores from annealing logs.")
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
        help="Specific iteration to show. If omitted, shows the last iteration.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    log_path = Path(args.log)
    if not log_path.exists():
        print(f"log not found: {log_path}")
        return

    selected = None
    with log_path.open("r", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            if args.target_iter is None:
                selected = row
            elif row.get("iter") == args.target_iter:
                selected = row
                break

    if selected is None:
        if args.target_iter is None:
            print("log is empty")
        else:
            print(f"iteration {args.target_iter} not found in log")
        return

    print("iter:", selected["iter"])
    print("temp:", selected["temp"])
    print("best_energy:", selected["best_energy"])
    print("serotype scores:")
    for serotype, details in sorted(selected["serotype_scores"].items()):
        print(
            f"  {serotype}: score={details['score']:.6f}, "
            f"chamfer={details['chamfer']:.6f}, "
            f"landmark_l2={details['landmark_l2']:.6f}, "
            f"weighted_score={details['weighted_score']:.6f}, "
            f"weight={details['weight']:.2f}"
        )


if __name__ == "__main__":
    main()
