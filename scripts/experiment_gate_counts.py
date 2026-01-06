#!/usr/bin/env python3
import argparse
import csv
import json
import os
import subprocess
from dataclasses import dataclass
from statistics import mean


@dataclass
class Agg:
    cx: list
    t: list
    tdg: list
    s: list
    sdg: list
    x: list
    z: list

    def __init__(self):
        self.cx, self.t, self.tdg, self.s, self.sdg, self.x, self.z = [], [], [], [], [], [], []

    def add(self, d):
        self.cx.append(d["cx"])
        self.t.append(d["t"] + d["tdg"])
        self.s.append(d["s"] + d["sdg"])
        self.x.append(d["x"])
        self.z.append(d["z"])

    def avg(self):
        return {
            "cx": mean(self.cx) if self.cx else 0.0,
            "t": mean(self.t) if self.t else 0.0,
            "s": mean(self.s) if self.s else 0.0,
            "x": mean(self.x) if self.x else 0.0,
            "z": mean(self.z) if self.z else 0.0,
            "clifford_sxz": mean([s + x + z for s, x, z in zip(self.s, self.x, self.z)]) if self.s else 0.0,
        }


def parse_int_list(s: str):
    out = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", required=True, help="comma-separated n values, e.g. 3,4,5")
    ap.add_argument("--card", required=True, help="comma-separated cardinalities, e.g. 2,4,8")
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--eps", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=1, help="base seed")
    ap.add_argument("--bin", default="./build/tools/prepare_state")
    ap.add_argument("--out", default="", help="write CSV to this path")
    args = ap.parse_args()

    ns = parse_int_list(args.n)
    cs = parse_int_list(args.card)
    if not ns or not cs:
        raise SystemExit("empty --n or --card")

    if not os.path.exists(args.bin):
        raise SystemExit(f"missing binary: {args.bin} (build first)")

    rows = []
    for n in ns:
        for card in cs:
            agg_prep = Agg()
            agg_ct = Agg()
            for t in range(args.trials):
                seed = args.seed + 100000 * n + 1000 * card + t
                cmd = [args.bin, "-n", str(n), "-c", str(card), "-s", str(seed), "-e", str(args.eps), "--json"]
                out = subprocess.check_output(cmd, text=True).strip()
                j = json.loads(out)
                agg_prep.add(j["prep"])
                agg_ct.add(j["ct"])

            prep = agg_prep.avg()
            ct = agg_ct.avg()
            rows.append(
                {
                    "n": n,
                    "cardinality": card,
                    "trials": args.trials,
                    "eps": args.eps,
                    "prep_cx_avg": prep["cx"],
                    "prep_t_avg": prep["t"],
                    "prep_clifford_sxz_avg": prep["clifford_sxz"],
                    "ct_cx_avg": ct["cx"],
                    "ct_t_avg": ct["t"],
                    "ct_clifford_sxz_avg": ct["clifford_sxz"],
                }
            )

    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    for r in rows:
        print(
            f"n={r['n']} card={r['cardinality']} trials={r['trials']} eps={r['eps']}"
            f" | prep: CX={r['prep_cx_avg']:.2f} T={r['prep_t_avg']:.2f} SXZ={r['prep_clifford_sxz_avg']:.2f}"
            f" | ct: CX={r['ct_cx_avg']:.2f} T={r['ct_t_avg']:.2f} SXZ={r['ct_clifford_sxz_avg']:.2f}"
        )


if __name__ == "__main__":
    main()


