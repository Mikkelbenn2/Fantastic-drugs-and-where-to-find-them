#!/usr/bin/env python3
"""Scores a finished control run against setup/controls/expected_bgcs.tsv and setup/controls/control_panel.tsv.
   python3 setup/controls/check_controls.py RESULTS/_combined/summary [OUT_DIR]
Writes three tables to OUT_DIR (default: the summary folder) and prints them:
   controls_recovery.tsv       positives: per tool, how many of the strain's MIBiG BGCs were found
   controls_expected_bgcs.tsv  every expected MIBiG BGC and which tools found it
   controls_loci.tsv           every control: loci called, novelty labels, best score, PASS/CHECK
A known BGC counts as found when it is among a predicted BGC's MIBiG matches (MIBiG_all_hits) with
>= MIN_SIM % similarity, not only the top hit: neighbouring clusters merged into one region count too."""
import csv, os, re, sys
from collections import defaultdict

summ = sys.argv[1]
out_dir = sys.argv[2] if len(sys.argv) > 2 else summ
here = os.path.dirname(os.path.abspath(__file__))
TOOLS = ["antismash", "gecco", "deepbgc", "bgc_prophet"]
MIN_SIM = 30   # % of the MIBiG cluster's genes hit (= the pipeline's "related" threshold)


def acc(x):
    """MIBiG accession without version: antiSMASH/MIBiG 4 data may write BGC0000176.5 for BGC0000176"""
    return re.sub(r"\.\d+$", "", (x or "").strip())


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.isfile(path) else []


def write(name, header, table):
    path = os.path.join(out_dir, name)
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in table:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"\n== {path}")
    print("\t".join(header))
    for r in table:
        print("\t".join(map(str, r)))


master_file = os.path.join(summ, "master_summary.tsv")
if not os.path.isfile(master_file):
    sys.exit(f"ERROR: {master_file} not found. The pipeline's combined/summary step has not finished: "
             "check slurm_logs/combined_*.out and _combined/logs/summary.log in that results folder")
os.makedirs(out_dir, exist_ok=True)

panel = {r["isolate"]: r for r in rows(os.path.join(here, "control_panel.tsv"))}
expected = defaultdict(list)
for r in rows(os.path.join(here, "expected_bgcs.tsv")):
    if r["control"] == "positive":
        r["mibig_accession"] = acc(r["mibig_accession"])
        expected[r["isolate"]].append(r)
found = {t: defaultdict(set) for t in TOOLS}
for t in TOOLS:
    for r in rows(os.path.join(summ, "per_tool", f"{t}_bgcs.tsv")):
        found[t][r["isolate"]].add(acc(r.get("MIBiG_accession")))
        for h in (r.get("MIBiG_all_hits") or "").split(";"):     # older results have no such column
            a, _, sim = h.rpartition(":")
            if a and sim.replace(".", "", 1).isdigit() and float(sim) >= MIN_SIM:
                found[t][r["isolate"]].add(acc(a))
master = rows(master_file)
# genome_summary.tsv lists every isolate the pipeline finished, also those with no BGC at all
in_run = {r["isolate"] for r in rows(os.path.join(summ, "genome_summary.tsv"))} \
    | {r["isolate"] for r in master} | {i for t in TOOLS for i in found[t]}

recovery, per_bgc = [], []
for iso in sorted(expected):
    exp = {e["mibig_accession"] for e in expected[iso]}
    hits = {t: exp & found[t][iso] for t in TOOLS}
    anyt = set().union(*hits.values())
    recovery.append([iso, iso in in_run, len(exp)] + [len(hits[t]) for t in TOOLS]
                    + [len(anyt), f"{100 * len(anyt) / len(exp):.0f}", ",".join(sorted(exp - anyt)) or "-"])
    for e in expected[iso]:
        by = [t for t in TOOLS if e["mibig_accession"] in hits[t]]
        per_bgc.append([iso, e["mibig_accession"], e["compound"], e["mibig_class"], e["bioactivity"],
                        ",".join(by) or "none"])
write("controls_recovery.tsv",
      ["isolate", "in_run", "expected"] + TOOLS + ["any_tool", "any_tool_pct", "missed_by_all"], recovery)
write("controls_expected_bgcs.tsv",
      ["isolate", "mibig_accession", "compound", "mibig_class", "bioactivity", "found_by"], per_bgc)

loci = []
for iso in sorted(set(panel) | {r["isolate"] for r in master}):
    m = [r for r in master if r["isolate"] == iso]
    count = lambda k: sum(r["novelty"] == k for r in m)
    strong = sum(r["novelty"] == "novel_candidate" and int(r["priority_score"]) >= 5 for r in m)
    best = max((int(r["priority_score"]) for r in m), default=0)
    role = panel.get(iso, {}).get("control", "not_a_control")
    if iso not in in_run:
        verdict = "NOT_IN_RUN"
    elif role == "negative":            # a negative passes when nothing is put forward as a novel antibiotic
        verdict = "PASS" if strong == 0 else "CHECK"
    elif role == "positive":            # a positive passes when its known BGCs are not called novel
        verdict = "PASS" if count("known") + count("related") > 0 else "CHECK"
    else:
        verdict = "-"
    loci.append([iso, role, len(m), count("known"), count("related"), count("novel_candidate"), count("novel_other"), count("low_support"),
                 strong, best, verdict])
write("controls_loci.tsv", ["isolate", "control", "loci", "known", "related", "novel_candidate", "novel_other", "low_support",
                            "novel_candidate_score>=5", "best_score", "verdict"], loci)
