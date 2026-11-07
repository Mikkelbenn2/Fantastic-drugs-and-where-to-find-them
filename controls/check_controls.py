#!/usr/bin/env python3
"""Scores a finished control run against controls/expected_bgcs.tsv.
   python3 controls/check_controls.py RESULTS/_combined/summary [controls/expected_bgcs.tsv]
Positive controls: per tool, how many of the strain's MIBiG BGCs were found (a BGC of that tool
whose MIBiG_accession is the expected one). Negative controls: how many loci the tools call and
how many of them the pipeline would wrongly put forward as novel antibiotic candidates."""
import csv, os, sys
from collections import defaultdict

summ = sys.argv[1]
expected_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "expected_bgcs.tsv")
TOOLS = ["antismash", "gecco", "deepbgc", "bgc_prophet"]


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.isfile(path) else []


expected = defaultdict(set)
for r in rows(expected_file):
    if r["control"] == "positive":
        expected[r["isolate"]].add(r["mibig_accession"])
found = {t: defaultdict(set) for t in TOOLS}
for t in TOOLS:
    for r in rows(os.path.join(summ, "per_tool", f"{t}_bgcs.tsv")):
        found[t][r["isolate"]].add(r.get("MIBiG_accession", ""))
master = rows(os.path.join(summ, "master_summary.tsv"))
isolates = sorted({r["isolate"] for r in master} | set(expected))

print("POSITIVE CONTROLS: expected MIBiG BGCs recovered (KnownClusterBlast-style top hit)")
print("isolate\texpected\t" + "\t".join(TOOLS) + "\tany_tool\tmissed_by_all")
for iso in sorted(expected):
    exp = expected[iso]
    per = [len(exp & found[t][iso]) for t in TOOLS]
    anyt = exp & set().union(*(found[t][iso] for t in TOOLS))
    print(f"{iso}\t{len(exp)}\t" + "\t".join(map(str, per)) + f"\t{len(anyt)}\t{','.join(sorted(exp - anyt)) or '-'}")

print("\nALL CONTROLS: loci in master_summary.tsv")
print("isolate\tloci\tknown\trelated\tnovel_candidate\tnovel_antibiotic_class_score>=5\tbest_score")
for iso in isolates:
    m = [r for r in master if r["isolate"] == iso]
    nov = [r for r in m if r["novelty"] == "novel_candidate"]
    strong = [r for r in nov if int(r["priority_score"]) >= 5]
    best = max((int(r["priority_score"]) for r in m), default=0)
    count = lambda k: sum(r["novelty"] == k for r in m)
    print(f"{iso}\t{len(m)}\t{count('known')}\t{count('related')}\t{len(nov)}\t{len(strong)}\t{best}")
