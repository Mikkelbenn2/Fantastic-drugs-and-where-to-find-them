#!/usr/bin/env python3
"""Performance of the pipeline on the controls: confusion matrices, sensitivity, specificity, accuracy, ROC/AUC.
   python3 controls/evaluate_controls.py RESULTS_CONTROLS/_combined/summary [OUT_DIR]
Run after the pipeline (and normally after check_controls.py). Needs only python3.

Two questions, each with its own truth:
1. DETECTION: are the characterised BGCs found?  Truth = controls/expected_bgcs.tsv (MIBiG BGCs of the strain).
   Unit = expected BGC. Found = among a predicted BGC's MIBiG matches with >= MIN_SIM % (same rule as
   check_controls.py). Gives sensitivity (recall) per tool. Over-calling is shown as loci per Mb in the negatives
   (their real BGCs, e.g. enterobactin, are unknown in number, so a true specificity cannot be computed here).
2. ANTIBIOTIC RANKING: does the pipeline's evidence point to antibiotic BGCs?  Unit = BGC locus of a control.
   Positive = locus matching an expected BGC with antimicrobial = yes (named antibiotic in the control paper, or
   MIBiG bioactivity antibacterial/antifungal). Negative = every locus of a negative control, and loci matching only
   non-antimicrobial expected BGCs (siderophores, pigments, terpenes ...). Loci of positive controls that match no
   expected BGC are left out (cryptic BGCs: truth unknown). Antimicrobial expected BGCs missed by every tool count
   as positives with the lowest score (false negatives).
   Score = "novelty-blind" priority score = priority_score minus the novelty points (novel_candidate 3, related 1):
   in the controls everything is known, so this asks "would the evidence have ranked it high if it were new?".
   Predicted positive = the pipeline's own rule for a strong novel candidate, without the novelty part:
   supported (antiSMASH, or >= 2 tools) AND antibiotic-relevant type AND blind score >= BLIND_CUTOFF
   (2: a novel candidate needs 3 + 2 = score 5 to be "strong"). Older results without the column
   antibiotic_type fall back to the blind score alone.
3. GENOME LEVEL: positive vs negative strains; predicted producer = >= 1 novel_candidate with score >= 5
   (the pipeline's own output). Only 12 genomes, so read it as a sanity check.
Writes eval_confusion.tsv, eval_metrics.tsv, eval_auc.tsv, eval_roc.tsv, eval_detection.tsv, eval_loci_truth.tsv."""
import csv, math, os, re, sys
from collections import defaultdict

summ = sys.argv[1]
out_dir = sys.argv[2] if len(sys.argv) > 2 else summ
here = os.path.dirname(os.path.abspath(__file__))
TOOLS = ["antismash", "gecco", "deepbgc", "bgc_prophet"]
MIN_SIM = 30          # same as check_controls.py
BLIND_CUTOFF = 2      # novelty-blind score at which a NEW locus would reach priority score 5
MIN_TOOLS = 2         # the pipeline's SUPPORT_MIN_TOOLS (read from settings.tsv when present)
NOVELTY_PTS = {"novel_candidate": 3, "related": 1}


def rows(*p):
    p = os.path.join(*p)
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.isfile(p) else []


def acc(x):
    return re.sub(r"\.\d+$", "", (x or "").strip())


def num(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def wilson(k, n, z=1.96):
    """95 % confidence interval of a proportion (Wilson score), good for small n"""
    if n == 0:
        return (None, None)
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def auc(pos, neg):
    """ROC AUC = probability that a random positive scores higher than a random negative (ties = 1/2),
       with the Hanley & McNeil (1982) standard error for a 95 % interval"""
    if not pos or not neg:
        return None, None, None
    a = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    se = math.sqrt((a * (1 - a) + (len(pos) - 1) * (q1 - a * a) + (len(neg) - 1) * (q2 - a * a)) / (len(pos) * len(neg)))
    return a, max(0.0, a - 1.96 * se), min(1.0, a + 1.96 * se)


def roc(pos, neg):
    pts = [(0.0, 0.0, "inf")]
    for t in sorted(set(pos + neg), reverse=True):
        pts.append((sum(n >= t for n in neg) / len(neg), sum(p >= t for p in pos) / len(pos), t))
    return pts


def metrics(tp, fp, fn, tn):
    def r(k, n):
        lo, hi = wilson(k, n)
        return (k / n if n else None, lo, hi)
    n = tp + fp + fn + tn
    sens, spec, prec = r(tp, tp + fn), r(tn, tn + fp), r(tp, tp + fp)
    m = {"sensitivity (recall, TPR)": sens, "specificity (TNR)": spec, "precision (PPV)": prec,
         "negative predictive value": r(tn, tn + fn), "accuracy": r(tp + tn, n)}
    if sens[0] is not None and spec[0] is not None:
        m["balanced accuracy"] = ((sens[0] + spec[0]) / 2, None, None)
    if prec[0] and sens[0]:
        m["F1 score"] = (2 * prec[0] * sens[0] / (prec[0] + sens[0]), None, None)
    d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    m["Matthews correlation (MCC)"] = ((tp * tn - fp * fn) / d if d else None, None, None)
    return m


def write(name, header, table):
    path = os.path.join(out_dir, name)
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for t in table:
            f.write("\t".join("NA" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v)) for v in t) + "\n")
    return path


master = rows(summ, "master_summary.tsv")
MIN_TOOLS = int(next((r["value"] for r in rows(summ, "settings.tsv") if r["setting"] == "min_tools"), MIN_TOOLS))
if not master:
    sys.exit(f"ERROR: no master_summary.tsv in {summ}")
os.makedirs(out_dir, exist_ok=True)
panel = {r["isolate"]: r["control"] for r in rows(here, "control_panel.tsv")}
expected = defaultdict(dict)                  # isolate -> {accession: antimicrobial yes/no}
for r in rows(here, "expected_bgcs.tsv"):
    if r["control"] == "positive":
        expected[r["isolate"]][acc(r["mibig_accession"])] = r.get("antimicrobial", "no") == "yes"
genome_mb = {r["isolate"]: num(r.get("genome_length"), 0) / 1e6 for r in rows(summ, "genome_summary.tsv")}
in_run = set(genome_mb) | {r["isolate"] for r in master}

# MIBiG matches of every tool call (>= MIN_SIM %), per locus and per tool
locus_hits, tool_hits, tool_loci = defaultdict(set), {t: defaultdict(set) for t in TOOLS}, defaultdict(int)
for t in TOOLS:
    for r in rows(summ, "per_tool", f"{t}_bgcs.tsv"):
        iso = r["isolate"]
        tool_loci[(t, iso)] += 1
        hits = {acc(r.get("MIBiG_accession"))} if num(r.get("MIBiG_similarity_pct"), 0) >= MIN_SIM else set()
        for h in (r.get("MIBiG_all_hits") or "").split(";"):
            a, _, s = h.rpartition(":")
            if a and num(s, 0) >= MIN_SIM:
                hits.add(acc(a))
        hits &= set(expected.get(iso, {}))
        tool_hits[t][iso] |= hits
        if r.get("locus_id"):
            locus_hits[r["locus_id"]] |= hits

# ---- 1. detection ----
det = []
for t in TOOLS + ["any_tool"]:
    k = n = 0
    for iso, exp in expected.items():
        if iso not in in_run:
            continue
        found = set().union(*(tool_hits[x][iso] for x in TOOLS)) if t == "any_tool" else tool_hits[t][iso]
        k += len(found & set(exp)); n += len(exp)
    lo, hi = wilson(k, n)
    neg = [i for i in in_run if panel.get(i) == "negative"]
    mb = sum(genome_mb.get(i, 0) for i in neg)
    calls = sum(tool_loci[(x, i)] for x in (TOOLS if t == "any_tool" else [t]) for i in neg) if t != "any_tool" else \
        sum(1 for r in master if r["isolate"] in neg)
    det.append([t, n, k, k / n if n else None, lo, hi, calls, calls / mb if mb else None])

# ---- 2. antibiotic ranking, one row per locus ----
truth = []
matched = defaultdict(set)
for r in master:
    iso, role = r["isolate"], panel.get(r["isolate"])
    hits = locus_hits.get(r["locus_id"], set())
    matched[iso] |= hits
    if role == "negative":
        label = "negative"
    elif role == "positive" and hits:
        label = "positive" if any(expected[iso][h] for h in hits) else "negative"
    else:
        continue                                  # cryptic locus in a producer: truth unknown
    blind = num(r["priority_score"], 0) - NOVELTY_PTS.get(r["novelty"], 0)
    truth.append({"isolate": iso, "locus_id": r["locus_id"], "truth": label, "matched_expected": ",".join(sorted(hits)) or "-",
                  "novelty": r["novelty"], "priority_score": r["priority_score"], "blind_score": blind,
                  "n_tools": num(r.get("n_tools"), 1), "antibiotic_type": 1 if r.get("antibiotic_type") == "yes" else 0,
                  "supported": "antismash" in r.get("tools", "") or num(r.get("n_tools"), 1) >= MIN_TOOLS,
                  "deepbgc_activity": num(r.get("deepbgc_max_antibacterial_or_antifungal"), 0)})
n_cryptic = sum(1 for r in master if panel.get(r["isolate"]) == "positive" and not locus_hits.get(r["locus_id"]))
for iso, exp in expected.items():                 # antimicrobial BGCs no tool found = false negatives
    if iso in in_run:
        for a, amr in exp.items():
            if amr and a not in matched[iso]:
                truth.append({"isolate": iso, "locus_id": "missed:" + a, "truth": "positive", "matched_expected": a,
                              "novelty": "not_detected", "priority_score": "NA", "blind_score": -99,
                              "n_tools": 0, "antibiotic_type": 0, "supported": False, "deepbgc_activity": 0})
P = [x for x in truth if x["truth"] == "positive"]
N = [x for x in truth if x["truth"] == "negative"]
has_type = any("antibiotic_type" in r for r in master)
pred = lambda x: x["blind_score"] >= BLIND_CUTOFF and (not has_type or (x["supported"] and x["antibiotic_type"]))
tp, fn = sum(map(pred, P)), sum(not pred(x) for x in P)
fp, tn = sum(map(pred, N)), sum(not pred(x) for x in N)

# ---- 3. genome level ----
strong = lambda i: any(r["isolate"] == i and r["novelty"] == "novel_candidate" and num(r["priority_score"], 0) >= 5 for r in master)
gp = [i for i in in_run if panel.get(i) == "positive"]
gn = [i for i in in_run if panel.get(i) == "negative"]
g = (sum(map(strong, gp)), sum(map(strong, gn)), sum(not strong(i) for i in gp), sum(not strong(i) for i in gn))

confusion = [["locus: antibiotic BGC (supported, antibiotic type, blind score >= %d)" % BLIND_CUTOFF if has_type
              else "locus: antibiotic BGC (blind score >= %d)" % BLIND_CUTOFF, tp, fp, fn, tn, n_cryptic],
             ["genome: producer (>= 1 novel candidate, score >= 5)", g[0], g[1], g[2], g[3], 0]]
met = []
for name, c in (("locus", (tp, fp, fn, tn)), ("genome", g)):
    for k, (v, lo, hi) in metrics(*c).items():
        met.append([name, k, v, lo, hi])
features = {"blind priority score": "blind_score", "number of tools": "n_tools",
            "DeepBGC antibacterial/antifungal probability": "deepbgc_activity"}
if any(r.get("antibiotic_type") for r in master):
    features["antibiotic-relevant type (yes/no)"] = "antibiotic_type"
aucs, rocs = [], []
for lab, k in features.items():
    a, lo, hi = auc([x[k] for x in P], [x[k] for x in N])
    aucs.append([lab, len(P), len(N), a, lo, hi])
    if P and N:
        rocs += [[lab, t, f, s] for f, s, t in roc([x[k] for x in P], [x[k] for x in N])]

paths = [
    write("eval_detection.tsv", ["tool", "expected_bgcs", "found", "sensitivity", "ci95_low", "ci95_high",
                                 "loci_in_negatives", "loci_per_Mb_in_negatives"], det),
    write("eval_confusion.tsv", ["level", "TP", "FP", "FN", "TN", "left_out_unknown"], confusion),
    write("eval_metrics.tsv", ["level", "metric", "value", "ci95_low", "ci95_high"], met),
    write("eval_auc.tsv", ["score", "positives", "negatives", "AUC", "ci95_low", "ci95_high"], aucs),
    write("eval_roc.tsv", ["score", "threshold", "FPR", "TPR"], rocs),
    write("eval_loci_truth.tsv", list(truth[0].keys()) if truth else ["isolate"], [list(x.values()) for x in truth]),
]
print("== Detection (expected MIBiG BGCs found)")
for d in det:
    print(f"  {d[0]:12s} {d[2]}/{d[1]}  sensitivity {d[3] or 0:.2f} (95 % CI {d[4] or 0:.2f}-{d[5] or 0:.2f})  "
          f"loci/Mb in negatives {d[7] if d[7] is None else round(d[7], 2)}")
print(f"== Antibiotic ranking: {len(P)} positive, {len(N)} negative loci ({n_cryptic} cryptic loci left out)")
print(f"  confusion ({confusion[0][0]}): TP {tp}  FP {fp}  FN {fn}  TN {tn}")
for m in met:
    if m[0] == "locus" and m[2] is not None:
        print(f"  {m[1]:28s} {m[2]:.2f}" + (f" (95 % CI {m[3]:.2f}-{m[4]:.2f})" if m[3] is not None else ""))
for a in aucs:
    print(f"  AUC {a[0]:45s} {a[3] if a[3] is None else round(a[3], 2)}")
print(f"== Genome level: TP {g[0]} FP {g[1]} FN {g[2]} TN {g[3]}")
print("Tables:", *paths, sep="\n  ")


if __name__ == "__main__" and os.environ.get("EVAL_SELFTEST"):
    a, lo, hi = auc([3, 2, 1], [0, 0, 1])
    assert abs(a - (8 + 0.5) / 9) < 1e-9 and lo < a < hi
    assert metrics(8, 2, 2, 8)["accuracy"][0] == 0.8 and wilson(0, 10)[0] == 0.0
    print("self-test passed")
