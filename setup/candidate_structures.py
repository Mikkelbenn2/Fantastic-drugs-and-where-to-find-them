#!/usr/bin/env python3
"""Structure predictions for the best novel BGC candidates, taken from antiSMASH (no new tools needed).

   python3 setup/candidate_structures.py [RESULTS_DIR] [--min-score 5] [--keep-edge]

Picks loci from RESULTS_DIR/_combined/summary/master_summary.tsv with novelty = novel_candidate and
priority score >= --min-score, not cut by a contig edge (incomplete clusters give wrong structures).
For each antiSMASH region it reads <isolate>/08_antismash/<region>.gbk and collects what antiSMASH
predicts about the product:
  NRPS/PKS : the polymer (monomer order) and SMILES of the predicted core, A-domain substrates
  RiPPs    : precursor class, leader and core peptide, molecular weight / monoisotopic mass
Writes RESULTS_DIR/_combined/structures/candidate_structures.tsv and copies the region GenBank files
to RESULTS_DIR/_combined/structures/gbk/ (input for PRISM or the antiSMASH web viewer).
Loci found only by the machine-learning tools have no antiSMASH region and are listed without a
prediction. Predictions are hypotheses (no tailoring steps for most classes), not confirmed structures.
Needs python3 only."""
import argparse, csv, os, re, shutil, sys

KEYS = ("SMILES", "polymer", "core_sequence", "leader_sequence", "predicted_class",
        "molecular_weight", "monoisotopic_mass", "peptide")


def gbk_features(path):
    """GenBank features -> list of (type, {qualifier: [values]}); stdlib only, multi-line values joined"""
    feats, cur, key, inside = [], None, None, False
    for line in open(path):
        if line.startswith("FEATURES"):
            inside = True; continue
        if not inside:
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            break
        if len(line) > 5 and line[5] != " ":                    # new feature: type at column 6
            cur = (line[5:21].strip(), {}); feats.append(cur); key = None; continue
        text = line[21:].rstrip("\n")
        m = re.match(r'/(\w+)(?:=(.*))?$', text)
        if m and cur:
            key = m.group(1); cur[1].setdefault(key, []).append((m.group(2) or "").strip('"'))
        elif key and cur:                                        # continuation of the last value
            sep = "" if key in ("translation", "SMILES", "core_sequence", "leader_sequence") else " "
            cur[1][key][-1] = (cur[1][key][-1] + sep + text.strip()).strip('"')
    return feats


def predictions(path):
    out = {k: [] for k in KEYS}
    out["A_domain_substrates"] = []
    for ftype, q in gbk_features(path):
        for k in KEYS:
            out[k] += [v for v in q.get(k, []) if v]
        if ftype == "aSDomain":                                  # A domains: consensus substrate
            out["A_domain_substrates"] += [v.split(":", 1)[1].strip() for v in q.get("specificity", [])
                                           if v.lower().startswith("consensus")]
    return {k: " | ".join(dict.fromkeys(v)) or "" for k, v in out.items()}   # unique, in order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default="3_results")
    ap.add_argument("--min-score", type=float, default=5)
    ap.add_argument("--keep-edge", action="store_true", help="also keep loci cut by a contig edge")
    a = ap.parse_args()
    master = os.path.join(a.results, "_combined", "summary", "master_summary.tsv")
    if not os.path.isfile(master):
        sys.exit(f"ERROR: {master} not found (run the pipeline's combined step first)")
    out_dir = os.path.join(a.results, "_combined", "structures")
    os.makedirs(os.path.join(out_dir, "gbk"), exist_ok=True)

    def score(r):
        try:
            return float(r["priority_score"])
        except (KeyError, ValueError):
            return float("-inf")

    rows = []
    for r in csv.DictReader(open(master), delimiter="\t"):
        if r.get("novelty") != "novel_candidate" or score(r) < a.min_score:
            continue
        if r.get("contig_edge") == "yes" and not a.keep_edge:
            continue
        base = {k: r.get(k, "") for k in ("rank", "priority_score", "isolate", "locus_id", "contig", "start", "end",
                                          "length_kb", "tools", "antismash_type", "closest_MIBiG_compound",
                                          "MIBiG_similarity_pct", "contig_edge")}
        regions = [x for x in (r.get("antismash_regions") or "").split(";") if x and x != "NA"]
        if not regions:
            rows.append({**base, "region": "", "note": "no antiSMASH region (machine-learning tools only)"})
        for reg in regions:
            gbk = os.path.join(a.results, r["isolate"], "08_antismash", reg + ".gbk")
            if not os.path.isfile(gbk):
                rows.append({**base, "region": reg, "note": f"missing {gbk}"}); continue
            shutil.copy(gbk, os.path.join(out_dir, "gbk", f"{r['isolate']}__{reg}.gbk"))
            p = predictions(gbk)
            note = "" if any(p.values()) else "antiSMASH gives no structure prediction for this class"
            rows.append({**base, "region": reg, **p, "note": note})

    cols = list(dict.fromkeys(k for r in rows for k in r))
    path = os.path.join(out_dir, "candidate_structures.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval=""); w.writeheader(); w.writerows(rows)
    n_pred = sum(1 for r in rows if r.get("region") and not r.get("note"))
    print(f"{len(rows)} candidate regions ({n_pred} with an antiSMASH structure prediction) -> {path}\n"
          f"GenBank files for PRISM / antiSMASH viewer: {os.path.join(out_dir, 'gbk')}")


if __name__ == "__main__":
    main()
