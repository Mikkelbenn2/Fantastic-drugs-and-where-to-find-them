#!/usr/bin/env python3
"""Rebuilds expected_bgcs.tsv (the truth set of the controls) from the MIBiG 4.0.1 release.
   python3 setup/controls/make_expected_bgcs.py [mibig_json_4.0.tar.gz] [--check]
Without a file it downloads the JSON release from Zenodo (record 14835872, ~10 MB).
--check: do not write, only compare with the existing expected_bgcs.tsv (exit 1 if different).

Rules (written in the README and the report):
1. MIBiG entry status "active" (retired and pending entries are left out).
2. Source organism = the control strain, matched by its strain-level NCBI TaxID (TAXID below), so a BGC
   characterised in another strain of the same species does not count, and old species names
   (e.g. Pseudomonas fluorescens Pf-5) cannot slip through.
3. antimicrobial = "yes" when the compound is a named antibiotic in the strain's paper
   (control_panel.tsv, expected_antibiotics) or MIBiG lists antibacterial/antifungal/antimicrobial/anti-Gram activity,
   except the CURATED entries, which were checked in the literature."""
import csv, io, json, os, re, sys, tarfile, urllib.request

here = os.path.dirname(os.path.abspath(__file__))
URL = "https://zenodo.org/api/records/14835872/files/mibig_json_4.0.tar.gz/content"
TAXID = {  # control isolate -> NCBI TaxID of exactly that strain (as used by MIBiG)
    "pos_Bvelezensis_FZB42": 326423, "pos_Scoelicolor_A3-2": 100226, "pos_Pprotegens_Pf-5": 220664,
    "pos_Svenezuelae_ATCC10712": 953739, "pos_Sgriseus_NBRC13350": 455632, "pos_Serythraea_NRRL2338": 405948,
    "pos_Ecoli_Nissle1917": 316435, "neg_Ecoli_K12_MG1655": 511145,
}
CURATED = {  # accession: (antimicrobial, evidence) - borderline cases checked in the literature
    "BGC0002152": ("yes", "curated: toxoflavin from Pf-5 inhibits plant-pathogenic bacteria (Philmus 2015, doi:10.1002/cbic.201500247)"),
    "BGC0000038": ("yes", "curated: the cpk cluster makes the antibiotic coelimycin A; P1 is its inactive derivative (Bednarz 2019, doi:10.1007/s00253-019-09975-w)"),
    "BGC0000349": ("no", "curated: putative siderophore, no antibacterial activity shown (Lazos 2010, doi:10.1016/j.chembiol.2010.01.011); MIBiG tag not followed"),
    "BGC0001043": ("no", "curated: new PTMs from a cryptic cluster, no antimicrobial activity reported (Luo 2013, doi:10.1038/ncomms3894)"),
    "BGC0002641": ("no", "curated: no antimicrobial activity of bacillothiazols found in MIBiG or a PubMed search"),
}
COLS = ["isolate", "control", "mibig_accession", "compound", "mibig_class", "bioactivity", "mibig_locus",
        "antimicrobial", "antimicrobial_evidence"]


def entries(path):
    data = open(path, "rb").read() if path else urllib.request.urlopen(URL, timeout=120).read()
    with tarfile.open(fileobj=io.BytesIO(data)) as t:
        for m in t:
            if m.name.endswith(".json"):
                yield json.load(t.extractfile(m))


def activity(b):
    """MIBiG 4 bioactivity: {"name": {"activity": "antibacterial"}} or {"name": "antibacterial"}"""
    n = b.get("name")
    return (n.get("activity", "") if isinstance(n, dict) else n or b.get("activity", "")) or ""


def paper_names(text):
    """'actinorhodin, calcium-dependent antibiotic (CDA)' -> ['actinorhodin', 'calcium-dependent antibiotic', 'cda']"""
    out = []
    for n in text.lower().split(","):
        n = n.strip()
        if not n or n.startswith("none"):
            continue
        out.append(n.split(" (")[0])
        if "(" in n:
            out.append(n.split("(")[-1].rstrip(")"))
    return out


def label(acc, isolate, compound, bioactivity, panel):
    if acc in CURATED:
        return CURATED[acc]
    comp = compound.lower()
    by_paper = any(n in comp for n in paper_names(panel.get(isolate, {}).get("expected_antibiotics", "")))
    by_mibig = bool(re.search(r"antibacterial|antifungal|antimicrobial|anti-gram", bioactivity.lower()))
    ev = "+".join(x for x, y in (("control paper", by_paper), ("MIBiG bioactivity", by_mibig)) if y) or "-"
    return ("yes" if by_paper or by_mibig else "no"), ev


def main():
    args = [a for a in sys.argv[1:] if a != "--check"]
    panel = {r["isolate"]: r for r in csv.DictReader(open(os.path.join(here, "control_panel.tsv")), delimiter="\t")}
    by_tax = {v: k for k, v in TAXID.items()}
    rows = []
    for d in entries(args[0] if args else None):
        iso = by_tax.get((d.get("taxonomy") or {}).get("ncbiTaxId"))
        if not iso or d.get("status") != "active":
            continue
        comps = [c.get("name", "") for c in d.get("compounds", [])]
        acts = sorted({activity(b) for c in d.get("compounds", []) for b in c.get("bioactivities", [])} - {""})
        compound = ";".join(comps) or "NA"
        bio = ",".join(acts) or "NA"
        r = {"isolate": iso, "control": panel[iso]["control"], "mibig_accession": d["accession"], "compound": compound,
             "mibig_class": "+".join(c["class"] for c in d.get("biosynthesis", {}).get("classes", [])) or "NA",
             "bioactivity": bio, "mibig_locus": (d.get("loci") or [{}])[0].get("accession", "NA")}
        r["antimicrobial"], r["antimicrobial_evidence"] = label(d["accession"], iso, compound, bio, panel)
        rows.append(r)
    order = list(panel)
    rows.sort(key=lambda r: (order.index(r["isolate"]), r["mibig_accession"]))
    pos = [r for r in rows if r["control"] == "positive"]
    print(f"{len(rows)} BGCs ({len(pos)} for positive controls; "
          f"{sum(r['antimicrobial'] == 'yes' for r in pos)} antimicrobial)", file=sys.stderr)

    out = os.path.join(here, "expected_bgcs.tsv")
    if "--check" in sys.argv:
        key = lambda rs: {tuple(r[c] for c in COLS) for r in rs}
        old = list(csv.DictReader(open(out), delimiter="\t"))
        diff = key(rows) ^ key(old)
        print("identical to expected_bgcs.tsv" if not diff else f"DIFFERENT: {sorted(diff)}", file=sys.stderr)
        sys.exit(1 if diff else 0)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"written {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
