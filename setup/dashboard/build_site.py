#!/usr/bin/env python3
"""Builds the shareable results dashboard: ONE html file with all data inside, opened in any browser.
   python3 setup/dashboard/build_site.py RESULTS/_combined/summary [--out FILE] [--project DIR]
Reads the summary tables of the pipeline (+ metadata.tsv, setup/controls/control_panel.tsv, software_versions.tsv)
and writes RESULTS/_combined/summary/dashboard.html. Needs only python3; nothing runs on the cluster
afterwards: copy the file and send it. plotly.js (charts) is put inside the file when found (R's plotly
package, e.g. setup/envs/dashboard), so it works offline; otherwise it is loaded from jsdelivr (needs internet)."""
import argparse, csv, datetime, glob, json, os, subprocess, sys
from html import escape

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("summary")
ap.add_argument("--out")
ap.add_argument("--project", help="folder with metadata.tsv (default: 3 levels above summary)")
ap.add_argument("--title", default="BGC screening of soil isolates")
ap.add_argument("--pipeline-version", default="")
ap.add_argument("--add", action="append", default=[], metavar="SUMMARY",
                help="another results/_combined/summary to merge in (e.g. split genomes); repeatable")
a = ap.parse_args()
summ = os.path.abspath(a.summary)
project = os.path.abspath(a.project or os.path.join(summ, "..", "..", ".."))
out = a.out or os.path.join(summ, "dashboard.html")


def rows(*parts):
    p = os.path.join(*parts)
    if not os.path.isfile(p) or os.path.getsize(p) == 0:
        return []
    with open(p, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if any(r.values())]


if not os.path.isfile(os.path.join(summ, "master_summary.tsv")):
    sys.exit(f"ERROR: no master_summary.tsv in {summ}: run the pipeline's combined/summary step first")
genomes = rows(summ, "genome_summary.tsv")
for extra in a.add:   # merged results (isolate, BGC, contig, read and resistance tables)
    genomes += rows(extra, "genome_summary.tsv")
names = {g["isolate"] for g in genomes}
# sample metadata: your metadata.tsv (first column = isolate) + the control panel, only isolates in this run
meta = {}
for src in (rows(project, "metadata.tsv"), rows(here, "..", "controls", "control_panel.tsv")):
    for r in src:
        iso = r.get("isolate") or next(iter(r.values()), "")
        if iso in names:
            meta.setdefault(iso, {}).update({k: v for k, v in r.items() if k != "isolate" and v})


def bgcquast(mode):
    """BGC-QUAST report.tsv files (metrics as rows, one column per tool or isolate) -> long rows"""
    out = []
    for p in sorted(glob.glob(os.path.join(summ, "..", "bgc_quast", mode, "*", "report.tsv"))):
        with open(p, newline="") as fh:
            t = list(csv.reader(fh, delimiter="\t"))
        if len(t) < 3:
            continue
        group = os.path.basename(os.path.dirname(p))      # isolate (compare_tools) or tool (compare_samples)
        for j, col in enumerate(t[0][1:], start=1):
            r = {"group": group, "label": col}
            r.update({m[0]: m[j] for m in t[2:] if len(m) > j and m[0] in KEEP or (len(m) > j and m[0].startswith("# BGCs ("))})
            out.append(r)
    return out


KEEP = {"# BGCs filtered by length", "# BGCs", "# unique BGCs", "Unique BGC rate", "Mean BGC length in bp",
        "Mean BGC length in genes", "Total BGC span in bp"}
def flye_graphs(summary, isolates):
    """Flye assembly_graph.gfa per isolate -> segments [name, length, depth, contig] and links [a, +/-, b, +/-]
       (sequences are skipped, so the dashboard stays small)"""
    out, res = {}, os.path.abspath(os.path.join(summary, "..", ".."))
    for iso in isolates:
        gfa = os.path.join(res, iso, "03_flye", "assembly_graph.gfa")
        if not os.path.isfile(gfa):
            continue
        ctg = {}
        for r in rows(res, iso, "03_flye", "assembly_info.txt"):   # graph_path "-5,*,12" -> edge_5, edge_12
            for e in (r.get("graph_path") or "").split(","):
                if e.strip("-+*"):
                    ctg.setdefault("edge_" + e.strip("-+"), r.get("#seq_name", ""))
        seg, link = [], set()
        with open(gfa) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "S":
                    tags = dict(t.split(":", 2)[0::2] for t in f[3:] if t.count(":") >= 2)
                    seg.append([f[1], int(tags.get("LN", len(f[2]))), float(tags.get("dp", tags.get("DP", 0)) or 0), ctg.get(f[1], "")])
                elif f[0] == "L":
                    link.add((f[1], f[2], f[3], f[4]))   # from, +/-, to, +/- (which ends touch)
        if len(seg) <= 3000:   # shortcut: larger graphs (fragmented assemblies) are left out; raise if needed
            out[iso] = {"s": seg, "l": sorted(link)}
    return out


versions = rows(summ, "software_versions.tsv") or rows(summ, "..", "..", "software_versions.tsv")

data = {
    "title": a.title, "pipeline_version": a.pipeline_version, "source": summ,
    "built": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    "genomes": genomes, "master": rows(summ, "master_summary.tsv"), "contigs": rows(summ, "contigs.tsv"),
    "reads": [[r["isolate"], int(r["length"]), float(r["meanQ"])] for r in rows(summ, "reads_sample.tsv")],
    "resistance": rows(summ, "resistance_genes.tsv"), "settings": rows(summ, "settings.tsv"),
    "versions": versions, "metadata": meta, "graphs": flye_graphs(summ, names),
    "bgcquast": {"tools": bgcquast("compare_tools"), "samples": bgcquast("compare_samples")},
    "controls": {"recovery": rows(summ, "controls_recovery.tsv"), "loci": rows(summ, "controls_loci.tsv"),
                 "expected": rows(summ, "controls_expected_bgcs.tsv"),
                 **{k: rows(summ, f"eval_{k}.tsv") for k in ("confusion", "metrics", "auc", "roc", "detection")}},
}


def plotly_js():
    rel = "plotly/htmlwidgets/lib/plotlyjs/plotly-latest.min.js"
    cand = [os.environ.get("BGC_PLOTLY_JS", "")] + glob.glob(os.path.join(project, "setup/envs/dashboard/lib/R/library", rel))
    try:
        cand.append(subprocess.run(["Rscript", "-e", 'cat(system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package="plotly"))'],
                                   capture_output=True, text=True, timeout=30).stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    for p in cand:
        if p and os.path.isfile(p):
            return "<script>" + open(p).read().replace("</script", "<\\/script") + "</script>", p
    return '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>', "jsdelivr (online)"


for extra in a.add:
    if not os.path.isfile(os.path.join(extra, "master_summary.tsv")):
        sys.exit(f"ERROR: no master_summary.tsv in {extra}")
    data["master"] += rows(extra, "master_summary.tsv"); data["contigs"] += rows(extra, "contigs.tsv")
    data["reads"] += [[r["isolate"], int(r["length"]), float(r["meanQ"])] for r in rows(extra, "reads_sample.tsv")]
    data["resistance"] += rows(extra, "resistance_genes.tsv")
    data["graphs"].update(flye_graphs(extra, {g["isolate"] for g in rows(extra, "genome_summary.tsv")}))
    data["source"] += " + " + os.path.abspath(extra)
if a.add:   # one ranking over all merged loci, same order as the summary step
    def f(r, k):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return 0.0
    data["master"].sort(key=lambda r: (-f(r, "priority_score"), -f(r, "n_tools"), -f(r, "length_kb")))
    for i, r in enumerate(data["master"], start=1):
        r["rank"] = str(i)

js, js_from = plotly_js()
html = open(os.path.join(here, "site_template.html")).read()
blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")   # safe inside <script>
html = html.replace("__TITLE__", escape(a.title)).replace("__DATA__", blob)
html = html.replace("__PLOTLY__", js)            # last: the plotly code must not be searched for the tokens
with open(out, "w") as fh:
    fh.write(html)
print(f"Dashboard: {out} ({os.path.getsize(out) / 1e6:.1f} MB, {len(genomes)} isolates; plotly.js from {js_from})")
