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
versions = rows(summ, "software_versions.tsv") or rows(summ, "..", "..", "software_versions.tsv")

data = {
    "title": a.title, "pipeline_version": a.pipeline_version, "source": summ,
    "built": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    "genomes": genomes, "master": rows(summ, "master_summary.tsv"), "contigs": rows(summ, "contigs.tsv"),
    "reads": [[r["isolate"], int(r["length"]), float(r["meanQ"])] for r in rows(summ, "reads_sample.tsv")],
    "resistance": rows(summ, "resistance_genes.tsv"), "settings": rows(summ, "settings.tsv"),
    "versions": versions, "metadata": meta,
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


js, js_from = plotly_js()
html = open(os.path.join(here, "site_template.html")).read()
blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")   # safe inside <script>
html = html.replace("__TITLE__", escape(a.title)).replace("__DATA__", blob)
html = html.replace("__PLOTLY__", js)            # last: the plotly code must not be searched for the tokens
with open(out, "w") as fh:
    fh.write(html)
print(f"Dashboard: {out} ({os.path.getsize(out) / 1e6:.1f} MB, {len(genomes)} isolates; plotly.js from {js_from})")
