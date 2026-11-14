#!/usr/bin/env bash
# Self-check of the summary code in bgc_pipeline.sh on tiny fake tool outputs.
# Run: bash setup/tests/test_summary.sh   (needs only python3; no conda envs or databases)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
awk '/^cat >"\$1" <<.PYTHON.$/{f=1;next} /^PYTHON$/{f=0} f' "$HERE/../bgc_pipeline.sh" >"$T/summary.py"

R="$T/results/iso1"; K="$T/kcb"
mkdir -p "$R"/{07_bakta,08_antismash,09_gecco,10_rgi_card,13_deepbgc,14_bgc_prophet,15_mibig} "$K"
# genes: g1-g3 = known NRPS cluster, g4-g5 = second cluster, g6 on another contig
printf '# Bakta\n#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\n' >"$R/07_bakta/iso1.tsv"
printf 'ctgA\tcds\t%s\t%s\t+\t%s\n' 1000 2000 g1 2100 3000 g2 3100 4000 g3 10000 11000 g4 11100 12000 g5 >>"$R/07_bakta/iso1.tsv"
printf 'ctgB\tcds\t500\t1500\t+\tg6\n' >>"$R/07_bakta/iso1.tsv"
printf 'BGC0000001\ttestomycin\tc1\tnrps\tP1;P2\tP1;P2\nBGC0000002\totherin\tc1\tt1pks\tP3;P4;P5;P6\tP3;P4;P5;P6\n' >"$K/clusters.txt"
# g4 -> P3 passes the cut-offs; g5 -> P4 has 40 % identity, below KnownClusterBlast's 45 %
printf '%s\t%s\t%s\t300\t330\t330\t1e-50\t500\n' g1 'BGC0000001|c1|1-10|+|P1|x|P1' 80 g2 'BGC0000001|c1|1-10|+|P2|x|P2' 80 \
  g4 'BGC0000002|c1|1-10|+|P3|x|P3' 50 g5 'BGC0000002|c1|1-10|+|P4|x|P4' 40 >"$R/15_mibig/iso1_vs_mibig.tsv"
cat >"$R/08_antismash/iso1.json" <<'EOF'
{"records": [{"id": "ctgA", "areas": [{"start": 899, "end": 4100, "products": ["NRPS"], "protoclusters": {"1": {"category": "NRPS"}}}],
  "modules": {"antismash.modules.clusterblast": {"knowncluster": {"results": [{"region_number": 1,
  "ranking": [[{"accession": "BGC0000001", "description": "testomycin"}, {"similarity": 100}]]}]}}}},
 {"id": "ctgC", "areas": [{"start": 0, "end": 5000, "products": ["NRP-metallophore", "NRPS"], "protoclusters": {"1": {"category": "NRPS"}}}]}]}
EOF
printf 'sequence_id\tcluster_id\tstart\tend\taverage_p\tmax_p\ttype\nctgA\tctgA_cluster_1\t950\t3900\t0.9\t0.99\tNRP\n' >"$R/09_gecco/iso1.clusters.tsv"
# GECCO call that only touches the end of the lone DeepBGC call on ctgB (51 bp of 1301): NOT agreement
printf 'ctgB\tctgB_cluster_1\t1350\t3000\t0.9\t0.99\tNRP\n' >>"$R/09_gecco/iso1.clusters.tsv"
printf 'sequence_id\tbgc_candidate_id\tnucl_start\tnucl_end\tdeepbgc_score\tproduct_class\tproduct_activity\tantibacterial\tantifungal\tcytotoxic\tinhibitor\n' >"$R/13_deepbgc/13_deepbgc.bgc.tsv"
printf 'ctgA\tctgA_9899-12100.1\t9899\t12100\t0.7\tPolyketide\tantibacterial\t0.8\t0.1\t0.2\t0.1\n' >>"$R/13_deepbgc/13_deepbgc.bgc.tsv"
# lone DeepBGC call on ctgB (no other tool, no MIBiG hit): must become low_support, not a candidate
printf 'ctgB\tctgB_99-1400.1\t99\t1400\t0.6\tPolyketide\tantibacterial\t0.9\t0.1\t0.1\t0.1\n' >>"$R/13_deepbgc/13_deepbgc.bgc.tsv"
# a 500 bp DeepBGC fragment: shorter than --min-bgc-length (1000 here), must be dropped before loci are made
printf 'ctgA\tctgA_19999-20500.1\t19999\t20500\t0.9\tPolyketide\tantibacterial\t0.9\t0.1\t0.1\t0.1\n' >>"$R/13_deepbgc/13_deepbgc.bgc.tsv"
cat >"$R/14_bgc_prophet/iso1_classified.csv" <<'EOF'
ID,sentence,labels,isBGC,TDsentence,TDlabels,probability
iso1_0,"['g4', 'g5']",Polyketide,Yes,g1 g2,0 0,0.910
iso1_1,"['g4', 'g5']",Polyketide,Yes,g1 g2,0 0,0.910
iso1_2,,NonBGC,No,g6,0,0.000
EOF
printf 'Contig\tStart\tStop\tBest_Hit_ARO\tDrug Class\tModel_type\nctgA_5\t10500\t10800\ttestR\ttetracycline antibiotic\tprotein homolog model\n' >"$R/10_rgi_card/iso1_rgi.txt"
# a mutation (variant) model in a housekeeping gene inside the same locus: must be left out (intrinsic resistome)
printf 'ctgA_6\t11000\t11500\tEF-Tu mutant\telfamycin antibiotic\tprotein variant model\n' >>"$R/10_rgi_card/iso1_rgi.txt"

# reads (MinKNOW header), Flye contigs and the Bakta report: for the dashboard's sequencing/assembly pages
mkdir -p "$R"/{01_nanoplot,03_flye}
printf '@r1 runid=abc123def456 start_time=2026-09-01T10:00:00Z flow_cell_id=PAW12345 barcode=barcode07 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\nACGT\n+\n5555\n@r2\nAC\n+\n++\n' | gzip -c >"$R/01_nanoplot/iso1.reads.fastq.gz"
printf '#seq_name\tlength\tcov.\tcirc.\nctgA\t12\t50\tY\nctgB\t4\t80\tN\n' >"$R/03_flye/assembly_info.txt"
printf '>ctgA desc\nACGTACGT\nACGT\n>ctgB\nACGT\n' >"$R/03_flye/iso1.fasta"
printf 'Database: v6.0, full\nCDSs: 6\ntRNAs: 2\n' >"$R/07_bakta/iso1.txt"
python3 "$T/summary.py" --results "$T/results" --samples iso1 --out "$T/out" --mibig "$K" --min-bgc-length 1000 2>/dev/null
python3 - "$T/out" <<'EOF'
import csv, os, sys
d = sys.argv[1]
rd = lambda p: list(csv.DictReader(open(os.path.join(d, p)), delimiter="\t"))
m = rd("master_summary.tsv")
assert len(m) == 4, m
top, known, lone, sidero = m
assert (top["tools"], top["novelty"], top["priority_score"]) == ("deepbgc;bgc_prophet", "novel_candidate", "9"), top
assert (top["closest_MIBiG"], top["MIBiG_similarity_pct"]) == ("BGC0000002", "25"), top        # 1 of 4 genes
assert top["resistance_gene_in_locus"].startswith("testR") and top["completeness"] == "NA"
assert (known["tools"], known["novelty"], known["priority_score"]) == ("antismash;gecco", "known", "3"), known
assert (lone["tools"], lone["novelty"], lone["priority_score"]) == ("gecco;deepbgc", "low_support", "1"), lone   # touching != agreeing
assert (top["antibiotic_type"], known["antibiotic_type"]) == ("yes", "yes")       # ML Polyketide / antiSMASH NRPS
# NRPS siderophore (antiSMASH NRP-metallophore): never antibiotic-relevant -> novel_other, no +2
assert (sidero["antismash_type"], sidero["antibiotic_type"], sidero["novelty"], sidero["priority_score"]) == \
    ("NRP-metallophore+NRPS", "no", "novel_other", "0"), sidero
assert known["MIBiG_from"].startswith("antismash (antiSMASH KnownClusterBlast)"), known
p = rd("per_tool/bgc_prophet_bgcs.tsv")
assert len(p) == 1 and (p[0]["start"], p[0]["end"]) == ("10000", "12000"), p                 # duplicate window removed
g = rd("per_tool/gecco_bgcs.tsv")[0]
assert (g["MIBiG_accession"], g["MIBiG_similarity_pct"], g["locus_id"]) == ("BGC0000001", "100", known["locus_id"]), g
assert g["MIBiG_all_hits"] == "BGC0000001:100", g
assert rd("per_tool/antismash_bgcs.tsv")[0]["MIBiG_all_hits"] == "BGC0000001:100"
assert rd("per_tool/deepbgc_bgcs.tsv")[0]["start"] == "9900"                                  # 0-based -> 1-based
gs = rd("genome_summary.tsv")[0]
assert gs["BGC_loci_all_tools"] == "4" and gs["short_calls_removed"] == "1" and gs["BGCs_deepbgc"] == "2", gs
assert (top["n_tools_agreeing"], known["n_tools_agreeing"], lone["n_tools_agreeing"]) == ("2", "2", "0")
pc = list(csv.DictReader(open(os.path.join(os.path.dirname(d), "results", "iso1", "14_bgc_prophet", "iso1_coords.tsv")), delimiter="\t"))
assert [(x["sequence_id"], x["start"], x["end"], x["type"]) for x in pc] == [("ctgA", "10000", "12000", "Polyketide")], pc
assert (gs["basecall_model"], gs["flow_cell_id"], gs["run_id"], gs["barcode"]) == \
    ("dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "PAW12345", "abc123def456", "barcode07"), gs
assert (gs["bakta_db"], gs["CDS"], gs["tRNA"]) == ("v6.0, full", "6", "2"), gs
r = rd("reads_sample.tsv")
assert [(x["length"], x["meanQ"]) for x in r] == [("4", "20.0"), ("2", "10.0")], r        # '5' = Q20, '+' = Q10
c = rd("contigs.tsv")
assert [(x["contig"], x["length"], x["depth"], x["circular"]) for x in c] == [("ctgA", "12", "50", "yes"), ("ctgB", "4", "80", "no")], c
rg = rd("resistance_genes.tsv")
assert [x["gene"] for x in rg] == ["testR"] and rg[0]["in_BGC_locus"] == top["locus_id"], rg     # EF-Tu variant left out
assert gs["resistance_genes"] == gs["resistance_genes_in_BGCs"] == "1", gs
print("summary self-check: all assertions passed")
EOF
# the shareable dashboard is built from the same tables and carries them as JSON
python3 "$HERE/../dashboard/build_site.py" "$T/out" --project "$T" --out "$T/dash.html" >/dev/null
python3 - "$T/dash.html" <<'EOF'
import json, re, sys
h = open(sys.argv[1]).read()
d = json.loads(re.search(r'<script id="data" type="application/json">(.*?)</script>', h, re.S).group(1))
assert len(d["master"]) == 4 and d["genomes"][0]["flow_cell_id"] == "PAW12345" and len(d["reads"]) == 2, d.keys()
assert "__DATA__" not in h and "__PLOTLY__" not in h
print("dashboard self-check: passed")
EOF
