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
  "ranking": [[{"accession": "BGC0000001", "description": "testomycin"}, {"similarity": 100}]]}]}}}}]}
EOF
printf 'sequence_id\tcluster_id\tstart\tend\taverage_p\tmax_p\ttype\nctgA\tctgA_cluster_1\t950\t3900\t0.9\t0.99\tNRP\n' >"$R/09_gecco/iso1.clusters.tsv"
printf 'sequence_id\tbgc_candidate_id\tnucl_start\tnucl_end\tdeepbgc_score\tproduct_class\tproduct_activity\tantibacterial\tantifungal\tcytotoxic\tinhibitor\n' >"$R/13_deepbgc/13_deepbgc.bgc.tsv"
printf 'ctgA\tctgA_9899-12100.1\t9899\t12100\t0.7\tPolyketide\tantibacterial\t0.8\t0.1\t0.2\t0.1\n' >>"$R/13_deepbgc/13_deepbgc.bgc.tsv"
# lone DeepBGC call on ctgB (no other tool, no MIBiG hit): must become low_support, not a candidate
printf 'ctgB\tctgB_99-1400.1\t99\t1400\t0.6\tPolyketide\tantibacterial\t0.9\t0.1\t0.1\t0.1\n' >>"$R/13_deepbgc/13_deepbgc.bgc.tsv"
cat >"$R/14_bgc_prophet/iso1_classified.csv" <<'EOF'
ID,sentence,labels,isBGC,TDsentence,TDlabels,probability
iso1_0,"['g4', 'g5']",Polyketide,Yes,g1 g2,0 0,0.910
iso1_1,"['g4', 'g5']",Polyketide,Yes,g1 g2,0 0,0.910
iso1_2,,NonBGC,No,g6,0,0.000
EOF
printf 'Contig\tStart\tStop\tBest_Hit_ARO\tDrug Class\nctgA_5\t10500\t10800\ttestR\ttetracycline antibiotic\n' >"$R/10_rgi_card/iso1_rgi.txt"

python3 "$T/summary.py" --results "$T/results" --samples iso1 --out "$T/out" --mibig "$K" 2>/dev/null
python3 - "$T/out" <<'EOF'
import csv, os, sys
d = sys.argv[1]
rd = lambda p: list(csv.DictReader(open(os.path.join(d, p)), delimiter="\t"))
m = rd("master_summary.tsv")
assert len(m) == 3, m
top, known, lone = m
assert (top["tools"], top["novelty"], top["priority_score"]) == ("deepbgc;bgc_prophet", "novel_candidate", "9"), top
assert (top["closest_MIBiG"], top["MIBiG_similarity_pct"]) == ("BGC0000002", "25"), top        # 1 of 4 genes
assert top["resistance_gene_in_locus"].startswith("testR") and top["completeness"] == "NA"
assert (known["tools"], known["novelty"], known["priority_score"]) == ("antismash;gecco", "known", "3"), known
assert (lone["tools"], lone["novelty"], lone["priority_score"]) == ("deepbgc", "low_support", "1"), lone
assert (top["antibiotic_type"], known["antibiotic_type"]) == ("yes", "yes")       # ML Polyketide / antiSMASH NRPS
assert known["MIBiG_from"].startswith("antismash (antiSMASH KnownClusterBlast)"), known
p = rd("per_tool/bgc_prophet_bgcs.tsv")
assert len(p) == 1 and (p[0]["start"], p[0]["end"]) == ("10000", "12000"), p                 # duplicate window removed
g = rd("per_tool/gecco_bgcs.tsv")[0]
assert (g["MIBiG_accession"], g["MIBiG_similarity_pct"], g["locus_id"]) == ("BGC0000001", "100", known["locus_id"]), g
assert g["MIBiG_all_hits"] == "BGC0000001:100", g
assert rd("per_tool/antismash_bgcs.tsv")[0]["MIBiG_all_hits"] == "BGC0000001:100"
assert rd("per_tool/deepbgc_bgcs.tsv")[0]["start"] == "9900"                                  # 0-based -> 1-based
assert rd("genome_summary.tsv")[0]["BGC_loci_all_tools"] == "3"
print("summary self-check: all assertions passed")
EOF
