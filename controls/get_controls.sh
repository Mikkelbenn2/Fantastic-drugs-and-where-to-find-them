#!/usr/bin/env bash
# Downloads the positive and negative controls (see controls/README.md; sources and papers in
# controls/control_panel.tsv). All real, published data:
#   genomes = the complete reference genome from NCBI (the pipeline uses a .fna as is: no assembly),
#   reads   = a real R10.4.1 ONT run from ENA (goes through read QC + assembly like an isolate).
#   bash controls/get_controls.sh OUT_DIR        e.g. ../controls_run/raw_data
# Needs internet: run on the login node. The genome files are the same as from
#   datasets download genome accession ACC --include genome
set -euo pipefail
OUT=${1:?usage: bash controls/get_controls.sh OUT_DIR}
mkdir -p "$OUT"

# isolate name              NCBI assembly     (positive = MIBiG BGCs expected; negative = none/known only)
GENOMES="
pos_Bvelezensis_FZB42       GCF_000015785.2
pos_Scoelicolor_A3-2        GCF_000203835.1
pos_Pprotegens_Pf-5         GCF_000012265.1
pos_Svenezuelae_ATCC10712   GCA_000253235.1
pos_Sgriseus_NBRC13350      GCF_000010605.1
pos_Serythraea_NRRL2338     GCF_000062885.1
pos_Ecoli_Nissle1917        GCF_000714595.1
neg_Buchnera_APS            GCF_000009605.1
neg_Rprowazekii_MadridE     GCF_000195735.1
neg_Cglutamicum_ATCC13032   GCF_000011325.1
neg_Ecoli_K12_MG1655_ref    GCF_000005845.2
"
# isolate name              ENA run           (E. coli K-12 MG1655 wild type, Kit14/R10.4.1, Dorado sup v5)
READS="
neg_Ecoli_K12_MG1655        SRR30993121
"

while read -r name acc; do
  [[ -n "$name" ]] || continue
  [[ -s "$OUT/$name.fna" ]] && { echo "$name: present"; continue; }
  echo "$name: downloading genome $acc"
  curl -fsSL -o "$OUT/$name.zip" \
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
  python3 -c 'import sys, zipfile; z = zipfile.ZipFile(sys.argv[1]); open(sys.argv[2], "wb").write(b"".join(z.read(n) for n in z.namelist() if n.endswith(".fna")))' \
    "$OUT/$name.zip" "$OUT/$name.fna.part"
  rm "$OUT/$name.zip"
  [[ -s "$OUT/$name.fna.part" ]] || { echo "$name: no genome FASTA in the download" >&2; exit 1; }
  mv "$OUT/$name.fna.part" "$OUT/$name.fna"
done <<<"$GENOMES"

while read -r name run; do
  [[ -n "$name" ]] || continue
  [[ -s "$OUT/$name.fastq.gz" ]] && { echo "$name: present"; continue; }
  url=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2 | cut -d';' -f1)
  echo "$name: downloading $run"
  curl -fSL -o "$OUT/$name.fastq.gz.part" "https://$url"
  mv "$OUT/$name.fastq.gz.part" "$OUT/$name.fastq.gz"
done <<<"$READS"
echo "Done: $(ls "$OUT"/*.fna "$OUT"/*.fastq.gz 2>/dev/null | wc -l) controls in $OUT"
