#!/usr/bin/env bash
# Downloads / simulates the positive and negative control read sets (see controls/README.md).
#   bash controls/get_controls.sh OUT_DIR        e.g. ../controls_run/raw_data
# Real reads come from ENA. Simulated reads: Badread (R10.4.1 "nanopore2023" model) on the
# published complete genome from NCBI. Needs internet: run on the login node.
set -euo pipefail
OUT=${1:?usage: bash controls/get_controls.sh OUT_DIR}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV="$HERE/../setup/envs/badread"
DEPTH=60x                      # coverage of the simulated read sets
mkdir -p "$OUT" "$OUT/../control_genomes"
GEN="$(cd "$OUT/../control_genomes" && pwd)"

CONDA=$(command -v mamba || command -v micromamba || command -v conda) || { echo "conda/mamba not found (module load miniforge?)" >&2; exit 1; }
[[ -d "$ENV/conda-meta" ]] || "$CONDA" create -y -p "$ENV" -c conda-forge -c bioconda badread >/dev/null

# isolate name           NCBI assembly     (positive = BGCs from MIBiG expected; negative = none/known only)
SIMULATED="
pos_Bvelezensis_FZB42       GCF_000015785.2
pos_Scoelicolor_A3-2        GCF_000203835.1
pos_Pprotegens_Pf-5         GCF_000012265.1
pos_Svenezuelae_ATCC10712   GCA_000253235.1
neg_Buchnera_APS            GCF_000009605.1
neg_Ecoli_K12_MG1655_sim    GCF_000005845.2
"
# isolate name           ENA run           (E. coli K-12 MG1655 wild type, Kit14/R10.4.1, Dorado sup v5)
REAL="
neg_Ecoli_K12_MG1655        SRR30993121
"

while read -r name acc; do
  [[ -n "$name" ]] || continue
  [[ -s "$OUT/$name.fastq.gz" ]] && { echo "$name: present"; continue; }
  fa="$GEN/$name.fna"
  if [[ ! -s "$fa" ]]; then
    echo "$name: downloading genome $acc"
    curl -fsSL -o "$GEN/$name.zip" \
      "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
    python3 -c 'import sys, zipfile; z = zipfile.ZipFile(sys.argv[1]); open(sys.argv[2], "wb").write(b"".join(z.read(n) for n in z.namelist() if n.endswith(".fna")))' \
      "$GEN/$name.zip" "$fa"
    rm "$GEN/$name.zip"
  fi
  echo "$name: simulating $DEPTH reads"
  "$CONDA" run -p "$ENV" badread simulate --reference "$fa" --quantity "$DEPTH" --seed 1 \
    --error_model nanopore2023 --qscore_model nanopore2023 | gzip >"$OUT/$name.fastq.gz.part"
  mv "$OUT/$name.fastq.gz.part" "$OUT/$name.fastq.gz"
done <<<"$SIMULATED"

while read -r name run; do
  [[ -n "$name" ]] || continue
  [[ -s "$OUT/$name.fastq.gz" ]] && { echo "$name: present"; continue; }
  url=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2 | cut -d';' -f1)
  echo "$name: downloading $run"
  curl -fSL -o "$OUT/$name.fastq.gz.part" "https://$url"
  mv "$OUT/$name.fastq.gz.part" "$OUT/$name.fastq.gz"
done <<<"$REAL"
echo "Done: $(ls "$OUT"/*.fastq.gz | wc -l) control read sets in $OUT"
