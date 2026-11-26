#!/bin/bash
#SBATCH --job-name=tree_isolates --cpus-per-task=16 --mem=32G --time=12:00:00
#SBATCH --output=/projects/students/bio-26-BT-7-1/Fantastic-drugs-and-where-to-find-them/3_results/slurm_logs/tree_isolates_%j.out
# Phylogenetic tree (GToTree, 74 bacterial single-copy genes) of all isolate genomes, the ISO1-2 split
# (3_results/split) and the ISO5-1 metagenome bins (3_results/meta), with the closest GTDB reference of each.
#  - A genome that was split is replaced by its parts: X is left out when X_<something> exists
#    (ISO1-2 -> ISO1-2_A, ISO1-2_B; ISO5-1 -> ISO5-1_bin.1, ...).
#  - Only genomes with >= 50 % completeness and <= 10 % contamination (CheckM2; MIMAG medium quality) are used.
# Run on BioCloud:  sbatch setup/tree_isolates.sh
# Output: 3_results/tree/tree/tree.tre and 3_results/tree/labels.tsv
#         -> copy to R_code/data/ as isolate_tree.tre and isolate_tree_labels.tsv
set -euo pipefail
cd /projects/students/bio-26-BT-7-1/Fantastic-drugs-and-where-to-find-them
OUT=3_results/tree; MIN_COMPL=50; MAX_CONTAM=10
# outgroup to root the tree: Thermotoga maritima MSB8 (phylum Thermotogota, not among the isolates)
OUTGROUP=GCF_000008545.1
OUTGROUP_TAX="d__Bacteria;p__Thermotogota;c__Thermotogae;o__Thermotogales;f__Thermotogaceae;g__Thermotoga;s__Thermotoga maritima"
rm -rf "$OUT"; mkdir -p "$OUT/genomes"

# 1. all finished genomes + their CheckM2 values, from every results folder that has a summary
for R in 3_results 3_results/split 3_results/meta; do
  S="$R/_combined/summary/genome_summary.tsv"; [[ -f "$S" ]] || { echo "skip $R (no summary)"; continue; }
  awk -F'\t' -v R="$R" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
                         {print R"\t"$h["isolate"]"\t"$h["completeness"]"\t"$h["contamination"]}' "$S"
done | sort -u -k2,2 > "$OUT/all.tsv"

# 2. keep: not replaced by split parts, and of sufficient quality; link the genome as <name>.fasta
#    (dots in names become "_", because GToTree cuts the label at the file extension)
: > "$OUT/genomes.txt"
while IFS=$'\t' read -r R name compl contam; do
  if cut -f2 "$OUT/all.tsv" | grep -q "^${name}_"; then echo "replaced by its parts: $name"; continue; fi
  if ! awk -v c="$compl" -v k="$contam" -v mc=$MIN_COMPL -v mk=$MAX_CONTAM 'BEGIN{exit !(c+0>=mc && k+0<=mk)}'; then
    echo "low quality, left out: $name ($compl % complete, $contam % contamination)"; continue; fi
  g=""; for f in "$R/$name/03_flye/$name.fasta" "$R/$name/03_genome/$name.fasta"; do [[ -s "$f" || -L "$f" ]] && { g=$f; break; }; done
  [[ -n "$g" ]] || { echo "no genome file: $name"; continue; }
  ln -sf "$PWD/$g" "$OUT/genomes/${name//./_}.fasta"; echo "$PWD/$OUT/genomes/${name//./_}.fasta" >> "$OUT/genomes.txt"
done < "$OUT/all.tsv"
echo "$(wc -l < "$OUT/genomes.txt") genomes in the tree"

# 3. closest GTDB reference of each genome + taxonomy labels (genome or accession -> GTDB taxonomy)
cat 3_results/_combined/gtdbtk/*.summary.tsv 3_results/{split,meta}/_combined/gtdbtk/*.summary.tsv 2>/dev/null | awk -F'\t' -v OFS='\t' '
  $1=="user_genome"{delete h; for(i=1;i<=NF;i++)h[$i]=i; next}
  { r = ("closest_genome_reference" in h) ? h["closest_genome_reference"] : h["fastani_reference"]
    t = ("closest_genome_taxonomy"  in h) ? h["closest_genome_taxonomy"]  : h["fastani_taxonomy"]
    n = $h["user_genome"]; gsub(/\./, "_", n); print n, $h["classification"]
    if (r && $r ~ /^GC[AF]_/) print $r, $t }' | sort -u > "$OUT/labels.tsv"
printf '%s\t%s\n' "$OUTGROUP" "$OUTGROUP_TAX" >> "$OUT/labels.tsv"
awk -F'\t' '$1 ~ /^GC[AF]_/{print $1}' "$OUT/labels.tsv" | sort -u > "$OUT/references.txt"

# 4. GToTree with the references + outgroup (needs internet to download them); without them if that fails
conda run --no-capture-output -p setup/envs/gtotree GToTree -f "$OUT/genomes.txt" -a "$OUT/references.txt" \
    -H Bacteria -j 16 -o "$OUT/tree" -F \
  || { echo "references could not be downloaded: tree of the isolates only"
       conda run --no-capture-output -p setup/envs/gtotree GToTree -f "$OUT/genomes.txt" -H Bacteria -j 16 -o "$OUT/tree" -F; }
ls -la "$OUT/tree/tree.tre"
