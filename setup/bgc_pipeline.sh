#!/usr/bin/env bash
# =============================================================================
#  bgc_pipeline.sh — ONT genomes of soil isolates -> screen for NOVEL antimicrobial BGCs
# =============================================================================
#
#  Author        Mikkel Bennetsen (GitHub: @Mikkelbenn2)
#                Aalborg University (AAU), 7th-semester project (P7), 2026
#  Repository    https://github.com/Mikkelbenn2/Fantastic-drugs-and-where-to-find-them
#  Version       3.4.0
#  Created       2026-10-06
#  Last updated  2026-11-06
#  License       MIT (see LICENSE)
#
#  AI assistance
#    This pipeline was developed with help from Claude (Anthropic).
#    It was used to draft and test the code, helped select and compare tools,
#    and helped write the documentation. The author decided the design, reviewed the code, and 
#    tweaked the parameters.
#
#  How to cite
#    Bennetsen M. (2026). bgc_pipeline.sh v3.4.0 — a Nanopore pipeline for
#    screening soil isolates for novel biosynthetic gene clusters.
#    https://github.com/Mikkelbenn2/Fantastic-drugs-and-where-to-find-them
#    Please also cite the individual programs it runs (see README.md).
#
#  Changelog
#    3.4.0  2026-11-07  mmlong2 removed (it belongs to metagenome_pipeline.sh): reads -> Filtlong ->
#                       Flye is the only assembly route. A finished genome (.fasta/.fna/.fa, may be .gz) in raw_data/ is used as is:
#                       read QC and assembly are skipped (03_genome). Used for the controls,
#                       which are now the real published NCBI genomes instead of simulated reads.
#    3.3.0  2026-11-06  Two more BGC finders: DeepBGC (13_deepbgc) and BGC-Prophet (14_bgc_prophet).
#                       Every tool's BGCs are compared with MIBiG the same way (15_mibig: DIAMOND
#                       against antiSMASH's KnownClusterBlast data, same cut-offs). New summaries:
#                       per_tool/<tool>_bgcs.tsv and master_summary.tsv (all tools merged per
#                       locus + QC + taxonomy + ranking). No polishing any more (Medaka removed):
#                       R10.4.1 + Dorado sup reads are accurate enough on their own.
#    3.2.1  2026-10-26  Brings back fixes from the 3.1.x line that 3.2.0 lost: GTDB-Tk >= 2.7
#                       without --skip_ani_screen (gone in 2.7), RGI threads as -n. BiG-SCAPE env
#                       pinned to Python 3.11 + SQLAlchemy 2.0.2.
#    3.2.0  2026-10-14  Setup links all BioCloud databases into setup/databases/biocloud/, uses a
#                       shared antiSMASH database if it fits, and only downloads what is missing;
#                       prints which database each program uses.
#    3.1.0  2026-10-14  ASSEMBLER setting: "mmlong2" (default) runs mmlong2-lite for assembly +
#                       Medaka polishing + curation in one step; "flye" keeps the old route.
#    3.0.6  2026-10-12  QUAST env pinned to Python 3.11 (QUAST needs distutils, gone in 3.12+).
#    3.0.5  2026-10-12  Bakta: skips databases too old for Bakta >= 1.10 (no "type" key in
#                       version.json) and downloads its own copy in setup instead.
#    3.0.4  2026-10-12  Medaka: MEDAKA_MODEL setting; falls back to the bacterial R10.4.1 model
#                       when --bacteria fails (reads from SRA have no model info in headers).
#    3.0.3  2026-10-12  Comment listing the AAU BioCloud database folders and GTDB releases.
#    3.0.2  2026-10-12  Uses the shared AAU server databases (/databases/GTDB, /databases/bakta,
#                       /databases/checkm2) by default; finds the database folder inside them
#    3.0.1  2026-10-12  Setup: SHARED_*_DB settings to reuse server databases; finds GTDB-Tk
#                       data inside the conda env; removes half-downloaded Bakta files
#    3.0.0  2026-10-12  Single self-contained script; GToTree, GECCO, BiG-SCAPE and
#                       ranked novelty summary; Slurm job array; README + flowchart
#    2.0.0  2026-10-06  Multi-file version (config + sample sheet)
#    1.0.0  2026-09-11  First version: NanoPlot -> Flye -> antiSMASH -> RGI -> CARD
#
#  Usage
#    bash setup/bgc_pipeline.sh setup | list | run [--dry-run] | submit | help
#    Full instructions: README.md
# =============================================================================

# ============================== SETTINGS =====================================
# Change these to fit your data and server. Everything below SETTINGS can stay.

THREADS=16                     # CPUs per program (inside Slurm: taken from the job)

# Reads -> genome: Filtlong -> Flye (no polishing). mmlong2 is only used by metagenome_pipeline.sh.

# Reads. "nano-hq" = R10.4.1 + Dorado hac/sup (Q15-20+). Use "nano-raw" for old R9.4.1 data.
# Without polishing, R10.4.1 reads should be basecalled with Dorado "sup" (Q20+); older/noisier
# reads give frameshifts in long NRPS/PKS genes.
FLYE_READ_MODE="nano-hq"
GENOME_SIZE=""                 # e.g. "8m" for Streptomyces; empty = Flye estimates it
FILTLONG_MIN_LEN=1000          # drop reads shorter than this (bp)
FILTLONG_KEEP_PERCENT=95       # keep the best 95 % of bases
FILTLONG_TARGET_BASES=""       # optional cap, e.g. 800000000 (= 100x of an 8 Mb genome)

BAKTA_DB_TYPE="light"          # "light" (~1.5 GB) or "full" (~32 GB download; only if the internet is fast)

# Databases that ALREADY EXIST on the server (shared copies), so "setup" does not have to
# download them. The script searches inside each folder for the actual database, so pointing
# at the top folder is enough. If a path does not exist, the script falls back to its own
# downloaded copy. On another server: change the paths, or set "" to let "setup" download it.
#
# ---- AAU BioCloud server ----------------------------------------------------------------
# If you are on the AAU BioCloud server, these database folders already exist and can be used
# (the defaults below already point to them):
#   Bakta     /databases/bakta   (NOTE: the 2022 database here is too old for Bakta >= 1.10;
#                                  the script skips it and "setup" downloads its own copy)
#   CheckM2   /databases/checkm2
#             /databases/checkm2/CheckM2_database          (the .dmnd file is in here)
#   GTDB      /databases/GTDB
#             /databases/GTDB/gtdb_data
#             /databases/GTDB/gtdbtk_packages               (one folder per GTDB release):
#               GTDB-TK_release207_2022_05_10   -> GTDB-Tk 2.1.x
#               GTDB-TK_release214_2023_04_27   -> GTDB-Tk 2.3.x
#               GTDB-TK_release220_2024_04_18   -> GTDB-Tk 2.4.x
#               GTDB-TK_release226_2025_04_11   -> GTDB-Tk 2.5.x
#               GTDB-TK_release232_2026_04_15   -> GTDB-Tk 2.7.x  (used by this pipeline)
# The GTDB release MUST match the GTDB-Tk version, otherwise GTDB-Tk stops with an error
# (version pairs from the GTDB-Tk docs: ecogenomics.github.io/GTDBTk).
# Check your version with:  conda run -p setup/envs/gtdbtk gtdbtk --version
# -------------------------------------------------------------------------------------------
SHARED_GTDBTK_DB="/databases/GTDB/gtdbtk_packages/GTDB-TK_release232_2026_04_15"
                               # GTDB-Tk 2.7.x needs release 232; an older GTDB-Tk needs an older release
SHARED_BAKTA_DB="/databases/bakta"
SHARED_CHECKM2_DB="/databases/checkm2"
SHARED_ANTISMASH_DB="/databases/AntiSMASH_databases"
                               # only used if a version there fits our antiSMASH (8.x); the
                               # AntiSMASH-5.1.0 folder does not, so "setup" downloads its own
# Not on BioCloud, so "setup" always downloads it: CARD (for RGI). BioCloud has ResFinder and
# deepARG, but RGI needs CARD.

# "setup" puts a link to EVERY folder in this server database folder into
# setup/databases/biocloud/, so all BioCloud databases are easy to find next to the pipeline's own.
# Empty "" = not on BioCloud, make no links.
SHARED_DB_ROOT="/databases"

# Optional programs (1 = run, 0 = skip)
RUN_GECCO=1                    # machine-learning BGC detection
RUN_DEEPBGC=1                  # deep-learning BGC detection + product class/activity prediction
RUN_BGCPROPHET=1               # protein-language-model BGC detection + product class
RUN_RGI=1                     # resistance genes (CARD)
RUN_GTDBTK=1                   # species ID — needs ~110 GB database + ~64-100 GB RAM
RUN_GTOTREE=1                  # phylogeny of the isolates (+ closest GTDB references)
RUN_BIGSCAPE=1                 # BGC families compared with MIBiG

GTOTREE_HMM="Bacteria"         # marker set; "Actinobacteria" if all isolates are actinomycetes

# antiSMASH extra analyses. --cb-knownclusters (MIBiG comparison) is ESSENTIAL for novelty.
ANTISMASH_EXTRA="--cb-knownclusters --cb-subclusters --cc-mibig --asf --rre --tfbs --pfam2go --clusterhmmer --tigrfam --smcog-trees"

DEEPBGC_SCORE=0.5              # DeepBGC: minimum BGC score per gene (its default)
BGCPROPHET_THRESHOLD=0.5       # BGC-Prophet: gene score cut-off (value used in its README)
BGCPROPHET_CLASSIFY_T=0.5      # BGC-Prophet: class score cut-off (README value)
BGCPROPHET_DEVICE="cpu"        # "cuda" if the job has a GPU (much faster)
BGCPROPHET_MODEL_URL="https://github.com/HUST-NingKang-Lab/BGC-Prophet/files/12733164/model.tar.gz"

BIGSCAPE_MIBIG_VERSION="4.0"
BIGSCAPE_CUTOFF="0.30"
KNOWN_SIMILARITY=80            # KnownClusterBlast % at/above which a BGC counts as "known"
RELATED_SIMILARITY=30          # below this (or no hit) = "novel_candidate"

# Slurm (only used by "submit")
SLURM_PARTITION=""             # empty = cluster default
SLURM_ACCOUNT=""
SLURM_SAMPLE_CPUS=16;  SLURM_SAMPLE_MEM="64G";   SLURM_SAMPLE_TIME="1-00:00:00"
SLURM_COMBINED_CPUS=32; SLURM_COMBINED_MEM="128G"; SLURM_COMBINED_TIME="1-00:00:00"

# Program versions (minimum). Exact installed versions are logged in results/software_versions.tsv
declare -A ENV_PACKAGES=(
  [reads]="nanoplot>=1.42 filtlong>=0.2.1 seqkit>=2.8"
  [flye]="flye>=2.9.5"
  [quast]="quast>=5.2 python=3.11"   # QUAST still imports distutils, which Python 3.12+ removed
  [checkm2]="checkm2>=1.0.2"
  [bakta]="bakta>=1.10"
  [antismash]="antismash>=8.0"
  [gecco]="gecco>=0.9.10"
  [deepbgc]="deepbgc>=0.1.31 python=3.7"   # needs TensorFlow 1.x, which only exists up to Python 3.7
  [bgcprophet]="bgc-prophet>=0.1.2"
  [rgi]="rgi>=6.0"
  [gtdbtk]="gtdbtk>=2.7"         # 2.7.x/2.8.x use GTDB release 232
  [gtotree]="gtotree>=1.8"
  [bigscape]="bigscape>=2.0 python=3.11 sqlalchemy=2.0.2"
)
# =========================== END OF SETTINGS =================================

PIPELINE_VERSION="3.4.0"       # keep in sync with the header at the top
PIPELINE_UPDATED="2026-11-07"


set -uo pipefail
# -u: unset variable = error (catches typos). pipefail: a pipe fails if any part fails.
# No -e on purpose: one failing isolate must not stop the others.

# ---------- folders (all relative to this script) ----------------------------
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SETUP_DIR="$(dirname "$SCRIPT_PATH")"
PROJECT_DIR="$(dirname "$SETUP_DIR")"
RAW_DIR="$PROJECT_DIR/raw_data"
RESULTS_DIR="$PROJECT_DIR/results"
ENV_DIR="$SETUP_DIR/envs"
DB_DIR="$SETUP_DIR/databases"
COMBINED_DIR="$RESULTS_DIR/_combined"

ANTISMASH_DB="$DB_DIR/antismash"
BAKTA_DB_PARENT="$DB_DIR/bakta"
CHECKM2_DB="$DB_DIR/checkm2"
CARD_DB="$DB_DIR/card"
GTDBTK_DB="$DB_DIR/gtdbtk"
DEEPBGC_DB="$DB_DIR/deepbgc"           # DeepBGC models + its Pfam (DEEPBGC_DOWNLOADS_DIR)
BGCPROPHET_DB="$DB_DIR/bgc_prophet"    # BGC-Prophet models + ESM-2 weights (TORCH_HOME)

# ---------- locating the databases (shared server copy first, then our own) ----------
# antiSMASH: setup tests the shared copies with "antismash --check-prereqs" (it checks that the
# database versions match the installed antiSMASH) and writes the one to use into this file.
ANTISMASH_DB_CHOICE="$DB_DIR/antismash_db_in_use.txt"
antismash_db_path() {
  local d
  if [[ -s "$ANTISMASH_DB_CHOICE" ]]; then
    d=$(head -n 1 "$ANTISMASH_DB_CHOICE"); [[ -d "$d" ]] && { echo "$d"; return; }
  fi
  echo "$ANTISMASH_DB"
}
antismash_db_ok() {   # true if antiSMASH accepts this database folder
  [[ -d "$1/clusterblast" ]] || return 1
  "$CONDA_RUN" run -p "$ENV_DIR/antismash" antismash --check-prereqs --databases "$1" >>"$MAIN_LOG" 2>&1
}

# find_db_dir TOP NAME : the folder below TOP that contains NAME (a file or folder that
# marks the database), e.g. bakta.db or taxonomy/. Prints nothing if not found.
find_db_dir() {
  [[ -n "$1" && -e "$1" ]] || return 0
  local hit
  hit=$(find -L "$1" -maxdepth 4 -name "$2" 2>/dev/null | sort | head -n 1)
  [[ -n "$hit" ]] && dirname "$hit"
  return 0
}

# Bakta: the folder that contains bakta.db. Only databases that fit Bakta >= 1.10 are used
# (schema v5/v6: their version.json has a "type" key). Older ones, e.g. the 2022 database
# (schema v4), stop Bakta with "Missed key 'type'", so they are skipped.
# Among the usable ones: full ("db") before light ("db-light"), then the newest folder.
bakta_db_ok() { [[ -f "$1/bakta.db" && -f "$1/version.json" ]] && grep -q '"type"' "$1/version.json"; }
bakta_db_path() {
  local top d best
  for top in "$SHARED_BAKTA_DB" "$BAKTA_DB_PARENT"; do
    [[ -n "$top" && -e "$top" ]] || continue
    best=$(find -L "$top" -maxdepth 4 -name bakta.db 2>/dev/null | while read -r d; do
             d=$(dirname "$d"); bakta_db_ok "$d" || continue
             printf '%s\t%s\n' "$([[ $(basename "$d") == db ]] && echo 0 || echo 1)" "$d"
           done | sort -k1,1n -k2,2r | head -n 1 | cut -f2)
    [[ -n "$best" ]] && { echo "$best"; return; }
  done
  echo "$BAKTA_DB_PARENT/db-$BAKTA_DB_TYPE"
}

# GTDB-Tk: the folder that contains taxonomy/. Also checks inside the conda env, because
# newer GTDB-Tk versions' download-db.sh unpacks there instead of the folder it is given.
gtdbtk_db_path() {
  local d top
  for top in "$SHARED_GTDBTK_DB" "$GTDBTK_DB" "$ENV_DIR/gtdbtk/share"; do
    d=$(find_db_dir "$top" taxonomy)
    [[ -n "$d" ]] && { echo "$d"; return; }
  done
  echo "$GTDBTK_DB"
}

# CheckM2: the .dmnd file
checkm2_dmnd() {
  local top hit
  for top in "$SHARED_CHECKM2_DB" "$CHECKM2_DB"; do
    [[ -n "$top" && -e "$top" ]] || continue
    hit=$(find -L "$top" -maxdepth 4 -name '*.dmnd' 2>/dev/null | sort | head -n 1)
    [[ -n "$hit" ]] && { echo "$hit"; return; }
  done
}

# MIBiG data that antiSMASH's KnownClusterBlast uses (proteins.dmnd + clusters.txt, newest
# version folder). All tools' BGCs are compared against this, so they share one MIBiG version.
kcb_dir() {
  local hit
  hit=$(find -L "$(antismash_db_path)/knownclusterblast" -maxdepth 2 -name proteins.dmnd 2>/dev/null | sort -V | tail -n 1)
  [[ -n "$hit" ]] && dirname "$hit"
  return 0
}
# BGC-Prophet: the folder with annotator.pt + classifier.pt (its README spells it "annoator.pt")
bgcprophet_model_dir() {
  local hit
  hit=$(find -L "$BGCPROPHET_DB" -maxdepth 4 -name 'anno*tor.pt' 2>/dev/null | sort | head -n 1)
  [[ -n "$hit" ]] && dirname "$hit"
  return 0
}

[[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && THREADS=$SLURM_CPUS_PER_TASK
DRY_RUN=0
FORCE_STEPS=()
MAIN_LOG=/dev/null
STEP_LOG=/dev/null

# ---------- small helpers ----------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MAIN_LOG" >&2; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR: $*"; exit 1; }

find_conda() {
  # conda, mamba and micromamba all understand: <tool> run -p <env> <command>
  CONDA_RUN=""; CONDA_CREATE=""
  for c in conda mamba micromamba; do command -v "$c" &>/dev/null && { CONDA_RUN=$c; break; }; done
  for c in mamba micromamba conda; do command -v "$c" &>/dev/null && { CONDA_CREATE=$c; break; }; done
  if [[ -z "$CONDA_RUN" ]]; then
    [[ $DRY_RUN -eq 1 ]] && { CONDA_RUN=conda; return; }
    die "conda/mamba not found. Try 'module load miniforge' (or install Miniforge) first."
  fi
}

# run_in_env ENV command args... : run a program inside its conda env, output -> step log
run_in_env() {
  local env=$1; shift
  if [[ $DRY_RUN -eq 1 ]]; then echo "    [dry-run] ($env) $(printf '%q ' "$@")" >&2; return 0; fi
  echo "### $(date '+%F %T') ($env) $(printf '%q ' "$@")" >>"$STEP_LOG"
  "$CONDA_RUN" run --no-capture-output -p "$ENV_DIR/$env" "$@" >>"$STEP_LOG" 2>&1
}

# run_plain command args... : same, for commands that need no conda env
run_plain() {
  if [[ $DRY_RUN -eq 1 ]]; then echo "    [dry-run] $(printf '%q ' "$@")" >&2; return 0; fi
  echo "### $(date '+%F %T') $(printf '%q ' "$@")" >>"$STEP_LOG"
  "$@" >>"$STEP_LOG" 2>&1
}

# check that a program produced its main output file
need() {
  [[ $DRY_RUN -eq 1 ]] && return 0
  [[ -s "$1" ]] || { log "    missing or empty: $1"; return 1; }
}

is_forced() { local f; for f in "${FORCE_STEPS[@]:-}"; do [[ "$f" == all || "$f" == "$1" ]] && return 0; done; return 1; }

# run_step STEP OUTPUT_FOLDER LOG_FOLDER FUNCTION
#   skips the step if LOG_FOLDER/STEP.done exists (= finished on an earlier run),
#   otherwise empties OUTPUT_FOLDER, runs FUNCTION and writes .done only on success.
run_step() {
  local step=$1 out=$2 logs=$3 func=$4 label=${5:-}
  if [[ -f "$logs/$step.done" ]] && ! is_forced "$step"; then
    log "  [$label] $step: done earlier, skipping"; return 0
  fi
  log "  [$label] $step: running"
  if [[ $DRY_RUN -eq 0 ]]; then
    rm -rf "$out"; mkdir -p "$out" "$logs"; rm -f "$logs/$step.done"
  fi
  STEP_LOG="$logs/$step.log"; [[ $DRY_RUN -eq 0 ]] && : >"$STEP_LOG"
  local t0=$SECONDS
  if "$func"; then
    [[ $DRY_RUN -eq 0 ]] && date '+%F %T' >"$logs/$step.done"
    log "  [$label] $step: ok ($(( (SECONDS - t0) / 60 )) min)"
  else
    log "  [$label] $step: FAILED — see $STEP_LOG"; return 1
  fi
}

# ---------- find the isolates in raw_data/ ------------------------------------
# Each FASTQ file = one isolate (name = file name without .fastq.gz etc.)
# Each folder containing FASTQ files = one isolate (name = folder name)
# Each FASTA file (.fasta/.fna/.fa, may be .gz) = a finished genome: no read QC or assembly
find_samples() {
  [[ -d "$RAW_DIR" ]] || die "No raw_data folder found at $RAW_DIR"
  SAMPLES=(); SAMPLE_SRC=()
  local p name
  while IFS= read -r p; do
    name=$(basename "$p")
    if [[ -d "$p" ]]; then
      compgen -G "$p/*.f*q*" >/dev/null || continue           # folder without FASTQ: ignore
    else
      [[ "$name" =~ \.(fastq|fq|fasta|fna|fa)(\.gz)?$ ]] || continue
      name=${name%.gz}; name=${name%.*}
    fi
    name=${name//[^A-Za-z0-9._-]/_}                            # no spaces/odd characters
    SAMPLES+=("$name"); SAMPLE_SRC+=("$p")
  done < <(find "$RAW_DIR" -mindepth 1 -maxdepth 1 \( -type f -o -type d -o -type l \) | sort)
  [[ ${#SAMPLES[@]} -gt 0 ]] || die "No FASTQ files, FASTQ folders or genome FASTA files found in $RAW_DIR"
}

# Where the files of one isolate live (one place, used by every step)
set_sample() {
  S=$1; SRC=$2
  SD="$RESULTS_DIR/$S"; SL="$SD/logs"
  READS="$SD/01_nanoplot/$S.reads.fastq.gz"   # all reads joined
  FILT="$SD/02_filtlong/$S.filtered.fastq.gz"
  DRAFT="$SD/03_flye/assembly.fasta"
  GENOME="$SD/03_flye/$S.fasta"            # final genome used by all later steps
  if is_genome_input; then DRAFT="$SD/03_genome/$S.fasta"; GENOME="$DRAFT"; fi
  GBFF="$SD/07_bakta/$S.gbff"
}
is_genome_input() { [[ -f "$SRC" && "$SRC" =~ \.(fasta|fna|fa)(\.gz)?$ ]]; }

# 03 genome — a finished genome was given instead of reads: copy it (unzipped) as the genome.
step_genome() {
  [[ $DRY_RUN -eq 1 ]] && return 0
  if [[ "$SRC" == *.gz ]]; then gzip -dc "$SRC" >"$GENOME"; else cp "$SRC" "$GENOME"; fi
  need "$GENOME"
}

# =============================================================================
#  PER-ISOLATE STEPS
# =============================================================================

# Join the reads of one isolate into one file (MinKNOW writes many small files per barcode).
# gzip files may simply be concatenated with cat.
join_reads() {
  [[ $DRY_RUN -eq 1 ]] && return 0
  mkdir -p "$(dirname "$READS")"
  if [[ -d "$SRC" ]]; then
    : >"$READS"
    local f
    for f in "$SRC"/*.f*q*; do
      case "$f" in
        *.gz)         cat "$f" >>"$READS" ;;
        *.fastq|*.fq) gzip -c "$f" >>"$READS" ;;
      esac
    done
  elif [[ "$SRC" == *.gz ]]; then ln -sf "$SRC" "$READS"
  else gzip -c "$SRC" >"$READS"
  fi
  need "$READS"
}

# 01 NanoPlot — read length + quality. Check: read N50 (>10 kb is good), mean quality
# (matches FLYE_READ_MODE?), total bases / genome size = coverage (aim >= 40x).
step_nanoplot() {
  join_reads || return 1
  run_in_env reads NanoPlot --fastq "$READS" --outdir "$SD/01_nanoplot" \
    --threads "$THREADS" --tsv_stats --loglength
}

# 02 Filtlong — remove short and low-quality reads so Flye gets cleaner input.
step_filtlong() {
  [[ $DRY_RUN -eq 1 || -s "$READS" ]] || join_reads || return 1
  local args=(--min_length "$FILTLONG_MIN_LEN" --keep_percent "$FILTLONG_KEEP_PERCENT")
  [[ -n "$FILTLONG_TARGET_BASES" ]] && args+=(--target_bases "$FILTLONG_TARGET_BASES")
  # Filtlong writes to the screen, so a pipe to gzip is needed -> run it via bash -c
  run_in_env reads bash -c 'set -o pipefail; out=$1; shift; filtlong "$@" | gzip -c > "$out"' \
    _ "$FILT" "${args[@]}" "$READS" || return 1
  need "$FILT" || return 1
  run_in_env reads bash -c 'seqkit stats -a -T -j "$1" "$2" "$3" > "$4"' \
    _ "$THREADS" "$READS" "$FILT" "$SD/02_filtlong/read_stats.tsv"
}

# 03 Flye — de novo assembly. Ideal: one chromosome contig (LINEAR in Streptomyces,
# that is normal) + plasmids. assembly_info.txt: length, depth and circular Y/N per contig.
step_flye() {
  local args=("--$FLYE_READ_MODE" "$FILT" --out-dir "$SD/03_flye" --threads "$THREADS")
  [[ -n "$GENOME_SIZE" ]] && args+=(--genome-size "$GENOME_SIZE")
  run_in_env flye flye "${args[@]}" || return 1
  need "$DRAFT" || return 1
  run_plain cp "$DRAFT" "$GENOME" || return 1     # no polishing: the Flye assembly is the genome
  need "$GENOME"
}

# 05 QUAST — contiguity: number of contigs, total length, N50, GC% (draft = final: no polishing).
step_quast() {
  run_in_env quast quast.py -o "$SD/05_quast" -t "$THREADS" \
    --labels "draft,final" "$DRAFT" "$GENOME"
}

# 06 CheckM2 — completeness (>95 % good) and contamination (<5 % good; high = mixed culture?)
step_checkm2() {
  local dmnd
  dmnd=$(checkm2_dmnd)
  [[ -n "$dmnd" || $DRY_RUN -eq 1 ]] || { log "    CheckM2 database missing — run: setup"; return 1; }
  run_in_env checkm2 checkm2 predict --threads "$THREADS" --input "$GENOME" \
    --output-directory "$SD/06_checkm2" --database_path "${dmnd:-CheckM2.dmnd}" --force
}

# 07 Bakta — finds and names all genes (successor of Prokka, which is no longer maintained).
# Its GenBank file (.gbff) is the input for antiSMASH, so all tools use the same genes.
step_bakta() {
  run_in_env bakta bakta --db "$(bakta_db_path)" --output "$SD/07_bakta" \
    --prefix "$S" --threads "$THREADS" --keep-contig-headers --skip-plot --force "$GENOME" \
    || return 1
  need "$GBFF"
}

# 08 antiSMASH — rule-based BGC detection (NRPS, PKS, RiPP, terpene, ...) and comparison
# with MIBiG (characterised BGCs). Open 08_antismash/index.html in a browser.
# If antiSMASH does not accept the Bakta file, it retries on the FASTA with Prodigal.
step_antismash() {
  local out="$SD/08_antismash"
  # shellcheck disable=SC2206
  local extra=($ANTISMASH_EXTRA)
  local common=(--cpus "$THREADS" --taxon bacteria --databases "$(antismash_db_path)"
                --output-basename "$S" --output-dir "$out" "${extra[@]}")
  if ! run_in_env antismash antismash "${common[@]}" --genefinding-tool none "$GBFF"; then
    warn "antiSMASH could not use the Bakta GenBank for $S — retrying on FASTA + Prodigal"
    [[ $DRY_RUN -eq 0 ]] && rm -rf "$out"
    run_in_env antismash antismash "${common[@]}" --genefinding-tool prodigal "$GENOME" || return 1
  fi
  need "$out/$S.json"
}

# 09 GECCO — machine-learning BGC detection (no hand-written rules). Found by both tools =
# well supported; found only by GECCO = possibly an unusual BGC worth a manual look.
step_gecco() {
  run_in_env gecco gecco run --genome "$GENOME" --output-dir "$SD/09_gecco" --jobs "$THREADS"
}

# 10 RGI/CARD — known resistance genes. A resistance gene INSIDE a BGC is often the
# producer's self-resistance gene -> hints that the product is an antibiotic + its target.
step_rgi() {
  local out="$SD/10_rgi_card"
  run_plain ln -sfn "$CARD_DB/localDB" "$out/localDB" || return 1     # rgi --local looks in ./localDB
  ( [[ $DRY_RUN -eq 1 ]] || cd "$out" || exit 1
    run_in_env rgi rgi main --input_sequence "$GENOME" --output_file "$out/${S}_rgi" \
      --input_type contig --alignment_tool DIAMOND -n "$THREADS" --local --clean )
}

# 13 DeepBGC — deep-learning (BiLSTM on Pfam domains) BGC detection; also predicts the product
# class and ACTIVITY (antibacterial, antifungal, cytotoxic, inhibitor) of each BGC. Finds its own
# genes (Prodigal) on the genome FASTA. Known to call more (and more fragmented) BGCs than
# antiSMASH, so a DeepBGC-only hit is weak evidence on its own.
step_deepbgc() {
  local out="$SD/13_deepbgc"
  run_in_env deepbgc bash -c 'export DEEPBGC_DOWNLOADS_DIR="$1"; shift; deepbgc "$@"' _ "$DEEPBGC_DB" \
    pipeline "$GENOME" -o "$out" --score "$DEEPBGC_SCORE" || return 1
  need "$out/13_deepbgc.bgc.tsv"      # DeepBGC names its files after the output folder
}

# 14 BGC-Prophet — transformer on ESM-2 protein embeddings; finds BGC genes and the product class.
# It reads proteins, so it gets Bakta's proteins (same genes as antiSMASH). Writes <isolate>_classified.csv.
step_bgcprophet() {
  local out="$SD/14_bgc_prophet" m
  m=$(bgcprophet_model_dir)
  [[ -n "$m" || $DRY_RUN -eq 1 ]] || { log "    BGC-Prophet models missing — run: setup"; return 1; }
  run_plain mkdir -p "$out/input" || return 1
  run_plain ln -sf "$SD/07_bakta/$S.faa" "$out/input/$S.faa" || return 1
  run_in_env bgcprophet bash -c 'export TORCH_HOME="$1"; cd "$2" && shift 2 && bgc_prophet "$@"' _ \
    "$BGCPROPHET_DB/torch" "$out" pipeline --genomesDir "$out/input" --outputPath "$out" --name "$S" \
    --modelPath "${m:-models}/annotator.pt" --classifierPath "${m:-models}/classifier.pt" \
    --threshold "$BGCPROPHET_THRESHOLD" --max_gap 3 --min_count 2 --classify_t "$BGCPROPHET_CLASSIFY_T" \
    --device "$BGCPROPHET_DEVICE" --threads "$THREADS" || return 1
  need "$out/${S}_classified.csv"
}

# 15 MIBiG comparison for EVERY tool: all Bakta proteins vs. the MIBiG proteins of antiSMASH's
# KnownClusterBlast (DIAMOND, as antiSMASH does). The summary then scores each predicted BGC
# like KnownClusterBlast (% of a MIBiG cluster's genes hit at > 45 % identity, > 40 % coverage).
step_mibig() {
  local d
  d=$(kcb_dir)
  [[ -n "$d" || $DRY_RUN -eq 1 ]] || { log "    no knownclusterblast/*/proteins.dmnd in $(antismash_db_path)"; return 1; }
  run_in_env antismash diamond blastp --query "$SD/07_bakta/$S.faa" --db "${d:-kcb}/proteins.dmnd" \
    --out "$SD/15_mibig/${S}_vs_mibig.tsv" --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --evalue 1e-05 --max-target-seqs 10000 --threads "$THREADS"
}

process_sample() {
  set_sample "$1" "$2"
  log "=== Isolate $S ==="
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$SD" "$SL"
  if is_genome_input; then
    run_step 03_genome    "$SD/03_genome"    "$SL" step_genome          "$S" || return 1
  else
    run_step 01_nanoplot  "$SD/01_nanoplot"  "$SL" step_nanoplot          "$S" || warn "NanoPlot failed for $S (continuing)"
    run_step 02_filtlong  "$SD/02_filtlong"  "$SL" step_filtlong          "$S" || return 1
    run_step 03_flye      "$SD/03_flye"      "$SL" step_flye              "$S" || return 1   # no polishing: no 04 folder
  fi
  run_step 05_quast     "$SD/05_quast"     "$SL" step_quast             "$S" || warn "QUAST failed for $S (continuing)"
  run_step 06_checkm2   "$SD/06_checkm2"   "$SL" step_checkm2           "$S" || warn "CheckM2 failed for $S (continuing)"
  run_step 07_bakta     "$SD/07_bakta"     "$SL" step_bakta             "$S" || return 1
  run_step 08_antismash "$SD/08_antismash" "$SL" step_antismash         "$S" || return 1
  [[ $RUN_GECCO -eq 1 ]] && { run_step 09_gecco    "$SD/09_gecco"    "$SL" step_gecco "$S" || warn "GECCO failed for $S (continuing)"; }
  [[ $RUN_RGI   -eq 1 ]] && { run_step 10_rgi_card "$SD/10_rgi_card" "$SL" step_rgi   "$S" || warn "RGI failed for $S (continuing)"; }
  [[ $RUN_DEEPBGC    -eq 1 ]] && { run_step 13_deepbgc     "$SD/13_deepbgc"     "$SL" step_deepbgc    "$S" || warn "DeepBGC failed for $S (continuing)"; }
  [[ $RUN_BGCPROPHET -eq 1 ]] && { run_step 14_bgc_prophet "$SD/14_bgc_prophet" "$SL" step_bgcprophet "$S" || warn "BGC-Prophet failed for $S (continuing)"; }
  run_step 15_mibig "$SD/15_mibig" "$SL" step_mibig "$S" || warn "MIBiG comparison failed for $S (continuing; only antiSMASH BGCs get a MIBiG hit)"
  log "=== Isolate $S finished ==="
}

# =============================================================================
#  COMBINED STEPS (all isolates together)
# =============================================================================
CL="$COMBINED_DIR/logs"

finished_genomes() {   # prints "name<TAB>genome path" for isolates whose genome exists
  local i
  for i in "${!SAMPLES[@]}"; do
    set_sample "${SAMPLES[$i]}" "${SAMPLE_SRC[$i]}"
    [[ $DRY_RUN -eq 1 || -s "$GENOME" ]] && printf '%s\t%s\n' "$S" "$GENOME"
  done
}

# GTDB-Tk — species from the genome (much finer than 16S, especially in Streptomyces).
# No species match ("s__" empty) = possibly a NEW species. Runs once for all isolates.
# GTDB-Tk >= 2.7 first compares every genome to all GTDB references with skani (the "ANI
# screen"); --skip_ani_screen no longer exists there.
step_gtdbtk() {
  local out="$COMBINED_DIR/gtdbtk" name g
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$out/genomes"
  while IFS=$'\t' read -r name g; do run_plain ln -sf "$g" "$out/genomes/$name.fasta"; done < <(finished_genomes)
  run_in_env gtdbtk bash -c 'export GTDBTK_DATA_PATH="$1"; shift; gtdbtk "$@"' _ "$(gtdbtk_db_path)" \
    classify_wf --genome_dir "$out/genomes" --extension fasta --out_dir "$out" \
    --cpus "$THREADS"
}

# GToTree — phylogenetic tree from single-copy marker genes, with the closest GTDB reference
# genome of each isolate added for context (downloaded from NCBI; if that is not possible,
# the tree is made from the isolates only).
step_gtotree() {
  local out="$COMBINED_DIR/gtotree" name g
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$out"
  : >"$out/genomes.txt" 2>/dev/null || true
  while IFS=$'\t' read -r name g; do
    [[ $DRY_RUN -eq 0 ]] && { ln -sf "$g" "$out/$name.fasta"; echo "$out/$name.fasta" >>"$out/genomes.txt"; }
  done < <(finished_genomes)
  local refs="$out/references.txt"
  if [[ $DRY_RUN -eq 0 ]]; then
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="closest_genome_reference"||$i=="fastani_reference") c=i; next}
                c && $c ~ /^GC[AF]_/ {print $c}' "$COMBINED_DIR"/gtdbtk/*.summary.tsv 2>/dev/null \
      | sort -u >"$refs"
  fi
  if [[ $DRY_RUN -eq 1 || -s "$refs" ]]; then
    run_in_env gtotree GToTree -f "$out/genomes.txt" -a "$refs" -H "$GTOTREE_HMM" \
      -j "$THREADS" -o "$out/tree" -F && return 0
    warn "GToTree with reference genomes failed (no internet on this node?) — isolates only"
  fi
  local n; n=$(wc -l <"$out/genomes.txt" 2>/dev/null || echo 0)
  [[ $DRY_RUN -eq 1 || $n -ge 3 ]] || { log "    fewer than 3 genomes: no tree possible"; return 0; }
  run_in_env gtotree GToTree -f "$out/genomes.txt" -H "$GTOTREE_HMM" -j "$THREADS" -o "$out/tree" -F
}

# BiG-SCAPE — groups all BGCs from all isolates + MIBiG into gene-cluster families.
# A family WITH a MIBiG member = resembles a known pathway; WITHOUT = no characterised relative.
step_bigscape() {
  local out="$COMBINED_DIR/bigscape" name g f n=0
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$out/input"
  while IFS=$'\t' read -r name g; do
    [[ $DRY_RUN -eq 1 ]] && continue
    for f in "$RESULTS_DIR/$name/08_antismash/"*.region*.gbk; do
      [[ -e "$f" ]] || continue
      cp "$f" "$out/input/${name}__$(basename "$f")"; n=$((n + 1))   # sample prefix avoids name clashes
    done
  done < <(finished_genomes)
  log "    $n BGC regions collected"
  [[ $DRY_RUN -eq 1 || $n -gt 0 ]] || return 0
  local pfam
  pfam=$(find -L "$(antismash_db_path)" -name 'Pfam-A.hmm' 2>/dev/null | sort | tail -n 1)   # reuse antiSMASH's Pfam
  [[ -n "$pfam" || $DRY_RUN -eq 1 ]] || { log "    Pfam-A.hmm not found in $(antismash_db_path)"; return 1; }
  run_in_env bigscape bigscape cluster -i "$out/input" -o "$out/output" -p "${pfam:-Pfam-A.hmm}" \
    --mibig-version "$BIGSCAPE_MIBIG_VERSION" --gcf-cutoffs "$BIGSCAPE_CUTOFF" \
    --include-singletons --label P7 -c "$THREADS"
}

# Summary — collects everything into ranked tables (Python code: write_summary_py below)
step_summary() {
  local names
  names=$(finished_genomes | cut -f1 | paste -sd, -)
  if [[ $DRY_RUN -eq 1 ]]; then echo "    [dry-run] python3 <summary> $names" >&2; return 0; fi
  write_summary_py "$CL/summary.py"
  python3 "$CL/summary.py" --results "$RESULTS_DIR" --samples "$names" \
    --out "$COMBINED_DIR/summary" --cutoff "$BIGSCAPE_CUTOFF" --mibig "$(kcb_dir)" \
    --known "$KNOWN_SIMILARITY" --related "$RELATED_SIMILARITY" >>"$STEP_LOG" 2>&1
}

process_combined() {
  log "=== Combined steps (all isolates) ==="
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$COMBINED_DIR" "$CL"
  [[ $RUN_GTDBTK   -eq 1 ]] && { run_step gtdbtk   "$COMBINED_DIR/gtdbtk"   "$CL" step_gtdbtk   combined || warn "GTDB-Tk failed (continuing)"; }
  [[ $RUN_GTOTREE  -eq 1 ]] && { run_step gtotree  "$COMBINED_DIR/gtotree"  "$CL" step_gtotree  combined || warn "GToTree failed (continuing)"; }
  [[ $RUN_BIGSCAPE -eq 1 ]] && { run_step bigscape "$COMBINED_DIR/bigscape" "$CL" step_bigscape combined || warn "BiG-SCAPE failed (continuing)"; }
  FORCE_STEPS+=(summary)                       # always rebuild the summary
  run_step summary "$COMBINED_DIR/summary" "$CL" step_summary combined || return 1
  log "Results: $COMBINED_DIR/summary/master_summary.tsv (all tools merged per BGC locus, ranked)"
  log "         $COMBINED_DIR/summary/per_tool/ (one table per BGC finder, with MIBiG hits)"
}

# =============================================================================
#  SETUP: programs + databases
# =============================================================================
do_setup() {
  mkdir -p "$ENV_DIR" "$DB_DIR"
  MAIN_LOG="$SETUP_DIR/setup.log"; STEP_LOG="$MAIN_LOG"
  log "Installing programs into $ENV_DIR (one conda environment each; skips existing ones)"
  local env fail=0
  for env in "${!ENV_PACKAGES[@]}"; do
    if [[ -d "$ENV_DIR/$env/conda-meta" ]]; then log "  $env: already installed"; continue; fi
    log "  $env: installing ${ENV_PACKAGES[$env]}"
    # shellcheck disable=SC2086
    "$CONDA_CREATE" create -y -p "$ENV_DIR/$env" -c conda-forge -c bioconda ${ENV_PACKAGES[$env]} \
      >>"$MAIN_LOG" 2>&1 || { warn "  $env: install FAILED (see $MAIN_LOG)"; fail=1; }
  done

  # Links to every BioCloud database (only for finding them; the pipeline reads the SHARED_* paths)
  if [[ -n "$SHARED_DB_ROOT" && -d "$SHARED_DB_ROOT" ]]; then
    local link_dir="$DB_DIR/biocloud" d n=0
    mkdir -p "$link_dir"
    for d in "$SHARED_DB_ROOT"/*; do
      [[ -e "$link_dir/$(basename "$d")" ]] || { ln -s "$d" "$link_dir/" && n=$((n + 1)); }
    done
    log "  BioCloud databases: linked in $link_dir ($n new links)"
  fi

  log "Databases: using shared server copies where they fit, downloading the rest into $DB_DIR"
  # antiSMASH: try each shared version folder (newest first), then our own copy, else download
  local as_db=""
  if [[ -n "$SHARED_ANTISMASH_DB" && -d "$SHARED_ANTISMASH_DB" ]]; then
    while read -r d; do
      if antismash_db_ok "$d"; then as_db=$d; break
      else log "  antiSMASH: shared $d does not fit antiSMASH $("$CONDA_RUN" run -p "$ENV_DIR/antismash" antismash --version 2>/dev/null | awk '{print $2}') — skipped"; fi
    done < <(find -L "$SHARED_ANTISMASH_DB" -maxdepth 3 -type d -name clusterblast \
               -printf '%h\n' 2>/dev/null | sort -rV)
  fi
  if [[ -z "$as_db" ]]; then
    if [[ ! -d "$ANTISMASH_DB/clusterblast" ]]; then
      log "  antiSMASH: downloading our own (~10 GB)"
      mkdir -p "$ANTISMASH_DB"
      run_in_env antismash download-antismash-databases --database-dir "$ANTISMASH_DB" \
        || { warn "antiSMASH databases FAILED"; fail=1; }
    fi
    as_db=$ANTISMASH_DB
  fi
  echo "$as_db" >"$ANTISMASH_DB_CHOICE"
  log "  antiSMASH: using $as_db"

  if bakta_db_ok "$(bakta_db_path)"; then log "  Bakta: present ($(bakta_db_path))"
  else
    [[ -n "$SHARED_BAKTA_DB" ]] && warn "  no Bakta database for Bakta >= 1.10 in $SHARED_BAKTA_DB (too old or missing) — downloading our own"
    mkdir -p "$BAKTA_DB_PARENT"
    # a half-finished download from an earlier attempt fails the MD5 check -> remove it first
    rm -f "$BAKTA_DB_PARENT"/*.tar.xz "$BAKTA_DB_PARENT"/*.tar.gz
    run_in_env bakta bakta_db download --output "$BAKTA_DB_PARENT" --type "$BAKTA_DB_TYPE" \
      || { warn "Bakta database FAILED (slow or blocked connection to zenodo.org?)"; fail=1; }
  fi

  if [[ -n "$(checkm2_dmnd)" ]]; then log "  CheckM2: present ($(checkm2_dmnd))"
  else
    [[ -n "$SHARED_CHECKM2_DB" ]] && warn "  no CheckM2 database found in $SHARED_CHECKM2_DB — downloading our own"
    mkdir -p "$CHECKM2_DB"
    run_in_env checkm2 checkm2 database --download --path "$CHECKM2_DB" \
      || { warn "CheckM2 database FAILED (slow or blocked connection to zenodo.org?)"; fail=1; }
  fi

  if [[ ! -d "$CARD_DB/localDB" ]]; then
    mkdir -p "$CARD_DB"
    ( cd "$CARD_DB" && curl -fsSL https://card.mcmaster.ca/latest/data -o card_data.tar.bz2 \
      && tar -xjf card_data.tar.bz2 && run_in_env rgi rgi load --card_json "$CARD_DB/card.json" --local ) \
      || { warn "CARD database FAILED"; fail=1; }
  else log "  CARD: present"; fi

  if [[ $RUN_GTDBTK -eq 1 ]]; then
    if [[ -d "$(gtdbtk_db_path)/taxonomy" ]]; then log "  GTDB-Tk: present ($(gtdbtk_db_path))"
    else
      [[ -n "$SHARED_GTDBTK_DB" ]] && warn "  no GTDB-Tk database found in $SHARED_GTDBTK_DB — downloading our own"
      mkdir -p "$GTDBTK_DB"
      run_in_env gtdbtk bash -c 'export GTDBTK_DATA_PATH="$1"; download-db.sh "$1"' _ "$GTDBTK_DB" \
        || { warn "GTDB-Tk database FAILED (no connection to data.gtdb.ecogenomic.org?) — set SHARED_GTDBTK_DB or RUN_GTDBTK=0"; fail=1; }
    fi
  fi

  # DeepBGC: models + its own Pfam 31 (~0.5 GB) into our folder (default would be ~/.local/share)
  if [[ $RUN_DEEPBGC -eq 1 ]]; then
    if [[ -n "$(find -L "$DEEPBGC_DB" -name '*.pkl' 2>/dev/null | head -n 1)" ]]; then log "  DeepBGC: present ($DEEPBGC_DB)"
    else
      mkdir -p "$DEEPBGC_DB"
      run_in_env deepbgc bash -c 'export DEEPBGC_DOWNLOADS_DIR="$1"; deepbgc download' _ "$DEEPBGC_DB" \
        || { warn "DeepBGC models FAILED"; fail=1; }
    fi
  fi

  # BGC-Prophet: its two models (GitHub) + the ESM-2 8M protein model it builds on (~30 MB,
  # fetched by torch on first use, so done here because compute nodes may have no internet)
  if [[ $RUN_BGCPROPHET -eq 1 ]]; then
    if [[ -n "$(bgcprophet_model_dir)" ]]; then log "  BGC-Prophet models: present ($(bgcprophet_model_dir))"
    else
      mkdir -p "$BGCPROPHET_DB"
      ( cd "$BGCPROPHET_DB" && curl -fsSL "$BGCPROPHET_MODEL_URL" -o model.tar.gz && tar -xzf model.tar.gz ) \
        && [[ -n "$(bgcprophet_model_dir)" ]] || { warn "BGC-Prophet models FAILED"; fail=1; }
    fi
    if [[ -n "$(find -L "$BGCPROPHET_DB/torch" -name 'esm2_t6_8M_UR50D*' 2>/dev/null | head -n 1)" ]]; then log "  ESM-2 (BGC-Prophet): present"
    else
      run_in_env bgcprophet bash -c 'export TORCH_HOME="$1"; python -c "import esm; esm.pretrained.esm2_t6_8M_UR50D()"' \
        _ "$BGCPROPHET_DB/torch" || { warn "ESM-2 download for BGC-Prophet FAILED"; fail=1; }
    fi
  fi

  log "Databases in use:"
  log "  antiSMASH  $(antismash_db_path)"
  log "  MIBiG      $(kcb_dir)   (antiSMASH KnownClusterBlast data, used for every tool)"
  log "  Bakta      $(bakta_db_path)"
  log "  CheckM2    $(checkm2_dmnd)"
  log "  CARD       $CARD_DB/localDB"
  [[ $RUN_GTDBTK -eq 1 ]] && log "  GTDB-Tk    $(gtdbtk_db_path)"
  [[ $RUN_DEEPBGC -eq 1 ]] && log "  DeepBGC    $DEEPBGC_DB"
  [[ $RUN_BGCPROPHET -eq 1 ]] && log "  BGC-Prophet $(bgcprophet_model_dir)"
  if [[ $fail -eq 0 ]]; then
    log "Setup complete. Put reads in $RAW_DIR, then: bash $SCRIPT_PATH run --dry-run"
  else
    log "Setup finished WITH ERRORS — see $MAIN_LOG, fix, and run setup again (finished parts are skipped)."
  fi
  return $fail
}

write_versions() {
  local out="$RESULTS_DIR/software_versions.tsv" env pkg ver
  printf 'environment\tpackage\tversion\n' >"$out"
  for env in "${!ENV_PACKAGES[@]}"; do
    for pkg in ${ENV_PACKAGES[$env]}; do
      pkg=${pkg%%[<>=]*}
      ver=$("$CONDA_RUN" list -p "$ENV_DIR/$env" 2>/dev/null | awk -v p="$pkg" '$1==p{print $2; exit}')
      printf '%s\t%s\t%s\n' "$env" "$pkg" "${ver:-not installed}" >>"$out"
    done
  done
}

# =============================================================================
#  SLURM: one array job (one task per isolate) + one combined job afterwards
# =============================================================================
do_submit() {
  command -v sbatch >/dev/null || die "sbatch not found — run this on the Slurm login node"
  local n=${#SAMPLES[@]} opts=() jid cid
  mkdir -p "$RESULTS_DIR/slurm_logs"
  [[ -n "$SLURM_PARTITION" ]] && opts+=(--partition "$SLURM_PARTITION")
  [[ -n "$SLURM_ACCOUNT"   ]] && opts+=(--account "$SLURM_ACCOUNT")
  jid=$(sbatch --parsable "${opts[@]}" --job-name bgc_isolates --array "1-$n" \
        --cpus-per-task "$SLURM_SAMPLE_CPUS" --mem "$SLURM_SAMPLE_MEM" --time "$SLURM_SAMPLE_TIME" \
        --output "$RESULTS_DIR/slurm_logs/isolate_%A_%a.out" \
        --wrap "bash '$SCRIPT_PATH' run --array-task --no-combined") || die "sbatch failed"
  cid=$(sbatch --parsable "${opts[@]}" --job-name bgc_combined --dependency "afterany:$jid" \
        --cpus-per-task "$SLURM_COMBINED_CPUS" --mem "$SLURM_COMBINED_MEM" --time "$SLURM_COMBINED_TIME" \
        --output "$RESULTS_DIR/slurm_logs/combined_%j.out" \
        --wrap "bash '$SCRIPT_PATH' combined") || die "sbatch failed"
  echo "Submitted: isolate array job $jid ($n isolates) and combined job $cid (starts after)."
  echo "Watch with: squeue --me      Logs: $RESULTS_DIR/slurm_logs/"
}

# =============================================================================
#  SUMMARY CODE (Python). Written to results/_combined/logs/summary.py and run by the
#  summary step. Bash only stores this text (quoted heredoc = nothing is expanded).
# =============================================================================
write_summary_py() {
cat >"$1" <<'PYTHON'
import argparse, csv, glob, json, os, re, sys
from collections import defaultdict

NA = "NA"
ANTIBIOTIC_CLASSES = {"NRPS", "PKS", "RiPP"}


def read_tsv(path):
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        lines = [l for l in fh if l.strip()]
    return list(csv.DictReader(lines, delimiter="\t"))


def pick(row, *names):
    low = {(k or "").lower().strip(): v for k, v in row.items()}
    for n in names:
        v = low.get(n.lower())
        if v not in (None, ""):
            return v
    return NA


def overlap(a1, a2, b1, b2):
    return a1 <= b2 and b1 <= a2


def genome_info(sd, s):
    g = {"isolate": s}
    rs = os.path.join(sd, "02_filtlong", "read_stats.tsv")
    if not os.path.isfile(rs):
        rs = os.path.join(sd, "03_mmlong2", "read_stats.tsv")   # results made before v3.4.0
    for r in read_tsv(rs):
        tag = "filtered" if "filtered" in pick(r, "file") else "raw"
        g[f"{tag}_reads"] = pick(r, "num_seqs"); g[f"{tag}_bases"] = pick(r, "sum_len")
        g[f"{tag}_read_N50"] = pick(r, "N50"); g[f"{tag}_meanQ"] = pick(r, "AvgQual")
    # "polished" = column name in results made before v3.3.0
    q = {r.get("Assembly"): r.get("final", r.get("polished", NA)) for r in read_tsv(os.path.join(sd, "05_quast", "report.tsv"))}
    g.update({"genome_length": q.get("Total length", NA), "contigs": q.get("# contigs", NA),
              "N50": q.get("N50", NA), "GC_percent": q.get("GC (%)", NA)})
    contigs = []
    info = os.path.join(sd, "03_flye", "assembly_info.txt")
    if not os.path.isfile(info):
        info = os.path.join(sd, "03_mmlong2", "assembly_info.txt")   # results made before v3.4.0
    if os.path.isfile(info):
        for line in open(info):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            contigs.append((int(f[1]), float(f[2]), f[3] == "Y"))
    g["circular_contigs"] = sum(c[2] for c in contigs) if contigs else NA
    g["mean_depth"] = round(sum(c[0] * c[1] for c in contigs) / sum(c[0] for c in contigs), 1) if contigs else NA
    ck = read_tsv(os.path.join(sd, "06_checkm2", "quality_report.tsv"))
    g["completeness"] = pick(ck[0], "Completeness") if ck else NA
    g["contamination"] = pick(ck[0], "Contamination") if ck else NA
    g["CDS"] = NA
    txt = os.path.join(sd, "07_bakta", f"{s}.txt")
    if os.path.isfile(txt):
        for line in open(txt):
            m = re.match(r"\s*CDSs?:\s*(\d+)", line)
            if m:
                g["CDS"] = m.group(1)
    return g


def antismash_regions(sd, s):
    d = os.path.join(sd, "08_antismash")
    path = os.path.join(d, f"{s}.json")
    if not os.path.isfile(path):
        return []
    data = json.load(open(path))
    edge = {os.path.basename(p)[:-4]: '/contig_edge="True"' in open(p).read()
            for p in glob.glob(os.path.join(d, "*.region*.gbk"))}
    out = []
    for rec in data.get("records", []):
        rid = rec.get("id", "?")
        kcb = {}
        try:
            for res in rec["modules"]["antismash.modules.clusterblast"]["knowncluster"]["results"]:
                if res.get("ranking"):
                    ref, score = res["ranking"][0]
                    kcb[res["region_number"]] = (ref.get("accession", NA), ref.get("description", NA),
                                                 score.get("similarity", NA))
        except (KeyError, TypeError, ValueError):
            pass
        for i, a in enumerate(rec.get("areas", []), start=1):
            cats = sorted({p.get("category", "") for p in a.get("protoclusters", {}).values() if p.get("category")})
            region = f"{rid}.region{i:03d}"
            hit = kcb.get(i, (NA, NA, NA))
            out.append({"isolate": s, "contig": rid, "region": region,
                        "start": int(a["start"]) + 1, "end": int(a["end"]),
                        "length_kb": round((int(a["end"]) - int(a["start"])) / 1000, 1),
                        "bgc_type": "+".join(a.get("products", [])) or NA,
                        "bgc_class": "+".join(cats) or NA,
                        "contig_edge": "yes" if edge.get(region) else "no",
                        "closest_MIBiG": hit[0], "closest_MIBiG_compound": hit[1],
                        "MIBiG_similarity_pct": hit[2]})
    return out


def gecco(sd):
    p = sorted(glob.glob(os.path.join(sd, "09_gecco", "*.clusters.tsv")))
    out = []
    for r in read_tsv(p[0] if p else None):
        try:
            out.append({"contig": pick(r, "sequence_id"), "id": pick(r, "cluster_id", "bgc_id"),
                        "start": int(pick(r, "start")), "end": int(pick(r, "end")), "type": pick(r, "type"),
                        "average_p": pick(r, "average_p"), "max_p": pick(r, "max_p")})
        except ValueError:
            pass
    return out


def deepbgc(sd):
    p = sorted(glob.glob(os.path.join(sd, "13_deepbgc", "*.bgc.tsv")))
    out = []
    for r in read_tsv(p[0] if p else None):
        try:
            out.append({"contig": pick(r, "sequence_id"), "id": pick(r, "bgc_candidate_id"),
                        "start": int(pick(r, "nucl_start")) + 1, "end": int(pick(r, "nucl_end")),   # 0-based start
                        "deepbgc_score": pick(r, "deepbgc_score"), "deepbgc_class": pick(r, "product_class"),
                        "deepbgc_activity": pick(r, "product_activity"),
                        **{f"deepbgc_{k}": pick(r, k) for k in ("antibacterial", "antifungal", "cytotoxic", "inhibitor")}})
        except ValueError:
            pass
    return out


def bakta_genes(sd, s):
    """locus tag -> (contig, start, end) from Bakta's .tsv (columns: contig, type, start, stop, strand, locus tag)"""
    genes = {}
    p = os.path.join(sd, "07_bakta", f"{s}.tsv")
    if os.path.isfile(p):
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 6 or not f[5]:
                continue
            try:
                genes[f[5]] = (f[0], int(f[2]), int(f[3]))
            except ValueError:
                pass
    return genes


def bgc_prophet(sd, s, genes):
    """BGC-Prophet gives each BGC as a list of gene IDs (Bakta locus tags) -> turn into coordinates"""
    p = os.path.join(sd, "14_bgc_prophet", f"{s}_classified.csv")
    if not os.path.isfile(p):
        return []
    out, seen = [], set()
    for r in csv.DictReader(open(p, newline="")):
        if r.get("isBGC") != "Yes":
            continue
        loc = [genes[g] for g in re.findall(r"[^\s\[\]',\"]+", r.get("sentence") or "") if g in genes]
        if not loc:
            continue
        contig = max({c for c, _, _ in loc}, key=[c for c, _, _ in loc].count)   # the 128-gene windows can cross contigs
        key = (contig, min(a for c, a, _ in loc if c == contig), max(b for c, _, b in loc if c == contig))
        if key in seen:                                     # the last window overlaps the one before
            continue
        seen.add(key)
        out.append({"contig": key[0], "id": r.get("ID", NA), "start": key[1], "end": key[2],
                    "bgc_prophet_class": r.get("labels", NA), "bgc_prophet_probability": r.get("probability", NA)})
    return out


MIN_IDENT, MIN_COVER = 45, 40   # same cut-offs as antiSMASH KnownClusterBlast (exclusive lower bounds)


def mibig_refs(d):
    """MIBiG cluster ('BGC0000001|c1') -> (compound, set of its protein IDs), from KnownClusterBlast clusters.txt"""
    refs = {}
    p = os.path.join(d, "clusters.txt") if d else ""
    if os.path.isfile(p):
        for line in open(p, encoding="utf-8"):
            t = line.rstrip("\n").split("\t")
            if len(t) >= 6:
                refs[f"{t[0]}|{t[2]}"] = (t[1], {x for x in t[5].split(";") if x})
    return refs


def mibig_hits(sd, s):
    """Bakta locus tag -> set of (MIBiG cluster, MIBiG protein) it hits (DIAMOND, KnownClusterBlast cut-offs)"""
    hits = defaultdict(set)
    p = os.path.join(sd, "15_mibig", f"{s}_vs_mibig.tsv")
    if os.path.isfile(p):
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, subj, ident, alen, qlen = f[0], f[1].split("|"), float(f[2]), float(f[3]), float(f[4])
            if len(subj) < 5 or ident <= MIN_IDENT or 100 * alen / qlen <= MIN_COVER:
                continue
            shift = 1 if subj[4].startswith("linearised") else 0   # as antiSMASH parses these headers
            hits[q].add((f"{subj[0]}|{subj[1]}", subj[4 + shift]))
    return hits


def best_mibig(contig, start, end, genes, hits, refs):
    """KnownClusterBlast-style best MIBiG match for one BGC: % of the MIBiG cluster's genes hit"""
    per_ref = defaultdict(set)
    for tag, (c, a, b) in genes.items():
        if c == contig and overlap(a, b, start, end):
            for ref, prot in hits.get(tag, ()):
                per_ref[ref].add(prot)
    best = None
    for ref, prots in per_ref.items():
        if ref not in refs or not refs[ref][1]:
            continue
        sim = int(100 * len(prots & refs[ref][1]) / len(refs[ref][1]))
        if best is None or (sim, len(prots)) > (best[2], best[3]):
            best = (ref.split("|")[0], refs[ref][0], sim, len(prots))
    return best


def merge_loci(bgcs):
    """group predictions of all tools that overlap on the same contig into one locus"""
    loci = []
    for b in sorted(bgcs, key=lambda x: (x["contig"], x["start"])):
        if loci and loci[-1]["contig"] == b["contig"] and b["start"] <= loci[-1]["end"]:
            loci[-1]["end"] = max(loci[-1]["end"], b["end"]); loci[-1]["members"].append(b)
        else:
            loci.append({"contig": b["contig"], "start": b["start"], "end": b["end"], "members": [b]})
    return loci


TOOLS = ["antismash", "gecco", "deepbgc", "bgc_prophet"]
QC_TAX = ["genome_length", "contigs", "N50", "GC_percent", "mean_depth", "circular_contigs", "completeness",
          "contamination", "gtdb_classification", "closest_reference", "closest_ANI", "possible_new_species"]
ANTIBIOTIC_WORDS = {"NRPS", "PKS", "NRP", "Polyketide", "RiPP"}   # antiSMASH classes + DeepBGC/GECCO/BGC-Prophet classes


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def tool_calls(s, regions, gec, dbg, pro, genes, hits, refs, mibig_ran):
    """one row per BGC per tool, all with the same MIBiG columns"""
    calls = []
    for tool, rows in zip(TOOLS, (regions, gec, dbg, pro)):
        for r in rows:
            c = {"isolate": s, "tool": tool, "bgc_id": r.get("region", r.get("id")), "contig": r["contig"],
                 "start": r["start"], "end": r["end"], "length_kb": round((r["end"] - r["start"] + 1) / 1000, 1)}
            if tool == "antismash":
                c.update(antismash_type=r["bgc_type"], antismash_class=r["bgc_class"], contig_edge=r["contig_edge"])
            elif tool == "gecco":
                c.update(gecco_type=r["type"], gecco_average_p=r["average_p"], gecco_max_p=r["max_p"])
            else:
                c.update({k: v for k, v in r.items() if k.startswith(tool + "_")})
            if tool == "antismash" and r["closest_MIBiG"] != NA:      # antiSMASH's own KnownClusterBlast hit
                mib = (r["closest_MIBiG"], r["closest_MIBiG_compound"], r["MIBiG_similarity_pct"], "antiSMASH KnownClusterBlast")
            else:
                b = best_mibig(c["contig"], c["start"], c["end"], genes, hits, refs)
                mib = (b[0], b[1], b[2], "DIAMOND vs KnownClusterBlast data") if b else \
                      (NA, NA, NA, "no MIBiG hit" if mibig_ran else "not_run")
            c.update(zip(("MIBiG_accession", "MIBiG_compound", "MIBiG_similarity_pct", "MIBiG_method"), mib))
            calls.append(c)
    return calls


def master_rows(s, g, calls, res, fam, members, bs_ran, known, related):
    """merge overlapping predictions of all tools into loci, score and describe each locus"""
    rows = []
    for n, L in enumerate(merge_loci(calls), start=1):
        m, lid = L["members"], f"{s}_locus{n:03d}"
        for c in m:
            c["locus_id"] = lid
        tools = sorted({c["tool"] for c in m}, key=TOOLS.index)

        def vals(k):
            return ";".join(sorted({str(c[k]) for c in m if c.get(k) not in (None, "", NA)})) or NA

        best = max((c for c in m if fnum(c["MIBiG_similarity_pct"]) is not None),
                   key=lambda c: fnum(c["MIBiG_similarity_pct"]), default=None)
        sim = fnum(best["MIBiG_similarity_pct"]) if best else 0.0
        rh = [h for h in res if h["contig"] == L["contig"] and overlap(h["start"], h["end"], L["start"], L["end"])]
        keys = [fam.get(f"{s}__{c['bgc_id']}") for c in m if c["tool"] == "antismash"]
        keys = [k for k in keys if k]
        fam_mibig = sorted({x for k in keys for x in members[k] if x.startswith("BGC")})
        novelty = "known" if sim >= known else ("related" if sim >= related or fam_mibig else "novel_candidate")
        words = set(re.split(r"[+;,\s]+", " ".join(vals(k) for k in
                    ("antismash_class", "gecco_type", "deepbgc_class", "bgc_prophet_class"))))
        activity = max([fnum(c.get(k)) or 0 for c in m for k in ("deepbgc_antibacterial", "deepbgc_antifungal")], default=0)
        edge = any(c.get("contig_edge") == "yes" for c in m)
        # ranking: novelty (3/1/0) + antibiotic-type class (2) + each extra tool that agrees (1)
        #          + resistance gene inside (2) + DeepBGC antibacterial/antifungal >= 0.5 (1) - cut off at contig edge (1)
        score = {"novel_candidate": 3, "related": 1, "known": 0}[novelty] + (2 if words & ANTIBIOTIC_WORDS else 0) \
            + (len(tools) - 1) + (2 if rh else 0) + (1 if activity >= 0.5 else 0) - (1 if edge else 0)
        rows.append({"priority_score": score, "novelty": novelty, "isolate": s, "locus_id": lid,
                     "contig": L["contig"], "start": L["start"], "end": L["end"],
                     "length_kb": round((L["end"] - L["start"] + 1) / 1000, 1),
                     "n_tools": len(tools), "tools": ";".join(tools),
                     "antismash_regions": ";".join(c["bgc_id"] for c in m if c["tool"] == "antismash") or NA,
                     **{k: vals(k) for k in ("antismash_type", "antismash_class", "gecco_type", "deepbgc_class",
                                             "deepbgc_activity", "bgc_prophet_class")},
                     "deepbgc_max_antibacterial_or_antifungal": round(activity, 3) if "deepbgc" in tools else NA,
                     "contig_edge": "yes" if edge else ("no" if "antismash" in tools else NA),
                     "closest_MIBiG": best["MIBiG_accession"] if best else NA,
                     "closest_MIBiG_compound": best["MIBiG_compound"] if best else NA,
                     "MIBiG_similarity_pct": best["MIBiG_similarity_pct"] if best else NA,
                     "MIBiG_from": f"{best['tool']} ({best['MIBiG_method']})" if best else NA,
                     "resistance_gene_in_locus": ";".join(f"{h['gene']} ({h['drug']})" for h in rh) or "none",
                     "bigscape_family": ";".join(sorted({k[1] for k in keys})) or (NA if bs_ran else "not_run"),
                     "family_MIBiG_members": ";".join(fam_mibig) or ("none" if keys else NA),
                     **{k: g.get(k, NA) for k in QC_TAX}})
    return rows


def rgi(sd, s):
    out = []
    for r in read_tsv(os.path.join(sd, "10_rgi_card", f"{s}_rgi.txt")):
        try:
            out.append({"contig": re.sub(r"_\d+$", "", pick(r, "Contig")), "start": int(pick(r, "Start")),
                        "end": int(pick(r, "Stop")), "gene": pick(r, "Best_Hit_ARO"), "drug": pick(r, "Drug Class")})
        except ValueError:
            pass
    return out


def gtdb(comb):
    t = {}
    for p in glob.glob(os.path.join(comb, "gtdbtk", "*.summary.tsv")):
        for r in read_tsv(p):
            cls = pick(r, "classification")
            t[pick(r, "user_genome")] = {
                "gtdb_classification": cls,
                "closest_reference": pick(r, "closest_genome_reference", "fastani_reference"),
                "closest_ANI": pick(r, "closest_genome_ani", "fastani_ani"),
                "possible_new_species": "yes" if cls.endswith("s__") else ("no" if cls != NA else NA)}
    return t


def bigscape(comb, cutoff):
    fam, members, ran = {}, defaultdict(set), False
    for p in glob.glob(os.path.join(comb, "bigscape", "output", "**", "*clustering_c*.tsv"), recursive=True):
        m = re.search(r"clustering_c([0-9.]+?)\.tsv$", os.path.basename(p))
        if not m or abs(float(m.group(1)) - float(cutoff)) > 1e-6:
            continue
        ran = True
        for r in read_tsv(p):
            g, f = pick(r, "GBK", "BGC", "Record"), pick(r, "Family", "GCF")
            if NA in (g, f):
                continue
            g = os.path.basename(g)
            g = g[:-4] if g.endswith(".gbk") else g
            key = (os.path.basename(p), f)
            fam.setdefault(g, key); members[key].add(g)
    return fam, members, ran


def write(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols or ["no_results"], delimiter="\t", restval=NA)
        w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--cutoff", default="0.30")
    ap.add_argument("--known", type=float, default=80); ap.add_argument("--related", type=float, default=30)
    ap.add_argument("--mibig", default="", help="antiSMASH knownclusterblast/<version> folder")
    a = ap.parse_args()
    comb = os.path.join(a.results, "_combined")
    taxa = gtdb(comb)
    fam, members, bs_ran = bigscape(comb, a.cutoff)
    refs = mibig_refs(a.mibig)
    genomes, bgcs, calls, master = [], [], [], []

    for s in [x for x in a.samples.split(",") if x]:
        sd = os.path.join(a.results, s)
        regions, gec, res = antismash_regions(sd, s), gecco(sd), rgi(sd, s)
        genes = bakta_genes(sd, s)
        dbg, pro = deepbgc(sd), bgc_prophet(sd, s, genes)
        g = genome_info(sd, s)
        g.update(taxa.get(s, {}))
        g.update({"BGCs_antismash": len(regions), "BGCs_gecco": len(gec), "BGCs_deepbgc": len(dbg),
                  "BGCs_bgc_prophet": len(pro), "resistance_genes": len(res)})
        genomes.append(g)

        mibig_ran = os.path.isfile(os.path.join(sd, "15_mibig", f"{s}_vs_mibig.tsv"))
        tc = tool_calls(s, regions, gec, dbg, pro, genes, mibig_hits(sd, s), refs, mibig_ran)
        loci = master_rows(s, g, tc, res, fam, members, bs_ran, a.known, a.related)
        g["BGC_loci_all_tools"] = len(loci)
        calls += tc; master += loci

        used = set()
        for r in regions:
            gh = [c for c in gec if c["contig"] == r["contig"] and overlap(c["start"], c["end"], r["start"], r["end"])]
            used.update(id(c) for c in gh)
            rh = [h for h in res if h["contig"] == r["contig"] and overlap(h["start"], h["end"], r["start"], r["end"])]
            key = fam.get(f"{s}__{r['region']}")
            if key:
                mem = members[key]
                mibig = sorted(m for m in mem if m.startswith("BGC"))
                fam_id, fam_mibig = key[1], ";".join(mibig) or "none"
                fam_iso = len({m.split("__")[0] for m in mem if "__" in m})
            else:
                mibig, fam_id, fam_mibig, fam_iso = [], (NA if bs_ran else "not_run"), NA, NA
            try:
                sim = float(r["MIBiG_similarity_pct"])
            except (TypeError, ValueError):
                sim = 0.0
            novelty = "known" if sim >= a.known else ("related" if sim >= a.related or mibig else "novel_candidate")
            score = {"novel_candidate": 3, "related": 1, "known": 0}[novelty]
            score += 2 if set(r["bgc_class"].split("+")) & ANTIBIOTIC_CLASSES else 0
            score += 1 if gh else 0
            score += 2 if rh else 0
            score -= 1 if r["contig_edge"] == "yes" else 0
            bgcs.append({"priority_score": score, "novelty": novelty, **r,
                         "gecco_agrees": "yes" if gh else "no",
                         "resistance_gene_in_BGC": ";".join(f"{h['gene']} ({h['drug']})" for h in rh) or "none",
                         "bigscape_family": fam_id, "family_MIBiG_members": fam_mibig,
                         "family_in_n_isolates": fam_iso})
        for c in gec:
            if id(c) not in used:
                bgcs.append({"priority_score": 1, "novelty": "gecco_only_check_manually", "isolate": s,
                             "contig": c["contig"], "region": f"GECCO:{c['id']}", "start": c["start"],
                             "end": c["end"], "length_kb": round((c["end"] - c["start"]) / 1000, 1),
                             "bgc_type": "not_found_by_antismash", "bgc_class": c["type"],
                             "gecco_agrees": "only_gecco"})

    bgcs.sort(key=lambda r: (-r["priority_score"], -float(r["length_kb"])))
    bgcs = [{"rank": i, **r} for i, r in enumerate(bgcs, start=1)]
    write(os.path.join(a.out, "genome_summary.tsv"), genomes)
    write(os.path.join(a.out, "bgc_summary.tsv"), bgcs)            # antiSMASH-based ranking (as before v3.3.0)
    for t in TOOLS:
        write(os.path.join(a.out, "per_tool", f"{t}_bgcs.tsv"), [c for c in calls if c["tool"] == t])
    master.sort(key=lambda r: (-r["priority_score"], -r["n_tools"], -r["length_kb"]))
    write(os.path.join(a.out, "master_summary.tsv"), [{"rank": i, **r} for i, r in enumerate(master, start=1)])

    # per-isolate copies: results/<isolate>/11_gtdbtk and 12_novelty
    for g in genomes:
        s = g["isolate"]
        write(os.path.join(a.results, s, "12_novelty", f"{s}_bgc_novelty.tsv"),
              [b for b in bgcs if b["isolate"] == s])
        write(os.path.join(a.results, s, "11_gtdbtk", f"{s}_taxonomy.tsv"),
              [{k: g.get(k, NA) for k in ("isolate", "gtdb_classification", "closest_reference",
                                          "closest_ANI", "possible_new_species")}])
    n_new = sum(b["novelty"] == "novel_candidate" for b in master)
    print(f"{len(genomes)} genomes, {len(calls)} BGC calls from {len(TOOLS)} tools = {len(master)} loci, "
          f"{n_new} novel candidates -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
PYTHON
}

# =============================================================================
#  MAIN
# =============================================================================
usage() {   # prints the comment block at the top of this file (up to SETTINGS)
  sed -n '3,/^# =* SETTINGS/p' "$SCRIPT_PATH" | sed '$d' | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

CMD=${1:-help}; shift || true
ONLY=""; ARRAY_TASK=0; DO_COMBINED=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) ONLY=$2; shift 2 ;;
    --force) IFS=',' read -r -a FORCE_STEPS <<<"$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --no-combined) DO_COMBINED=0; shift ;;
    --array-task) ARRAY_TASK=1; shift ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

case "$CMD" in
  setup)
    find_conda; do_setup ;;
  list)
    find_samples
    for i in "${!SAMPLES[@]}"; do printf '%3d  %-25s %s\n' "$((i + 1))" "${SAMPLES[$i]}" "${SAMPLE_SRC[$i]}"; done ;;
  run|combined)
    find_conda; find_samples
    [[ $DRY_RUN -eq 0 ]] && mkdir -p "$RESULTS_DIR/logs"
    MAIN_LOG="$RESULTS_DIR/logs/pipeline_$(date +%Y%m%d_%H%M%S)_$$.log"
    [[ $DRY_RUN -eq 1 ]] && MAIN_LOG=/dev/null
    log "bgc_pipeline.sh v$PIPELINE_VERSION ($PIPELINE_UPDATED) — Mikkel Bennetsen, AAU — github.com/Mikkelbenn2/Fantastic-drugs-and-where-to-find-them"
    log "Project: $PROJECT_DIR | isolates: ${#SAMPLES[@]} | threads: $THREADS | dry run: $DRY_RUN"
    [[ $DRY_RUN -eq 0 ]] && write_versions
    FAILED=()
    if [[ "$CMD" == run ]]; then
      if [[ $ARRAY_TASK -eq 1 ]]; then
        idx=$(( ${SLURM_ARRAY_TASK_ID:?not inside a Slurm array} - 1 ))
        ONLY=${SAMPLES[$idx]}
      fi
      found=0
      for i in "${!SAMPLES[@]}"; do
        [[ -n "$ONLY" && "${SAMPLES[$i]}" != "$ONLY" ]] && continue
        found=1
        process_sample "${SAMPLES[$i]}" "${SAMPLE_SRC[$i]}" || FAILED+=("${SAMPLES[$i]}")
      done
      [[ $found -eq 1 ]] || die "Isolate '$ONLY' not found in raw_data (see: list)"
    fi
    if [[ "$CMD" == combined || ( $DO_COMBINED -eq 1 && -z "$ONLY" ) ]]; then
      process_combined || FAILED+=("combined")
    fi
    if [[ ${#FAILED[@]} -gt 0 ]]; then log "Finished with failures: ${FAILED[*]} (check results/<isolate>/logs/)"; exit 1; fi
    log "All done." ;;
  submit)
    find_samples; do_submit ;;
  help|-h|--help) usage 0 ;;
  *) echo "Unknown command: $CMD" >&2; usage 1 ;;
esac
exit 0
