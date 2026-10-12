# =============================================================================
#  bgc_pipeline.sh — ONT genomes of soil isolates -> screen for NOVEL antimicrobial BGCs
#  P7 project (AAU). One self-contained script: paste it into  setup/bgc_pipeline.sh
# =============================================================================

# ============================== SETTINGS =====================================
# Change these to fit your data and server. Everything below SETTINGS can stay.

THREADS=16                     # CPUs per program (inside Slurm: taken from the job)

# Reads. "nano-hq" = R10.4.1 + Dorado hac/sup (Q15-20+). Use "nano-raw" for old R9.4.1 data.
FLYE_READ_MODE="nano-hq"
GENOME_SIZE=""                 # e.g. "8m" for Streptomyces; empty = Flye estimates it
FILTLONG_MIN_LEN=1000          # drop reads shorter than this (bp)
FILTLONG_KEEP_PERCENT=95       # keep the best 95 % of bases
FILTLONG_TARGET_BASES=""       # optional cap, e.g. 800000000 (= 100x of an 8 Mb genome)

POLISH=1                       # 1 = polish with Medaka (--bacteria model), 0 = skip
BAKTA_DB_TYPE="light"          # "light" (~1.5 GB) or "full" (~40 GB, slightly better names)

# Optional programs (1 = run, 0 = skip)
RUN_GECCO=1                    # machine-learning BGC detection
RUN_RGI=1                      # resistance genes (CARD)
RUN_GTDBTK=1                   # species ID — needs ~110 GB database + ~64-100 GB RAM
RUN_GTOTREE=1                  # phylogeny of the isolates (+ closest GTDB references)
RUN_BIGSCAPE=1                 # BGC families compared with MIBiG

GTOTREE_HMM="Bacteria"         # marker set; "Actinobacteria" if all isolates are actinomycetes

# antiSMASH extra analyses. --cb-knownclusters (MIBiG comparison) is ESSENTIAL for novelty.
ANTISMASH_EXTRA="--cb-knownclusters --cb-subclusters --cc-mibig --asf --rre --tfbs --pfam2go --clusterhmmer --tigrfam --smcog-trees"

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
  [medaka]="medaka>=2.0 samtools>=1.20"
  [quast]="quast>=5.2"
  [checkm2]="checkm2>=1.0.2"
  [bakta]="bakta>=1.10"
  [antismash]="antismash>=8.0"
  [gecco]="gecco>=0.9.10"
  [rgi]="rgi>=6.0"
  [gtdbtk]="gtdbtk>=2.4"
  [gtotree]="gtotree>=1.8"
  [bigscape]="bigscape>=2.0"
)
# =========================== END OF SETTINGS =================================


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
find_samples() {
  [[ -d "$RAW_DIR" ]] || die "No raw_data folder found at $RAW_DIR"
  SAMPLES=(); SAMPLE_SRC=()
  local p name
  while IFS= read -r p; do
    name=$(basename "$p")
    if [[ -d "$p" ]]; then
      compgen -G "$p/*.f*q*" >/dev/null || continue           # folder without FASTQ: ignore
    else
      [[ "$name" =~ \.(fastq|fq)(\.gz)?$ ]] || continue
      name=${name%.gz}; name=${name%.fastq}; name=${name%.fq}
    fi
    name=${name//[^A-Za-z0-9._-]/_}                            # no spaces/odd characters
    SAMPLES+=("$name"); SAMPLE_SRC+=("$p")
  done < <(find "$RAW_DIR" -mindepth 1 -maxdepth 1 \( -type f -o -type d -o -type l \) | sort)
  [[ ${#SAMPLES[@]} -gt 0 ]] || die "No FASTQ files or FASTQ folders found in $RAW_DIR"
}

# Where the files of one isolate live (one place, used by every step)
set_sample() {
  S=$1; SRC=$2
  SD="$RESULTS_DIR/$S"; SL="$SD/logs"
  READS="$SD/01_nanoplot/$S.reads.fastq.gz"   # all reads joined
  FILT="$SD/02_filtlong/$S.filtered.fastq.gz"
  DRAFT="$SD/03_flye/assembly.fasta"
  GENOME="$SD/04_medaka/$S.fasta"          # final genome used by all later steps
  GBFF="$SD/07_bakta/$S.gbff"
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

# 02 Filtlong — remove short and low-quality reads so Flye and Medaka get cleaner input.
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
  need "$DRAFT"
}

# 04 Medaka — polishing. Fixes small indels; one indel in a gene = frameshift = broken
# NRPS/PKS gene. --bacteria selects the model trained on native (methylated) bacterial DNA.
step_medaka() {
  if [[ "$POLISH" -eq 1 ]]; then
    run_in_env medaka medaka_consensus -i "$FILT" -d "$DRAFT" -o "$SD/04_medaka" \
      -t "$THREADS" --bacteria || return 1
    run_plain cp "$SD/04_medaka/consensus.fasta" "$GENOME" || return 1
  else
    run_plain cp "$DRAFT" "$GENOME" || return 1
  fi
  need "$GENOME"
}

# 05 QUAST — contiguity: number of contigs, total length, N50, GC%; draft vs. polished.
step_quast() {
  run_in_env quast quast.py -o "$SD/05_quast" -t "$THREADS" \
    --labels "flye_draft,polished" "$DRAFT" "$GENOME"
}

# 06 CheckM2 — completeness (>95 % good) and contamination (<5 % good; high = mixed culture?)
step_checkm2() {
  local dmnd
  dmnd=$(find "$CHECKM2_DB" -name '*.dmnd' 2>/dev/null | head -n 1)
  [[ -n "$dmnd" || $DRY_RUN -eq 1 ]] || { log "    CheckM2 database missing — run: setup"; return 1; }
  run_in_env checkm2 checkm2 predict --threads "$THREADS" --input "$GENOME" \
    --output-directory "$SD/06_checkm2" --database_path "${dmnd:-CheckM2.dmnd}" --force
}

# 07 Bakta — finds and names all genes (successor of Prokka, which is no longer maintained).
# Its GenBank file (.gbff) is the input for antiSMASH, so all tools use the same genes.
step_bakta() {
  run_in_env bakta bakta --db "$BAKTA_DB_PARENT/db-$BAKTA_DB_TYPE" --output "$SD/07_bakta" \
    --prefix "$S" --threads "$THREADS" --keep-contig-headers --skip-plot --force "$GENOME" \
    || run_in_env bakta bakta --db "$BAKTA_DB_PARENT/db" --output "$SD/07_bakta" \
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
  local common=(--cpus "$THREADS" --taxon bacteria --databases "$ANTISMASH_DB"
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
      --input_type contig --alignment_tool DIAMOND --num_threads "$THREADS" --local --clean )
}

process_sample() {
  set_sample "$1" "$2"
  log "=== Isolate $S ==="
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$SD" "$SL"
  run_step 01_nanoplot  "$SD/01_nanoplot"  "$SL" step_nanoplot          "$S" || warn "NanoPlot failed for $S (continuing)"
  run_step 02_filtlong  "$SD/02_filtlong"  "$SL" step_filtlong          "$S" || return 1
  run_step 03_flye      "$SD/03_flye"      "$SL" step_flye              "$S" || return 1
  run_step 04_medaka    "$SD/04_medaka"    "$SL" step_medaka            "$S" || return 1
  run_step 05_quast     "$SD/05_quast"     "$SL" step_quast             "$S" || warn "QUAST failed for $S (continuing)"
  run_step 06_checkm2   "$SD/06_checkm2"   "$SL" step_checkm2           "$S" || warn "CheckM2 failed for $S (continuing)"
  run_step 07_bakta     "$SD/07_bakta"     "$SL" step_bakta             "$S" || return 1
  run_step 08_antismash "$SD/08_antismash" "$SL" step_antismash         "$S" || return 1
  [[ $RUN_GECCO -eq 1 ]] && { run_step 09_gecco    "$SD/09_gecco"    "$SL" step_gecco "$S" || warn "GECCO failed for $S (continuing)"; }
  [[ $RUN_RGI   -eq 1 ]] && { run_step 10_rgi_card "$SD/10_rgi_card" "$SL" step_rgi   "$S" || warn "RGI failed for $S (continuing)"; }
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
step_gtdbtk() {
  local out="$COMBINED_DIR/gtdbtk" name g
  [[ $DRY_RUN -eq 0 ]] && mkdir -p "$out/genomes"
  while IFS=$'\t' read -r name g; do run_plain ln -sf "$g" "$out/genomes/$name.fasta"; done < <(finished_genomes)
  run_in_env gtdbtk bash -c 'export GTDBTK_DATA_PATH="$1"; shift; gtdbtk "$@"' _ "$GTDBTK_DB" \
    classify_wf --genome_dir "$out/genomes" --extension fasta --out_dir "$out" \
    --cpus "$THREADS" --skip_ani_screen
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
  pfam=$(find "$ANTISMASH_DB" -name 'Pfam-A.hmm' 2>/dev/null | sort | tail -n 1)   # reuse antiSMASH's Pfam
  [[ -n "$pfam" || $DRY_RUN -eq 1 ]] || { log "    Pfam-A.hmm not found in $ANTISMASH_DB"; return 1; }
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
    --out "$COMBINED_DIR/summary" --cutoff "$BIGSCAPE_CUTOFF" \
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
  log "Results: $COMBINED_DIR/summary/bgc_summary.tsv (all BGCs, ranked)"
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

  log "Downloading databases into $DB_DIR (skips existing ones; GTDB-Tk takes hours)"
  if [[ ! -d "$ANTISMASH_DB/clusterblast" ]]; then
    mkdir -p "$ANTISMASH_DB"
    run_in_env antismash download-antismash-databases --database-dir "$ANTISMASH_DB" \
      || { warn "antiSMASH databases FAILED"; fail=1; }
  else log "  antiSMASH: present"; fi

  if [[ ! -d "$BAKTA_DB_PARENT/db-$BAKTA_DB_TYPE" && ! -d "$BAKTA_DB_PARENT/db" ]]; then
    mkdir -p "$BAKTA_DB_PARENT"
    run_in_env bakta bakta_db download --output "$BAKTA_DB_PARENT" --type "$BAKTA_DB_TYPE" \
      || { warn "Bakta database FAILED"; fail=1; }
  else log "  Bakta: present"; fi

  if ! find "$CHECKM2_DB" -name '*.dmnd' 2>/dev/null | grep -q .; then
    mkdir -p "$CHECKM2_DB"
    run_in_env checkm2 checkm2 database --download --path "$CHECKM2_DB" \
      || { warn "CheckM2 database FAILED"; fail=1; }
  else log "  CheckM2: present"; fi

  if [[ ! -d "$CARD_DB/localDB" ]]; then
    mkdir -p "$CARD_DB"
    ( cd "$CARD_DB" && curl -fsSL https://card.mcmaster.ca/latest/data -o card_data.tar.bz2 \
      && tar -xjf card_data.tar.bz2 && run_in_env rgi rgi load --card_json "$CARD_DB/card.json" --local ) \
      || { warn "CARD database FAILED"; fail=1; }
  else log "  CARD: present"; fi

  if [[ $RUN_GTDBTK -eq 1 ]]; then
    if [[ ! -d "$GTDBTK_DB/taxonomy" ]]; then
      mkdir -p "$GTDBTK_DB"
      run_in_env gtdbtk bash -c 'export GTDBTK_DATA_PATH="$1"; download-db.sh "$1"' _ "$GTDBTK_DB" \
        || { warn "GTDB-Tk database FAILED"; fail=1; }
    else log "  GTDB-Tk: present"; fi
  fi

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
import argparse, csv, glob, json, os, re, shutil, sys
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
    for r in read_tsv(os.path.join(sd, "02_filtlong", "read_stats.tsv")):
        tag = "filtered" if "filtered" in pick(r, "file") else "raw"
        g[f"{tag}_reads"] = pick(r, "num_seqs"); g[f"{tag}_bases"] = pick(r, "sum_len")
        g[f"{tag}_read_N50"] = pick(r, "N50"); g[f"{tag}_meanQ"] = pick(r, "AvgQual")
    q = {r.get("Assembly"): r.get("polished", NA) for r in read_tsv(os.path.join(sd, "05_quast", "report.tsv"))}
    g.update({"genome_length": q.get("Total length", NA), "contigs": q.get("# contigs", NA),
              "N50": q.get("N50", NA), "GC_percent": q.get("GC (%)", NA)})
    contigs = []
    info = os.path.join(sd, "03_flye", "assembly_info.txt")
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
                        "start": int(pick(r, "start")), "end": int(pick(r, "end")), "type": pick(r, "type")})
        except ValueError:
            pass
    return out


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
    a = ap.parse_args()
    comb = os.path.join(a.results, "_combined")
    taxa = gtdb(comb)
    fam, members, bs_ran = bigscape(comb, a.cutoff)
    genomes, bgcs = [], []

    for s in [x for x in a.samples.split(",") if x]:
        sd = os.path.join(a.results, s)
        regions, gec, res = antismash_regions(sd, s), gecco(sd), rgi(sd, s)
        g = genome_info(sd, s)
        g.update(taxa.get(s, {}))
        g.update({"BGCs_antismash": len(regions), "BGCs_gecco": len(gec), "resistance_genes": len(res)})
        genomes.append(g)

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
    write(os.path.join(a.out, "bgc_summary.tsv"), bgcs)

    # per-isolate copies: results/<isolate>/11_gtdbtk and 12_novelty
    for g in genomes:
        s = g["isolate"]
        write(os.path.join(a.results, s, "12_novelty", f"{s}_bgc_novelty.tsv"),
              [b for b in bgcs if b["isolate"] == s])
        write(os.path.join(a.results, s, "11_gtdbtk", f"{s}_taxonomy.tsv"),
              [{k: g.get(k, NA) for k in ("isolate", "gtdb_classification", "closest_reference",
                                          "closest_ANI", "possible_new_species")}])
    n_new = sum(b["novelty"] == "novel_candidate" for b in bgcs)
    print(f"{len(genomes)} genomes, {len(bgcs)} BGCs, {n_new} novel candidates -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
PYTHON
}

# =============================================================================
#  MAIN
# =============================================================================
usage() { sed -n '3,61p' "$SCRIPT_PATH" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

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
