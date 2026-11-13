#!/usr/bin/env bash
# One-off tidy-up of the BioCloud project folder for pipeline v3.8 (numbered folders, all code in setup/).
#   bash setup/tidy_project.sh        shows what would be moved (changes nothing)
#   bash setup/tidy_project.sh --do   moves it
# Only MOVES (mv -n: never overwrites); nothing is deleted. Old and unused things go to 0_junk_yard/,
# delete that folder yourself when you are sure. Refuses to run while bgc_ Slurm jobs are running.
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DO=0; [[ "${1:-}" == "--do" ]] && DO=1
cd "$P"

if command -v squeue >/dev/null && [[ -n "$(squeue --me -h -o '%j' 2>/dev/null | grep '^bgc_' || true)" ]]; then
  echo "Pipeline jobs are running (squeue --me). Wait until they are done, then run this again." >&2; exit 1
fi

J=0_junk_yard
move() {   # move SRC DEST: only if SRC exists and DEST does not
  local src=$1 dst=$2
  [[ -e "$src" || -L "$src" ]] || return 0
  if [[ -e "$dst" ]]; then echo "  skip  $src  ($dst already exists)"; return 0; fi
  echo "  move  $src  ->  $dst"
  if [[ $DO -eq 1 ]]; then mkdir -p "$(dirname "$dst")" && mv -n "$src" "$dst"; fi
}

echo "Project: $P"; [[ $DO -eq 1 ]] || echo "(dry run: nothing is changed; add --do to move)"
echo "== junk yard"
move "0. junk_yard" "$J"
[[ $DO -eq 1 ]] || J="0. junk_yard"            # dry run: the rename above has not happened
move "$J/metagenome_pipeline.sh"  "0_junk_yard/metagenome_pipeline/metagenome_pipeline.sh"
move "$J/setup_metagenomes.log"   "0_junk_yard/metagenome_pipeline/setup_metagenomes.log"
move "$J/controls_run"            "0_junk_yard/old_runs/controls_run"
move "0. raw_data_mags"           "0_junk_yard/metagenome_pipeline/raw_data_mags"
move "0. results_metagenomes"     "0_junk_yard/metagenome_pipeline/results_metagenomes"
move "tests"                      "0_junk_yard/old_files/tests"            # old copy; the real tests are setup/tests
move "controls_summary.tar.gz"    "0_junk_yard/old_files/controls_summary.tar.gz"
move "gbk_to_table.py"            "0_junk_yard/old_files/gbk_to_table.py"
move "setup/bgc_pipeline_v3.1.0_backup.sh" "0_junk_yard/old_files/bgc_pipeline_v3.1.0_backup.sh"
move "Rplots.pdf"                 "0_junk_yard/old_files/Rplots.pdf"
# v3.8.1: controls/ and dashboard/ live in setup/. Once the new copies are there, the old ones are not used.
[[ -d setup/controls ]] && move "controls" "0_junk_yard/old_files/controls_before_v3.8.1"
[[ -d setup/dashboard ]] && move "dashboard" "0_junk_yard/old_files/dashboard_before_v3.8.1"
echo "== numbered data and result folders"
move "raw_data"           "1_raw_data"
move "raw_data_controls"  "2_raw_data_controls"
move "results"            "3_results"
move "results_controls"   "4_results_controls"

[[ $DO -eq 1 ]] || { echo "Nothing changed. To do it:  bash setup/tidy_project.sh --do"; exit 0; }

# Links inside the results that point to the old folder names (e.g. reads, Bakta proteins): repoint them
echo "== repointing links to the old folder names"
n=0
while IFS= read -r l; do
  t=$(readlink "$l"); new=$t
  for pair in "raw_data_controls/:2_raw_data_controls/" "raw_data/:1_raw_data/" \
              "results_controls/:4_results_controls/" "results/:3_results/"; do
    old=${pair%%:*}; nw=${pair##*:}
    [[ "$new" == "$P/$old"* ]] && { new="$P/$nw${new#"$P/$old"}"; break; }
  done
  [[ "$new" != "$t" ]] && { ln -sfn "$new" "$l"; n=$((n + 1)); }
done < <(find 3_results 4_results_controls -type l 2>/dev/null)
echo "  $n links updated"

echo "== dashboards"
bash setup/bgc_pipeline.sh dashboard || true
echo "Done. Main folder now:"; ls -1 "$P"
