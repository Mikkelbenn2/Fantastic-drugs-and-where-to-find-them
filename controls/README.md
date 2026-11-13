# Positive and negative controls

Real, published data only. They check that `setup/bgc_pipeline.sh` finds what it should (positives) and does not invent novel antibiotic candidates (negatives).

- **Positives** are strains whose antibiotics are proven in the literature and whose BGCs are characterised in MIBiG 4.0.
- **Negatives** are strains with no known antibiotic production.
- All controls except one are the **complete reference genome** from NCBI. The exception is one real R10.4.1 ONT read set (*E. coli* K-12), which also tests read QC and assembly.

Every accession, genome paper and paper showing (non-)production is in **`control_panel.tsv`**, with DOIs. The expected MIBiG BGCs per positive are in **`expected_bgcs.tsv`**.

## The panel

### Positives: proven antibiotic producers

| Isolate name | Strain | Data | Antibiotics (evidence) | MIBiG BGCs |
|---|---|---|---|---|
| `pos_Bvelezensis_FZB42` | *Bacillus velezensis* FZB42 | GCF_000015785.2 | difficidin, bacilysin, bacillaene, macrolactin, lipopeptides ([Chen 2007](https://doi.org/10.1038/nbt1325); [Chen 2009](https://doi.org/10.1016/j.jbiotec.2008.10.011)) | 11 |
| `pos_Scoelicolor_A3-2` | *Streptomyces coelicolor* A3(2) | GCF_000203835.1 (chromosome + SCP1 + SCP2) | actinorhodin, undecylprodigiosin, CDA, methylenomycin ([Bentley 2002](https://doi.org/10.1038/417141a)) | 17 |
| `pos_Pprotegens_Pf-5` | *Pseudomonas protegens* Pf-5 | GCF_000012265.1 | pyoluteorin, DAPG, pyrrolnitrin, orfamide, rhizoxin ([Paulsen 2005](https://doi.org/10.1038/nbt1110); [Loper 2012](https://doi.org/10.1371/journal.pgen.1002784)) | 4 |
| `pos_Svenezuelae_ATCC10712` | *Streptomyces venezuelae* ATCC 10712 | GCA_000253235.1 | chloramphenicol ([Ehrlich 1947](https://doi.org/10.1126/science.106.2757.417)); genome: [Pullan 2011](https://doi.org/10.1186/1471-2164-12-175) | 5 |
| `pos_Sgriseus_NBRC13350` | *Streptomyces griseus* NBRC 13350 | GCF_000010605.1 | streptomycin ([Schatz, Bugie & Waksman 1944](https://doi.org/10.1097/01.blo.0000175887.98112.fe)); genome: [Ohnishi 2008](https://doi.org/10.1128/JB.00204-08) | 14 |
| `pos_Serythraea_NRRL2338` | *Saccharopolyspora erythraea* NRRL 2338 | GCF_000062885.1 | erythromycin ([Oliynyk 2007](https://doi.org/10.1038/nbt1297)) | 5 |
| `pos_Ecoli_Nissle1917` | *Escherichia coli* Nissle 1917 | GCF_000714595.1 | microcins M and H47 ([Patzer 2003](https://doi.org/10.1099/mic.0.26396-0)); genome: [Reister 2014](https://doi.org/10.1016/j.jbiotec.2014.07.442) | 1 |

### Negatives: no known antibiotic production

| Isolate name | Strain | Data | Why it is a negative |
|---|---|---|---|
| `neg_Ecoli_K12_MG1655` | *E. coli* K-12 MG1655 | **real reads**: ENA SRR30993121 (Kit14/R10.4.1, Dorado sup v5, ~508 Mb ≈ 110x; BioProject PRJNA1172918) | Lab strain; its only MIBiG entry is the siderophore enterobactin (BGC0002476), which must come out as `known` ([Blattner 1997](https://doi.org/10.1126/science.277.5331.1453)) |
| `neg_Ecoli_K12_MG1655_ref` | *E. coli* K-12 MG1655 | GCF_000005845.2 | Technical control: the published genome of the same strain as the reads, to see what assembly changes |
| `neg_Buchnera_APS` | *Buchnera aphidicola* APS | GCF_000009605.1 | Endosymbiont with a 0.66 Mb genome and almost no secondary metabolism ([Shigenobu 2000](https://doi.org/10.1038/35024074)) |
| `neg_Rprowazekii_MadridE` | *Rickettsia prowazekii* Madrid E | GCF_000195735.1 | Obligate intracellular pathogen, 1.1 Mb, relies on its host's metabolism ([Andersson 1998](https://doi.org/10.1038/24094)) |
| `neg_Cglutamicum_ATCC13032` | *Corynebacterium glutamicum* ATCC 13032 | GCF_000011325.1 | High-GC actinobacterium used industrially for amino acids, with no antibiotic described ([Ikeda & Nakagawa 2003](https://doi.org/10.1007/s00253-003-1328-1)). It tests whether the ML tools call false BGCs in high-GC genomes like *Streptomyces* |

**Why these strains?**

- **Positives** cover the main producer groups: *Bacillus*, *Pseudomonas*, *Streptomyces* (3 species), a non-*Streptomyces* actinomycete (*Saccharopolyspora*) and *E. coli*. Each makes an antibiotic proven in a paper, and most have several MIBiG BGCs from exactly that strain.
- **Nissle vs K-12** is the hardest test: the same species, but one makes small antibacterial microcins and the other doesn't.
- **Negatives** cover low-GC reduced genomes (*Buchnera*, *Rickettsia*), a free-living Gram-negative (*E. coli* K-12) and a free-living high-GC actinobacterium (*C. glutamicum*).
- ***Mycoplasma* was left out** on purpose. It uses genetic code 4, but the pipeline assumes code 11, so its genes would be called wrongly.

**Why genomes and not reads for most of them?**

- A `.fna` file in `1_raw_data/` is used as a finished genome: the pipeline skips read QC and assembly and runs everything from QUAST onwards. The BGC part is then tested on exactly the published sequence.
- In ENA (searched 2026-11-08), none of these strains has good R10.4.1 ONT reads. For example, the *S. venezuelae* ATCC 10712 run SRR33220699 has a median quality of about Q5. The other ONT runs are MinION runs of other strains or of mutants.
- The read and assembly part is tested by the real *E. coli* K-12 run (median Q ≈ 23). Comparing it with `neg_Ecoli_K12_MG1655_ref` shows what assembly changes.

**Expected BGCs.** `expected_bgcs.tsv` is taken directly from the MIBiG 4.0.1 JSON release (Zenodo record 14835872). It uses only `active` entries whose organism is exactly that strain.

## How to run

The controls run with the normal pipeline script, from the main project folder. They get their own input folder `2_raw_data_controls/` and output folder `4_results_controls/`, so they never mix with your isolates in `1_raw_data/` and `3_results/`. They use the same installed programs and databases.

```bash
cd /projects/students/bio-26-BT-7-1/Fantastic-drugs-and-where-to-find-them
bash setup/bgc_pipeline.sh controls --dry-run     # downloads the controls if missing (~10 min), shows what will run
bash setup/bgc_pipeline.sh controls               # downloads if missing, then submits the Slurm jobs
# when "squeue --me" shows no bgc_ jobs any more:
bash setup/bgc_pipeline.sh check-controls
```

The checker writes three tables into `4_results_controls/_combined/summary/` and also prints them:

| File | What it shows |
|---|---|
| `controls_recovery.tsv` | Per positive: how many expected MIBiG BGCs each tool found, and which were missed by all |
| `controls_expected_bgcs.tsv` | Every expected MIBiG BGC, and which tools found it |
| `controls_loci.tsv` | Every control: loci, novelty labels, best score, and `PASS`/`CHECK`/`NOT_IN_RUN` |

**Figures for the report** (R with ggplot2, dplyr, tidyr, readr, forcats, patchwork), after `check-controls`:

```bash
Rscript controls/plot_controls.R 4_results_controls/_combined/summary
```

PDF + PNG in `4_results_controls/_combined/summary/figures/`:

| File | Content |
|---|---|
| `fig_controls_overview` | A: known BGCs recovered per positive control and tool; B: BGC loci per control by novelty, with the verdict; C: priority score of every locus (dashed line = 5) |
| `fig_recovery` | Panel A alone, for the main text |
| `fig_expected_bgcs` | Every expected MIBiG BGC and which tool found it (supplementary) |

On your own computer: copy the summary folder first (see "Explore the results on your own computer" in the main README), then `Rscript controls/plot_controls.R controls_summary`. In RStudio: open `controls/plot_controls.R`, set `summ_dir` and `controls_dir` at the top to your folders, and press Source.

If it stops with "master_summary.tsv not found", the pipeline's last (combined) job has not finished or failed. Look in `4_results_controls/slurm_logs/combined_*.out`.

## Performance metrics: confusion matrix, sensitivity, specificity, ROC/AUC

`check-controls` also runs `controls/evaluate_controls.py`. It writes `eval_*.tsv` next to the other tables, the dashboard's Controls page shows them, and `plot_controls.R` draws them as `fig_evaluation` (A: ROC curves, B: confusion matrix, C: sensitivity per tool). It answers two questions, each with its own truth:

| Question | Unit | Truth | Metrics |
|---|---|---|---|
| **1. Detection:** are the known BGCs found? | expected MIBiG BGC of a positive control | `expected_bgcs.tsv` | Sensitivity (recall) per tool and for any tool, with 95 % CI (Wilson). Over-calling as loci per Mb in the negatives |
| **2. Antibiotic ranking:** does the evidence rank antibiotic BGCs high? | BGC locus of a control | **positive** = matches an expected BGC with `antimicrobial = yes`; **negative** = every locus of a negative control, and loci matching only non-antimicrobial BGCs (siderophores, pigments, terpenes …) | Confusion matrix; sensitivity, specificity, precision, NPV, accuracy, balanced accuracy, F1, MCC; ROC curve and AUC (with 95 % CI, Hanley & McNeil) |
| **3. Genome level** | strain | positive vs negative control | 2 × 2 table: predicted producer = at least 1 novel candidate with score ≥ 5. Only 12 genomes, so treat it as a sanity check |

How question 2 works:
- **Truth column.** `antimicrobial` in `expected_bgcs.tsv` is `yes` when the compound is a named antibiotic in the control's paper (`control_panel.tsv`) or MIBiG lists it as antibacterial/antifungal. Five borderline BGCs were checked in the literature and set by hand: toxoflavin and coelimycin are `yes`; erythrochelin, the SGR PTMs and the bacillothiazols are `no`. The column `antimicrobial_evidence` gives the source, with a DOI for the hand-set ones.
- **Score.** In the controls every BGC is known, so the priority score would always give it 0 novelty points. The evaluation therefore uses the **novelty-blind score**: the priority score without the novelty points. It asks: "would the evidence (type, tool agreement, resistance gene, DeepBGC activity) have ranked this BGC high if it were new?"
- **Flagged** uses the pipeline's own rule for a strong candidate, minus the novelty part. The locus must be supported (antiSMASH, or ≥ 2 tools), be of an antibiotic-relevant type (siderophores excluded), and have a blind score ≥ 2, so that a new locus would reach priority score ≥ 5. Results made before v3.6.0 have no `antibiotic_type` column; for those, only the blind score is used.
- **Left out.** Loci of producers that match no expected BGC are left out. They may be real, uncharacterised clusters, so their truth is unknown.
- **Misses.** Antimicrobial BGCs that no tool found count as false negatives.

**Limits.** The truth only covers characterised BGCs. A specificity for detection itself cannot be computed, because even the negative strains have real BGCs (e.g. enterobactin). The AUCs are based on about 25 positive loci, so their confidence intervals are wide.

## What counts as a pass

| Check | Expected |
|---|---|
| Positives, antiSMASH | Most `expected_bgcs.tsv` entries recovered as the closest MIBiG hit, above all the large NRPS/PKS clusters and the named antibiotics (streptomycin, erythromycin, chloramphenicol, actinorhodin, difficidin, pyoluteorin) |
| Positives, other tools | Fewer than antiSMASH is normal. A miss by every tool (`missed_by_all`) for a large NRPS/PKS cluster is a real problem |
| Positives, `novelty` | Expected clusters labelled `known` or `related`, **not** `novel_candidate` |
| `pos_Ecoli_Nissle1917` vs `neg_Ecoli_K12_MG1655_ref` | Nissle has loci that K-12 lacks (microcin, colibactin, yersiniabactin). Missing the small microcin cluster is a known weak spot of BGC tools, so write it down rather than treating it as a failure |
| `neg_Ecoli_K12_MG1655` (reads) | Enterobactin locus `known`; no `novel_candidate` with score ≥ 5 |
| `neg_Buchnera_APS`, `neg_Rprowazekii_MadridE` | No loci at all, or only DeepBGC-only/BGC-Prophet-only loci with a low score |
| `neg_Cglutamicum_ATCC13032` | At most a few `known`/terpene loci; no `novel_candidate` with score ≥ 5. Many ML-only loci here would show that those tools over-call in high-GC genomes |
| Real reads vs published genome, *E. coli* | Same loci with both. If they differ, the assembly from reads is the cause |
| QC | CheckM2 completeness ≥ 95 % and contamination ≤ 5 % for the free-living strains. *Buchnera* and *Rickettsia* can score lower because of their reduced genomes, which is expected. QUAST genome length equal to the reference |

**Counting rule.** A known BGC counts as found when it is among a predicted BGC's MIBiG matches with ≥ 30 % similarity (column `MIBiG_all_hits`), not only the best match, so clusters that antiSMASH merges into one region still count. Target: ≥ 90 % with any tool.

**Limits.**
- **Absence of evidence:** a negative means "no antibiotic known". Proving that a strain makes none is not possible.
