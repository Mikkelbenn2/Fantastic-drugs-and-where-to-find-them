# Positive and negative controls

Six genomes with **published, complete reference sequences** and BGCs **characterised in MIBiG 4.0**, used to check that `setup/bgc_pipeline.sh` finds what it should (positives) and does not invent novel antibiotic candidates (negatives).

## The panel

| Isolate name (in `raw_data/`) | Control | Reads | Source | Why |
|---|---|---|---|---|
| `pos_Bvelezensis_FZB42` | positive | simulated R10.4.1, 60x | NCBI GCF_000015785.2 (complete, 3.9 Mb) | Best-studied *Bacillus* producer: 11 MIBiG BGCs, e.g. surfactin, bacillomycin D, fengycin, difficidin, macrolactin, bacillaene, plantazolicin, amylocyclicin |
| `pos_Scoelicolor_A3-2` | positive | simulated R10.4.1, 60x | NCBI GCF_000203835.1 (complete, 9.1 Mb, linear chromosome + SCP1/SCP2) | Model *Streptomyces*: 17 MIBiG BGCs, e.g. actinorhodin, undecylprodigiosin, CDA, methylenomycin (on SCP1) |
| `pos_Pprotegens_Pf-5` | positive | simulated R10.4.1, 60x | NCBI GCF_000012265.1 (complete, 7.1 Mb) | Gram-negative biocontrol strain: pyoluteorin, orfamide, pyoverdine, toxoflavin |
| `pos_Svenezuelae_ATCC10712` | positive | simulated R10.4.1, 60x | NCBI GCA_000253235.1 (FR845719, 8.2 Mb) | Chloramphenicol, jadomycin, venezuelin, watasemycin/thiazostatin, gaburedins |
| `neg_Ecoli_K12_MG1655` | negative | **real** ONT, Kit14 (R10.4.1), Dorado sup v5, ~508 Mb (~110x) | ENA **SRR30993121** (BioProject PRJNA1172918, wild type) | Lab strain with no antibiotic BGC; only enterobactin (BGC0002476), which must come out as `known` |
| `neg_Buchnera_APS` | negative | simulated R10.4.1, 60x | NCBI GCF_000009605.1 (0.66 Mb) | Reduced endosymbiont genome with essentially no specialised metabolism: anything called is a false positive |
| `neg_Ecoli_K12_MG1655_sim` | technical | simulated R10.4.1, 60x | NCBI GCF_000005845.2 | Same genome as the real *E. coli* reads, so you can see whether simulated reads behave like real ones |

The expected BGCs per positive strain are in `expected_bgcs.tsv`. That table was taken directly from the MIBiG 4.0.1 JSON release (Zenodo record 14835872), using only `active` entries whose organism is exactly that strain.

**Why simulated reads for the positives?** I searched ENA for ONT reads of these producer strains. None had good-quality R10.4.1 data. The *S. venezuelae* ATCC 10712 wild-type run (SRR33220699) has a median read quality of about Q5, and other producer runs were about Q12–14 (older chemistry). This pipeline no longer polishes, so those reads would give frameshifted genes. That would test the reads, not the pipeline. Badread's `nanopore2023` model is trained on real R10.4.1 reads (Wick, github.com/rrwick/Badread). The real *E. coli* run (median Q ≈ 23, mean read length ≈ 11 kb) and its simulated twin show how close the simulation is.

## How to run

Run the controls in a **separate copy** of the repository, so they do not mix with your isolates. Re-use the installed programs and databases:

```bash
git clone https://github.com/Mikkelbenn2/Fantastic-drugs-and-where-to-find-them controls_run
ln -s "$PWD/Fantastic-drugs-and-where-to-find-them/setup/envs" controls_run/setup/envs
ln -s "$PWD/Fantastic-drugs-and-where-to-find-them/setup/databases" controls_run/setup/databases
bash controls_run/controls/get_controls.sh controls_run/raw_data      # needs internet; Badread takes ~0.5-2 h per genome
bash controls_run/setup/bgc_pipeline.sh submit
# when all jobs are done:
python3 controls_run/controls/check_controls.py controls_run/results/_combined/summary
```

## What counts as a pass

| Check | Expected |
|---|---|
| Positives, antiSMASH | Most `expected_bgcs.tsv` entries recovered as the closest MIBiG hit, above all the large NRPS/PKS clusters (they are MIBiG reference entries from these very strains) |
| Positives, other tools | Fewer than antiSMASH is normal: GECCO, DeepBGC and BGC-Prophet are compared with MIBiG by the pipeline's DIAMOND step. A miss by every tool (`missed_by_all`) for a large NRPS/PKS cluster is a real problem |
| Positives, `novelty` | Expected clusters labelled `known` or `related`, **not** `novel_candidate` |
| `neg_Ecoli_K12_MG1655` | Enterobactin locus `known`; no `novel_candidate` with score ≥ 5 |
| `neg_Buchnera_APS` | No loci at all, or only DeepBGC-only/BGC-Prophet-only loci with a low score |
| Real vs simulated *E. coli* | Same loci with both. If they differ, the simulated positives are less trustworthy |
| QC | CheckM2 completeness ≥ 95 %, contamination ≤ 5 % for every control; QUAST genome length close to the reference |

**Limits.** The checker counts a BGC as found only when it is the *top* MIBiG hit of a predicted BGC (as antiSMASH reports it). Neighbouring clusters that antiSMASH merges into one region, such as the bacillaene and other *Bacillus* lipopeptide regions, can therefore show as missed. Look those up in `per_tool/antismash_bgcs.tsv` before calling them misses.
