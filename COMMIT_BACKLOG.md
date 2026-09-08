# Commit Backlog

This file tracks the history of commits made by the AI agent, providing a high-level summary of the changes and the reasoning behind them.

## [2026-09-08] - `curation.isoquant_args`: Pass-Through Options for IsoQuant
- **Branch**: `dev-IsoQuant`
- **Goal**: Let users customise IsoQuant without adding a pipeline key per edge case. First use: Iso-Seq FLNC reads have their poly(A) tails removed by `isoseq refine`; IsoQuant (3.10 and 3.13.1 checked) builds novel mono-exon transcripts only from reads with a detected tail and uses the tail as the only strand evidence for unspliced reads, so on Arabidopsis 3 of 15,080 models were mono-exon and Tier 1 contained no single-exon gene. `--polya_trimmed all` (artificial poly(A) at the read 3′ end from the mapped strand; valid for 5′→3′-oriented reads only) is the fix; `--polya_requirement never` would not have helped (it only relaxes reference mono-exon isoforms).
- **Summary**:
    - New key `curation.isoquant_args` (default `""`), appended verbatim to the `isoquant` call in `rule run_isoquant`.
    - `scripts/input_check.py`: `validate_isoquant_args()` rejects options the rule already sets (`--reference/-r`, `--data_type/-d`, `--prefix/-p`, `--threads/-t`, `-o/--output`, input flags); enforced both pre-flight and at Snakemake parse time via the shared module.
    - Documented in `config.yaml` and `README.md` with the FLNC example `"--polya_trimmed all --stranded forward"` and the orientation caveat; logged in the configuration summary. Tests: +8 (60 total).
    - Design decision (Pablo): no dedicated key for `--polya_trimmed`; users add it through `isoquant_args`.

## [2026-09-08] - Split-Mode Augustus Rules: QoS and Runtime From the Same Tier
- **Branch**: `dev-IsoQuant`
- **Goal**: Fix the human full run failing at submission with `QOSMaxWallDurationPerJobLimit`.
- **Summary**: `rules/split_augustus.smk` (`ed_augustus_per_chromosome`, `ab_augustus_per_chromosome`) requested `--qos` from the `small` tier (`short`, 1-day limit) but `runtime` from the `big` tier. With `big.time: 24h` this sat exactly at the limit; the human config raised `big.time` to 48h for IsoQuant and every Augustus job was rejected. Both rules now take QoS and runtime from `big`; CPU/memory stay on `small` (Augustus is single-threaded, ~2.5 GB). Audited the other rules: no other qos/runtime tier mix.

## [2026-09-07] - Miniprot Alignment Rule on the Medium Resource Tier
- **Branch**: `dev-IsoQuant`
- **Goal**: Prepare the human full run: `align_proteins_miniprot` used the `small` tier (2 CPUs / 8 GB), enough for Arabidopsis but not for indexing a 3 Gb genome.
- **Summary**: `rules/evidence_driven.smk`: `align_proteins_miniprot` now uses `resources.medium` (cpus, mem, time, qos) and `threads` accordingly.

## [2026-09-07] - BUSCO Gene Ids Unified (mixed training now really mixed) and Deterministic Augustus Training (roadmap 1.7b, 1.4)
- **Branch**: `dev-IsoQuant`
- **Goal**: Make BUSCO genes actually reach the training set, and make `etraining` reproducible from an identical training set.
- **Summary**:
    - Diagnosis: BUSCO's miniprot GFF keys records by `ID=MP######`/`Parent=MP######`, while `busco_complete_aa.py` names proteins by BUSCO id (`10052at3699`), which is the prefix of the GFF `Target=` attribute. `filter_busco_overlaps.py` and `assemble_training_gff.py` keyed genes by the MP id, so the FAA subset was always empty (0 of 727 on Arabidopsis) and no BUSCO gene ever matched a CD-HIT representative: every `mixed` run trained on SQANTI genes only, and `busco_only` would have produced an empty training set through the same path.
    - Fix: both scripts derive the gene key from `Target=` (`busco_gene_key`), keep an ID→key map so Parent-only children (`stop_codon`) follow their mRNA (`resolve_gene_key`), and fall back to `gene_id`/`ID=`/`Parent=` for other GFFs. `filter_and_write_faa` raises if retained genes exist but no FAA sequence matches. Real-data check: 3,199 BUSCO genes parsed, 727 retained, 727 sequences written.
    - 1.4 follow-up: `export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0` before `gff2gbSmallDNA.pl`, `new_species.pl`, `filterGenes.pl`. Perl's per-process hash randomisation swapped 9 of 4,696 GenBank loci between two runs of the same training set on Arabidopsis and changed the `.pbl` model files in the 3rd decimal.
    - Tests: +4 (real-shaped ids through filter and assembly, Parent-only children, GFF/FAA mismatch error); 52 total.
    - Consequence: the next `mixed` run will be the first with BUSCO genes in the training set (≈727 + SQANTI up to `test_size`); the trained model will differ from all previous runs.

## [2026-09-07] - Miniprot Paralog Recovery Becomes the Default
- **Branch**: `dev-IsoQuant`
- **Goal**: Adopt the validated Miniprot settings as the pipeline default (roadmap 1.8 / 1.6 validation).
- **Summary**:
    - Two Arabidopsis reruns (Garnatxa jobs 3038924, 3039034) vs the custom-filter baseline, CDS level: locus Sn/Sp 55.9/74.2 → 57.5/75.8 (`--outs=0.5`) → 57.7/76.2 (+ fixed hint converter); transcript 39.3/67.9 → 40.6/69.0; missed reference loci 7,198 → 7,072. Tier 2 genes with protein support 28 → 1,459; fully supported multi-exon genes 84% exact vs 38% unhinted.
    - `prediction.miniprot_args` default changed from `"-N 5 -p 0.6"` to `"-N 30 -p 0.6 --outs=0.5"` in `rules/setup/functions.smk`, `config.yaml`, `README.md`.
    - Register of all Arabidopsis variants and their metrics: `genomes_benchmark/Arabidopsis_tests.md` on Garnatxa (gitignored benchmark area).

## [2026-09-07] - Tier Resolver: Nested-Interval-Safe Overlap Search, Log Output (roadmap 1.5)
- **Branch**: `dev-IsoQuant`
- **Goal**: Fix the overlap search in `scripts/resolve_transcript_tiers.py`, which could let Augustus genes nested inside long SQANTI genes through as duplicated loci, and make the rule log useful.
- **Summary**:
    - Bug: intervals were sorted by start but binary-searched over their *ends*, which are not sorted when a long gene contains a shorter one; the search then jumped past the containing gene and reported no overlap. Same construction in the monoexon hint-support check. Confirmed on the Arabidopsis run (about ten leaked genes; small there, larger on genomes with nested/overlapping loci).
    - Fix: `build_interval_index()` (sorted starts + running maximum of ends, built once per chrom/strand) and `index_overlaps()`; `check_locus_overlap` / `check_hint_support` kept as wrappers. Identical answers on non-nested data.
    - `resolve_tiers()` validates `filter_mode`, writes its summary to the rule's log file (was empty: `script:` rules do not capture stdout) and returns the counts; new `main_snakemake()` entry point.
    - `tests/test_resolve_transcript_tiers.py`: +8 tests (nested locus, first-longer-than-second, opposite strand, nested hints, 3,000-query brute-force cross-check, end-to-end nested discard with log assertion, invalid filter_mode, real Snakemake entry point).

## [2026-09-07] - Miniprot Hints: Unique Groups and Evidence-Gated Start/Stop (roadmap 1.6)
- **Branch**: `dev-IsoQuant`
- **Goal**: Stop the protein-hint converter from asserting gene boundaries Miniprot never claimed, and from merging paralog loci into one Augustus hint group. Both defects were dormant while Miniprot wrote almost no secondary alignments (`--outs` default 0.99) and became live with `--outs=0.5`: in the experiment run, Tier 2 genes with *all* introns P-supported were 84.8% exact, those with only *some* supported 25.9% (mostly class j).
- **Summary**:
    - `scripts/miniprot_to_hints.py`: `grp=` is now the Miniprot alignment `ID` (unique per alignment) instead of the protein name; the protein is kept as `target=`. `start` is emitted only when the aligned protein range begins at residue 1 (`Target` attribute); `stop` only when Miniprot emitted a `stop_codon` feature, at that feature's coordinates (Miniprot CDS features include the stop codon, so this coincides with the previous last-codon position for complete alignments). `Rank` and `Identity` parsed and kept. Missing/empty input now warns. Counts of suppressed start/stop hints reported.
    - No cross-check against TD2 ORFs: `lrRNA` start/stop hints (from TD2 via `generate_hints.py`) are untouched; the change only removes `P` boundary hints for partial alignments.
    - Real-data check (Arabidopsis, 18,992 alignments): CDSpart 134,065 and intron 115,073 unchanged; start 18,992 → 17,707; stop 18,992 → 16,277; hint groups 12,447 → 18,992.
    - `tests/test_miniprot_to_hints.py` rewritten (8 tests): complete +/- strand, partial secondary alignment (no start/stop), N-terminal-complete without stop codon (start only), two alignments of one protein (distinct groups), numeric score column, empty input, Snakemake entry point.
    - `prediction.miniprot_args` default left at `"-N 5 -p 0.6"` until the rerun with the fixed converter is validated; then `-N 30 -p 0.6 --outs=0.5` becomes the default.

## [2026-09-07] - Configurable Miniprot Options for Self-Protein Hints
- **Branch**: `dev-IsoQuant`
- **Goal**: Make the Miniprot call configurable so the paralog-recovery experiment (`--outs`) can be run from a config file without editing rules.
- **Summary**:
    - Diagnosis on the Arabidopsis run: Miniprot 0.18 only writes secondary alignments scoring >= `--outs` x best (default 0.99), so `-N 5 -p 0.6` produced 12,620 alignments for 12,447 proteins and only 0.2% of `src=P` hints fell outside Tier 1 loci; 10,866 of 10,894 Tier 2 transcripts had no hint support (47% exact CDS match vs 86-90% where P hints exist).
    - New config key `prediction.miniprot_args` (default `"-N 5 -p 0.6"`, i.e. unchanged behaviour) consumed by `rule align_proteins_miniprot`; documented in `config.yaml`, `README.md`, and the configuration log summary.

## [2026-09-07] - Idempotent Augustus Species Handling and Deterministic Training Set
- **Branch**: `dev-IsoQuant`
- **Commit**: (pending)
- **Goal**: Roadmap 1.3 and 1.4: reruns behave like first runs, and the training set does not depend on filesystem listing order.
- **Summary**:
    - `rule new_species`: fails fast with a clear message if `$AUGUSTUS_CONFIG_PATH/species` is unset or not writable (checked inside the Augustus env, where the variable is defined); logs where the species is created; removed the always-`None` `augustus_dir` param.
    - `scripts/modify_SC_freq.py`: rewritten. Stop-codon parameters are matched by name and their value column replaced regardless of current content (the old literal `0.33`/`0.34` substitution was a silent no-op on rerun); raises if the three codon frequencies or the three parameter lines are not found; writes a copy of the modified `<species>_parameters.cfg` into the training directory; CLI and Snakemake entry points.
    - `rule modify_stop_codon_freq`: new `params_cfg` output alongside the `.done` sentinel.
    - `scripts/busco_complete_aa.py`: sorted directory listing; `rule concatenate_gff`: sorted concatenation. Both keep CD-HIT input order, and therefore cluster representatives and the training set, identical between runs.
    - `tests/test_modify_SC_freq.py` (7 tests, incl. rerun and non-default template values).

## [2026-09-07] - Wrapper Exit Code and Configuration Pre-flight Hardening
- **Branch**: `dev-IsoQuant`
- **Commit**: (pending)
- **Goal**: Roadmap 1.1 and 1.2: make failures visible to Slurm chaining and reject invalid configurations before Snakemake starts.
- **Summary**:
    - `sqanti_evidence`: returns Snakemake's exit code (`sys.exit`), and adds its own directory to `sys.path` so it can be invoked from any working directory.
    - `scripts/input_check.py`: enum validation for `training.mode`, `prediction.mode`, `prediction.filter_mode`, `curation.mode`, `curation.data_type`; `curation.mode: ab_initio` requires `training.mode: busco_only` (explanatory error); `training.lineage` required for `mixed`/`busco_only`; optional files (`hint_config`, `hint_weights`, `filter_rules`, `reference_gtf`, `prediction_genome`) must exist when set; `TypeError` handled in numeric checks; `prediction.species` required; writability check of the tools directory (and of `AUGUSTUS_CONFIG_PATH` if exported) for the Augustus species model.
    - `rules/setup/functions.smk`: imports `ALLOWED_VALUES` and `validate_modes` from `scripts/input_check.py` (single source of truth for enums and the `ab_initio`/`busco_only` rule) so `snakemake -n` reports the explanation instead of a `CyclicGraphException`; `user_gtf` required in `user` mode; `data_type` validated.
    - `tests/test_input_check.py` (15 tests) and `tests/dryrun/configs/busco_only_split_abinitio.yaml` fixture.

## [2026-09-07] - README Rewrite for the Tier-Resolver Architecture
- **Branch**: `dev-IsoQuant`
- **Commit**: (pending)
- **Goal**: Make the only version-controlled document describe the current pipeline (roadmap item 0.7).
- **Summary**:
    - Rewrote `README.md`: dual-track (Tier 1 / Tier 2) overview, correct installation (`conda env create -f env.yaml`), 6-block configuration tables matching `config.yaml`, IsoQuant input expectations, 8-step workflow including tier resolution, output tree, testing section, tool citations.
    - Renamed the conda environment in `env.yaml` from `wapisimo` to `sqanti_evidence`.

## [2026-09-07] - Config Template, Validator, Logger and DAG Regression Test (`f27d35f`)
- **Branch**: `dev-IsoQuant`
- **Commit**: `f27d35f`
- **Goal**: Facewash items 0.3–0.6 of the stabilisation roadmap: dead configuration removed, Garnatxa resource convention adopted, working dry-run test, rule tidy-up.
- **Summary**:
    - `config.yaml`: `hint_weights` points to the existing default cfg; personal `toolsdir` replaced by a placeholder; `prediction.filter_mode` and `curation.data_type` documented; unused `small_bigMem` tier removed.
    - `rules/setup/functions.smk`: resource defaults use `mem: "8GB"` / `time: "2h"`; `mem_mb`/`time_min` rewrite path deleted; missing resource keys filled per tier; `filter_mode` validation aligned with `resolve_transcript_tiers.py` (`strict`/`medium`/`none`, default `medium`).
    - `rules/setup/logging_setup.smk`: fixed wrong keys (`miniprot_id`, `flanking`) and extended the configuration summary.
    - `tests/dryrun/`: fixtures, six configs (3 training × 2 prediction modes) and `run_dryruns.sh` as the DAG regression test; outputs gitignored.
    - `rules/split_augustus.smk`: `ed_augusuts_per_chromosome` renamed to `ed_augustus_per_chromosome`; the file is now included once from `snakefile` instead of from both `ab_initio.smk` and `evidence_driven.smk`; scope comment added to `ab_initio.smk`.

## [2026-09-07] - Remove Orphaned Iso-Seq, TSEBRA and Legacy Files (`dc36bb7`)
- **Branch**: `dev-IsoQuant`
- **Commit**: `dc36bb7`
- **Goal**: Facewash item 0.2: delete code no rule references.
- **Summary**:
    - Deleted `rules/isoseq.smk`; TSEBRA leftovers (`scripts/genome_anno.py`, `scripts/evidence.py`, `scripts/filter_monoexons.py`, `envs/tsebra.yaml`, `envs/tsebra.cfg`); Iso-Seq era files (`scripts/fastq2bam.py`, `scripts/create_tama_filelist.py`, `scripts/rename_lima_output.py`, `envs/tama.yaml`, `envs/isoseq.yaml`, `envs/minimap2.yaml`, `envs/pacbio_mock.bam`); superseded training scripts (`scripts/concatenate_GFF.py`, `scripts/generate_subset.py`, `scripts/filterClassification.R`); old-schema `.dummy_config.yaml`.
    - Moved `scripts/download_sm_genomes.sh` to `genomes_benchmark/bin/` (outside version control).
    - Removed the dangling `prediction.tsebra_config` default from `config.yaml` and `functions.smk`.
    - Kept `envs/extrinsic.M.RM.PB.cfg` (valid low-bonus preset declaring `lrRNA` and `P`).

## [2026-09-03] - Integrate Snakemake-Native Hierarchical Tier Resolver
- **Branch**: `dev-IsoQuant`
- **Commit**: `7c678b7`
- **Goal**: Replace TSEBRA with native Snakemake hierarchical tier resolution (Tier 1 SQANTI pass-through + Tier 2 Augustus gap-filling), eliminating duplicate transcripts and precision collapse.
- **Summary**:
    - Created `scripts/resolve_transcript_tiers.py` supporting both Snakemake's native `script:` directive and standalone CLI execution.
    - Replaced `rule run_tsebra` in `rules/evidence_driven.smk` with `rule resolve_transcript_tiers`.
    - Updated `rules/quality_control.smk` to connect `resolved_prediction.gtf` to `agat_cleaning` and `gaqet2_setup`.
    - Added comprehensive unit tests in `tests/test_resolve_transcript_tiers.py` verifying Tier 1 preservation, overlapping Augustus discard, gap-filler retention, and monoexon noise filtering.
    - Updated `.gitignore` to ignore Python and Pytest caches.

## [2026-09-03] - Integrate TSEBRA Ensembling & Miniprot Paralog Recovery
- **Branch**: `dev-IsoQuant`
- **Goal**: Implement TSEBRA transcript selection to preserve curated SQANTI3 models without Augustus mutation, and Miniprot self-protein hinting for unexpressed paralog recovery.
- **Summary**:
    - Created `envs/tsebra.yaml` and `envs/tsebra.cfg` for automated transcript arbitration.
    - Implemented `rule run_tsebra` in `rules/evidence_driven.smk` combining SQANTI3 curated GTF, Augustus predictions, and unified hints.
    - Updated `rules/quality_control.smk` to feed `tsebra_prediction.gtf` directly into `agat_cleaning` and `gaqet2_setup`.
    - Implemented `scripts/miniprot_to_hints.py` to convert spliced protein alignments into `src=P` Augustus hints (with frame phase 0, 1, 2).
    - Added `align_proteins_miniprot`, `convert_miniprot_hints`, and `combine_evidence_hints` in `rules/evidence_driven.smk`.
    - Added source `P` and weights across all `envs/extrinsic*.cfg` files.
    - Added comprehensive unit tests in `tests/test_miniprot_to_hints.py`.

## [2026-07-15] - Resolve GffCompare Dot Suffix Output Naming Bug
- **Branch**: `dev-GffCompareFix`
- **Goal**: Prevent missing stats file errors during evaluation when a sample name contains dots.
- **Summary**:
    - Defined `sample_qc = sample.replace('.', '-')` in both `rules/quality_control.smk` and `snakefile` to safely isolate name formatting.
    - Updated `gffcompare_eval` and `gffcompare_cds_eval` to use `sample_qc` for output paths and `-o` prefix arguments.
    - Updated target outputs in `rule all` to check for `sample_qc` based stats files.

## [2026-07-14] - Integrate IsoQuant with Input Autodetection and Cleanup IsoSeq Traces
- **Branch**: `dev-IsoQuant`
- **Goal**: Implement IsoQuant rules and dynamic input flag selection, simplifying directories and functions.
- **Summary**:
    - Created `envs/isoquant.yaml` for conda environment.
    - Updated `rules/setup/directories.smk` to rename and remove old IsoSeq paths, retaining only `isoquant`.
    - Updated `rules/setup/functions.smk` to add `get_isoquant_input_flag` autodetection helper and default config for `isoquant.data_type`.
    - Rewrote `rules/transcript_modelling.smk` to only include `run_isoquant`, using autodetection.
    - Updated `rules/evidence_driven.smk` to take `isoforms` directly from IsoQuant's output.

## [2026-06-04] - Revert Transcript Preprocessing to Align-only Mode
- **Branch**: `dev-IsoSeq`
- **Goal**: Revert preprocessing pipeline to start directly from alignment, bypass read clustering.
- **Summary**:
    - Removed `cluster` rule from `rules/transcript_modelling.smk`.
    - Restored `get_pbmm2_input` function as the direct source for the `reads` parameter in the `mapping_reads_pbmm2` rule.

## [2026-06-04] - Configure Iso-Seq Read Clustering and Resource Parameterization
- **Branch**: `dev-hint_testing`
- **Goal**: Integrate consensus transcript clustering into the preprocessing phase and configure production resource properties.
- **Summary**:
    - Integrated `isoseq cluster2` into `rules/transcript_modelling.smk` under a new `cluster` rule.
    - Updated `mapping_reads_pbmm2` to utilize the clustered outputs (`{sample}.cluster.bam`) instead of raw input files.
    - Parameterized Slurm resource limits in `config.yaml` to use human-readable GB and hour units (e.g., `8GB`, `2h`) instead of `mem_mb`/`time_min`.
    - Renamed source key `PB` to `lrRNA` in `envs/extrinsic.hints_conf.Ale.cfg`.
    - Cleaned up `localrules` definition in `snakefile` to only register `all`.

## [2026-06-02] - Resource Specification for CDS Subset Rules
- **Branch**: `dev-hint_testing`
- **Goal**: Enable Slurm cluster execution for CDS subsetting rules by specifying memory, CPU, and time resources.
- **Summary**:
    - Added explicit `resources` definitions to `subset_reference_cds` and `subset_prediction_cds` rules in `rules/quality_control.smk`.

## [2026-06-02] - QC Output Path Refactoring
- **Branch**: `dev-hint_testing`
- **Goal**: Restructure quality control directories and targets to use gffcompare.
- **Summary**:
    - Updated `rules/setup/directories.smk` to replace directory `qc_agat` with `qc_gffcompare` and cleaned up unused directories (`qc_omark`, `qc_busco`).
    - Standardized rule output paths in `rules/quality_control.smk` to use the new `qc_gffcompare` directory.
    - Updated `rule all` in `snakefile` to track `gffcompare` stats.

## [2026-05-27] - Intron Hint Filtering by CDS Coordinates
- **Branch**: `dev-hint_testing`
- **Goal**: Prevent out-of-bounds intron predictions when exon hints are disabled.
- **Summary**:
    - Modified `scripts/generate_hints.py` to filter intron hints, writing them only if they fall within the minimum and maximum CDS boundaries of the transcript when exon-related hints are disabled.

## [2026-05-26] - Integrate CDS-only Comparison into Pipeline
- **Branch**: `dev-hint_testing`
- **Goal**: Automate CDS-only benchmark evaluations and prevent gffcompare dot-suffix output naming bugs.
- **Summary**:
    - Added `scripts/subset_cds.py` to subset annotations to coding sequences only.
    - Updated `rules/quality_control.smk` to add `subset_reference_cds`, `subset_prediction_cds`, and `gffcompare_cds_eval` rules, generating `{sample}_cds.stats`.
    - Added the `-T` flag to both `gffcompare_eval` and `gffcompare_cds_eval` to avoid intermediate file generation.
    - Modified `snakefile` to add `{sample}_cds.stats` to rule `all` inputs.

## [2026-05-21] - Refactor Monoexon Filtering Parameter Handling
- **Branch**: `dev-hint_testing`
- **Goal**: Improve robustness and explicit configuration handling for monoexon filtering.
- **Summary**:
    - Updated `rules/evidence_driven.smk` to pass `filter_mode` explicitly to the `filter_monoexons` script via rule parameters.
    - Added validation in `rules/setup/functions.smk` to ensure `prediction.filter_mode` is set to a valid value ("monoexon" or "all").
    - Simplified `scripts/filter_monoexons.py` by retrieving the filter mode from `snakemake.params` instead of the global configuration object.

## [2026-05-21] - Rule Refactoring and Alternative Splicing Support
- **Branch**: `dev-hint_testing`
- **Goal**: Improve code maintainability and enhance Augustus prediction accuracy using long-read evidence.
- **Summary**:
    - Performed a major style refactoring of `rules/ab_initio.smk` and `rules/evidence_driven.smk` to improve consistency and readability.
    - Enabled `--alternatives-from-evidence=true` in the `augustus_hints` rule (`rules/evidence_driven.smk`) to better leverage long-read evidence for predicting alternative splicing isoforms.
    - Cleaned up `config.yaml` by removing the redundant `training.skip` parameter.

## [2026-05-21] - Fix Default Filter Rules Path Resolution
- **Branch**: `dev-hint_testing`
- **Goal**: Fix a bug where the default filter rules path was resolving incorrectly in a Snakemake context.
- **Summary**:
    - Updated `rules/setup/functions.smk` to use `workflow.basedir` for resolving the `envs` directory, ensuring consistent path resolution across different execution environments (local and HPC).

## [2026-05-21] - Default SQANTI Filtering Rules
- **Branch**: `dev-hint_testing`
- **Goal**: Improve user experience by providing sensible defaults for transcript curation.
- **Summary**:
    - Modified `rules/setup/functions.smk` to automatically default to `envs/filter_rules.json` if `curation.filter_rules` is missing or points to a non-existent file.

## [2026-05-21] - Documentation and HPC Synchronization
- **Branch**: `dev-hint_testing`
- **Goal**: Improve project traceability and synchronize local and remote states.
- **Summary**:
    - Updated `GEMINI.md` to include the project's path on the Garnatxa HPC cluster.
    - Synchronized `rules/evidence_driven.smk` with the fix applied on Garnatxa (setting `ref_genome` to `config.project.genome`).

## [2026-05-20] - Stable Version Update
- **Branch**: `dev-hint_testing`
- **Commit**: `6b959e2`
- **Summary**:
    - Refined SQANTI3 filtering rules in `envs/filter_rules.json` (added `CDS_length` and updated TTS checks).
    - Simplified `rules/ab_initio.smk` by removing redundant temporary file definitions.
    - Updated `rules/evidence_driven.smk` to invoke `sqanti3_qc.py` and `sqanti3_filter.py` more idiomatically, removing hardcoded paths and manual `LD_LIBRARY_PATH` exports.
- **Goal**: Ensure pipeline robustness when users provide absolute output paths in the configuration.

## [2026-05-20] - Benchmarking Infrastructure Refactoring
- **Branch**: `dev-hint_testing`
- **Goal**: Enable the Human Chr19 benchmarking strategy by decoupling training and prediction genomes and integrating evaluation tools.
- **Summary**:
    - Introduced `prediction_genome` in `config.yaml` to allow genome-wide training with subset-specific prediction.
    - Refactored `split_augustus.smk` to respect the `augustus.config` parameter, enabling weight-sensitivity testing.
    - Integrated `gffcompare` into the workflow via `rules/quality_control.smk` (using a dedicated `envs/gffcompare.yaml`) for automated F1 score calculation against reference GTFs.
    - Updated `rule all` to dynamically include evaluation outputs when a reference GTF is provided.
## [2026-05-20] - Configuration Schema Transition (Revamp)
- **Status**: COMPLETED
- **Goal**: Reorganize `config.yaml` into a logical hierarchy to improve modularity and support benchmarking.
- **Implementation**:
    - Refactored `rules/setup/functions.smk` to validate the new 5-block structure (`project`, `training`, `prediction`, `curation`, `evaluation`) and handle default path resolutions.
    - Updated `scripts/input_check.py` to perform strict pre-launch validation on hierarchical keys.
    - Systematically updated all rules (`ab_initio.smk`, `evidence_driven.smk`, `split_augustus.smk`, `quality_control.smk`, `transcript_modelling.smk`) and the `snakefile` to use new configuration accessors.
- **Mapping Reference**:
    | Old Accessor | New Accessor |
    | :--- | :--- |
    | `config.required.genome` | `config.project.genome` |
    | `config.required.input` | `config.project.input` |
    | `config.required.outdir` | `config.project.outdir` |
    | `config.required.toolsdir` | `config.project.toolsdir` |
    | `config.required.prediction_genome` | `config.project.prediction_genome` |
    | `config.required.log_level` | `config.project.log_level` |
    | `config.ab_initio.lineage` | `config.training.lineage` |
    | `config.ab_initio.miniprot_threshold` | `config.training.miniprot_threshold` |
    | `config.ab_initio.flanking_region` | `config.training.flanking_region` |
    | `config.ab_initio.test_size` | `config.training.test_size` |
    | `config.ab_initio.busco_downloads` | `config.training.busco_downloads` |
    | `config.augustus.species_name` | `config.prediction.species` |
    | `config.augustus.mode` | `config.prediction.mode` |
    | `config.augustus.utr` | `config.prediction.utr` |
    | `config.augustus.config` | `config.prediction.hint_weights` |
    | `config.augustus.hint_config` | `config.prediction.hint_config` |
    | `config.sqanti.gtf_type` | `config.curation.mode` |
    | `config.sqanti.reference_gtf` | `config.curation.user_gtf` |
    | `config.sqanti.json_rules` | `config.curation.filter_rules` |
    | `config.augustus.reference_gtf` | `config.evaluation.reference_gtf` |
    | `config.qc.omark_db` | `config.evaluation.omark_db` |
    | `config.qc.omark_taxid` | `config.evaluation.omark_taxid` |
- **Note**: This mapping serves as the blueprint for the manual cleanup and final rule refactoring.

## [2026-05-19] - Initial Pre-trial State
- **Branch**: `dev-hint_testing`
- **Commit**: `111b5a5`
- **Summary**: Initial setup of the `dev-hint_testing` branch to prepare for Human Chr19 benchmarking.

## [2026-05-19] - Pre-Refactoring Save
- **Branch**: `dev-augustus_config`
- **Commit**: `54d3503`
- **Summary**: Saved the state of the `dev-augustus_config` branch before initiating the major refactoring of the hint generation logic.
- **Goal**: Ensure a recovery point exists before changing the core hint extraction and weighting mechanisms.

## [2026-05-08] - Enhanced Monoexon Filtering
- **Branch**: `dev-augustus_config`
- **Commit**: `dd21c7e`
- **Summary**: Added configurable options to filter monoexon-only genes or apply filtering across all genes to improve annotation specificity.

## [2026-05-08] - Merge Placebo Annotation
- **Branch**: `dev-augustus_config`
- **Commit**: `bb7eb52`
- **Summary**: Merged `dev-placebo_annotation` into the main configuration branch. Included removal of protection from GAQET plots and SQ3 reports to allow for easier troubleshooting.

## [2026-04-29] - Augustus Hint Configuration
- **Branch**: `dev-augustus_config`
- **Commit**: `47e268b`
- **Summary**: Initial changes to how Augustus hints are configured, moving towards a more flexible setup.

## [2026-04-29] - Local Rule Optimization
- **Branch**: `dev-augustus_config`
- **Commit**: `0ee01b9`
- **Summary**: Set `gaqet2` as a local rule to prevent unnecessary job overhead on clusters.

## [2026-04-29] - Input and IsoSeq Improvements
- **Branch**: `dev-augustus_config`
- **Commit**: `49ccc3a`
- **Summary**: Changed default input format to BAM and integrated PacBio conversion for IsoSeq workflows.

## [2026-04-22] - Placebo GTF Development
- **Branch**: `dev-augustus_config`
- **Commit**: `69de755`
- **Summary**: Initial development of the placebo GTF generation logic to support benchmarking.

## [2026-04-22] - GAQET2 Integration and Garnatxa Fixes
- **Branch**: `dev`
- **Commit**: `ceba269`, `db5041c`
- **Summary**: Merged GAQET2 development and applied specific fixes for the Garnatxa HPC cluster.

## [2026-04-22] - Hint Generation Refactoring
- **Branch**: `dev-gaqet2`
- **Commit**: `d5401f2`, `c5891be`
- **Summary**: Introduced a new Python-based hint generation script and initial statement for automated hint creation.

## [2026-04-15] - Resource Management
- **Branch**: `dev-gaqet2`
- **Commit**: `c343535`, `71d00b3`
- **Summary**: Modifications to resource usage and temporary file handling to improve pipeline robustness.

## [2026-04-15] - Core Fixes for Workflow Stability
- **Branch**: `dev-gaqet2`
- **Commit**: `db8dee1`
- **Summary**: Critical changes to ensure the pipeline runs correctly after recent updates.

## [2026-04-08] - SQANTI Update
- **Branch**: `dev-gaqet2`
- **Commit**: `6b756d9`
- **Summary**: Updated the pipeline to maintain compatibility with the new SQANTI release.
