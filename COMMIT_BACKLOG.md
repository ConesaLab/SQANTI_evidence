# Commit Backlog

This file tracks the history of commits made by the AI agent, providing a high-level summary of the changes and the reasoning behind them.

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
