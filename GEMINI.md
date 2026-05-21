# SQANTI-evidence: Long-Read Evidence-Driven Structural Annotation Pipeline

SQANTI-evidence is a **Snakemake**-based workflow designed for high-quality structural genome annotation using long-read RNA-seq data (PacBio or Oxford Nanopore). It implements an evidence-driven strategy that combines long-read transcriptomes (collapsed at the transcript level) with gene prediction tools like **Augustus** and curation tools like **SQANTI3**.

## Project Overview

- **Purpose**: Automates the integration of long-read transcriptomics into genome annotation workflows.
- **Key Technologies**: Snakemake, Python, Conda, Augustus, SQANTI3, BUSCO, Minimap2, GAQET2.
- **Architecture**: 
  - **HPC Path (Garnatxa)**: `/home/patienza/oscars/LR_annotation`
  - **Orchestration**: Managed by a central `snakefile` and a Python wrapper `sqanti_evidence`.
  - **Modular Rules**: Pipeline logic is split into functional modules under `rules/` (e.g., `ab_initio.smk`, `evidence_driven.smk`).
  - **Environment Management**: Conda environments for each step are defined in `envs/`.
  - **Auxiliary Scripts**: Custom Python and R scripts in `scripts/` handle data transformation and specialized filtering.

## Getting Started

### Prerequisites
- Conda or Mamba installed.
- Git.

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/pabloati/SQANTI_evidence.git
   cd SQANTI_evidence
   ```
2. Create the base Snakemake environment:
   ```bash
   conda env create -f env.yaml
   conda activate snakemake
   ```

## Usage

The pipeline is typically invoked via the `sqanti_evidence` wrapper script, which handles input validation and Snakemake execution.

### Basic Command
```bash
./sqanti_evidence --config your_config.yaml --jobs 8
```

### Options
- `--config`, `-c`: Path to the YAML configuration file (Required).
- `--dryrun`, `-n`: Perform a dry run to see what rules will be executed.
- `--slurm`, `-s`: Execute the pipeline on a Slurm-managed HPC cluster.
- `--jobs`, `-j`: Number of parallel jobs/cores.
- `--unlock`, `-u`: Unlock the working directory if a previous run was interrupted.

### Configuration
Use `config.yaml` as a template. Key sections include:
- `required`: Paths to genome FASTA, input RNA-seq (FASTQ/BAM), and output/tools directories.
- `augustus`: Parameters for gene prediction (species name, UTR, prediction mode).
- `sqanti`: Rules and reference GTF for transcript curation.
- `ab_initio`: Training parameters for Augustus (BUSCO lineage, flanking regions).
- `resources`: Resource allocation (CPU, memory, time) for different rule categories.

## Project Structure

- `rules/`: Snakemake rule definitions.
  - `setup/`: Pipeline initialization, logging, and tool installation logic.
- `scripts/`: Custom logic for GTF manipulation, filtering, and tool wrappers.
- `envs/`: YAML files defining Conda environments for various tools.
- `config.yaml`: Main configuration template.
- `snakefile`: Primary Snakemake workflow definition.
- `sqanti_evidence`: Python wrapper for simplified pipeline execution.

## Development Conventions

- **Snakemake First**: All workflow logic must be implemented as Snakemake rules or called via Snakemake.
- **Conda Isolation**: Every rule requiring external tools should specify a `conda:` environment file from `envs/`.
- **Modular Rules**: Keep rules focused; use the `include:` directive in the main `snakefile` to pull in modules.
- **Input Validation**: Use `scripts/input_check.py` and `rules/setup/functions.smk` to validate user configuration before execution.
- **Logging**: Use the centralized logging setup defined in `rules/setup/logging_setup.smk`.
- **Commit Backlog**: Maintain a `COMMIT_BACKLOG.md` file in the project root. ALWAYS update the commit backlog before performing a commit. After every commit, ensure the entry accurately reflects the branch, the commit hash, and the strategic goal of the work.
