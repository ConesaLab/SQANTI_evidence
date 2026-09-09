# SQANTI-evidence

**Long-Read Evidence-Driven Structural Annotation Pipeline**

SQANTI-evidence is a Snakemake workflow that produces a structural genome annotation from a genome assembly and long-read RNA-seq data (PacBio Iso-Seq / HiFi or Oxford Nanopore). It follows the findings of [Paniagua et al. 2025](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124), where the best-performing strategy was an evidence-driven gene prediction guided by a long-read transcriptome curated at the transcript level.

The pipeline combines two tracks:

1. **Curated empirical track (Tier 1)**: long reads are assembled into transcript models with **IsoQuant**, curated with **SQANTI3**, and the curated transcripts are written to the final annotation **unchanged**.
2. **Evidence-guided prediction track (Tier 2)**: **Augustus** is trained on the species' own dominant complete ORFs (optionally mixed with BUSCO orthologs) and guided by hints derived from the curated transcripts and from Miniprot self-protein alignments. Its predictions are added only where no Tier 1 transcript exists, to recover unexpressed genes.

---

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Configuration](#configuration)
4. [Input data](#input-data)
5. [Pipeline workflow](#pipeline-workflow)
6. [Output](#output)
7. [Testing](#testing)
8. [License & citations](#license--citations)

---

## Installation

Clone the repository and create the Snakemake environment from `env.yaml`. Every other tool (IsoQuant, SQANTI3, BUSCO, Augustus, AGAT, GAQET2, GffCompare, ...) is installed automatically by Snakemake into isolated conda environments on first use.

```bash
git clone https://github.com/pabloati/SQANTI_evidence.git
cd SQANTI_evidence
conda env create -f env.yaml      # creates the environment "sqanti_evidence"
conda activate sqanti_evidence
```

Tool environments and databases are stored under the directory given in `project.toolsdir`. Use the same `toolsdir` across runs so environments are built once and reused.

---

## Usage

Run the pipeline through the `sqanti_evidence` wrapper **from the repository root**. It validates the configuration and launches Snakemake.

```bash
./sqanti_evidence --config my_config.yaml --jobs 8            # local execution
./sqanti_evidence --config my_config.yaml --slurm --jobs 20   # Slurm cluster
./sqanti_evidence --config my_config.yaml --dryrun            # show the DAG only
```

<details> <summary>All options</summary>

```bash
usage: sqanti_evidence [-h] --config CONFIG [--dryrun] [--unlock] [--rerun_incomplete] [--slurm] [--jobs JOBS] [--extras EXTRAS] [--slurm_extras SLURM_EXTRAS]

options:
  -h, --help            show this help message and exit
  --config, -c CONFIG   Config file
  --dryrun, -n          Dry run
  --unlock, -u          Unlock
  --rerun_incomplete, -R
                        Rerun incomplete
  --slurm, -s           Run on slurm
  --jobs, -j JOBS       Number of jobs
  --extras, -e EXTRAS   Extra arguments passed to snakemake
  --slurm_extras, -se SLURM_EXTRAS
                        Extra arguments for the slurm executor
```

</details>

---

## Configuration

All inputs and parameters are given in a YAML file. Use `config.yaml` in this repository as a template. The file has six blocks.

### `project`

| Key | Description |
| :--- | :--- |
| `genome` | Reference genome FASTA to annotate. **Required.** A softmasked assembly is recommended. |
| `input` | Long reads in FASTQ, FASTA or unaligned BAM format. **Required.** See [Input data](#input-data). |
| `outdir` | Output directory (default `SQANTI_evidence_results`). |
| `toolsdir` | Directory where conda environments and databases are installed. **Required.** Use an absolute path. |
| `prediction_genome` | Genome used for the Augustus prediction step. Defaults to `genome`. |
| `log_level` | Verbosity of the pipeline log (`INFO` by default). |

### `training`

Controls how the Augustus species model is trained.

| Key | Description |
| :--- | :--- |
| `mode` | `mixed` (default): BUSCO single-copy orthologs plus dominant SQANTI3 ORFs. `sqanti_only`: SQANTI3 ORFs only (skips the BUSCO run). `busco_only`: BUSCO orthologs only. |
| `lineage` | BUSCO lineage dataset (e.g. `brassicales_odb12`). Required for `mixed` and `busco_only`. |
| `miniprot_threshold` | Minimum identity for a BUSCO/Miniprot gene model to be kept (default `0.95`). |
| `flanking_region` | Flanking bases added around each training gene, and minimum distance between BUSCO and SQANTI training genes (default `1000`). |
| `test_size` | Maximum number of genes in the training set (default `5000`). |

### `prediction`

| Key | Description |
| :--- | :--- |
| `species` | Name of the Augustus species model that will be created and trained. **Required.** |
| `mode` | `split` (default): one Augustus job per chromosome/scaffold. `full`: a single Augustus job. Results are identical; `split` is faster on a cluster. |
| `utr` | Boolean, reserved for UTR-aware prediction. |
| `filter_mode` | Filter for single-exon Augustus genes added as Tier 2: `strict` (must be supported by a hint), `medium` (default; at least 300 bp or hint-supported), `none`. |
| `hint_config` | TSV switching hint types on/off (default `envs/hint_config.tsv`). |
| `hint_weights` | Augustus extrinsic configuration with the bonus/malus per hint source (default `envs/extrinsic.hints_weights_default.cfg`). |
| `miniprot_args` | Options passed to Miniprot when aligning the dominant proteins back to the genome for `src=P` hints (default `-N 30 -p 0.6 --outs=0.5`). Miniprot writes secondary (paralog) alignments only above `--outs`; its own default of 0.99 suppresses them. |

### `curation`

| Key | Description |
| :--- | :--- |
| `mode` | Reference annotation given to SQANTI3. `placebo` (default): a dummy annotation, so that classification is reference-independent. `user`: the GTF in `user_gtf`. `ab_initio`: an ab initio Augustus prediction, only possible together with `training.mode: busco_only`. |
| `user_gtf` | Reference GTF when `mode` is `user`. |
| `data_type` | Sequencing technology passed to IsoQuant: `pacbio` (CCS/HiFi, default), `nanopore`, or `assembly`. |
| `isoquant_args` | Extra options appended verbatim to the IsoQuant command (default empty). Options the pipeline sets itself, such as `--reference`, `--data_type`, `--prefix`, `--threads`, `-o` and the input flag, are rejected. Typical use: Iso-Seq FLNC reads have their poly(A) tails removed by `isoseq refine`, and IsoQuant needs a tail to build novel single-exon transcripts and to strand unspliced reads; for these oriented reads set `"--polya_trimmed all --stranded forward"`. `--polya_trimmed all` assumes reads oriented 5'→3', so do not use it for raw ONT cDNA. |
| `filter_rules` | SQANTI3 rules-filter JSON (default `envs/filter_rules.json`). The default rules avoid structural categories so that filtering does not depend on the reference used. |

### `evaluation`

| Key | Description |
| :--- | :--- |
| `reference_gtf` | Optional gold-standard annotation. When set, GffCompare sensitivity/precision statistics are produced. |
| `omark_db`, `omark_taxid` | OMArk database name and NCBI taxid used by GAQET2. |

### `resources`

Four tiers, `small`, `medium`, `big` and `busco`, each with `cpus`, `mem`, `time` and `qos`. Memory and time use Slurm-style human-readable values (`mem: "8GB"`, `time: "2h"`). Any missing key falls back to the defaults in `rules/setup/functions.smk`.

```yaml
resources:
  small:  {cpus: 2,  mem: "8GB",  time: "2h",  qos: "short"}
  medium: {cpus: 8,  mem: "12GB", time: "10h", qos: "short"}
  big:    {cpus: 10, mem: "20GB", time: "24h", qos: "long"}
  busco:  {cpus: 30, mem: "60GB", time: "72h", qos: "long"}
```

---

## Input data

`project.input` is passed to IsoQuant, which aligns the reads itself with minimap2. The format is detected from the file extension:

| Extension | Handled as |
| :--- | :--- |
| `.fastq`, `.fq`, `.fasta`, `.fa` (optionally `.gz`) | Raw reads |
| `.bam` | **Unaligned** reads (e.g. PacBio FLNC or clustered BAM) |

Set `curation.data_type` to match the technology (`pacbio` or `nanopore`). Already-aligned BAM files are not supported yet.

---

## Pipeline workflow

1. **Transcript reconstruction.** IsoQuant aligns the reads to the genome and builds transcript models with full-length read counts.
2. **Curation.** SQANTI3 classifies the models against the reference chosen in `curation.mode` and the rules filter removes artefacts (intra-priming, RT-switching, non-canonical junctions, short or non-coding models). The surviving transcripts are the **Tier 1** set.
3. **Training set assembly.** One dominant complete ORF per gene is extracted from the curated transcripts. In `mixed`/`busco_only` mode BUSCO single-copy genes are added after removing those that overlap or lie within `flanking_region` of a SQANTI gene. Proteins are clustered with CD-HIT at 80% identity to remove redundancy, and the set is capped at `test_size` genes.
4. **Augustus training.** A new species model is created and trained with `etraining`; genes that fail training are removed and the model is retrained. Stop-codon frequencies are set from the training data.
5. **Hint generation.** Two hint sources are combined: `lrRNA` hints (introns, CDS, start and stop codons) from the curated transcripts, and `P` hints from Miniprot spliced alignments of the dominant proteins back to the genome, which point Augustus to unexpressed paralogs.
6. **Prediction.** Augustus runs with the trained model and the combined hints, per chromosome in `split` mode.
7. **Tier resolution.** Every Tier 1 transcript is written unchanged. An Augustus gene is added as **Tier 2** only if it does not overlap any Tier 1 gene on the same strand; single-exon Tier 2 genes are additionally filtered according to `prediction.filter_mode`.
8. **Standardisation and QC.** AGAT converts the result into a clean GFF3, GAQET2 evaluates structural quality against the reads, and, if `evaluation.reference_gtf` is set, GffCompare computes sensitivity and precision at transcript and CDS level.

---

## Output

Inside `project.outdir`:

```
outdir/
├── isoquant/<sample>/                  IsoQuant transcript models and counts
├── ab_initio/
│   ├── busco/                          BUSCO run (mixed / busco_only)
│   └── augustus/
│       ├── model/                      training set: sqanti_dominant.gff/.faa, busco_*.gff/.faa,
│       │                               training_genes.gff, cdhit.lst
│       └── training/                   etraining logs, bad.lst, filtered.gb, SC_freq.txt
├── evidence_driven/
│   ├── sqanti/                         <species>_classification.txt, <species>.filtered.gtf, <species>.filtered.regrouped.gtf (IsoQuant gene grouping restored), ...
│   ├── hints/                          <species>.rna.hints.gff, <species>.protein.hints.gff,
│   │                                   <species>.hints.gff (combined)
│   ├── augustus/                       Augustus_prediction.gff, resolved_prediction.gtf (Tier 1 + Tier 2)
│   └── Final_clean_prediction.gff      final annotation (AGAT-standardised GFF3)
├── quality_control/
│   ├── gaqet2/                         <sample>_GAQET.stats.tsv, <sample>_GAQET.plot.png
│   └── gffcompare/                     <sample>.stats, <sample>_cds.stats (if reference_gtf is set)
└── logs/                               one log per rule plus the pipeline log
```

The trained Augustus species model is written to `$AUGUSTUS_CONFIG_PATH/species/<prediction.species>`.

---

## Testing

Unit tests for the custom scripts and a DAG regression test are provided:

```bash
pytest tests/ -v                      # unit tests
bash tests/dryrun/run_dryruns.sh      # dry-run the workflow for every training/prediction mode
```

No example dataset ships with the repository yet.

---

## License & citations

License: to be defined.

If you use SQANTI-evidence, please cite this repository and the tools it relies on:

* IsoQuant: Prjibelski, A. D., et al. (2023). *Nature Biotechnology*, 41(7), 939–948.
* SQANTI3: Pardo-Palacios, F. J., et al. (2024). *Nature Methods*, 21(5), 789–798.
* Augustus: Stanke, M., & Morgenstern, B. (2005). *Nucleic Acids Research*, 33, W465–W467.
* BUSCO: Manni, M., et al. (2021). *Molecular Biology and Evolution*, 38(10), 4647–4654.
* Miniprot: Li, H. (2023). *Bioinformatics*, 39(1), btad014.
* CD-HIT: Fu, L., et al. (2012). *Bioinformatics*, 28(23), 3150–3152.
* AGAT: Jacques, C., et al. (2023). *Zenodo*, doi:10.5281/zenodo.3552717.
* GAQET2: https://github.com/compbio-gang/GAQET2
* GffCompare: Pertea, G., & Pertea, M. (2020). *F1000Research*, 9, 304.
* Snakemake: Mölder, F., et al. (2021). *F1000Research*, 10, 33.
