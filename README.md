# SQANTI-evidence

**Long-Read Evidence-Driven Structural Annotation Pipeline**  
SQANTI-evidence is a workflow that integrates all the steps needed to generate a structural annotation of a genome using long-read RNA-seq data. The strategy chosen is based on the findings of Paniagua et al. 2024, where the best combiation was to use an evidence-driven gene prediction strategy with the long-read transcriptome collapsed at transcript level. It mainly relies on Augustus for gene prediciton and SQANTI3 for the long-read transcriptome curation.  The pipeline has been developed using Snakemake. 

---

## Table of Contents

1. [Overview](#overview)  
2. [Features](#features)  
3. [Requirements](#requirements)  
4. [Installation](#installation)  
5. [Usage](#usage)  
6. [Configuration](#configuration)  
7. [Pipeline Workflow](#pipeline-workflow)  
8. [Scripts & Rules](#scripts--rules)  
9. [Output](#output)  
10. [Examples](#examples)  
11. [License & Citations](#license--citations)  
12. [Contact / Support](#contact--support)

---

## Overview

This repository implements a **Snakemake** pipeline (with auxiliary scripts) to generate structural genome annotations guided by long-read sequencing data (e.g. PacBio, Oxford Nanopore)
It aims to produce high-quality annotations by combining transcript evidence from long reads with conventional annotation strategies. The main structure of the pipeline and use of the long-read transcriptomics is derived from [this paper](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124).

---


## Installation

In order to install SQANTI-evidence, you only need to clone the GitHub repository and create a conda environment with snakemake. In order to do so, you can use the `env.yaml` file you will find in the repository.

```bash
   git clone https://github.com/pabloati/SQANTI_evidence.git
   cd SQANTI_evidence
   conda create -n snakemake
```
The rest of the needed packages and dependencies installation will be handled automatically by the pipeline. You only need to provide the path to the directory where you wish to install the tools. It is recommended to use the same directory for the installation of tools for different runs, since SQANTI-evidence will autodetect them and only do the installation process once. If this directory is changed, the conda environments and the non-conda packages will be installed again. You can control this by using the `toolsdir` option in the configuration file. 

---

## Usage

The main script to run the pipeline is `sqanti_evidence`. In order to run the script you need to include a valid configuration file. To do so, take as template the `config.yaml` file you will find in the directory. Extra options can be added to `sqanti_evidence` to control snakemake's behaviour.

<details> <summary>See the options </summary>

```bash
usage: sqanti_evidence [-h] --config CONFIG [--dryrun] [--unlock] [--rerun_incomplete] [--slurm] [--jobs JOBS] [--extras EXTRAS] [--slurm_extras SLURM_EXTRAS]

Run snakemake

options:
  -h, --help            show this help message and exit
  --config, -c CONFIG   Config file
  --dryrun, -n          Dry run
  --unlock, -u          Unlock
  --rerun_incomplete, -R
                        Rerun incomplete
  --slurm, -s           Run on slurm
  --jobs, -j JOBS       Number of jobs
  --extras, -e EXTRAS   Extra arguments
  --slurm_extras, -se SLURM_EXTRAS
                        Extra arguments for slurm
```

</details>

---

### Configuration file

In SQANTI evidence, all of the inputs are parameters are defined via the `config.yaml`file. As the name implies, the file follows a YAML format, which can be broken down into the following sections:

1. Required paramters: Fundamental to run SQANTI-evidence
   - genome: Path to the reference genome to be annotated in fasta format
   - input: Path to the long reads, in fastq or fasta format. For more information see the [preprocessing section](#Preprocessing)
   - outdir: Desired output directory for the run
   - toolsdir: directory where the software will be installed and the different databases downloaded. It is recommended to use a full path rather than relative
   - log_level: The amount of log to display by SQANTI-evidence
2. SQANTI parameters: Needed to run SQANTI
   - json_rules: Path to the json file with the rules to filter transcripts. A template is provided with the best rules for a human transcriptome.

2. Augustus Parameters: Parameters to control the gene prediction variables
   - prediction: To choose between `ab_initio` or `evidence_driven`. (Recheck this thing how it works)
   - reference_gtf: In the event of selecting an evidence_driven prediciton, the user must specify the reference annotation for SQANTI3 to work. 
   - species_name: Desired name for the species model to be stored in augustus.
   - utr: Boolean value to control if Untrasnlated Regions will be included in the annotation
   - mode: `split` or `full` (perhaps change to a boolean called split). Determines if the genome will be chunked for Augustus to run faster. This won't impact the final result, only speed up the analysis if you have more cores available
3. QC: Parameters needed to control the QC
   - omark_db: Name of the database used by OMARK
   - omark_taxid: NCBI taxid of the species, for OMARK to run its analysis
4. ab_intio parameters: These parameters control the training of the augustus gene model
   - lineage: Lineage used in BUSCO to find the set of core genes. These will also be used in the quality control step
   - miniprot_threshold: Minimum score to consider a miniprot hit of the agustus set of genes as valid
   - flanking_region: Number of extra nucleotides added to the ends of each gene for augustus training
   - test_size: Maximum number of genes used for agustus training
5. Resources: Parameters used to control the resources used by each step of the pipeline. They are divided into 3 categories:
   - big
   - small
   - medium
Busco is a exception. These parameters can be changed to suit your needs, and are of high importance when running the pipeline in a HPC cluster. 


### Input preprocessing

The long-reads transcriptomics data has to be given to SQANTI-evidence at a certain point of pre-processing. This will depend on the sequecing technology used:

1. PacBio: IsoSeq data should follow the [IsoSeq3](https://isoseq.how/) until before mapping (running cluster2). If reads are at any other point of the preprocessing pipeline, the results are not guaranteed to be optimal or fail at any point.
2. NanoPore: Only alskkdjfñaslk has to be done
---

## Pipeline Workflow

The initial step is to find the core gene set of the organism using busco. It uses minimmap mode, as it is faster and than using the augustus mode and produces the same results when the non-exact matches are discarded. This gene set is then expaned to include the flanking region determined in the configuration. After agusutus is trained, in the case of not having a reference annotation, the SQANTI-evidence will produce an ab-initio annotation with the trained agusutus model for SQANTI3 to run. 

Then, the long reads are mapped to the genome and collapsed into transcripts using the ISoSeq3 approach. These transcripts are then run through SQANTI3, and cleaned using the filter module, with the rules strategy. The rules are tailored so the structural categories are not used, since these would bias the final result. Most of the classifciation colyumns that are taken into account for filter are independen from the reference used. 

The final set of transcripts is then used as evidence to guide the final gene prediciton by augustus. The quality control is perofm using [GAQET2](https://github.com/compbio-gang/GAQET2).

---

## Output

The output will be saved in the directory stated in `outdir`. The structure of the output files will follow the following tree structure:


---

## Examples

*(You may want to include a small example or test dataset to demonstrate pipeline execution. If you have one, mention it here. E.g.):*

* `example/` — folder with toy genome + reads, config, and expected outputs

* Usage:

  ```bash
  cd example
  snakemake -j 4
  ```

* Compare output GFF3 with expected reference.

If you don’t have an example, you could add one in future to help users.

---

## License & Citations

State your license (e.g. MIT, GPL, etc.) here.
Also include citations to relevant tools or papers used in this pipeline.

```text
MIT License
(c) 2025 Genomics of Gene Expression Lab(or your name)
```

Please cite this repository as:

