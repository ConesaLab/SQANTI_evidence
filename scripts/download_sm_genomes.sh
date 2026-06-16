#!/bin/bash
set -e

# Target data directory (relative to genomes_benchmark root)
DATA_DIR="data"

if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory '$DATA_DIR' not found. Please run this script from the genomes_benchmark directory."
    exit 1
fi

echo "Starting download of softmasked genomes from Ensembl..."

# 1. Mus musculus (Ensembl release 115)
echo "Downloading Mus musculus..."
mkdir -p "$DATA_DIR/mus_musculus"
wget -O "$DATA_DIR/mus_musculus/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/mus_musculus/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz"

# 2. Saccharomyces cerevisiae (Ensembl release 115)
echo "Downloading Saccharomyces cerevisiae..."
mkdir -p "$DATA_DIR/sac_cerevisiae"
wget -O "$DATA_DIR/sac_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/sac_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"

# 3. Arabidopsis thaliana (Ensembl Plants release 62)
echo "Downloading Arabidopsis thaliana..."
mkdir -p "$DATA_DIR/arabidopsis"
wget -O "$DATA_DIR/arabidopsis/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/ensembl/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/arabidopsis/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz"

# 4. Caenorhabditis elegans (Ensembl release 115)
echo "Downloading Caenorhabditis elegans..."
mkdir -p "$DATA_DIR/c_elegans"
wget -O "$DATA_DIR/c_elegans/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/c_elegans/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz"

# 5. Danio rerio (Ensembl release 115)
echo "Downloading Danio rerio..."
mkdir -p "$DATA_DIR/danio_rerio"
wget -O "$DATA_DIR/danio_rerio/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/danio_rerio/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz"

# 6. Drosophila melanogaster (Ensembl release 115)
echo "Downloading Drosophila melanogaster..."
mkdir -p "$DATA_DIR/drosophila"
wget -O "$DATA_DIR/drosophila/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/drosophila/Drosophila_melanogaster.BDGP6.54.dna_sm.toplevel.fa.gz"

# 7. Gallus gallus (Ensembl release 115)
echo "Downloading Gallus gallus..."
mkdir -p "$DATA_DIR/gallus_gallus"
wget -O "$DATA_DIR/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz" \
    "https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz"
gunzip -f "$DATA_DIR/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz"

echo "All downloads and extractions completed successfully!"
