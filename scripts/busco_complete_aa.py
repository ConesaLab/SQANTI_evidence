'''
Author: Alejandro Paniagua
Adapted by: Pablo Atienza
Scripts to obtain the protein sequences in fasta foramt from all the BUSCO 
genes that where found in single copy in the genome
'''

import os
from Bio import SeqIO
import argparse



def main():
    path = os.path.join(snakemake.input[0],f"run_{snakemake.params.lineage}",
                        "busco_sequences",f"{snakemake.params.gene_type}_copy_busco_sequences")
    outfile = snakemake.output[0]
    # The BUSCO output is a individual fasta for each gene that was found in the 
    # genome. All of this fasta are stored in the same path. 
    # sorted: os.listdir order is filesystem-dependent, and the FASTA order decides which
    # sequence CD-HIT picks as cluster representative -> keeps the training set reproducible
    archivos = sorted(os.listdir(path))
    l_records = []
    for archivo in archivos:
        # Protein sequences fasta end with faa while nucleotide sequences 
        # end with fasta. In this case the protein sequences are obtained
        if archivo.endswith("faa"):
            # The header of the fasta file is not informative, while the name
            # of the file is, because it is the orthoDB id of the sequence.
            # The header of the fasta is replace by the orthoDB id
            record = SeqIO.read(path + "/" + archivo, "fasta")
            record.id = archivo.split(".")[0]
            record.description = archivo.split(".")[0]
            l_records.append(record)
    SeqIO.write(l_records, outfile, "fasta")

if __name__=="__main__":
    main()