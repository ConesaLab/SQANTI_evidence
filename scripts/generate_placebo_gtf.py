import sys
from Bio import SeqIO

def generate_sqanti_dummy_gtf(fasta_path, output_path):
    """
    Streams a FASTA file using Biopython and generates a dummy GTF file.
    Uses hardcoded coordinates since SQANTI only needs the chromosome IDs.
    """
    try:
        with open(output_path, 'w') as out:
            out.write("##gff-version 3\n")
            
            # Stream the fasta file to save memory and time
            for record in SeqIO.parse(fasta_path, "fasta"):
                seqid = record.id
                
                # Hardcoded dummy coordinates
                start, end = 1, 1000
                exon1_start, exon1_end = 1, 100
                exon2_start, exon2_end = 900, 1000
                source = "sqanti_faker"
                strand = "+"
                
                # Write Gene
                out.write(f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"gene_{seqid}\"\n")
                
                # Write Transcript (mRNA)
                out.write(f"{seqid}\t{source}\ttranscript\t{start}\t{end}\t.\t{strand}\t.\ttranscript_id \"transcript_{seqid}\";gene_id \"gene_{seqid}\"\n")
                
                # Write Exons
                out.write(f"{seqid}\t{source}\texon\t{exon1_start}\t{exon1_end}\t.\t{strand}\t.\ttranscript_id \"transcript_{seqid}\";gene_id \"gene_{seqid}\"\n")
                out.write(f"{seqid}\t{source}\texon\t{exon2_start}\t{exon2_end}\t.\t{strand}\t.\ttranscript_id \"transcript_{seqid}\";gene_id \"gene_{seqid}\"\n")
    
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)
        
# Snakemake integration using unnamed input/output
fasta_in = snakemake.input[0]
gff_out = snakemake.output[0]

generate_sqanti_dummy_gtf(fasta_in, gff_out)