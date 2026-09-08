import argparse
import re

def main():
    first_header = True
    info_line = False
    gene_name = None  # set on each 'gene' line; features before the first gene line are left untouched
    with open(snakemake.output[0], "w") as outf:
        # Open the Augustus output file
        with open(snakemake.input[0]) as f:
            # Iterate throught the lines
            for line in f:
                if line.startswith("#"):
                    if line.startswith("# start gene"):
                        info_line = True
                    elif line.startswith("# end gene"):
                        info_line = False
                    if first_header or info_line:
                        outf.write(line)
                        if line.startswith("# admissible"):
                            first_header = False
                else:
                    columns = line.split("\t")
                    chr = columns[0]
                    # Add the chromosome to the gene id and trasncript id 
                    if columns[2] == "gene":
                        gene_name = columns[8].strip('\n')
                        columns[8] = f"{chr}_{columns[8]}"
                    elif columns[2] == "transcript":
                        columns[8] = f"{chr}_{columns[8]}"
                    elif gene_name is not None:
                        columns[8] = columns[8].replace(f'\"{gene_name}',f'\"{chr}_{gene_name}')
                        
                    outf.write("\t".join(columns))



if __name__ == "__main__":
    main()  