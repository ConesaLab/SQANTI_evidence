import argparse

def process_transcript(t_id, chrom, strand, exons, cds, out_file, cds_dict):
    """Processes a single transcript and writes its hints immediately to save memory."""
    if not exons:
        return
        
    hint_attrs = f"grp={t_id};pri=1;src=lrRNA"
    
    # Sort coordinates just in case they appear out of order in the GFF
    exons.sort()
    
    # 1. Generate Exon hints
    for i, (ex_start, ex_end) in enumerate(exons):
        hint_type = "exonpart" if i == 0 or i == len(exons) - 1 else "exon"
        out_file.write(f"{chrom}\tHints\t{hint_type}\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t{hint_attrs}\n")
        
    # 2. Generate Intron hints
    for i in range(len(exons) - 1):
        intron_start = exons[i][1] + 1
        intron_end = exons[i+1][0] - 1
        if intron_start <= intron_end:
            out_file.write(f"{chrom}\tHints\tintron\t{intron_start}\t{intron_end}\t.\t{strand}\t.\t{hint_attrs}\n")
            
    # 3. Generate Start and Stop hints
    if cds:
        cds.sort()
        min_cds = cds[0][0]
        max_cds = cds[-1][1]
        
        # Fast dictionary lookup (defaults to 'internal' if missing)
        cds_type = cds_dict.get(t_id, 'internal')
        has_start = cds_type in ['complete', '3prime_partial']
        has_stop = cds_type in ['complete', '5prime_partial']
        
        if strand == '+':
            if has_start:
                out_file.write(f"{chrom}\tHints\tstart\t{min_cds}\t{min_cds+2}\t.\t+\t0\t{hint_attrs}\n")
            if has_stop:
                out_file.write(f"{chrom}\tHints\tstop\t{max_cds-2}\t{max_cds}\t.\t+\t0\t{hint_attrs}\n")
        elif strand == '-':
            if has_start:
                out_file.write(f"{chrom}\tHints\tstart\t{max_cds-2}\t{max_cds}\t.\t-\t0\t{hint_attrs}\n")
            if has_stop:
                out_file.write(f"{chrom}\tHints\tstop\t{min_cds}\t{min_cds+2}\t.\t-\t0\t{hint_attrs}\n")

def main():
    gff_file = snakemake.input.gtf
    table_file = snakemake.input.classification
    out_file_path = snakemake.output[0]
    sep = '\t'

    print("1. Loading CDS metadata into memory...")
    cds_dict = {}
    
    # Read table sequentially without loading full columns into memory
    with open(table_file, 'r') as f:
        header = f.readline().strip('\n').split(sep)
        
        # Hardcoded check assuming inputs are always the same
        iso_idx = header.index('isoform')
        cds_idx = header.index('CDS_type')
        
        for line in f:
            parts = line.strip('\n').split(sep)
            if len(parts) > max(iso_idx, cds_idx):
                cds_dict[parts[iso_idx]] = parts[cds_idx].strip().lower()

    print("2. Processing GFF and generating hints (Streaming mode)...")
    current_t_id = None
    chrom = strand = None
    exons = []
    cds = []

    # Stream the GFF line by line and immediately write outputs
    with open(gff_file, 'r') as f, open(out_file_path, 'w') as out:
        for line in f:
            if line[0] == '#' or not line.strip():
                continue
                
            parts = line.strip('\n').split('\t')
            if len(parts) < 9:
                continue
                
            # Ultra-fast string splitting (replaces slow regex)
            try:
                t_id = parts[8].split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue
            
            # When we detect a new transcript_id, process the PREVIOUS transcript
            if t_id != current_t_id:
                if current_t_id is not None:
                    process_transcript(current_t_id, chrom, strand, exons, cds, out, cds_dict)
                
                # Reset variables for the new transcript
                current_t_id = t_id
                chrom = parts[0]
                strand = parts[6]
                exons.clear()
                cds.clear()
                
            # Collect coordinates for the CURRENT transcript
            feature = parts[2]
            if feature == 'exon':
                exons.append((int(parts[3]), int(parts[4])))
            elif feature == 'CDS':
                cds.append((int(parts[3]), int(parts[4])))

        # Process the very last transcript in the file once the loop finishes
        if current_t_id is not None:
            process_transcript(current_t_id, chrom, strand, exons, cds, out, cds_dict)

    print(f"Done! Hints successfully saved to {out_file_path}")

if __name__ == '__main__':
    main()
