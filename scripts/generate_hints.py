import argparse
import os

def process_transcript(t_id, chrom, strand, exons, cds, out_file, cds_dict, config):
    """Processes a single transcript and writes its hints based on config."""
    if not exons and not cds:
        return
        
    hint_attrs = f"grp={t_id};pri=1;src=lrRNA"
    
    # Sort coordinates
    exons.sort()
    cds.sort()

    # Get CDS type from classification
    cds_type = cds_dict.get(t_id, 'internal')
    has_start = cds_type in ['complete', '3prime_partial']
    has_stop = cds_type in ['complete', '5prime_partial']

    # 1. Exon hints
    if config.get('exon', False) or config.get('exonpart', False):
        for i, (ex_start, ex_end) in enumerate(exons):
            if config.get('exonpart', False) and (i == 0 or i == len(exons) - 1):
                out_file.write(f"{chrom}\tHints\texonpart\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t{hint_attrs}\n")
            elif config.get('exon', False):
                out_file.write(f"{chrom}\tHints\texon\t{ex_start}\t{ex_end}\t.\t{strand}\t.\t{hint_attrs}\n")
        
    # 2. Intron hints
    if config.get('intron', False):
        for i in range(len(exons) - 1):
            intron_start = exons[i][1] + 1
            intron_end = exons[i+1][0] - 1
            if intron_start <= intron_end:
                out_file.write(f"{chrom}\tHints\tintron\t{intron_start}\t{intron_end}\t.\t{strand}\t.\t{hint_attrs}\n")
            
    # 3. CDS hints
    if cds and config.get('CDS', False):
        for i, (c_start, c_end) in enumerate(cds):
            h_type = "CDS"
            # Logic for CDSpart if missing start/stop
            if strand == '+':
                # Start is at the beginning (min_cds), Stop is at the end (max_cds)
                if (i == 0 and not has_start) or (i == len(cds) - 1 and not has_stop):
                    h_type = "CDSpart"
            else: # strand == '-'
                # Start is at the end (max_cds), Stop is at the beginning (min_cds)
                if (i == len(cds) - 1 and not has_start) or (i == 0 and not has_stop):
                    h_type = "CDSpart"
            
            out_file.write(f"{chrom}\tHints\t{h_type}\t{c_start}\t{c_end}\t.\t{strand}\t.\t{hint_attrs}\n")

    # 4. Start and Stop hints
    if cds:
        min_cds = cds[0][0]
        max_cds = cds[-1][1]
        
        if strand == '+':
            if has_start and config.get('start', False):
                out_file.write(f"{chrom}\tHints\tstart\t{min_cds}\t{min_cds+2}\t.\t+\t0\t{hint_attrs}\n")
            if has_stop and config.get('stop', False):
                out_file.write(f"{chrom}\tHints\tstop\t{max_cds-2}\t{max_cds}\t.\t+\t0\t{hint_attrs}\n")
        elif strand == '-':
            if has_start and config.get('start', False):
                out_file.write(f"{chrom}\tHints\tstart\t{max_cds-2}\t{max_cds}\t.\t-\t0\t{hint_attrs}\n")
            if has_stop and config.get('stop', False):
                out_file.write(f"{chrom}\tHints\tstop\t{min_cds}\t{min_cds+2}\t.\t-\t0\t{hint_attrs}\n")

def main():
    gff_file = snakemake.input.gtf
    table_file = snakemake.input.classification
    config_file = snakemake.input.hint_config
    out_file_path = snakemake.output[0]
    sep = '\t'

    # Load hint configuration
    hint_config = {}
    with open(config_file, 'r') as f:
        next(f) # skip header
        for line in f:
            feat, enabled = line.strip().split('\t')
            hint_config[feat] = enabled.lower() == 'true'

    print("1. Loading CDS metadata into memory...")
    cds_dict = {}
    with open(table_file, 'r') as f:
        header = f.readline().strip('\n').split(sep)
        iso_idx = header.index('isoform')
        cds_idx = header.index('CDS_type')
        for line in f:
            parts = line.strip('\n').split(sep)
            if len(parts) > max(iso_idx, cds_idx):
                cds_dict[parts[iso_idx]] = parts[cds_idx].strip().lower()

    print("2. Processing GFF and generating hints...")
    current_t_id = None
    chrom = strand = None
    exons = []
    cds = []

    with open(gff_file, 'r') as f, open(out_file_path, 'w') as out:
        for line in f:
            if line[0] == '#' or not line.strip():
                continue
            parts = line.strip('\n').split('\t')
            if len(parts) < 9:
                continue
            try:
                t_id = parts[8].split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue
            
            if t_id != current_t_id:
                if current_t_id is not None:
                    process_transcript(current_t_id, chrom, strand, exons, cds, out, cds_dict, hint_config)
                current_t_id = t_id
                chrom = parts[0]
                strand = parts[6]
                exons.clear()
                cds.clear()
                
            feature = parts[2]
            if feature == 'exon':
                exons.append((int(parts[3]), int(parts[4])))
            elif feature == 'CDS':
                cds.append((int(parts[3]), int(parts[4])))

        if current_t_id is not None:
            process_transcript(current_t_id, chrom, strand, exons, cds, out, cds_dict, hint_config)

if __name__ == '__main__':
    main()
