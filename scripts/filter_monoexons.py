#!/usr/bin/env python3
# ==============================================================
# Filter for monoexon genes based on TSEBRA's logic
# ==============================================================
import argparse
import sys
import os

def main():
    # Import local TSEBRA classes
    from genome_anno import Anno
    from evidence import Evidence

    gtf_file = snakemake.input.gff
    hint_input = snakemake.input.hints
    hintfiles = hint_input if isinstance(hint_input, list) else [hint_input]
    out = snakemake.output.gff
    
    # Extract filter mode from configuration
    filter_mode = 'monoexon'
    try:
        if 'augustus' in snakemake.config and 'filter_mode' in snakemake.config['augustus']:
            filter_mode = snakemake.config['augustus']['filter_mode']
        elif hasattr(snakemake.config, 'augustus') and hasattr(snakemake.config.augustus, 'filter_mode'):
            filter_mode = snakemake.config.augustus.filter_mode
    except Exception:
        pass

    quiet = False

    if not quiet:
        sys.stderr.write(f'### FILTER MODE: {filter_mode.upper()}\n')
        sys.stderr.write(f'### READING GENE PREDICTION: [{gtf_file}]\n')
    
    # Load the GTF
    anno = Anno(gtf_file, 'anno1')
    anno.addGtf()
    anno.norm_tx_format()

    # Load extrinsic evidence (hints)
    evi = Evidence()
    for h in hintfiles:
        if not quiet:
            sys.stderr.write(f'### READING EXTRINSIC EVIDENCE: [{h}]\n')
        evi.add_hintfile(h)

    if not quiet:
        sys.stderr.write('### FILTERING UNSUPPORTED TRANSCRIPTS\n')

    keep_txs = {}
    filtered_out_count = 0

    # Iterate over all transcripts
    for tx_id, tx in anno.transcripts.items():
        # Check if the transcript has introns
        has_introns = 'intron' in tx.transcript_lines and len(tx.transcript_lines['intron']) > 0
        
        if not has_introns:
            # Monoexonic transcript, check for hint support
            supported = False
            
            # Search for support in start codons
            for sc_line in tx.transcript_lines.get('start_codon', []):
                # sc_line[3] is start, sc_line[4] is end
                if evi.get_hint(tx.chr, sc_line[3], sc_line[4], 'start', tx.strand):
                    supported = True
                    break
                    
            # If not supported by start, check stop codons
            if not supported:
                for sc_line in tx.transcript_lines.get('stop_codon', []):
                    if evi.get_hint(tx.chr, sc_line[3], sc_line[4], 'stop', tx.strand):
                        supported = True
                        break
                        
            if supported:
                keep_txs[tx_id] = tx
            else:
                filtered_out_count += 1
        else:
            if filter_mode == 'all':
                # Multi-exon transcript, check if at least one intron is supported
                supported = False
                for intron_line in tx.transcript_lines.get('intron', []):
                    if evi.get_hint(tx.chr, intron_line[3], intron_line[4], 'intron', tx.strand):
                        supported = True
                        break
                
                if supported:
                    keep_txs[tx_id] = tx
                else:
                    filtered_out_count += 1
            else:
                # Keep all multi-exon transcripts based on "monoexon" mode
                keep_txs[tx_id] = tx

    if not quiet:
        sys.stderr.write(f'### FILTERED OUT {filtered_out_count} TRANSCRIPTS\n')
        sys.stderr.write('### WRITING RESULTS\n')

    # Create a new annotation with filtered transcripts
    filtered_anno = Anno('', 'filtered_annotation')
    filtered_anno.transcripts = keep_txs  # Assign directly to avoid a bug in add_transcripts
    filtered_anno.find_genes()
    filtered_anno.write_anno(out)

    if not quiet:
        sys.stderr.write('### FINISHED\n\n')
        sys.stderr.write(f'### Filtered prediction is available at {out}.\n')

if __name__ == '__main__':
    main()
