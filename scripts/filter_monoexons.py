#!/usr/bin/env python3
# ==============================================================
# Filter for monoexon genes based on TSEBRA's logic
# ==============================================================
import argparse
import sys
import os

def main():
    # Importar las clases locales de TSEBRA
    from genome_anno import Anno
    from evidence import Evidence

    gtf_file = snakemake.input.gff
    hint_input = snakemake.input.hints
    hintfiles = hint_input if isinstance(hint_input, list) else [hint_input]
    out = snakemake.output.gff
    
    # Extraer parámetro de modo de filtrado de la configuración
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
        sys.stderr.write(f'### MODO DE FILTRADO: {filter_mode.upper()}\n')
        sys.stderr.write(f'### LEYENDO PREDICCIÓN DE GENES: [{gtf_file}]\n')
    
    # Cargar el GTF
    anno = Anno(gtf_file, 'anno1')
    anno.addGtf()
    anno.norm_tx_format()

    # Cargar las evidencias (hints)
    evi = Evidence()
    for h in hintfiles:
        if not quiet:
            sys.stderr.write(f'### LEYENDO EVIDENCIA EXTRÍNSECA: [{h}]\n')
        evi.add_hintfile(h)

    if not quiet:
        sys.stderr.write('### FILTRANDO TRANSCRITOS MONOEXÓNICOS SIN SOPORTE\n')

    keep_txs = {}
    filtered_out_count = 0

    # Iterar sobre todos los transcritos
    for tx_id, tx in anno.transcripts.items():
        # Comprobar si el transcrito tiene intrones
        has_introns = 'intron' in tx.transcript_lines and len(tx.transcript_lines['intron']) > 0
        
        if not has_introns:
            # Es un transcrito monoexónico, comprobamos si tiene soporte de hints
            supported = False
            
            # Buscar soporte en los codones de inicio (start)
            for sc_line in tx.transcript_lines.get('start_codon', []):
                # sc_line[3] es start, sc_line[4] es end
                if evi.get_hint(tx.chr, sc_line[3], sc_line[4], 'start', tx.strand):
                    supported = True
                    break
                    
            # Si no está soportado por start, buscar en los codones de parada (stop)
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
                # Transcrito multiexónico, comprobamos si tiene al menos un intrón soportado
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
                # Conservamos todos los multiexónicos según el modo "monoexon"
                keep_txs[tx_id] = tx

    if not quiet:
        sys.stderr.write(f'### SE FILTRARON {filtered_out_count} TRANSCRITOS MONOEXÓNICOS\n')
        sys.stderr.write('### ESCRIBIENDO RESULTADOS\n')

    # Crear una nueva anotación con los transcritos filtrados
    filtered_anno = Anno('', 'filtered_annotation')
    filtered_anno.transcripts = keep_txs  # Asignamos directamente para evitar un bug en add_transcripts
    filtered_anno.find_genes()
    filtered_anno.write_anno(out)

    if not quiet:
        sys.stderr.write('### FINALIZADO\n\n')
        sys.stderr.write(f'### La predicción filtrada se encuentra en {out}.\n')

if __name__ == '__main__':
    main()
