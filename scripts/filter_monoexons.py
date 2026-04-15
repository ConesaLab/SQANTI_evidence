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

    args = parseCmd()

    gtf_file = args.gtf
    hintfiles = args.hintfiles.split(',') if args.hintfiles else []
    out = args.out
    quiet = args.quiet

    if not quiet:
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
            # Transcrito multiexónico, lo conservamos automáticamente
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

def parseCmd():
    parser = argparse.ArgumentParser(description='Filtra genes monoexónicos que no están respaldados por hints.')
    parser.add_argument('-g', '--gtf', type=str, required=True,
        help='Archivo de predicción de genes en formato GTF (solo uno).')
    parser.add_argument('-e', '--hintfiles', type=str, required=True,
        help='Lista separada por comas de archivos con evidencias en GFF.')
    parser.add_argument('-o', '--out', type=str, required=True,
        help='Archivo de salida para la predicción filtrada en GTF.')
    parser.add_argument('-q', '--quiet', action='store_true',
        help='Modo silencioso (no imprime logs en stderr).')
    return parser.parse_args()

if __name__ == '__main__':
    main()
