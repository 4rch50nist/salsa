#!/bin/sh
set -e

ARTIC_URL="${ARTIC_URL:-https://zenodo.org/records/4486528/files}"
SCHEME_URL="https://raw.githubusercontent.com/artic-network/primer-schemes/e6ddb7c4a21a65e1e4ae3c21129f9d08c2cac12f/nCoV-2019/V3"
GENOME_READS="${GENOME_READS:-1800000}"
DORADO="${DORADO:-dorado}"
MINIMAP2="${MINIMAP2:-minimap2}"

FASTQ=data/fastq.tar.gz
SUMMARY=data/sequencing_summary_FAN53258_cabfcaea.txt

fetch() {
  [ -f "$2" ] && return 0
  curl -L --fail -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}

mkdir -p data runs

fetch "$ARTIC_URL/fastq.tar.gz?download=1" "$FASTQ"
fetch "$ARTIC_URL/sequencing_summary_FAN53258_cabfcaea.txt?download=1" "$SUMMARY"
fetch "$SCHEME_URL/nCoV-2019.primer.bed" data/nCoV-2019.primer.bed
fetch "$SCHEME_URL/nCoV-2019.reference.fasta" data/nCoV-2019.reference.fasta

gcc -O2 -o runs/hazard hazard.c -lz

if [ ! -f runs/genome.ops ]; then
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $c["passes_filtering"] == "TRUE" { print $c["start_time"] "\t" $c["read_id"] }' "$SUMMARY" \
    | sort -t "$(printf '\t')" -k1,1g | head -n "$GENOME_READS" | cut -f2 > runs/genome_ids.txt
  tar -xzOf "$FASTQ" | awk 'NR == FNR { keep[$1]; next } FNR % 4 == 1 { id = substr($1, 2); p = (id in keep) } p' runs/genome_ids.txt - > runs/genome_reads.fastq
  "$DORADO" demux --kit-name EXP-NBD196 --emit-fastq --emit-summary --output-dir runs/demux runs/genome_reads.fastq
  "$MINIMAP2" -x map-ont --secondary=no -t 4 data/nCoV-2019.reference.fasta runs/genome_reads.fastq > runs/genome.paf
  python3 genome_ops.py runs/genome_ids.txt runs/demux/barcoding_summary.txt runs/genome.paf data/nCoV-2019.primer.bed runs/genome.ops
  mv genome_meta.txt runs/
fi

cd runs
./hazard genome.ops
python3 ../make_figure.py
cp genome_fig.tex ..
