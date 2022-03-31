# rbh.pl
Finds reciprocal-best-hit single-copy matches between two all-by-all BLAST searches

It will find the reciprocal-best single-copy genes between two sets of BLAST outputs in format 7. If FASTA files are provided it will generate FASTA files with the sequence pairs found. It will also generate a list of gene names. A threshold for e-values can be provided using the `-filter` argument.

## Installation

    mkdir ~/bin/build
    cd ~/bin/build
    git clone https://github.com/lskatz/rbh
    cd ~/bin
    ln -v ~/bin/build/rbh/rbh.pl ~/bin/

## Example

    for i in 1.fasta 2.fasta; do
      makeblastdb -dbtype nucl -in $i
    done;
    blastn -db 1.fasta -query 2.fasta -outfmt 6 -max_target_seqs 10000 > 1.tsv
    blastn -db 2.fasta -query 1.fasta -outfmt 6 -max_target_seqs 10000 > 2.tsv

    perl rbh.pl --blast1 1.tsv --blast2 2.tsv --filter 0.05 > map.tsv

## Usage

Run with `-h` for more details.

    perl rbh.pl -blast1 <blast_output_tab_1>
              -blast2 <blast_output_tab_2>
              -filter <e-value>

