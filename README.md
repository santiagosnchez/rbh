# rbh.pl
Finds reciprocal-best-hit single-copy matches between two al-by-all BLAST searches

It will find the reciprocal-best single-copy genes between two sets of BLAST outputs in format 7. If FASTA files are provided it will generate FASTA files with the sequence pairs found. It will also generate a list of gene names. A threshold for e-values can be provided using the `-filter` argument.

## Installation

    git clone https://github.com/santiagosnchez/rbh
    cd rbh
    chmod +x rbh.pl
    sudo cp rbh.pl /usr/local/bin

## Running the script

Run with `-h` for more details.

    perl rbh.pl -h

    USEAGE:
    
    perl rbh.pl -blast1 <blast_output_tab_1> \
                -blast2 <blast_output_tab_2> \
                -fasta1 <file_1.fasta> \
                -fasta2 <file_2.fasta> \
                -outdir <pairs> \
                -outfile <outfile> \ 
                -filter <e-value> \

Example:

    perl rbh.pl -bast1 sp1_BLAST.txt -bast2 sp2_BLAST.txt -fasta1 sp1_seqs.fasta -fasta2 sp2_seqs.fasta -outdir my_rbh -outfile myout_rbh.txt -filter 0

