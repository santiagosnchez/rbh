#!/usr/bin/perl

my $usage = "\nUSEAGE:\nperl rbh.pl -blast1 <blast_output_tab_1> -blast2 <blast_output_tab_2>
-fasta1 <file_1.fasta> -fasta2 <file_2.fasta> -outdir <pairs> -outfile <outfile> -filter <e-value>\n";

my $blast1;
my $blast2;
my $fasta1;
my $fasta2;
my $outdir;
my $outfile;
my $tab1;
my $tab2;
my $filter;

if ($ARGV[0] eq "-h"){
	die "$usage\n";
} else {
	for (my $i=0; $i<scalar(@ARGV); ++$i){
		if ($ARGV[$i] eq "-blast1"){
			$blast1 = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-blast2"){
			$blast2 = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-fasta1"){
			$fasta1 = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-fasta2"){
			$fasta2 = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-outdir"){
			$outdir = $ARGV[$i+1];
			$outdir =~ s/\///;
		}
		if ($ARGV[$i] eq "-outfile"){
			$outfile = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-filter"){
			$filter = $ARGV[$i+1];
		}
	}
}

open(BLAST1, "<$blast1");
open(BLAST2, "<$blast2");

if (length($filter) != 0){
	print "BLAST hits above the evalue $filter will be ignored\n";
}

while(<BLAST1>){
	next if ($_ =~ m/^# [A-Z0].+/);
	if (length($filter) != 0){
		if (!/^#/){
			@line = split /\t/, $_;
			next if ($line[10] > $filter);
		}
	}
	$tab1 .= $_;
	print "reading $. lines for $blast1      \r";
}
print "\n";
close BLAST1;

my @hits1 = split "#", $tab1;

while(<BLAST2>){
	next if ($_ =~ m/^# [A-Z0].+/);
	if (length($filter) != 0){
		if (!/^#/){
			@line = split /\t/, $_;
			next if ($line[10] > $filter);
		}
	}
	$tab2 .= $_;
	print "reading $. lines for $blast2      \r";
}
print "\n";
close BLAST2;

my @hits2 = split "#", $tab2;

my @best1 = &getBest(@hits1);
my @best2 = &getBest(@hits2);

my %hashBest1 = &hashBest(@best1);
my %hashBest2 = &hashBest(@best2);

my @recip1=();
my @recip2=();

print "finding reciprocal-best single-copy genes...  ";

if (scalar(@best1) > scalar(@best2)){
	for (my $i=0; $i<scalar(@best2); ++$i){
		@line = split /\t/, $best2[$i];
		if (($hashBest1{$line[1]} eq $line[0]) and ($hashBest2{$line[0]} eq $line[1])){
			push @recip1, $line[1];
			push @recip2, $line[0];
		}
	}
} else {
	for (my $i=0; $i<scalar(@best1); ++$i){
                @line = split /\t/, $best1[$i];
                if (($hashBest2{$line[1]} eq $line[0]) and ($hashBest1{$line[0]} eq $line[1])){
                        push @recip1, $line[1];
                        push @recip2, $line[0];
                }
        }
}

print "done\n";

if (scalar(@recip1) == scalar(@recip2)){
	open(OUTFILE1, ">$outfile" . "_list1");
	open(OUTFILE2, ">$outfile" . "_list2");
	foreach(@recip1){
		#print "$_\n";
		@line = split "\t", $_;
		print OUTFILE1 $line[0] . "\n";
	}
	close OUTFILE1;
	foreach(@recip2){
		#print "$_\n";
		@line = split "\t", $_;
		print OUTFILE2 $line[0] . "\n";
	}
	close OUTFILE2;
} else {
	die "...  Error: number of reciprocal best hits not equal\n";
}

open(FASTA1, "<$fasta1");
open(FASTA2, "<$fasta2");

my $paste_fasta1;
my $paste_fasta2;

while(<FASTA1>){
	$paste_fasta1 .= $_;
}
close FASTA1;

while(<FASTA2>){
	$paste_fasta2 .= $_;
}
close FASTA2;

my @fasta_rec1 = split ">", $paste_fasta1;
my @fasta_rec2 = split ">", $paste_fasta2;

shift @fasta_rec1;
shift @fasta_rec2;

my %hash_fasta1 = &hashFasta(@fasta_rec1);
my %hash_fasta2 = &hashFasta(@fasta_rec2);

for (my $i=0; $i<scalar(@recip1); ++$i){
	open(OUT, ">" . $outdir . "\/" . "$outfile" . "_$i");
	print OUT ">" . $recip1[$i] . "__1" . "\n" . $hash_fasta1{$recip1[$i]} . ">" . $recip2[$i] . "__2" . "\n" . $hash_fasta2{$recip2[$i]};
	close OUT;
}

print "#############\n";
print scalar(@recip1) . " reciprocal-best single-copy genes found\n";
print "fasta files with sequences can be found in the following directory: $outdir\n";
print "a list of gene labels can be found in $outfile list 1 and 2\n";
print "#############\n";

# Subroutines

sub hashFasta {
	my %hash=();
	for (my $i=0; $i<scalar(@_); ++$i){
		@line = split /\n/, $_[$i];
		#print ">" . $line[0] . "\n";
		my $seq;
		for (my $j=1; $j<scalar(@line); ++$j){
			chomp($line[$j]);
			$seq .= $line[$j];
		}
		$hash{$line[0]} = "$seq\n";
	}
	return %hash;
}

sub hashBest {
	my %hash=();
	for (my $i=0; $i<scalar(@_); ++$i){
		@line = split /\t/, $_[$i];
		$hash{$line[0]} = $line[1];
	}
	return %hash;
}

sub getBest {
	shift @_;
	foreach (@_){
		$_ =~ s/^ //;
	}
	my @list=();
	for (my $i=0; $i<scalar(@_); ++$i){
		@bhit = split "\n", $_[$i];
		shift @bhit;
		my $match=0;
		my @query=();
		my @subject=();
		for (my $j=0; $j<scalar(@bhit); ++$j){
			@line = split "\t", $bhit[$j];
			push @query, $line[0];
			push @subject, $line[1];
		}
		%seenQ=();
		@sameQ = grep { $seenQ{$_} ++ } @query;
		%seenS=();
		@sameS = grep { $seenS{$_} ++ } @subject;
		if (scalar(@sameQ) == scalar(@sameS)){
			push @list, @bhit[$#bhit] . "\n";
		}
	}
	return @list;
}

