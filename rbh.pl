#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;

my $VERSION = "0.1";

local $0 = basename $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{

  my $settings = {};
  GetOptions($settings, qw(help blast1=s blast2=s filter=f)) or die $!;

  usage() if($$settings{help});

  # Mandatory options
  for my $opt (qw(blast1 blast2 filter)){
    if(!defined $$settings{$opt}){
      print "ERROR $opt not set\n";
      usage();
    }
  }

  # Hashes of best hits
  my $hits1 = readTabBlast($$settings{blast1}, $settings);
  my $hits2 = readTabBlast($$settings{blast2}, $settings);

  while(my($query,$hit) = each(%$hits1)){
    logmsg "$query => $hit";
    if($$hits2{$hit} eq $query){
      logmsg " AND $hit => $query";
      print join("\t", $query, $hit)."\n";
    }
  }

  return 0;
}

sub usage{
  my $usage = "\nrbh.pl v$VERSION\nUSAGE:\n
  perl rbh.pl -blast1 <blast_output_tab_1> \
              -blast2 <blast_output_tab_2> \
              -filter <e-value>\n\n";
  print $usage;
  exit 0;
}

sub readTabBlast{
  my($infile, $settings) = @_;

  my %bestHits = ();

  # Sort the results by name, then bitscore, then evalue
  open(my $fh, "sort -k1,1r -k12,12nr -k11,11g '$infile' | ") or die "ERROR: could not read $infile: $!";

  while(<$fh>){
    chomp;
    my($query, $hit, $alnPerc, $alnLen, $mismatches, $gaps, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
    # skip comment lines which are usually headers
    next if ($_ =~ m/^#/);
    next if ($evalue > $$settings{filter});

    # The tabbed file is already sorted by best hit and so
    # if it was already defined in the hash then skip! 
    # Subsequential hits are worse than the first.
    if(!defined($bestHits{$query})){
      $bestHits{$query} = $hit;
    }
  }
  close $fh;

  return \%bestHits;
}

