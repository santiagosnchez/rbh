#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/dirname/;
use Data::Dumper;

use Test::More tests => 1;

my $thisDir = dirname($0);

# Prioritize this rbh in the path
$ENV{PATH}="$thisDir/..:$ENV{PATH}";

subtest 'best blast hits' => sub{
  my $logfile = "$thisDir/map.log";
  my $mapfile = "$thisDir/map.tsv";
  system("rbh.pl -blast1 blast1.tsv -blast2 blast2.tsv -filter 0.05 > $mapfile 2> $logfile");
  if($?){
    my $log = `cat $mapfile`;
    BAIL_OUT("rbh.pl resulted in an error. Log was:\n$log");
  }

  my %expected = (
    SALM_18928 => "STMMW_45271",
    SALM_18910 => "STMMW_45081",
    SALM_15703 => "STMMW_11291",
  );

  my %obs;
  open(my $fh, $mapfile) or BAIL_OUT("ERROR could not read $mapfile: $!");
  while(<$fh>){
    chomp;
    my($query, $target) = split /\t/;
    $obs{$query} = $target;
    $obs{$target} = $query;
  }
  close $fh;
  unlink($mapfile);
  unlink($logfile);

  while(my($query, $target) = each(%expected)){
    is($expected{$query}, $obs{$query}, "Check query $query matches to target $target");
  }
}
