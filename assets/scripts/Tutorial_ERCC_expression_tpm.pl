#!/usr/bin/perl

use strict;
use warnings;

use IO::File;

my $ercc_file = 'ERCC_Controls_Analysis.txt';
my $tpm_file = 'gene_tpm_all_samples.tsv';
my $ercc_tpm_file = 'ercc_tpm.tsv';

my $ercc_fh = IO::File->new($ercc_file,'r');
unless ($ercc_fh) { die('Failed to find file: '. $ercc_file) }

my %ercc_data;
while (my $ercc_line = $ercc_fh->getline) {
    chomp($ercc_line);
    if ($ercc_line =~ /^Re/) { next; }
    #my ($resort,$id,$subgroup,$mix1,$mix2,$fold_change,$log2)
    my @ercc_entry = split("\t",$ercc_line);
    $ercc_data{$ercc_entry[1]} = \@ercc_entry;
}

my @labels = qw/HBR_Rep1 HBR_Rep2 HBR_Rep3 UHR_Rep1 UHR_Rep2 UHR_Rep3/;

my $tpm_fh = IO::File->new($tpm_file,'r');
unless ($tpm_fh) { die('Failed to find file: '. $tpm_file); }

my $ercc_tpm_fh = IO::File->new($ercc_tpm_file,'w');
unless ($ercc_tpm_fh) { die('Failed to open file: '. $ercc_tpm_file); }

my %tpm_data;
print $ercc_tpm_fh "ID\tSubgroup\tLabel\tMix\tConcentration\tTPM\n";
while (my $tpm_line = $tpm_fh->getline) {
    chomp($tpm_line);
    #print "$tpm_line\n";
    my @tpm_entry = split('\t',$tpm_line);
    if ($ercc_data{$tpm_entry[0]}) {
        my $id = $tpm_entry[0];
        my $subgroup = $ercc_data{$id}->[2];
	#print "$id\t$subgroup\n";
 	for (my $i = 0; $i < scalar(@labels); $i++) {
            my $tpm = $tpm_entry[$i+1];
            my $label = $labels[$i];
            my $conc;
            my $mix;
            if ($label =~ /UHR/) {
                $mix = 1;
                $conc = $ercc_data{$id}->[3];
            } else {
                $mix = 2;
                $conc = $ercc_data{$id}->[4];
            }
            print $ercc_tpm_fh $id ."\t". $subgroup ."\t". $label ."\t". $mix ."\t". $conc ."\t". $tpm ."\n";
        }
    }
}


exit;
