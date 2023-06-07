#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;
use IO::File;
use File::Spec;
use File::Basename;

umask 002;

die "Provide SOMA workflow output dir" unless @ARGV and @ARGV == 1;
my $dir = $ARGV[0];

unless (-d $dir) {
    die "SOMA workflow output dir: $dir not existing";
}

my @metrics_names = qw(
    HAPLOTECT_SCORE
    INFORMATIVE_SITES
    TOTAL_GIGA_BASE 
    TOTAL_MEGA_READ 
    PCT_MAPPED_READ 
    PCT_PROPERLY_PAIRED_READ 
    MISMATCH_RATE_R1
    MISMATCH_RATE_R2
    PCT_Q30_BASE
    MEAN_INSERT_LENGTH
    AVG_ALIGN_COVERAGE
    PCT_ALIGN_READ_ON_TARGET
    PCT_ALIGNED_BASE_ON_TARGET
    PCT_TARGET_COVERAGE_GT_1500
    PCT_TARGET_COVERAGE_GT_100
);

my @headers = ('Case', @metrics_names);

my $out_file = $dir.'/QC_metrics.tsv';
my $out_fh = IO::File->new(">$out_file") or die "Failed to write to $out_file";
$out_fh->print(join "\t", @headers);
$out_fh->print("\n");

my @case_dirs = glob($dir."/H_QT*lib*");

for my $case_dir (@case_dirs) {
    my ($case_name) = basename($case_dir) =~ /^(H_\S+lib\d+)/;
    my $dragen_dir = File::Spec->join($case_dir, 'dragen');
    my ($mapping_metrics) = glob($dragen_dir."/*.mapping_metrics.csv");
    my ($target_metrics)  = glob($dragen_dir."/*.target_bed_coverage_metrics.csv");
    my ($haplotect_out)   = File::Spec->join($case_dir, $case_name.'.haplotect.txt');

    unless (-s $mapping_metrics and -s $target_metrics and -s $haplotect_out) {
        die "No dragen mapping and or target metrics and or haplotect out for $case_dir";
    }
   
    my %QC;
    my $str = `tail -1 $haplotect_out`;
    my @cols = split /\s+/, $str;
    $QC{HAPLOTECT_SCORE}   = $cols[6];
    $QC{INFORMATIVE_SITES} = $cols[2];

    my $map_fh = IO::File->new($mapping_metrics) or die "fail to open $mapping_metrics for read";
    while (my $line = $map_fh->getline) {
        last if $line =~ /^MAPPING\/ALIGNING\sPER\sRG,/;
        chomp $line;
        if ($line =~ /Total input reads,(\d+)/) {
            $QC{TOTAL_MEGA_READ} = sprintf("%.1f", $1/1000000);
        }
        elsif ($line =~ /Mapped reads,\d+,(\S+)/) {
            $QC{PCT_MAPPED_READ} = $1;
        }
        elsif ($line =~ /Properly paired reads,\d+,(\S+)/) {
            $QC{PCT_PROPERLY_PAIRED_READ} = $1;
        }
        elsif ($line =~ /Total bases,(\d+)/) {
            $QC{TOTAL_GIGA_BASE} = sprintf("%.1f", $1/1000000000);
        }
        elsif ($line =~ /Mismatched bases R1,\d+,(\S+)/) {
            $QC{MISMATCH_RATE_R1} = $1;
        }
        elsif ($line =~ /Mismatched bases R2,\d+,(\S+)/) {
            $QC{MISMATCH_RATE_R2} = $1;
        }
        elsif ($line =~ /Q30 bases,\d+,(\S+)/) {
            $QC{PCT_Q30_BASE} = $1;
        }
        elsif ($line =~ /Insert length: mean,(\S+)/) {
            $QC{MEAN_INSERT_LENGTH} = $1;
        }
    }
    $map_fh->close;

    my $target_fh = IO::File->new($target_metrics) or die "fail to open $target_metrics for read";
    while (my $l = $target_fh->getline) {
        chomp $l;
        if ($l =~ /Aligned bases in target region,\d+,(\S+)/) {
            $QC{PCT_ALIGNED_BASE_ON_TARGET} = $1;
        }
        elsif ($l =~ /Average alignment coverage over target region,(\S+)/) {
            $QC{AVG_ALIGN_COVERAGE} = $1;
        }
        elsif ($l =~ /PCT of target region with coverage \[1500x\: inf\),(\S+)/) {
            $QC{PCT_TARGET_COVERAGE_GT_1500} = $1;
        }
        elsif ($l =~ /PCT of target region with coverage \[\s100x\: inf\),(\S+)/) {
            $QC{PCT_TARGET_COVERAGE_GT_100} = $1;
        }
        elsif ($l =~ /Aligned reads in target region,\d+,(\S+)/) {
            $QC{PCT_ALIGN_READ_ON_TARGET} = $1;
        }
    }
    $target_fh->close;
   
    my @values = ($case_name);
    for my $metrics_name (@metrics_names) {
        push @values, $QC{$metrics_name};
    }
    $out_fh->print(join "\t", @values);
    $out_fh->print("\n");
}

$out_fh->close;
