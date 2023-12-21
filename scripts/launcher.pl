#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use lib "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/perl5/lib/perl5";
use Spreadsheet::Read;
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;
use Cwd 'abs_path';

die "Provide rundir, excel sample spreadsheet, batch name(transfer label) in order" unless @ARGV == 3;

my ($rundir, $sample_sheet, $batch_name) = @ARGV;
$sample_sheet = abs_path($sample_sheet);

die "$rundir is not valid" unless -d $rundir;
die "$sample_sheet is not valid" unless -s $sample_sheet;

my $dir = '/storage1/fs1/gtac-mgi/Active/CLE/assay/SOMA/batchdir';
my $git_dir = '/storage1/fs1/gtac-mgi/Active/CLE/assay/SOMA/process/git/cle-soma';

my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'Soma.wdl');
my $json_template = File::Spec->join($git_dir, 'Soma.json');

my $group  = '/cle/wdl/tcp';
my $queue  = 'gtac-mgi';
my $docker = 'mgibio/genome_perl_environment:compute1-38';

my $user_group = 'compute-gtac-mgi';

my $out_dir = File::Spec->join($dir, $batch_name);
unless (-d $out_dir) {
    unless (mkdir $out_dir) {
        die "Failed to make directory $out_dir";
    }
}

#parse sample spreadsheet
my $data = Spreadsheet::Read->new($sample_sheet);

#error-checking
my $qc = 'QC Metrics';
my $sheet = $data->sheet($qc);
die "$qc is not a valid sheet in $sample_sheet" unless $sheet;

my $flag = 0;
my %qc_samples;
for my $row ($sheet->rows()) {
    if ($row->[0] =~ /\sNUMBER/) {
        unless ($row->[0] =~ /^ACCESSION\sNUMBER/ and $row->[2] =~ /^SAMPLE\sID/) {
            die "No ACCESSION NUMBER and or SAMPLE ID as the first and third column header for $qc sheet";
        }
        $flag = 1;
    }
    else {
        my $acc = $row->[0];
        $acc =~ s/\s+//g;
        die "$acc not starting with G" unless $acc =~ /^G/;
        unless ($acc =~ /^G\-/) {
            unless ($acc =~ /^G\d+\-\d+$/) {
                die "$acc is not a valid accession number in $qc sheet";
            }
        }
        my $sample = $row->[2];
        $sample =~ s/\s+//g;
        if ($qc_samples{$sample}) {
            die "There are multiple $sample in $qc sheet";
        }
        else {
            $qc_samples{$sample} = 1;
        }
    }
}
die "There is no row for valid column headers in $qc sheet" unless $flag;

my $ss = 'Samplesheet';
$sheet = $data->sheet($ss);
die "$ss is not a valid sheet in $sample_sheet" unless $sheet;

my $ds_str;
my $si_str;
my $seq_id = 2900000000;

my @cases_excluded;
my %samples;
for my $row ($sheet->rows()) {
    next if $row->[0] =~ /Run|Lane/i;
    unless ($row->[0] =~ /\d+/) {
        die "Lane number is expected, Check $ss sheet";
    }
    my ($lane, $flowcell, $name, $index, $exception) = @$row;

    $name =~ s/\s+//g;
    if ($samples{$name}) {
        die "There are multiple $name in $ss sheet";
    }
    else {
        if ($qc_samples{$name}) {
            delete $qc_samples{$name};
        }
        else {
            die "$name in $ss sheet does not exist in $qc sheet";
        }
        $samples{$name} = 1;
    }

    my $lib = $name;
    $lib .= '-lib1' unless $lib =~ /lib/;

    my ($index1, $index2) = $index =~ /([ATGC]{10})\-([ATGC]{10})/;
    
    $exception = 'NONE' unless $exception;
    
    $ds_str .= join ',', $lane, $name, $name, '', $index1, $index2;
    $ds_str .= "\n";
    $si_str .= join "\t", $index1.'-'.$index2, $name, $seq_id, $flowcell, $lane, $lib, $name;
    $si_str .= "\n";

    $seq_id++;
}
if (%qc_samples) {
    my $sample_str = join ",", sort keys %qc_samples;
    die "$sample_str in $qc sheet but not in $ss sheet";
}

## DRAGEN sample sheet
my $dragen_ss  = File::Spec->join($out_dir, 'demux_sample_sheet.csv'); 
my $ss_fh = IO::File->new(">$dragen_ss") or die "Fail to write to $dragen_ss";
$ss_fh->print("[Settings]\n");
$ss_fh->print("AdapterBehavior,trim\n");
$ss_fh->print("AdapterRead1,AGATCGGAAGAGCACACGTCTGAAC\n");
$ss_fh->print("AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGA\n");
$ss_fh->print("OverrideCycles,Y151;I10U9;I10;Y151\n");
$ss_fh->print("[Data]\n");
$ss_fh->print("Lane,Sample_ID,Sample_Name,Sample_Project,index,index2\n");
$ss_fh->print($ds_str);
$ss_fh->close;

## Sample Index
my $si = File::Spec->join($out_dir, 'sample_index');
my $si_fh = IO::File->new(">$si") or die "Fail to write to $si";
$si_fh->print($si_str);
$si_fh->close;

## Input JSON
my $inputs = from_json(`cat $json_template`);
$inputs->{'Soma.OutputDir'}        = $out_dir;
$inputs->{'Soma.IlluminaDir'}      = $rundir;
$inputs->{'Soma.SampleSheet'}      = $si;
$inputs->{'Soma.XferLabel'}        = $batch_name;
$inputs->{'Soma.DemuxSampleSheet'} = $dragen_ss;
$inputs->{'Soma.InputSpreadSheet'} = $sample_sheet;

my $input_json = File::Spec->join($out_dir, 'Soma.json');
my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

$json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
$json_fh->close;

my $out_log = File::Spec->join($out_dir, 'out.log');
my $err_log = File::Spec->join($out_dir, 'err.log');

my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $input_json $wdl";

system $cmd;
#print $cmd."\n";
