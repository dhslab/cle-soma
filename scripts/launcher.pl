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
#my $json_template = File::Spec->join($git_dir, 'Soma.json');
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
my $sheet = $data->sheet(1);

my $ds_str;
my $si_str;
my $seq_id = 2900000000;

my @cases_excluded;

for my $row ($sheet->rows()) {
    next if $row->[0] =~ /Run|Lane/i;
    unless ($row->[0] =~ /\d+/) {
        die "Lane number is expected, Check sample sheet spreadsheet";
    }
    my ($lane, $flowcell, $name, $index, $exception) = @$row;

    $name =~ s/\s+//g;
    my $lib = $name;
    $lib .= '-lib1' unless $lib =~ /lib/;

    my ($index1, $index2) = $index =~ /([ATGC]{10})\-([ATGC]{10})/;
    my $fix_index2 = rev_comp($index2);
    
    $exception = 'NONE' unless $exception;
    
    $ds_str .= join ',', $lane, $name, $name, '', $index1, $fix_index2;
    $ds_str .= "\n";
    $si_str .= join "\t", $index1.'-'.$fix_index2, $name, $seq_id, $flowcell, $lane, $lib, $name;
    $si_str .= "\n";

    $seq_id++;
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

sub rev_comp {
    my $index = shift;
    my $revcomp = reverse $index;
    $revcomp =~ tr/ATGCatgc/TACGtacg/;

    return $revcomp;
}
