process PREPARE_SAMPLES {
    label 'process_single'
    container params.ubuntu_container

    input:
    path(sample_sheet)
    path(read1_list)
    path(read2_list)

    output:
    path('sample_sheet.txt'), emit: sample_sheet
    path('*.fastq.gz')      , emit: fastq
    path('versions.yml')    , emit: versions

    script:
    def previous_fastq_dir = params.previous_fastq_dir ?: ''
    """
    cp ${read1_list} 1.tmp.txt
    cp ${read2_list} 2.tmp.txt

    perl -e '
        use File::Basename;
        use Cwd "abs_path";
        open(R1,"1.tmp.txt"); my @r1 = <R1>; chomp @r1; close R1;
        open(R2,"2.tmp.txt"); my @r2 = <R2>; chomp @r2; close R2;
        open(SS,"${sample_sheet}");
        while(<SS>) {
            chomp;
            my @l = split("\\t", \$_);
            my \$s = \$l[1]."_";
            my \$orig_r1 = (grep /\$s/, @r1)[0];
            my \$orig_r2 = (grep /\$s/, @r2)[0];
            
            my \$base_r1 = basename(\$orig_r1);
            my \$base_r2 = basename(\$orig_r2);
            
            # Create local copies to avoid mutating shared files
            system("cp \$orig_r1 \$base_r1");
            system("cp \$orig_r2 \$base_r2");

            my \$prev_dir = "${previous_fastq_dir}";
            if (\$prev_dir) {
                my (\$n) = \$base_r1 =~/^(\\S+?)_/;
                my @p_r1 = glob(\$prev_dir."/".\$n."*_R1_001.fastq.gz");
                my @p_r2 = glob(\$prev_dir."/".\$n."*_R2_001.fastq.gz");
                unless (@p_r1 and @p_r1 == 1 and @p_r2 and @p_r2 == 1) {
                    die "fail to get previous R1 and or R2 for \$n";
                }
                my \$rc1 = system "cat \$p_r1[0] >> \$base_r1";
                my \$rc2 = system "cat \$p_r2[0] >> \$base_r2";
                unless (\$rc1 == 0 and \$rc2 == 0) {
                    die "R1 and or R2 cat failed for \$n";
                }
            }
            # Use absolute path of the local copy
            my \$abs_r1 = abs_path(\$base_r1);
            my \$abs_r2 = abs_path(\$base_r2);
            print join("\\t", @l, \$abs_r1, \$abs_r2),"\\n";
        }
        close SS;
    ' > sample_sheet.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print \$];')
    END_VERSIONS
    """
}
