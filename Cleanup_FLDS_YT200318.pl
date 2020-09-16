#!/usr/bin/perl -w

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print " -lst file (fastq-file-name) -seqDir (FASTQ DIRECTORY)\n";

    print "OPTION:
    -trimo  (1)/0 execute trimming of adaptor & low-quality
    -bowtie2  (1)/0 remove phiX
    -contami  index file for bowtie2
    -cut  (1)/0 execute trimming of primers
    -low  (1)/0 execute Low-complexity remove
    -dup  (T)/F execute duplication remove
    -lib (FLDS)/oligo/BL   Library type
    -ribo (T)/F  excute SORTMERNA\n";

    print "OUTPUT: ./cleanup/ LOG: ./logs/
    TRIMMOMATIC xxx_TP_R1.fq xxx_TP_R2.fq xxx_TU_R1.fq xxx_TU_R1.fq
    DECONTAMI   xxx_nocont_R1.fq xxx_nocont_R2.fq
    CUTADAPT    xxx_CP_R1.fq xxx_CP_R2.fq
    PRINSEQ     xxx_PP_R1.fq xxx_PP_R2.fq
    SORTMERNA   xxx_SP_R1.fq xxx_SP_R2.fq xx_SP_rrn_R1.fq xxx_SP_rrn_R2.fq
    \n";
    exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};
$|=1;

$HOME=$ENV{"HOME"};

## 0)  SETTING BIN & DB PATH
    $TRIMMOMATIC="${HOME}/biotools/local/assembly/Trimmomatic-0.36";
	$BOWTIE2="${HOME}/biotools/local/mapping/bowtie2-2.3.5.1";
    $CUTADAPT="${HOME}/biotools/rhel6/miniconda3/bin";
	$PRINSEQ="${HOME}/biotools/local/assembly/prinseq-lite-0.20.4";
    $PRINSEQP="${HOME}/biotools/local/assembly/prinseq_parallel";
#	$SORTMERNA="${HOME}/biotools/local/assembly/sortmerna-2.1-linux-64";
    $SORTMERNA="${HOME}/biotools/local/assembly/sortmerna-2.1b";
	$phiX_index="${HOME}/biotools/local/mapping/mappedIndex/PhiX_index/PhiX_IDX";
#	$MoCV1_index="${HOME}/biotools/local/mapping/mappedIndex/MoCV1_A_IDX/MoCV1_A_IDX";
#    $ADAPTORS_D="${HOME}/Desktop/work_CLC/bin";

    $OUT_DIR= exists $option->{-outdir} ? $option->{-outdir}->[0] : "cleanup";;


## 1)  CHECK FASTQ

$seq_dir= exists $option->{-seqDir} ? $option->{-seqDir}->[0] : "seq_data";

($sampinfo)= &create_sampinfo($option->{-lst}->[0]);
($check)= &check_fastq_files($sampinfo,$seq_dir);

    print "CHECK FASTQ $check\n";
	exit if ($check == 0);


## 2) MAKE FOLDER

if (-d "$OUT_DIR") { } else {system ("mkdir $OUT_DIR")};
if (-d "logs") { } else {system ("mkdir logs")};

## 3) RUN PIPELINE
$exec_trimo= exists $option->{-trimo} ? $option->{-trimo}->[0] : 1;
$exec_bowtie2= exists $option->{-bowtie2} ? $option->{-bowtie2}->[0] : 1;
$contami= exists $option->{-contami} ? $option->{-contami}->[0] : "F";
$exec_cutadapt= exists $option->{-cut} ? $option->{-cut}->[0] : 1;
$libtype= exists $option->{-lib} ? $option->{-lib}->[0] : "FLDS";
$duplicate= exists $option->{-dup} ? $option->{-dup}->[0] : "T";
$exec_low= exists $option->{-low} ? $option->{-low}->[0] : 1;
$rmribo= exists $option->{-ribo} ? $option->{-ribo}->[0] : "T";

($run)= &run_meta_cleanup($sampinfo,$contami,$duplicate,$libtype,$seq_dir,$rmribo);

## 4)  COUNT SUMMARY

($run_steps)= &summary_count_reads($sampinfo,$seq_dir);

print "RUN STEP: @$run_steps\n";
print "OUTCOUNT: summary_count_cleanup.txt\n";

#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====

sub run_meta_cleanup{
    my ($sampinfo,$contami,$duplicate,$libtype,$seq_dir,$rmribo)= @_;
    my ($name,$readfile1,$readfile2);
    
    $run= 0;
    
    foreach $name (@{$sampinfo->{samples}}){
        $readfile1= "${seq_dir}/$sampinfo->{pe}->{$name}->[0]";
        $readfile2= "${seq_dir}/$sampinfo->{pe}->{$name}->[1]";
        print "SAMPLE: $name\n";
        print "FASRQ:\t$readfile1\n\t$readfile2\n";
        
        #########    trimmomatic    #########
        print "NAME: $name\n";
        
        if ($exec_trimo == 1) {
            ($outfiles)= &run_trimmomatic($name,$readfile1,$readfile2);
            ($readfile1,$readfile2)= @$outfiles;
            print "FIN trimmomatic $name \n";
        }
        
        #########    bowtie2    #########
        if ($exec_bowtie2 == 1) {
            
            ($outfiles)= &run_rm_contami($name,$readfile1,$readfile2,$phiX_index,$contami);
            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN bowtie2 phiX & CONTAMI $name \n";
        }
        
        #########    cutadaptor    #########
        
        if ($exec_cutadapt == 1) {
            
            if ($libtype eq "FLDS") {
                ($outfiles)= &run_cutadapt_FLDS($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            } elsif ($libtype eq "oligo") {
                ($outfiles)= &run_cutadapt_oligo($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            }elsif ($libtype eq "BL") {
                ($outfiles)= &run_cutadapt_FLDS_BL($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            }
            
            
            print "FIN cutadaptor $name \n";
            
        }
        
        #########    DUST low complexicy    ###########
        if ($exec_low == 1) {
            ($outfiles)= &run_rm_low_complexity($name,$readfile1,$readfile2);
            
            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN $name DUST32\n";
            
        }
        
        
        #########    duplication    ###########
        
        if ($duplicate eq "T") {
            ($outfiles)= &run_rm_duplicate($name,$readfile1,$readfile2);
            
            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN $name duplication\n";
        }
        
        
        #########    SORTMERNA    ###########
        
        if ($rmribo eq 'T') {
            print "****** START SORTMERNA\n";
            
            ## PAIRED-END
            $readType= "pe";
            ($outfiles)= &run_sortmerna($name,$readType,$readfile1,$readfile2);
                        
            print "FIHISHED SORTMERNA\n";
            
        }
        
        
        
    }
    
    ###########   logs
    @flists=<logs_*>;
    system("mv logs_* ./logs/") if ( @flists > 0 );

    $run= 1;
    
    return($run);
}

#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


sub run_trimmomatic {
    my($name,$readfile1,$readfile2) = @_;
    my ($outfiles);
    
    system ("java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 24 -phred33 -trimlog logs_trimo_${name} $readfile1 $readfile2 ./${OUT_DIR}/${name}_TP_R1.fq ./${OUT_DIR}/${name}_OP_R1.fq ./${OUT_DIR}/${name}_TP_R2.fq ./${OUT_DIR}/${name}_OP_R2.fq ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50");
    ## TruSeq3-PE-1_NISHI_FIN2.fa
    
    $outfiles= ["./${OUT_DIR}/${name}_TP_R1.fq",
    "./${OUT_DIR}/${name}_TP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_contami {
    my($name,$readfile1,$readfile2,$phiX_index,$contami) = @_;
    
    #PhiX
    
    system ("${BOWTIE2}/bowtie2 -p 4 -N 1 --un-conc ./${OUT_DIR}/${name}_nocont.fq --al-conc ./${OUT_DIR}/${name}_phi.fq -x $phiX_index -1 $readfile1 -2 $readfile2 > /dev/null 2> logs_rmphi_${name}.txt");
    
    system ("mv ./${OUT_DIR}/${name}_nocont.1.fq ./${OUT_DIR}/${name}_nocont_R1.fq");
    system ("mv ./${OUT_DIR}/${name}_nocont.2.fq ./${OUT_DIR}/${name}_nocont_R2.fq");
    
    $outfiles= ["./${OUT_DIR}/${name}_nocont_R1.fq",
    "./${OUT_DIR}/${name}_nocont_R2.fq"];
    
    
    if ($contami eq "T") {
        $contami_ref_index=$contami;
        #MoCV1
        system ("${BOWTIE2}/bowtie2 -p 4 -N 1 --un-conc ./${OUT_DIR}/${name}_nocont2.fq --al-conc ./${OUT_DIR}/${name}_contami2.fq -x $contami_ref_index -1 ./${OUT_DIR}/${name}_nocont.1.fq -2 ./${OUT_DIR}/${name}_nocont.2.fq > /dev/null 2> logs_rmmov_${name}.txt");
        
        $outfiles= ["./${OUT_DIR}/${name}_nocont2.1.fq",
        "./${OUT_DIR}/${name}_nocont2.2.fq"];
    }
    
    return ($outfiles);
    
}

sub run_cutadapt_FLDS {
    my($name,$readfile1,$readfile2) = @_;
    
    $UPM="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG";
    $UPM_comp="CCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG";
    $UPM_short="CTAATACGACTCACTATAGGGC";
    $UPM_short_comp="GCCCTATAGTGAGTCGTATTAG";
    $U2= "GACGTAAGAACGTCGCACCA";
    $U2_comp="TGGTGCGACGTTCTTACGTC";
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 30 -n 5 -O 10 \
    -g $UPM -a $U2 -G $U2_comp -A $UPM_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 30 -n 5 -O 10 \
    -g $U2_comp -a $UPM_comp -G $UPM -A $U2 \
    -o out_cut2_R1.fq -p out_cut2_R2.fq out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 30 -n 5 -O 10 \
    -g $UPM_short -A $UPM_short_comp \
    -o out_cut3_R1.fq -p out_cut3_R2.fq out_cut2_R1.fq out_cut2_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD3: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 30 -n 5 -O 10 \
    -a $UPM_short_comp -G $UPM_short \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq out_cut3_R1.fq out_cut3_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD4: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq","out_cut2_R1.fq",
    "out_cut2_R2.fq","out_cut3_R1.fq","out_cut3_R2.fq",);
    
    print "FIN cutadaptor $name \n";
    
    $outfiles= ["./${OUT_DIR}/${name}_CP_R1.fq","./${OUT_DIR}/${name}_CP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_cutadapt_FLDS_BL {
    my($name,$readfile1,$readfile2) = @_;
    
    $UPM_ALL= "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG";
    $UPM_ALL_comp= "CCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG";
    $UPM_short= "CTAATACGACTCACTATAGGGC";
    $UPM_short_comp= "GCCCTATAGTGAGTCGTATTAG";
    $BL1= "CTGTAGGCACCATCAAT";
    $BL1_comp= "ATTGATGGTGCCTACAG";
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $UPM_ALL -a $BL1_comp -G $BL1 -A $UPM_ALL_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $BL1 -a $UPM_ALL_comp -G $UPM_ALL -A $BL1_comp \
    -o out_cut2_R1.fq -p out_cut2_R2.fq out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $UPM_short -A $UPM_short_comp \
    -o out_cut3_R1.fq -p out_cut3_R2.fq out_cut2_R1.fq out_cut2_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD3: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -a $UPM_short_comp -G $UPM_short \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq out_cut3_R1.fq out_cut3_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD4: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq","out_cut2_R1.fq",
    "out_cut2_R2.fq","out_cut3_R1.fq","out_cut3_R2.fq",);
    
    print "FIN cutadaptor $name \n";
    
    $outfiles= ["./${OUT_DIR}/${name}_CP_R1.fq","./${OUT_DIR}/${name}_CP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_cutadapt_oligo {
    my($name,$readfile1,$readfile2) = @_;
    
    
    open (LOG,">log_cutadaptor_${name}.txt");
    
    $SMARTer_oligo_II= "AAGCAGTGGTATCAACGCAGAGTACATGGG";
    $SMARTer_oligo_II_comp= "CCCATGTACTCTGCGTTGATACCACTGCTT";
    $SMARTer_oligo_II_modify= "AAGCAGTGGTATCAACGCAGAGTACT";
    $SMARTer_oligo_II_modify_comp= "AGTACTCTGCGTTGATACCACTGCTT";
    
    ##### plus #####
    print LOG "###### $readfile1 \n##SMARTer-oligo-II / SMARTer_oligo_II_modify_comp\n";
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $SMARTer_oligo_II -a $SMARTer_oligo_II_modify_comp -G $SMARTer_oligo_II_modify -A $SMARTer_oligo_II_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    ##### minus #####
    print LOG "###### $readfile1 \n##SMARTer_oligo_II_modify / SMARTer_oligo_II_comp\n";
    
    $cmd_cutadapt="cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $SMARTer_oligo_II_modify -a $SMARTer_oligo_II_comp -G $SMARTer_oligo_II -A $SMARTer_oligo_II_modify_comp \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print "CMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq");
    
    print "FIN cutadaptor $name \n";
    
    $outfiles= ["./${OUT_DIR}/${name}_CP_R1.fq","./${OUT_DIR}/${name}_CP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_low_complexity {
    my ($name,$readfile1,$readfile2) = @_;
    my ($cwd,$abs_readfile1,$abs_readfile2,$cmd_prinseq_low);
    
    use File::Spec;
    use Cwd;
    
    $cwd=Cwd::getcwd();
    
    $abs_readfile1= "${cwd}/$readfile1";
    $abs_readfile2= "${cwd}/$readfile2";
    
    $cmd_prinseq_low="$PRINSEQP/prinseq_parallel_rev.sh -verbose \
    -out_format 3 -log ${cwd}/logs_dust32_${name}.txt -lc_method dust -lc_threshold 32 \
    $abs_readfile1 $abs_readfile2 ${cwd}/${OUT_DIR}/${name}_prinseq PE 8";
    $cmd_prinseq_low=~ tr/\n//d;
    print "$cmd_prinseq_low\n";
    
    system ("$cmd_prinseq_low");
    system("mv ./${OUT_DIR}/${name}_prinseq_1.fastq ./${OUT_DIR}/${name}_LP_R1.fq");
    system("mv ./${OUT_DIR}/${name}_prinseq_2.fastq ./${OUT_DIR}/${name}_LP_R2.fq");
    system("rm -rf ${cwd}/${OUT_DIR}/split*");
    
    $outfiles= ["./${OUT_DIR}/${name}_LP_R1.fq","./${OUT_DIR}/${name}_LP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_duplicate {
    my($name,$readfile1,$readfile2) = @_;
    
    system ("${PRINSEQ}/prinseq-lite.pl -verbose -fastq $readfile1 -fastq2 $readfile2 -derep 2 -out_format 3 -log logs_dup_${name}.txt");
    
    system ("mv ./${OUT_DIR}/${name}_*_R1_prinseq_good_*.fastq ./${OUT_DIR}/${name}_PP_R1.fq");
    system ("mv ./${OUT_DIR}/${name}_*_R2_prinseq_good_*.fastq ./${OUT_DIR}/${name}_PP_R2.fq");
    
    system ("mv ./${OUT_DIR}/${name}_*_R1_prinseq_good_singletons_*.fastq ./${OUT_DIR}/${name}_PU_R1.fq");
    system ("mv ./${OUT_DIR}/${name}_*_R2_prinseq_good_singletons_*.fastq ./${OUT_DIR}/${name}_PU_R2.fq");
    
    system ("rm ./${OUT_DIR}/${name}_*_prinseq_bad_*.fastq");
    
    $outfiles= ["./${OUT_DIR}/${name}_PP_R1.fq", "./${OUT_DIR}/${name}_PP_R2.fq"];
    
    return ($outfiles);
    
    
}

sub run_rm_duplicate_P {
    my($name,$readfile1,$readfile2) = @_;
    my ($cwd,$abs_readfile1,$abs_readfile2,$cmd_prinseq_low);
    
    use File::Spec;
    use Cwd;
    
    $cwd=Cwd::getcwd();
    $abs_readfile1= "${cwd}/$readfile1";
    $abs_readfile2= "${cwd}/$readfile2";
    
    $cmd_prinseq_dup="$PRINSEQP/prinseq_parallel_rev.sh -verbose \
    -out_format 3 -log ${cwd}/logs_dup_${name}.txt -derep 2 \
    $abs_readfile1 $abs_readfile2 ${cwd}/${OUT_DIR}/${name}_prinseq PE 8";
    $cmd_prinseq_dup=~ tr/\n//d;
    print "$cmd_prinseq_dup\n";
    
    system ("$cmd_prinseq_dup");
    system("mv ./${OUT_DIR}/${name}_prinseq_1.fastq ./${OUT_DIR}/out_PE_R1.fq");
    system("mv ./${OUT_DIR}/${name}_prinseq_2.fastq ./${OUT_DIR}/out_PE_R2.fq");
    system("rm -rf ${cwd}/${OUT_DIR}/split*");
    
    system ("rm ./${OUT_DIR}/${name}_*_prinseq_bad_*.fastq");
    
    $outfiles= ["./${OUT_DIR}/out_PE_R1.fq","./${OUT_DIR}/out_PE_R2.fq"];
    
    return ($outfiles);
    
    
}

sub run_sortmerna {
    my($name,$readType,$readfile1,$readfile2) = @_;
    my ($infile,$cmd);
    print "$readfile1 $readfile2\n";
    
    if ($readType eq "pe") {
        $infile="./${OUT_DIR}/out_PE_R12.fq";
        system ("bash ${SORTMERNA}/scripts/merge-paired-reads.sh $readfile1 $readfile2 $infile");
    } elsif ($readType eq "se") {
        $infile= $readfile1;
    }
    
    $cmd="${SORTMERNA}/sortmerna  --ref \
${SORTMERNA}/rRNA_databases/silva-bac-16s-id90.fasta,${SORTMERNA}/index/silva-bac-16s-db:\
${SORTMERNA}/rRNA_databases/silva-bac-23s-id98.fasta,${SORTMERNA}/index/silva-bac-23s-db:\
${SORTMERNA}/rRNA_databases/silva-arc-16s-id95.fasta,${SORTMERNA}/index/silva-arc-16s-db:\
${SORTMERNA}/rRNA_databases/silva-arc-23s-id98.fasta,${SORTMERNA}/index/silva-arc-23s-db:\
${SORTMERNA}/rRNA_databases/silva-euk-18s-id95.fasta,${SORTMERNA}/index/silva-euk-18s-db:\
${SORTMERNA}/rRNA_databases/silva-euk-28s-id98.fasta,${SORTMERNA}/index/silva-euk-28s-db:\
${SORTMERNA}/rRNA_databases/rfam-5s-database-id98.fasta,${SORTMERNA}/index/rfam-5s-db:\
${SORTMERNA}/rRNA_databases/rfam-5.8s-database-id98.fasta,${SORTMERNA}/index/rfam-5.8s-db \
 --reads $infile --sam --num_alignments 1 --fastx \
 --aligned ./${OUT_DIR}/${name}_rrn --other ./${OUT_DIR}/${name}_norrn \
 --log -v -a 12";
        $cmd=~ s/\n//g;
    
    system ($cmd);
    unlink ("$infile");
    
    if ($readType eq "pe") {
        system ("bash ${SORTMERNA}/scripts/unmerge-paired-reads.sh ./${OUT_DIR}/${name}_norrn.fq ./${OUT_DIR}/${name}_SP_R1.fq ./${OUT_DIR}/${name}_SP_R2.fq");
        system ("bash ${SORTMERNA}/scripts/unmerge-paired-reads.sh ./${OUT_DIR}/${name}_rrn.fq ./${OUT_DIR}/${name}_rrn_R1.fq ./${OUT_DIR}/${name}_rrn_R2.fq");
        $outfiles= ["./${OUT_DIR}/${name}_SP_R1.fq",
        "./${OUT_DIR}/${name}_SP_R2.fq"];
        
        unlink ("./${OUT_DIR}/${name}_norrn.fq","./${OUT_DIR}/${name}_rrn.fq");
        
    } else {
        $outfiles= ["./${OUT_DIR}/${name}_norrn.fq",
        "./${OUT_DIR}/${name}_rrn.fq"];
    }
    
    return ($outfiles);
    
}


##==== === === === === === === === === === === === === === ===
sub summary_count_reads {
    my($sampinfo,$seq_dir)= @_;
    my ($ifo);
    $info= {};
    
    $info->{steps}= [qw(raw trim no_cont cutadapt low prin_pe sort_pe sort_se1 sort_se2)];
    
    
    foreach $name (@{$sampinfo->{samples}}) {
        
        $onefile_path= "${seq_dir}/$sampinfo->{pe}->{$name}->[0]";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{raw}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_TP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{trim}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_nocont_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{no_cont}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_CP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{cutadapt}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_LP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{low}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_PP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{prin_pe}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_SP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{sort_pe}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_SU_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{sort_se1}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_SU_R2.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{sort_se2}= $count;
        
    }
    
    
    open (CO,">summary_count_cleanup.txt");
    print CO (join "\t",("SAMPLE",@{$info->{steps}}))."\n";
    
    foreach $name (@{$sampinfo->{samples}}) {
        @oneres= ($name);
        
        foreach $step (@{$info->{steps}}) {
            push @oneres,$info->{sample}->{$name}->{$step};
        }
        
        print CO (join "\t",@oneres)."\n";
    }
    
    return ($info->{steps})
}



sub count_reads {
    my($readfile)=@_;
    my ($count);
    
    if ($readfile=~ /.gz$/) {
        $cmd="zcat $readfile".' | echo "$((`wc -l`/4))"';
        $count= `$cmd`;
    } else {
        $cmd="cat $readfile".' | echo "$((`wc -l`/4))"';
        $count= `$cmd`;
    }
    #print "$cmd\n";
    chomp ($count);
    
    return $count;
    
}



#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
sub create_sampinfo {
    my($samplist)= @_;
    
    open (LIS,"$samplist");
    $sampinfo= {};
    while(<LIS>){
        chomp;
        next if (/^(\#|\s*$)/);
            @item= split/\t/;
        if (@item == 2) {
            $name= $1 if ($item[0] =~ /^([^_]+)/);
            if ($item[1] =~ /^$name/) {
                push @{$sampinfo->{samples}},$name;
                $sampinfo->{pe}->{$name}=[@item];
                #print "@item\n";
            } else {
                $name= $1 if ($item[0] =~ /^(.+).(fq|fastq|fq.gz|fastq.gz)$/);
                push @{$sampinfo->{samples}},$name;
                $sampinfo->{pe}->{$name}=[@item];
            }
            $sampinfo->{sampnum}->{$_}++;
            
            if (exists $sampinfo->{name_check}->{$name}) {
                print "ERROR: SAMPEL NAME IS NOT UNIQUE NAME\n";
                exit;
            } else {
                $sampinfo->{name_check}->{$name}++;
            }
        } elsif (@item == 3) {
            $name= $item[2];
            push @{$sampinfo->{samples}},$name;
            $sampinfo->{pe}->{$name}=[@item[0,1]];
            $sampinfo->{sampnum}->{$_}++;
            if (exists $sampinfo->{name_check}->{$name}) {
                print "ERROR: SAMPEL NAME IS NOT UNIQUE NAME\n";
                exit;
            } else {
                $sampinfo->{name_check}->{$name}++;
            }
            
        } elsif (@item > 2) {
            print "ERROR: More than 2 items in a line\n";
            exit;
        }
        
    }
    
    return ($sampinfo);
    
}


sub check_fastq_files {
    my ($sampinfo,$seq_dir)= @_;
    
    my($check);
    $check= 1;
    
    if (-d $seq_dir) {
        foreach $name (@{$sampinfo->{samples}}) {
            
            $files = $sampinfo->{pe}->{$name};
            
            if (-f "${seq_dir}/$$files[0]") {
                
            } else {
                print "ERROR: DO NOT EXISTS $$files[0]\n";
                $check=0;
                last;
            }
            if (-f "${seq_dir}/$$files[1]") {
                
            } else {
                print "ERROR: DO NOT EXISTS $$files[1]\n";
                $check=0;
                last;
            }
        }
    } else {
        print "ERROR: DO NOT EXISTS \"$seq_dir\" DIRECTORY\n";
        $check=0;
    }
    
    return ($check);
}










