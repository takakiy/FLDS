#!/usr/bin/perl -w


## update: 2022.04.18


if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print " -lst file (fastq-file-name) -seqDir (FASTQ DIRECTORY)\n";

    print "OPTION:
    -trimo  (T)/F trimming illumina adaptor & low-quality (trimmomatic)
    -bowtie2  (T)/F removing phiX (bowtie2)
         -contami  T/(F) removing contami reads
         -idx  index file for bowtie2
    -cut  (T)/F trimming primers/adaptor (cutadapt)
         -lib (FLDS)/oligo/BL   Library type
    -low  (T)/F removing Low-complexity (prinseq++)
    -dup  (T)/F removing duplication (prinseq-lite.pl)
    -ribo (T)/F  removing ribosomal RNA (sortmerna)
    
    -outdir Output directory (cleanup)
        \n";

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

## 0)  SETTING BIN DIRECTORY & DB FILE PATH
    $TRIMMOMATIC="${HOME}/biotools/local/assembly/Trimmomatic-0.39";
	$BOWTIE2="${HOME}/biotools/local/mapping/bowtie2-2.3.5.1";
    $CUTADAPT="${HOME}/biotools/rhel6/miniconda3/bin";
    $PRINSEQ="${HOME}/biotools/local/assembly/prinseq-lite-0.20.4";

 ## prinseq++
    $PRINSEQ_PP="${HOME}/biotools/local/assembly/bin";

    $SMR="${HOME}/biotools/local/assembly/sortmerna-4.3.4-Linux/bin";
    $SMRDB="${HOME}/biotools/local/assembly/sortmerna-4.3.4-Linux/sortmerna/data/rRNA_databases";

	$phiX_index="${HOME}/biotools/local/mapping/mappedIndex/PhiX_index/PhiX_IDX";
#	$MoCV1_index="${HOME}/biotools/local/mapping/mappedIndex/MoCV1_A_IDX/MoCV1_A_IDX";

    $OUT_DIR= exists $option->{-outdir} ? $option->{-outdir}->[0] : "cleanup";;


## 1)  CHECK FASTQ

#$seq_dir= exists $option->{-seqDir} ? $option->{-seqDir}->[0] : "seq_data";

#($sampinfo)= &create_sampinfo($option->{-lst}->[0]);
($sampinfo)= &create_sampinfo_path($option->{-lst}->[0]);
($check)= &check_fastq_files($sampinfo);

    print "CHECK FASTQ $check\n";
	exit if ($check == 0);


## 2) MAKE FOLDER

if (-d "$OUT_DIR") { } else {system ("mkdir $OUT_DIR")};
if (-d "logs") { } else {system ("mkdir logs")};

## 3) RUN PIPELINE CONTROL
$execinfo= {};
$execinfo->{trimo}= exists $option->{-trimo} ? $option->{-trimo}->[0] : "T";
$execinfo->{bowtie2}= exists $option->{-bowtie2} ? $option->{-bowtie2}->[0] : "T";
$execinfo->{contami}= exists $option->{-contami} ? $option->{-contami}->[0] : "F";
$execinfo->{cut}= exists $option->{-cut} ? $option->{-cut}->[0] : "T";
$execinfo->{lib}= exists $option->{-lib} ? $option->{-lib}->[0] : "FLDS";
$execinfo->{dup}= exists $option->{-dup} ? $option->{-dup}->[0] : "T";
$execinfo->{low}= exists $option->{-low} ? $option->{-low}->[0] : "T";
$execinfo->{ribo}= exists $option->{-ribo} ? $option->{-ribo}->[0] : "T";
$execinfo->{idx}= exists $option->{-idx} ? $option->{-idx}->[0] : "F";


open (LOG2,">out_Cleanup_FLDS.logs.txt");

($run)= &run_FLDS_cleanup($sampinfo,$execinfo);

## 4)  COUNT SUMMARY

($run_steps)= &summary_count_reads($sampinfo,$seq_dir);

print "RUN STEP: @$run_steps\n";
print "OUTCOUNT: summary_count_cleanup.txt\n";

#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====

sub run_FLDS_cleanup{
    my ($sampinfo,$execinfo)= @_;
    my ($name,$readfile1,$readfile2);
    
    $run= 0;
    
    foreach $name (@{$sampinfo->{samples}}){
        $readfile1= "$sampinfo->{pe}->{$name}->[0]";
        $readfile2= "$sampinfo->{pe}->{$name}->[1]";
        print "SAMPLE: $name\n";
        print "FASRQ:\t$readfile1\n\t$readfile2\n";
        print LOG2 "## SAMPLE: $name\n";
        print LOG2 "## INPUT FASRQ:\t$readfile1 $readfile2\n";

        #########    trimmomatic    #########
        print "NAME: $name\n";
        
        if ($execinfo->{trimo} eq "T") {
            print LOG2 "\n## TRIMO INPUT:\t$readfile1 $readfile2\n";

            ($outfiles)= &run_trimmomatic($name,$readfile1,$readfile2);
            ($readfile1,$readfile2)= @$outfiles;
            print "FIN trimmomatic $name \n";
        }
#exit;
        
        #########    bowtie2    #########
        if ($execinfo->{bowtie2} eq "T") {
            print LOG2 "\n## RM_CONTAMI INPUT:\t$readfile1 $readfile2\n";
            ($outfiles)= &run_rm_contami($name,$readfile1,$readfile2,$phiX_index,$execinfo->{contami},$execinfo->{idx});
            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN bowtie2 phiX & CONTAMI $name \n";
        }
        
        #########    cutadaptor    #########
        
        if ($execinfo->{bowtie2} eq "T") {
            
            if ($execinfo->{lib} eq "FLDS") {
                print LOG2 "\n## CUTADAPT FLDS INPUT:\t$readfile1 $readfile2\n";
                ($outfiles)= &run_cutadapt_FLDS($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            } elsif ($execinfo->{lib} eq "oligo") {
                print LOG2 "\n## CUTADAPT oligo INPUT:\t$readfile1 $readfile2\n";
                ($outfiles)= &run_cutadapt_oligo($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            }elsif ($execinfo->{lib} eq "BL") {
                print LOG2 "\n## CUTADAPT BL INPUT:\t$readfile1 $readfile2\n";
                ($outfiles)= &run_cutadapt_FLDS_BL($name,$readfile1,$readfile2);
                ($readfile1,$readfile2)= @$outfiles;
            }
            
            
            print "FIN cutadaptor $name \n";
            
        }
         
        #########    DUST low complexicy    ###########
        if ($execinfo->{low} eq "T") {
            print LOG2 "\n## LOW_COMPLEXITY INPUT:\t$readfile1 $readfile2\n";
            ($outfiles)= &run_rm_low_complexity_PP($name,$readfile1,$readfile2);
            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN $name DUST32\n";
            
        }

        #########    duplication    ###########
        
        if ($execinfo->{dup} eq "T") {
            print LOG2 "\n## DUPLICATION INPUT:\t$readfile1 $readfile2\n";
     #      ($readfile1,$readfile2)=
     #       ("./${OUT_DIR}/${name}_LP_R1.fq","./${OUT_DIR}/${name}_LP_R2.fq");

#            ($outfiles)= &run_rm_duplicate_PP($name,$readfile1,$readfile2);
            
            ($outfiles)= &run_rm_duplicate($name,$readfile1,$readfile2);

            ($readfile1,$readfile2)= @$outfiles;
            
            print "FIN $name duplication\n";
        }
        
        
        #########    SORTMERNA    ###########
        
        if ($execinfo->{ribo} eq "T") {
            print LOG2 "\n## SORTMERNA INPUT:\t$readfile1 $readfile2\n";
            print "****** START SORTMERNA\n";
#            ($readfile1,$readfile2)=("./${OUT_DIR}/${name}_PP_R1.fq","./${OUT_DIR}/${name}_PP_R2.fq");
            
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
    
    $cmd_trimmomatic="java -jar ${TRIMMOMATIC}/trimmomatic-0.39.jar PE -threads 24 -phred33 -trimlog logs_trimo_${name} $readfile1 $readfile2 ./${OUT_DIR}/${name}_TP_R1.fq ./${OUT_DIR}/${name}_OP_R1.fq ./${OUT_DIR}/${name}_TP_R2.fq ./${OUT_DIR}/${name}_OP_R2.fq ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50";
    
    system ($cmd_trimmomatic);
    ## TruSeq3-PE-1_NISHI_FIN2.fa
    print LOG2 "\tCMD: $cmd_trimmomatic\n";
    
    $outfiles= ["./${OUT_DIR}/${name}_TP_R1.fq",
    "./${OUT_DIR}/${name}_TP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_contami {
    my($name,$readfile1,$readfile2,$phiX_index,$contami,$otidx) = @_;
    
    #PhiX
    
    $cmd_bowtie2="${BOWTIE2}/bowtie2 -p 4 -N 1 --un-conc ./${OUT_DIR}/${name}_nocont.fq --al-conc ./${OUT_DIR}/${name}_phi.fq -x $phiX_index -1 $readfile1 -2 $readfile2 > /dev/null 2> logs_rmphi_${name}.txt";
    
    system ($cmd_bowtie2);
    print LOG2 "\tCMD: $cmd_bowtie2\n";

    system ("mv ./${OUT_DIR}/${name}_nocont.1.fq ./${OUT_DIR}/${name}_nocont_R1.fq");
    system ("mv ./${OUT_DIR}/${name}_nocont.2.fq ./${OUT_DIR}/${name}_nocont_R2.fq");
    
    $outfiles= ["./${OUT_DIR}/${name}_nocont_R1.fq",
    "./${OUT_DIR}/${name}_nocont_R2.fq"];
    
    
    if ($contami eq "T") {
        $contami_ref_index=$otidx;
        #MoCV1
        $cmd_contami="${BOWTIE2}/bowtie2 -p 4 -N 1 --un-conc ./${OUT_DIR}/${name}_nocont2.fq --al-conc ./${OUT_DIR}/${name}_contami2.fq -x $contami_ref_index -1 ./${OUT_DIR}/${name}_nocont.1.fq -2 ./${OUT_DIR}/${name}_nocont.2.fq > /dev/null 2> logs_rmmov_${name}.txt";
        
        print LOG2 "\tCMD: $cmd_contami\n";
        system ($cmd_contami);
        
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
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 50 -n 5 -O 10 \
    -g $UPM -a $U2 -G $U2_comp -A $UPM_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 50 -n 5 -O 10 \
    -g $U2_comp -a $UPM_comp -G $UPM -A $U2 \
    -o out_cut2_R1.fq -p out_cut2_R2.fq out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 50 -n 5 -O 10 \
    -g $UPM_short -A $UPM_short_comp \
    -o out_cut3_R1.fq -p out_cut3_R2.fq out_cut2_R1.fq out_cut2_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD3: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 50 -n 5 -O 10 \
    -a $UPM_short_comp -G $UPM_short \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq out_cut3_R1.fq out_cut3_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD4: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq","out_cut2_R1.fq",
    "out_cut2_R2.fq","out_cut3_R1.fq","out_cut3_R2.fq",);
    
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
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $UPM_ALL -a $BL1_comp -G $BL1 -A $UPM_ALL_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $BL1 -a $UPM_ALL_comp -G $UPM_ALL -A $BL1_comp \
    -o out_cut2_R1.fq -p out_cut2_R2.fq out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $UPM_short -A $UPM_short_comp \
    -o out_cut3_R1.fq -p out_cut3_R2.fq out_cut2_R1.fq out_cut2_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD3: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -a $UPM_short_comp -G $UPM_short \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq out_cut3_R1.fq out_cut3_R2.fq >> logs_cutadapt_${name}.txt";
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD4: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq","out_cut2_R1.fq",
    "out_cut2_R2.fq","out_cut3_R1.fq","out_cut3_R2.fq",);
    
    print "FIN cutadaptor $name \n";
    
    $outfiles= ["./${OUT_DIR}/${name}_CP_R1.fq","./${OUT_DIR}/${name}_CP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_cutadapt_oligo {
    my($name,$readfile1,$readfile2) = @_;
    
        $SMARTer_oligo_II= "AAGCAGTGGTATCAACGCAGAGTACATGGG";
    $SMARTer_oligo_II_comp= "CCCATGTACTCTGCGTTGATACCACTGCTT";
    $SMARTer_oligo_II_modify= "AAGCAGTGGTATCAACGCAGAGTACT";
    $SMARTer_oligo_II_modify_comp= "AGTACTCTGCGTTGATACCACTGCTT";
    
    ##### plus #####
     
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 50 -n 5 -O 10 \
    -g $SMARTer_oligo_II -a $SMARTer_oligo_II_modify_comp \
    -G $SMARTer_oligo_II_modify -A $SMARTer_oligo_II_comp \
    -o out_cut1_R1.fq -p out_cut1_R2.fq \
    $readfile1 $readfile2 > logs_cutadapt_${name}.txt";
    
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD1: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    ##### minus #####
    
    $cmd_cutadapt="${CUTADAPT}/cutadapt -j 8 -e 0.1 --minimum-length 20 -n 5 -O 10 \
    -g $SMARTer_oligo_II_modify -a $SMARTer_oligo_II_comp \
    -G $SMARTer_oligo_II -A $SMARTer_oligo_II_modify_comp \
    -o ./${OUT_DIR}/${name}_CP_R1.fq -p ./${OUT_DIR}/${name}_CP_R2.fq \
    out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_${name}.txt";
    
    $cmd_cutadapt=~ tr/\n//d;
    print LOG2 "\tCMD2: $cmd_cutadapt\n";
    system($cmd_cutadapt);
    
    unlink ("out_cut1_R1.fq","out_cut1_R2.fq");
    
    print "FIN cutadaptor $name \n";
    
    $outfiles= ["./${OUT_DIR}/${name}_CP_R1.fq","./${OUT_DIR}/${name}_CP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_low_complexity_PP {
    my ($name,$readfile1,$readfile2) = @_;
    my ($cmd_prinseq_low);
    
    use File::Spec;
    use Cwd;
    
#    $cwd=Cwd::getcwd();
        
    $cmd_prinseq_low="$PRINSEQ_PP/prinseq++ -threads 8 -VERBOSE 2 \
    -fastq $readfile1 -fastq2 $readfile2 \
    -min_len 50 -lc_dust=0.5 \
    -out_format 0 \
    -out_good ./${OUT_DIR}/${name}_LP_R1.fq -out_good2 ./${OUT_DIR}/${name}_LP_R2.fq >  logs_low_${name}.txt";
    
    $cmd_prinseq_low=~ tr/\n//d;
    print LOG2 "\tCMD: $cmd_prinseq_low\n";
    system ("$cmd_prinseq_low");
    
    @rmfiles= glob("*_bad_out_R*.fastq *_single_out_R[12].fastq");
    unlink (@rmfiles);

    $outfiles= ["./${OUT_DIR}/${name}_LP_R1.fq","./${OUT_DIR}/${name}_LP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_duplicate_PP {
    my ($name,$readfile1,$readfile2) = @_;
    my ($cmd_prinseq_dupli);
    
    use File::Spec;
    use Cwd;
    
#    $cwd=Cwd::getcwd();
    
    $cmd_prinseq_dupli="$PRINSEQ_PP/prinseq++ -threads 8 -VERBOSE 2 \
    -fastq $readfile1 -fastq2 $readfile2 \
    -min_len 50 -derep \
    -out_format 0 \
    -out_good ./${OUT_DIR}/${name}_PP_R1.fq -out_good2 ./${OUT_DIR}/${name}_PP_R2.fq >  logs_dup_${name}.txt";

    $cmd_prinseq_dupli=~ tr/\n//d;
    print LOG2 "\tCMD: $cmd_prinseq_dupli\n";
    system ("$cmd_prinseq_dupli");
    
    @rmfiles= glob("*_bad_out_R*.fastq *_single_out_R[12].fastq");
    unlink (@rmfiles);
    
    $outfiles= ["./${OUT_DIR}/${name}_PP_R1.fq","./${OUT_DIR}/${name}_PP_R2.fq"];
    
    return ($outfiles);
    
}

sub run_rm_duplicate {
    my($name,$readfile1,$readfile2) = @_;
    
    $cmd_prinseq_dup="${PRINSEQ}/prinseq-lite.pl -verbose -fastq $readfile1 -fastq2 $readfile2 -derep 123 -out_format 3 -out_good ./${OUT_DIR}/${name}_PP -out_bad ./${OUT_DIR}/${name}_bad -log logs_dup_${name}.txt";
    
    print LOG2 "\tCMD: $cmd_prinseq_dup\n";
    
    system ("$cmd_prinseq_dup");
    
    system ("mv ./${OUT_DIR}/${name}_PP_1.fastq ./${OUT_DIR}/${name}_PP_R1.fq");
    system ("mv ./${OUT_DIR}/${name}_PP_2.fastq ./${OUT_DIR}/${name}_PP_R2.fq");

 #   system ("rm ./${OUT_DIR}/${name}_PP_1_singletons.fastq");
 #   system ("rm ./${OUT_DIR}/${name}_PP_2_singletons.fastq");

    @rmfiles= glob("./${OUT_DIR}/*_bad_[12].fastq");
    unlink (@rmfiles);

    $outfiles= ["./${OUT_DIR}/${name}_PP_R1.fq", "./${OUT_DIR}/${name}_PP_R2.fq"];
    
    return ($outfiles);
    
    
}

sub run_sortmerna {
    my($name,$readType,$readfile1,$readfile2) = @_;
    my ($cmd_sortmerna);
    print "$readfile1 $readfile2\n";

    $cmd_sortmerna="$SMR/sortmerna \
    --ref $SMRDB/silva-bac-16s-id90.fasta \
    --ref $SMRDB/silva-bac-23s-id98.fasta \
    --ref $SMRDB/silva-arc-16s-id95.fasta \
    --ref $SMRDB/silva-arc-23s-id98.fasta \
    --ref $SMRDB/silva-euk-18s-id95.fasta \
    --ref $SMRDB/silva-euk-28s-id98.fasta \
    --ref $SMRDB/rfam-5.8s-database-id98.fasta \
    --ref $SMRDB/rfam-5s-database-id98.fasta \
    --reads $readfile1 --reads $readfile2 --workdir ${OUT_DIR} \
    --paired_in true -sam -fastx -num_alignments 1 -v \
    --aligned ./${OUT_DIR}/${name}_rrn --other ./${OUT_DIR}/${name}_norrn --out2 true \
    --threads 12";
    $cmd_sortmerna=~ s/\n//g;
    print "CMD: $cmd_sortmerna\n";
    
    system ($cmd_sortmerna);

    system ("rm -rf ./${OUT_DIR}/kvdb");

    if ($readType eq "pe") {
        system ("mv ./${OUT_DIR}/${name}_norrn_fwd.fq ./${OUT_DIR}/${name}_SP_R1.fq");
        system ("mv ./${OUT_DIR}/${name}_norrn_rev.fq ./${OUT_DIR}/${name}_SP_R2.fq");
        $outfiles= ["./${OUT_DIR}/${name}_SP_R1.fq",
        "./${OUT_DIR}/${name}_SP_R2.fq"];

        
    } else {
        $outfiles= ["./${OUT_DIR}/${name}_norrn.fq",
        "./${OUT_DIR}/${name}_rrn.fq"];
    }
    
    return ($outfiles);
    
}

sub run_sortmerna_o {
    my($name,$readType,$readfile1,$readfile2) = @_;
    my ($infile,$cmd);
    print "$readfile1 $readfile2\n";
    
    if ($readType eq "pe") {
        $infile="./${OUT_DIR}/out_PE_R12.fq";
        system ("bash ${SORTMERNA}/scripts/merge-paired-reads.sh $readfile1 $readfile2 $infile");
    } elsif ($readType eq "se") {
        $infile= $readfile1;
    }
    
    $cmd_sortmerna="${SORTMERNA}/sortmerna  --ref \
${SORTMERNA}/rRNA_databases/silva-bac-16s-id90.fasta,${SORTMERNA}/index/silva-bac-16s-db:\
${SORTMERNA}/rRNA_databases/silva-bac-23s-id98.fasta,${SORTMERNA}/index/silva-bac-23s-db:\
${SORTMERNA}/rRNA_databases/silva-arc-16s-id95.fasta,${SORTMERNA}/index/silva-arc-16s-db:\
${SORTMERNA}/rRNA_databases/silva-arc-23s-id98.fasta,${SORTMERNA}/index/silva-arc-23s-db:\
${SORTMERNA}/rRNA_databases/silva-euk-18s-id95.fasta,${SORTMERNA}/index/silva-euk-18s-db:\
${SORTMERNA}/rRNA_databases/silva-euk-28s-id98.fasta,${SORTMERNA}/index/silva-euk-28s-db:\
${SORTMERNA}/rRNA_databases/rfam-5s-database-id98.fasta,${SORTMERNA}/index/rfam-5s-db:\
${SORTMERNA}/rRNA_databases/rfam-5.8s-database-id98.fasta,${SORTMERNA}/index/rfam-5.8s-db \
 --reads $infile --sam --num_alignments 1 --fastx  --paired_in \
 --aligned ./${OUT_DIR}/${name}_rrn --other ./${OUT_DIR}/${name}_norrn \
 --log -v -a 12";
        $cmd_sortmerna=~ s/\n//g;
    
    print LOG2 "\tCMD: $cmd_sortmerna\n";

    system ($cmd_sortmerna);
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
    
    $info->{steps}= [qw(raw trim no_cont cutadapt low dupli sortme)];
    
    
    foreach $name (@{$sampinfo->{samples}}) {
        
        $onefile_path= "$sampinfo->{pe}->{$name}->[0]";
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
        $info->{sample}->{$name}->{dupli}= $count;
        
        $onefile_path= "./${OUT_DIR}/${name}_SP_R1.fq";
        $count= -f $onefile_path ? &count_reads($onefile_path) : 0;
        $info->{sample}->{$name}->{sortme}= $count;
        
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
sub create_sampinfo_path {
    my($samplist)= @_;
    
    open (LIS,"$samplist");
    $sampinfo= {};
    while(<LIS>){
        chomp;
        next if (/^(\#|\s*$)/);
        @item= split/\t/;
        $name="xxx";
        $fname= (split/\//,$item[0])[-1];

        if (@item == 2) {
            
            $name= $1 if ($fname =~ /^([^_]+)/);
            if ($fname =~ /^$name/) {
                push @{$sampinfo->{samples}},$name;
                $sampinfo->{pe}->{$name}=[$item[0],$item[1]];
                #print "@item\n";
            } else {
                $name= $1 if ($fname =~ /^(.+).(fq|fastq|fq.gz|fastq.gz)$/);
                push @{$sampinfo->{samples}},$name;
                $sampinfo->{pe}->{$name}=[$item[0],$item[1]];
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
    my ($sampinfo)= @_;
    
    my($check);
    $check= 1;
    
    foreach $name (@{$sampinfo->{samples}}) {
            
        $files = $sampinfo->{pe}->{$name};
        
        if (-f "$$files[0]") {
                
        } else {
            print "ERROR: DO NOT EXISTS $$files[0]\n";
            $check=0;
            last;
        }
        if (-f "$$files[1]") {
                
        } else {
            print "ERROR: DO NOT EXISTS $$files[1]\n";
            $check=0;
            last;
        }
    }
        
    return ($check);
}










