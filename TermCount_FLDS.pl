#!/usr/bin/perl -w

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print " -ref reference -r1 read1 -r2 read2 -name xxx -mname xx -tmpdir yyy\n";

    print "OPTION:
      -cutadapt   0/(1)  1: execute cutadapt
      -mapupm 0/(1)  1: mapping UPM reads
      -mapu2 0/(1)  1: mapping U2 reads
    \n";

    exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};
$|=1;

$HOME=$ENV{"HOME"};

## 0)  SETTING BIN
	$BOWTIE2="${HOME}/biotools/local/mapping/bowtie2-2.3.5.1";
    $CUTADAPT="${HOME}/biotools/rhel6/miniconda3/bin";
	$SAMTOOLS="${HOME}/biotools/local/mapping/samtools-1.9";
    $COUNT_PIPE="${HOME}/Desktop/work_CLC/FLDS/bin/count_readTermSpecific_sam_PE.pl";
 
    $tmpdir= exists $option->{-tmpdir} ? $option->{-tmpdir}->[0] : "tmp";;


## 1)  CHECK FASTQ

if (exists $option->{-r1} && exists $option->{-r2}) {
    $read1= $option->{-r1}->[0];
    $read2= $option->{-r2}->[0];
} else {
    print "ERROR: PLEASE SPECIFY PAIRED-END READ\n";
    exit;
}
if (exists $option->{-ref}) {
    $ref= $option->{-ref}->[0];
} else {
    print "ERROR: PLEASE SPECIFY REFERENCE FASTA FILE\n";
    exit;
}

if (-d "$tmpdir") { } else {system ("mkdir $tmpdir")};
if (-d "logs") { } else {system ("mkdir logs")};

$name= exists $option->{-name} ? $option->{-name}->[0] : "xx";
$mname= exists $option->{-mname} ? $option->{-mname}->[0] : $name;

#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
    
    $exe_cutadapt= exists $option->{-cutadapt} ?  $option->{-cutadapt}->[0] : 1;
    $exe_mapping_upm= exists $option->{-mapupm} ?  $option->{-mapupm}->[0] : 1;
    $exe_mapping_u2= exists $option->{-mapu2} ?  $option->{-mapu2}->[0] : 1;

if ($exe_cutadapt == 1) {
        ($ok)= &run_cutadapt_term($name,$read1,$read2,$tmpdir);
    
        print "FIN cutadapt $name \n";
    }
#exit;
    if ($exe_mapping_upm == 1) {
        ($ok)= &mapping_bam_select_UPM($name,$ref,$tmpdir,$mname);
    
        print "FIN mapping_bam_select_UPM $mname \n";
    }
#exit;
    if ($exe_mapping_u2 == 1) {
        ($ok)= &mapping_bam_select_U2($name,$ref,$tmpdir,$mname);
    
        print "FIN mapping_bam_select_U2 $mname \n";
    }
    ## UPM
    $cmd_count_outlier= "$COUNT_PIPE -f ${tmpdir}/${mname}_mm.pls1.sort.sam -r ${tmpdir}/${mname}_mm.min1.sort.sam";
    system("$cmd_count_outlier");
    
    # RENAME
    system("mv out_countTerm_summary.txt out1_${mname}_countTerm_summary.txt");
    system("mv out_countTerm_list.txt out1_${mname}_countTerm_list.txt");
    system("mv out_countTerm_outlier.txt out1_${mname}_countTerm_outlier.txt");
    
    # U2
    $cmd_count_outlier= "$COUNT_PIPE -f ${tmpdir}/${mname}_mm.pls2.sort.sam -r ${tmpdir}/${mname}_mm.min2.sort.sam";
    system("$cmd_count_outlier");

    system("mv out_countTerm_summary.txt out2_${mname}_countTerm_summary.txt");
    system("mv out_countTerm_list.txt out2_${mname}_countTerm_list.txt");
    system("mv out_countTerm_outlier.txt out2_${mname}_countTerm_outlier.txt");

    
    
    ###########   logs
    @flists=<log_*>;
    system("mv log_* ./logs/") if ( @flists > 0 );

    print "OUTPUT UPM:
            out1_${mname}_countTerm_summary.txt
            out1_${mname}_countTerm_list.txt
            out1_${mname}_countTerm_outlier.txt\n";
    print "OUTPUT U2:
            out2_${mname}_countTerm_summary.txt
            out2_${mname}_countTerm_list.txt
            out2_${mname}_countTerm_outlier.txt\n";



#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====

sub run_cutadapt_term {
   my ($name,$read1,$read2,$tmpdir)= @_;
    
    $UPM="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG";
    $UPM_comp="CCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG";
    $U2="GACGTAAGAACGTCGCACCA";
    $U2_comp="TGGTGCGACGTTCTTACGTC";

    ###  STRAND1    UPM ====> (U2)  [UPM ---> U2 / UPM ---> ||  U2-C ---> UPM-C
      $cmd_cut111= "$CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 --discard-untrimmed --pair-filter both \
       -g $UPM \
       -o ${tmpdir}/${name}_cut110_R1.fq -p ${tmpdir}/${name}_cut110_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_110.txt; \
      $CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 \
       -g $UPM -g $U2_comp -a $U2 -a $UPM_comp \
       -G $U2_comp -G $UPM -A $UPM_comp -A $U2  \
       -o ${tmpdir}/${name}_cut111_R1.fq -p ${tmpdir}/${name}_cut111_R2.fq \
       ${tmpdir}/${name}_cut110_R1.fq ${tmpdir}/${name}_cut110_R2.fq > ${tmpdir}/log_cutadapt_${name}_111.txt";

    ###  STRAND2    (U2_C) ====> UPM-C
      $cmd_cut121= "$CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 --discard-untrimmed --pair-filter both \
       -G $UPM \
       -o ${tmpdir}/${name}_cut120_R1.fq -p ${tmpdir}/${name}_cut120_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_120.txt; \
      $CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 \
       -g $UPM -g $U2_comp -a $U2 -a $UPM_comp \
       -G $U2_comp -G $UPM -A $UPM_comp -A $U2  \
       -o ${tmpdir}/${name}_cut121_R1.fq -p ${tmpdir}/${name}_cut121_R2.fq \
       ${tmpdir}/${name}_cut120_R1.fq ${tmpdir}/${name}_cut120_R2.fq > ${tmpdir}/log_cutadapt_${name}_121.txt";

    ###  STRAND1    (UPM) ====> U2

      $cmd_cut211= "$CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 --discard-untrimmed --pair-filter both \
       -G $U2_comp \
       -o ${tmpdir}/${name}_cut210_R1.fq -p ${tmpdir}/${name}_cut210_R2.fq\
       $read1 $read2 > log_cutadapt_${name}_210.txt; \
      $CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 \
       -g $UPM -g $U2_comp -a $U2 -a $UPM_comp \
       -o ${tmpdir}/${name}_cut211_R1.fq -p ${tmpdir}/${name}_cut211_R2.fq \
       ${tmpdir}/${name}_cut210_R1.fq ${tmpdir}/${name}_cut210_R2.fq > ${tmpdir}/log_cutadapt_${name}_211.txt";

    ###  STRAND2    U2-C ===> (UPM-C)
      $cmd_cut221= "$CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 --discard-untrimmed --pair-filter both \
       -g $U2_comp \
       -o ${tmpdir}/${name}_cut220_R1.fq -p ${tmpdir}/${name}_cut220_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_220.txt; \
      $CUTADAPT/cutadapt -j 8 -e 0.15 --minimum-length 20 -n 5 -O 10 \
       -g $UPM -g $U2_comp -a $U2 -a $UPM_comp \
       -o ${tmpdir}/${name}_cut221_R1.fq -p ${tmpdir}/${name}_cut221_R2.fq \
       ${tmpdir}/${name}_cut220_R1.fq ${tmpdir}/${name}_cut220_R2.fq > ${tmpdir}/log_cutadapt_${name}_221.txt";
    
    
    $cmd_cut111=~ tr/\n//d;
    print "CMD: $cmd_cut111\n";
    system("$cmd_cut111");
    $cmd_cut121=~ tr/\n//d;
    system("$cmd_cut121");
    $cmd_cut211=~ tr/\n//d;
    system("$cmd_cut211");
    $cmd_cut221=~ tr/\n//d;
    system("$cmd_cut221");

    
    return (1);
    
    
}


sub mapping_bam_select_UPM {
    my ($name,$ref,$tmpdir,$mname)= @_;
    
    
    #**  Mapping with R1(UPM) PE & SELECT READS ** ** ** **
     $read1="${tmpdir}/${name}_cut111_R1.fq";
     $read2="${tmpdir}/${name}_cut111_R2.fq";
     $nameR1="${mname}_S111";

     ### SELECT R1 mapped, R2 mapped
     ## first in pair
     ##  SELECT READ MAPPED ON THE PLUS STRAND
     ## second in pair
     ## first in pair
     ## second in pair

    $cmd_mapping_select_S111= "$BOWTIE2/bowtie2-build -f ${ref} ${ref}; \
    $BOWTIE2/bowtie2 --rf --minins 100 --maxins 1000 -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR1}.sam;\
    $SAMTOOLS/samtools view -@ 12 -S -b ${tmpdir}/${nameR1}.sam > ${tmpdir}/${nameR1}.bam;\
    $SAMTOOLS/samtools view -u -h -f 1 -F 12 ${tmpdir}/${nameR1}.bam > ${tmpdir}/${nameR1}_mm.bam;\
    $SAMTOOLS/samtools view -b -f 64 -F 16 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.fwd1.bam;\
    $SAMTOOLS/samtools view -b -f 144 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.fwd2.bam;\
    $SAMTOOLS/samtools view -b -f 80 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.rev1.bam;\
    $SAMTOOLS/samtools view -b -f 128 -F 16 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.rev2.bam";
    
    $cmd_mapping_select_S111=~ tr/\n//d;
    system("$cmd_mapping_select_S111");
    
    #**  Mapping with R2(UPM) PE & SELECT READS ** ** ** **
    $read1="${tmpdir}/${name}_cut121_R1.fq";
    $read2="${tmpdir}/${name}_cut121_R2.fq";
    $nameR2="${mname}_S121";

     $cmd_mapping_select_S121= "$BOWTIE2/bowtie2 --rf --minins 100 --maxins 1000 -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR2}.sam; \
    $SAMTOOLS/samtools view -@ 12 -S -b ${tmpdir}/${nameR2}.sam > ${tmpdir}/${nameR2}.bam; \
    $SAMTOOLS/samtools view -u -h -f 1 -F 12 ${tmpdir}/${nameR2}.bam > ${tmpdir}/${nameR2}_mm.bam; \
    $SAMTOOLS/samtools view -b -f 64 -F 16 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.fwd1.bam; \
    $SAMTOOLS/samtools view -b -f 144 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.fwd2.bam; \
    $SAMTOOLS/samtools view -b -f 80 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.rev1.bam; \
    $SAMTOOLS/samtools view -b -f 128 -F 16 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.rev2.bam";

    $cmd_mapping_select_S121=~ tr/\n//d;
    system("$cmd_mapping_select_S121");

    #**  COLLECT READS ON THE PLUS / MINUS STRAND
     $cmd_select_pls_min="$SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_mm.pls1.bam ${tmpdir}/${nameR1}_mm.fwd1.bam ${tmpdir}/${nameR2}_mm.rev2.bam; \
    $SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_mm.min1.bam ${tmpdir}/${nameR1}_mm.rev1.bam ${tmpdir}/${nameR2}_mm.fwd2.bam";

    #**  SORT
     $cmd_sort_pls_min= "$SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_mm.pls1.sort.bam ${tmpdir}/${mname}_mm.pls1.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_mm.pls1.sort.bam; \
    $SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_mm.min1.sort.bam ${tmpdir}/${mname}_mm.min1.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_mm.min1.sort.bam";

    #**  BAM TO SAM
     $cmd_convert_bam2sam= "$SAMTOOLS/samtools view -h ${tmpdir}/${mname}_mm.pls1.sort.bam > ${tmpdir}/${mname}_mm.pls1.sort.sam; \
    $SAMTOOLS/samtools view -h ${tmpdir}/${mname}_mm.min1.sort.bam > ${tmpdir}/${mname}_mm.min1.sort.sam";

    $cmd_select_pls_min=~ tr/\n//d;
    system("$cmd_select_pls_min");
    $cmd_sort_pls_min=~ tr/\n//d;
    system("$cmd_sort_pls_min");
    $cmd_convert_bam2sam=~ tr/\n//d;
    system("$cmd_convert_bam2sam");

    #** ** ** ** ** ** ** ** ** ** ** **

    return (1);
    
}


sub mapping_bam_select_U2 {
    my ($name,$ref,$tmpdir,$mname)= @_;
    
    
    #**  Mapping with R1(UPM) PE & SELECT READS ** ** ** **
     $read1="${tmpdir}/${name}_cut211_R1.fq";
     $read2="${tmpdir}/${name}_cut211_R2.fq";
     $nameR1="${mname}_S211";

     ### SELECT R1 mapped, R2 mapped
     ## first in pair
     ##  SELECT READ MAPPED ON THE PLUS STRAND
     ## second in pair
     ## first in pair
     ## second in pair

    $cmd_mapping_select_S211= "$BOWTIE2/bowtie2-build -f ${ref} ${ref}; \
    $BOWTIE2/bowtie2 --rf --minins 100 --maxins 1000 -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR1}.sam;\
    $SAMTOOLS/samtools view -@ 12 -S -b ${tmpdir}/${nameR1}.sam > ${tmpdir}/${nameR1}.bam;\
    $SAMTOOLS/samtools view -u -h -f 1 -F 12 ${tmpdir}/${nameR1}.bam > ${tmpdir}/${nameR1}_mm.bam;\
    $SAMTOOLS/samtools view -b -f 64 -F 16 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.fwd1.bam;\
    $SAMTOOLS/samtools view -b -f 144 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.fwd2.bam;\
    $SAMTOOLS/samtools view -b -f 80 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.rev1.bam;\
    $SAMTOOLS/samtools view -b -f 128 -F 16 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.rev2.bam";
    
    $cmd_mapping_select_S211=~ tr/\n//d;
    system("$cmd_mapping_select_S211");
    
    #**  Mapping with R2(UPM) PE & SELECT READS ** ** ** **
    $read1="${tmpdir}/${name}_cut221_R1.fq";
    $read2="${tmpdir}/${name}_cut221_R2.fq";
    $nameR2="${mname}_S221";

     $cmd_mapping_select_S221= "$BOWTIE2/bowtie2 --rf --minins 100 --maxins 1000 -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR2}.sam; \
    $SAMTOOLS/samtools view -@ 12 -S -b ${tmpdir}/${nameR2}.sam > ${tmpdir}/${nameR2}.bam; \
    $SAMTOOLS/samtools view -u -h -f 1 -F 12 ${tmpdir}/${nameR2}.bam > ${tmpdir}/${nameR2}_mm.bam; \
    $SAMTOOLS/samtools view -b -f 64 -F 16 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.fwd1.bam; \
    $SAMTOOLS/samtools view -b -f 144 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.fwd2.bam; \
    $SAMTOOLS/samtools view -b -f 80 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.rev1.bam; \
    $SAMTOOLS/samtools view -b -f 128 -F 16 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.rev2.bam";

    $cmd_mapping_select_S221=~ tr/\n//d;
    system("$cmd_mapping_select_S221");

    #**  COLLECT READS ON THE PLUS / MINUS STRAND
     $cmd_select_pls_min="$SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_mm.pls2.bam ${tmpdir}/${nameR1}_mm.rev2.bam ${tmpdir}/${nameR2}_mm.fwd1.bam; \
    $SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_mm.min2.bam ${tmpdir}/${nameR1}_mm.fwd2.bam ${tmpdir}/${nameR2}_mm.rev1.bam";

    #**  SORT
     $cmd_sort_pls_min= "$SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_mm.pls2.sort.bam ${tmpdir}/${mname}_mm.pls2.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_mm.pls2.sort.bam; \
    $SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_mm.min2.sort.bam ${tmpdir}/${mname}_mm.min2.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_mm.min2.sort.bam";

    #**  BAM TO SAM
     $cmd_convert_bam2sam= "$SAMTOOLS/samtools view -h ${tmpdir}/${mname}_mm.pls2.sort.bam > ${tmpdir}/${mname}_mm.pls2.sort.sam; \
    $SAMTOOLS/samtools view -h ${tmpdir}/${mname}_mm.min2.sort.bam > ${tmpdir}/${mname}_mm.min2.sort.sam";

    $cmd_select_pls_min=~ tr/\n//d;
    system("$cmd_select_pls_min");
    $cmd_sort_pls_min=~ tr/\n//d;
    system("$cmd_sort_pls_min");
    $cmd_convert_bam2sam=~ tr/\n//d;
    system("$cmd_convert_bam2sam");

    #** ** ** ** ** ** ** ** ** ** ** **

    return (1);
    
}





