#!/usr/bin/perl -w


## version v5

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print " -ref reference -r1 read1 -r2 read2 -name xxx -mname xx -tmpdir yyy\n";

    print "OPTION:
      -name   name of producted fastq
      -mname   name of mapped bam file
      -cutadapt   0/(1)  1: execute cutadapt
      -mapupm 0/(1)  1: mapping UPM reads
      -mapu2 0/(1)  1: mapping U2 reads
    
      -bwtopt ' xx yy' Bowtie2 option:  --end-to-end --minins 100 --maxins 1000
    
    \n";

    exit;
}


use FindBin;
use File::Basename;

sub get_full_path {
    return $FindBin::Bin . '/' .  $FindBin::Script;
}

$scrpath= &get_full_path();
$scrdir= dirname $scrpath;

## print "START PATH: $scrpath $scrdir\n";

## exit;




foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};
$|=1;

$HOME=$ENV{"HOME"};

## 0)  SETTING BIN
	$BOWTIE2="{your path}/bowtie2-2.3.5.1";
    $CUTADAPT="{your path}/miniconda3/bin";
	$SAMTOOLS="{your path}/samtools-1.15.1";
    $COUNT_PIPE="${scrdir}/count_readTermSpecific_sam_PE.pl";
 
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
$bwtopt= exists $option->{-bwtopt} ? $option->{-bwtopt}->[0] : '';


#print "$bwtopt\n"; exit;

#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
    
$exe_cutadapt= exists $option->{-cutadapt} ?  $option->{-cutadapt}->[0] : 1;
$exe_mapping_upm= exists $option->{-mapupm} ?  $option->{-mapupm}->[0] : 1;
$exe_mapping_u2= exists $option->{-mapu2} ?  $option->{-mapu2}->[0] : 1;

if ($exe_cutadapt == 1) {
        ($ok)= &run_cutadapt_term($name,$read1,$read2,$tmpdir);
    
        print "FIN cutadapt $name $ok\n";
}

$exe_select=1;
    ## UPM
if ($exe_mapping_upm == 1) {
       ($map_upm_c111,$map_upm_c121)= &mapping_bowtie2_UPM($name,$ref,$tmpdir,$mname,$bwtopt)
} else {
       $map_upm_c111= "${tmpdir}/${mname}_S111.bam";
       $map_upm_c121= "${tmpdir}/${mname}_S121.bam";
}

#exit;

if ($exe_select == 1) {
        ($map_upm_pos,$map_upm_neg)= &mapping_bam_select_UPM($name,$ref,$tmpdir,$mname,$map_upm_c111,$map_upm_c121);
} else {
        $map_upm_pos= "${tmpdir}/${mname}_upm.pos.sort.sam";
        $map_upm_neg= "${tmpdir}/${mname}_upm.neg.sort.sam";
}

print "FIN mapping_bam_select_UPM $mname \n";
        
    ## U2
if ($exe_mapping_u2 == 1) {
     ($map_u2_c211,$map_u2_c221)= &mapping_bowtie2_U2($name,$ref,$tmpdir,$mname,$bwtopt)
        
} else {
        $map_u2_c211= "${tmpdir}/${mname}_S211.bam";
        $map_u2_c221= "${tmpdir}/${mname}_S221.bam";
}

if ($exe_select == 1) {
        ($map_u2_pos,$map_u2_neg)= &mapping_bam_select_U2($name,$ref,$tmpdir,$mname,$map_u2_c211,$map_u2_c221);
        
} else {
        $map_u2_pos= "${tmpdir}/${mname}_u2.pos.sort.sam";
        $map_u2_neg= "${tmpdir}/${mname}_u2.neg.sort.sam";
}

print "FIN mapping_bam_select_U2 $mname \n";


## EVALUATION OF STATICS

## UPM
$cmd_count_outlier= "$COUNT_PIPE -f $map_upm_pos -r $map_upm_neg -adapt upm";
system("$cmd_count_outlier");
    
# RENAME
system("mv out_countTerm_summary.txt out_${mname}_countTerm_summary.txt");
system("mv out_countTerm_list.txt out_${mname}_countTerm_list.txt");
system("mv out_countTerm_outlier.txt out_${mname}_countTerm_outlier.txt");
    
# U2
$cmd_count_outlier= "$COUNT_PIPE -f $map_u2_pos -r $map_u2_neg -adapt u2";
system("$cmd_count_outlier");

system("tail -n +2 out_countTerm_summary.txt >> out_${mname}_countTerm_summary.txt");
system("tail -n +2 out_countTerm_list.txt >> out_${mname}_countTerm_list.txt");
system("tail -n +2 out_countTerm_outlier.txt >> out_${mname}_countTerm_outlier.txt");

unlink ("out_countTerm_summary.txt","out_countTerm_list.txt","out_countTerm_outlier.txt");

    
###########   logs
@flists=<log_*>;
system("mv log_* ./logs/") if ( @flists > 0 );

print "OUTPUT UPM/U2:
    out_${mname}_countTerm_summary.txt
    out_${mname}_countTerm_list.txt
    out_${mname}_countTerm_outlier.txt\n";


#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
#==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====

sub run_cutadapt_term {
    my ($name,$read1,$read2,$tmpdir)= @_;
    my ($cut_options);
    
    $UPM="CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG";
    $UPM_comp="CCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTATTAG";
    $UPM_short="CTAATACGACTCACTATAGGGC";
    $UPM_short_comp="GCCCTATAGTGAGTCGTATTAG";
    $U2="GACGTAAGAACGTCGCACCA";
    $U2_comp="TGGTGCGACGTTCTTACGTC";
    
    $cut_options="-j 8 -e 0.1 --minimum-length 50 -n 5 -O 10";

    ###  STRAND1    UPM ====> (U2)  [UPM ---> U2 / UPM ---> ||  U2-C ---> UPM-C
    $cmd_cut110= "$CUTADAPT/cutadapt $cut_options --discard-untrimmed --pair-filter both \
       -g $UPM \
       -o ${tmpdir}/${name}_cut110_R1.fq -p ${tmpdir}/${name}_cut110_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_110.txt";

    ###  STRAND2    (U2_C) ====> UPM-C
    $cmd_cut120= "$CUTADAPT/cutadapt $cut_options --discard-untrimmed --pair-filter both \
       -G $UPM \
       -o ${tmpdir}/${name}_cut120_R1.fq -p ${tmpdir}/${name}_cut120_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_120.txt";

    ###  STRAND1    (UPM) ====> U2

    $cmd_cut210= "$CUTADAPT/cutadapt $cut_options --discard-untrimmed --pair-filter both \
       -G $U2_comp \
       -o ${tmpdir}/${name}_cut210_R1.fq -p ${tmpdir}/${name}_cut210_R2.fq\
       $read1 $read2 > log_cutadapt_${name}_210.txt";

    ###  STRAND2    U2-C ===> (UPM-C)
    $cmd_cut220= "$CUTADAPT/cutadapt $cut_options --discard-untrimmed --pair-filter both \
       -g $U2_comp \
       -o ${tmpdir}/${name}_cut220_R1.fq -p ${tmpdir}/${name}_cut220_R2.fq \
       $read1 $read2 > log_cutadapt_${name}_220.txt";
    
    
    $cmd_cut110=~ tr/\n//d;
    print "CMD: $cmd_cut110\n";
    system("$cmd_cut110");
    $cmd_cut120=~ tr/\n//d;
    system("$cmd_cut120");
    $cmd_cut210=~ tr/\n//d;
    system("$cmd_cut210");
    $cmd_cut220=~ tr/\n//d;
    system("$cmd_cut220");

    foreach $group (110,120,210,220) {
        
        $in1st_1f="${tmpdir}/${name}_cut${group}_R1.fq";
        $in1st_2f="${tmpdir}/${name}_cut${group}_R2.fq";
        $out2nd_1f=$in1st_1f;
        $out2nd_2f=$in1st_2f;
        $out2nd_1f=~ s/0_R1.fq/1_R1.fq/;
        $out2nd_2f=~ s/0_R2.fq/1_R2.fq/;

        $cmd_cutadapt1="${CUTADAPT}/cutadapt $cut_options \
        -g $UPM -a $U2 -G $U2_comp -A $UPM_comp \
        -o out_cut1_R1.fq -p out_cut1_R2.fq $in1st_1f $in1st_2f > logs_cutadapt_1st.txt";
        $cmd_cutadapt1=~ tr/\n//d;
        system($cmd_cutadapt1);
        
        $cmd_cutadapt2="${CUTADAPT}/cutadapt $cut_options \
        -g $U2_comp -a $UPM_comp -G $UPM -A $U2 \
        -o out_cut2_R1.fq -p out_cut2_R2.fq out_cut1_R1.fq out_cut1_R2.fq >> logs_cutadapt_1st.txt";
        $cmd_cutadapt2=~ tr/\n//d;
        system($cmd_cutadapt2);
        
        $cmd_cutadapt3="${CUTADAPT}/cutadapt $cut_options \
        -g $UPM_short -A $UPM_short_comp \
        -o out_cut3_R1.fq -p out_cut3_R2.fq out_cut2_R1.fq out_cut2_R2.fq >> logs_cutadapt_1st.txt";
        $cmd_cutadapt3=~ tr/\n//d;
        system($cmd_cutadapt3);
        
        $cmd_cutadapt4="${CUTADAPT}/cutadapt $cut_options \
        -a $UPM_short_comp -G $UPM_short \
        -o $out2nd_1f -p $out2nd_2f out_cut3_R1.fq out_cut3_R2.fq >> logs_cutadapt_1st.txt";
        $cmd_cutadapt4=~ tr/\n//d;
        system($cmd_cutadapt4);
        
        unlink("out_cut1_R1.fq","out_cut1_R2.fq");
        unlink("out_cut2_R1.fq","out_cut2_R2.fq");
        unlink("out_cut3_R1.fq","out_cut3_R2.fq");

    }
    
    
    return (1);
    
    
}

sub mapping_bowtie2_UPM {
    my ($name,$ref,$tmpdir,$mname,$bwtopt)= @_;
    
    print "##### START MAPPING UPM\n";
    
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
    $BOWTIE2/bowtie2 $bwtopt --fr -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR1}.sam; \
    $SAMTOOLS/samtools view -@ 12 -S -b -o ${tmpdir}/${nameR1}.bam ${tmpdir}/${nameR1}.sam";
    
    $cmd_mapping_select_S111=~ tr/\n//d;
    system("$cmd_mapping_select_S111");
    
    #**  Mapping with R2(UPM) PE & SELECT READS ** ** ** **
    $read1="${tmpdir}/${name}_cut121_R1.fq";
    $read2="${tmpdir}/${name}_cut121_R2.fq";
    $nameR2="${mname}_S121";

     $cmd_mapping_select_S121= "$BOWTIE2/bowtie2-build -f ${ref} ${ref}; \
     $BOWTIE2/bowtie2 $bwtopt --fr -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR2}.sam; \
     $SAMTOOLS/samtools view -@ 12 -S -b -o ${tmpdir}/${nameR2}.bam ${tmpdir}/${nameR2}.sam";

    $cmd_mapping_select_S121=~ tr/\n//d;
    system("$cmd_mapping_select_S121");

    print "$cmd_mapping_select_S121\n";
    
    #** ** ** ** ** ** ** ** ** ** ** **

    return ("${tmpdir}/${nameR1}.bam","${tmpdir}/${nameR2}.bam");
    
}


sub mapping_bam_select_UPM {
    my ($name,$ref,$tmpdir,$mname,$map_upm_c111,$map_upm_c121)= @_;
    
    print "##### START SELECT UPM\n";
    
    #**   SELECT READS with R1(UPM) ** ** ** **
    $nameR1="${mname}_S111";

    $cmd_mapping_select_S111= "$SAMTOOLS/samtools view -u -h -f 1 -F 12 $map_upm_c111 > ${tmpdir}/${nameR1}_mm.bam; \
    $SAMTOOLS/samtools view -b -f 64 -F 272 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.1st.pos.bam; \
    $SAMTOOLS/samtools view -b -f 80 -F 256 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.1st.neg.bam";
    
    $cmd_mapping_select_S111=~ tr/\n//d;
    system("$cmd_mapping_select_S111");
    
    #**   SELECT READS with R2(UPM) ** ** ** **

    $nameR2="${mname}_S121";

    $cmd_mapping_select_S121= "$SAMTOOLS/samtools view -u -h -f 1 -F 12 $map_upm_c121 > ${tmpdir}/${nameR2}_mm.bam; \
    $SAMTOOLS/samtools view -b -f 128 -F 272 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.2nd.pos.bam; \
    $SAMTOOLS/samtools view -b -f 144 -F 256 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.2nd.neg.bam";

    $cmd_mapping_select_S121=~ tr/\n//d;
    system("$cmd_mapping_select_S121");

    #**  COLLECT READS ON THE PLUS / MINUS STRAND
     $cmd_merge_pos_neg="$SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_upm.pos.bam ${tmpdir}/${nameR1}_mm.1st.pos.bam ${tmpdir}/${nameR2}_mm.2nd.pos.bam; \
    $SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_upm.neg.bam ${tmpdir}/${nameR1}_mm.1st.neg.bam ${tmpdir}/${nameR2}_mm.2nd.neg.bam";

    #**  SORT
     $cmd_sort_pos_neg= "$SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_upm.pos.sort.bam ${tmpdir}/${mname}_upm.pos.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_upm.pos.sort.bam; \
    $SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_upm.neg.sort.bam ${tmpdir}/${mname}_upm.neg.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_upm.neg.sort.bam";

    #**  BAM TO SAM
     $cmd_convert_bam2sam= "$SAMTOOLS/samtools view -h ${tmpdir}/${mname}_upm.pos.sort.bam > ${tmpdir}/${mname}_upm.pos.sort.sam; \
    $SAMTOOLS/samtools view -h ${tmpdir}/${mname}_upm.neg.sort.bam > ${tmpdir}/${mname}_upm.neg.sort.sam";

    $cmd_merge_pos_neg=~ tr/\n//d;
    system("$cmd_merge_pos_neg");
    $cmd_sort_pos_neg=~ tr/\n//d;
    system("$cmd_sort_pos_neg");
    $cmd_convert_bam2sam=~ tr/\n//d;
    system("$cmd_convert_bam2sam");

    #** ** ** ** ** ** ** ** ** ** ** **

    return ("${tmpdir}/${mname}_upm.pos.sort.sam","${tmpdir}/${mname}_upm.neg.sort.sam");
    
}

sub mapping_bowtie2_U2 {
    my ($name,$ref,$tmpdir,$mname,$bwtopt)= @_;
    
    print "##### START MAPPING U2\n";

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
    $BOWTIE2/bowtie2 $bwtopt --fr -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR1}.sam; \
    $SAMTOOLS/samtools view -@ 12 -S -b -o ${tmpdir}/${nameR1}.bam ${tmpdir}/${nameR1}.sam";
    
    $cmd_mapping_select_S211=~ tr/\n//d;
    system("$cmd_mapping_select_S211");
    
    #**  Mapping with R2(UPM) PE & SELECT READS ** ** ** **
    $read1="${tmpdir}/${name}_cut221_R1.fq";
    $read2="${tmpdir}/${name}_cut221_R2.fq";
    $nameR2="${mname}_S221";

    $cmd_mapping_select_S221= "$BOWTIE2/bowtie2-build -f ${ref} ${ref}; \
    $BOWTIE2/bowtie2 $bwtopt --fr -q -p 8 -x ${ref} -1 $read1 -2 $read2 -S ${tmpdir}/${nameR2}.sam; \
    $SAMTOOLS/samtools view -@ 12 -S -b -o ${tmpdir}/${nameR2}.bam ${tmpdir}/${nameR2}.sam";

    $cmd_mapping_select_S221=~ tr/\n//d;
    system("$cmd_mapping_select_S221");

    print "$cmd_mapping_select_S221\n";
    
    return ("${tmpdir}/${nameR1}.bam","${tmpdir}/${nameR2}.bam");

}



sub mapping_bam_select_U2 {
    my ($name,$ref,$tmpdir,$mname,$map_u2_c211,$map_u2_c221)= @_;
    
    print "##### START SELECT U2\n";

     ### SELECT R1 mapped, R2 mapped
     ## SENCOND IN PAIR  211
     ##    SELECT READ MAPPED ON NEGATIVE  -->   POS
     ##    SELECT READ MAPPED ON POSITIVE  -->   NEG
     ## FIRST IN PAIR   221
     ##    SELECT READ MAPPED ON NEGATIVE  -->   POS
     ##    SELECT READ MAPPED ON POSITIVE  -->   NEG

    #**   SELECT READS with R2(U2) ** ** ** **
    $nameR1="${mname}_S211";

    $cmd_mapping_select_S211= "$SAMTOOLS/samtools view -u -h -f 1 -F 12 $map_u2_c211 > ${tmpdir}/${nameR1}_mm.bam; \
    $SAMTOOLS/samtools view -b -f 128 -F 272 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.2nd.neg.bam; \
    $SAMTOOLS/samtools view -b -f 144 -F 256 ${tmpdir}/${nameR1}_mm.bam > ${tmpdir}/${nameR1}_mm.2nd.pos.bam";
    
    $cmd_mapping_select_S211=~ tr/\n//d;
    system("$cmd_mapping_select_S211");
    
    #**   SELECT READS with R1(U2) ** ** ** **
    $nameR2="${mname}_S221";

     $cmd_mapping_select_S221= "$SAMTOOLS/samtools view -u -h -f 1 -F 12 $map_u2_c221 > ${tmpdir}/${nameR2}_mm.bam; \
     $SAMTOOLS/samtools view -b -f 64 -F 272 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.1st.neg.bam; \
     $SAMTOOLS/samtools view -b -f 80 -F 256 ${tmpdir}/${nameR2}_mm.bam > ${tmpdir}/${nameR2}_mm.1st.pos.bam";

    $cmd_mapping_select_S221=~ tr/\n//d;
    system("$cmd_mapping_select_S221");

    #**  COLLECT READS ON THE PLUS / MINUS STRAND
    $cmd_merge_pos_neg="$SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_u2.pos.bam ${tmpdir}/${nameR1}_mm.2nd.pos.bam ${tmpdir}/${nameR2}_mm.1st.pos.bam; \
     $SAMTOOLS/samtools merge -f ${tmpdir}/${mname}_u2.neg.bam ${tmpdir}/${nameR1}_mm.2nd.neg.bam ${tmpdir}/${nameR2}_mm.1st.neg.bam";

    #**  SORT
     $cmd_sort_pos_neg= "$SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_u2.pos.sort.bam ${tmpdir}/${mname}_u2.pos.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_u2.pos.sort.bam; \
    $SAMTOOLS/samtools sort -@ 12 -o ${tmpdir}/${mname}_u2.neg.sort.bam ${tmpdir}/${mname}_u2.neg.bam; \
    $SAMTOOLS/samtools index ${tmpdir}/${mname}_u2.neg.sort.bam";

    #**  BAM TO SAM
     $cmd_convert_bam2sam= "$SAMTOOLS/samtools view -h ${tmpdir}/${mname}_u2.pos.sort.bam > ${tmpdir}/${mname}_u2.pos.sort.sam; \
    $SAMTOOLS/samtools view -h ${tmpdir}/${mname}_u2.neg.sort.bam > ${tmpdir}/${mname}_u2.neg.sort.sam";
    
    
    $cmd_merge_pos_neg=~ tr/\n//d;
    system("$cmd_merge_pos_neg");
    $cmd_sort_pos_neg=~ tr/\n//d;
    system("$cmd_sort_pos_neg");
    $cmd_convert_bam2sam=~ tr/\n//d;
    system("$cmd_convert_bam2sam");

    #** ** ** ** ** ** ** ** ** ** ** **

    return ("${tmpdir}/${mname}_u2.pos.sort.sam","${tmpdir}/${mname}_u2.neg.sort.sam");

}





