# FLDS Tools

Two Perl pipline scripts for a cleanup of FLDS raw reads and a determination of both termini of virus genome/segement.

## Script 1: Cleanup_FLDS.pl

Raw Illumina reads are sequentially processed.

### Manual

**Dependencies:**

	 * Trimmomatic ver. 0.39
	 * Bowtie 2 ver. 2.3.5.1
	 * Cutadapt ver. 4.1
	 * PRINSEQ ver. 0.20.4
     * PRINSEQ++ ver. 1.2
	 * SortMeRNA ver. 4.3.4

**Config:**

You need to correct for the pipline script "Cleanup_FLDS.pl" about paths of the above programs in your environment.

**Execute:**

Input: Paired-end FLDS reads (FASTQ format) in the "fastq" directory

Output: Clean FLDS reads (FASTQ format) in the "cleanup" directory


Step 1: Make "samplelist.txt" file.

```
rm samplelist.txt; \
 for i in ./fastq/*_S*_L001_R1_001.fastq.gz; do j=${i##./*/}; \
   k=${j%%_S*_L001_R1_001.fastq.gz}; echo "$k"; \
	echo -e "$j\t${j/_R1/_R2}\t$k" >> samplelist.txt; \
  done
```

Step 2: Run Cleanup_FLDS.pl.


```
Cleanup_FLDS.pl -lst samplelist.txt -seqDir fastq -lib FLDS -outdir cleanup
```

**Optional parameters:**
```
    * -lst samplelist.txt
    * -seqDir Directory path containg fastq files
    * -lib a kind of library, default FLDS
    * -outdir Directory path of output, default cleanup

    # Skipping a step in script
    * -trimo  T/F(execute/skip) trim adaptor sequences and low-quality sequences, default T
    * -bowtie2  T/F(execute/skip) remove the contamination of control library, default T
    * -cut  T/F(execute/skip) trim cDNA synthesis adaptors, default T
    * -low  T/F(execute/skip) exclude low-complexity reads, default T
    * -dup  T/F(execute/skip) exclude PCR duplicates, default T
    * -ribo T/F(execute/skip)  remove rRNA-derived reads, default T
```

**Output:**
 This pipline script generates the FASTQ files for each step in the cleanup processes.

```
* ./cleanup/"Sample name"_PP_R1.fq ./cleanup/"Sample name"_PP_R2.fq  ## after excluding PCR duplicates
* ./cleanup/"Sample name"_SP_R1.fq ./cleanup/"Sample name"_SP_R2.fq  ## after removing rRNA-derived reads
```



## Script 2: TermCount_FLDS.pl

Terminal ends of viral genome/segments are determined in the candidate sequence for virus genome/segment.

### Manual

**Dependencies:**

     * Bowtie2 ver. 2.3.5.1
     * Cutadapt ver. 4.1
     * samtools ver. 1.15.0
     * count_readTermSpecific_sam_PE.pl

**Config**
  You need to correct for the pipline script "TermCount_FLDS.pl" about paths of the above programs in your environment.

**Execute:**

Step 1: Run TermCount_FLDS.pl

Reference sequence: ref.fna (Fasta file containing the candidate sequence for virus genome/segment)

Paired-end reads: read_R1.fq read_R2.fq (Paired-end FLDS reads in FASTQ format after the trimmomatic's treatment)


```
TermCount_FLDS.pl 
   -ref ref.fna \
   -r1 read_R1.fq -r2 read_R2.fq -name "header of output" -tmpdir tmp -mname "header of sam/bam file"\
```

**Optional parameters:**
```
	 * -name header name used output files
     * -mname header name for the mapped file (sam/bam)
	 * -tmpdir temporary directory

	# Skipping a step in script
     * -cutadapt  1/0(execute/skip) trim cDNA synthesis adaptors, default 1
     * -mapu2  1/0(execute/skip) mapping using Bowtie2, default 1
```

**Outputs:**
```
     * out2_(out_header)_countTerm_summary.txt
     * out2_(out_header)_countTerm_list.txt
     * out2_(out_header)_countTerm_outlier.txt
```

Step 2: Graphic wirh R

```
 library(ggplot2)
 library(reshape2)

 dat<-read.table("out2_(out_header)_countTerm_list.txt",header=T,sep="\t")

## Select a reference sequence ("Target")
  target<-"Target"
  datS<-dat[dat[,1] == target,]

  cx<-c(0,datS[1,2])
  clen<-data.frame(cx,cy=c(0,0),type=c(3,3))

  g<-ggplot(data=datS, aes(x=loc, group=type, color=type)) +
     geom_bar(data=datS, aes(y=ifelse(..group.. == 2, datS$count, -1 * datS$count)), stat="identity",width=1) + ylab("count")   
　g<-g + geom_line(data=clen,aes(x=cx,y=cy),linetype="solid",color="black",size=0.7)

  plot(g)

```



