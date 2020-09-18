#### FLDS Tools

Two Perl scripts for a cleanup of FLDS raw reads and a determination of both termini of virus genome/segement.

=====

#### Script 1: Cleanup_FLDS_YT200318.pl

Raw Illumina reads are sequentially processed.

=====

## Manual
=====

**Dependencies:**

	 * Trimmomatic ver. 0.39
	 * Bowtie 2 ver. 2.3.5.1
	 * Cutadapt ver. 2.4
	 * PRINSEQ ver. 0.20.4
	 * SortMeRNA ver. 2.1b

**Runbook:**

Input: Paired-end FLDS reads in FASTQ format

Output: Clean FLDS reads in FASTQ format


Step 1: Make "sample_list.txt" file.

```
rm samplelist.txt; \
 for i in ./fastq/*_S*_L001_R1_001.fastq.gz; do j=${i##./*/}; \
   k=${j%%_S*_L001_R1_001.fastq.gz}; echo "$k"; \
	echo -e "$j\t${j/_R1/_R2}\t$k" >> samplelist.txt; \
  done
```

Step 2: Run Cleanup_FLDS_YT200318.pl.


```
Cleanup_FLDS_YT200318.pl -lst samplelist.txt -seqDir fastq -lib FLDS -outdir cleanup
```

**Optional parameters:**
```
    * -lst sample_list1.txt
    * -seqDir Directory path containg fastq files
    * -lib a kind of library, default FLDS
    * -outdir Directory path of output, default cleanup

    # Skipping a step in script
    * -trimo  1/0(execute/skip) trim adaptor sequences and low-quality sequences, default 1
    * -bowtie2  1/0(execute/skip) remove the contamination of control library, default 1
    * -cut  1/0(execute/skip) trim cDNA synthesis adaptors, default 1
    * -low  1/0(execute/skip) exclude low-complexity reads, default 1
    * -dup  T/F(execute/skip) exclude PCR duplicates, default T
    * -ribo T/F(execute/skip)  remove rRNA-derived reads, default T
```

**Output:**
 This script generates the output files for each step in the cleanup processes.

```
* ./cleanup/"Sample name"_PP_R1.fq ./cleanup/"Sample name"_PP_R2.fq  ## after excluding PCR duplicates
* ./cleanup/"Sample name"_SP_R1.fq ./cleanup/"Sample name"_PP_R2.fq  ## after removing rRNA-derived reads
```

#### Script 2: TermCount_FLDS.pl

Terminal ends of viral genome/segments are determined in the candidate sequence for virus genome/segment.

## Manual
=====

**Dependencies:**

     * Bowtie 2 ver. 2.3.5.1
     * Cutadapt ver. 2.4
     * samtools ver. 1.9
     * count_term_mapping_sam_v2.pl

**Runbook:**

Input: ref.fna (Fasta file containing the candidate sequence for virus genome/segment)
Input: read_R1.fq read_R2.fq (Paired-end FLDS reads in FASTQ format after the trimmomatic's treatment)

Output: Clean FLDS reads in FASTQ format


Step 1: Run TermCount_FLDS.pl

```
TermCount_FLDS.pl 
   -ref ref.fna \
   -r1 read_R1 -r2 read_R2.fq -name out_header -tmpdir tmp \
```

**Optional parameters:**
```
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



