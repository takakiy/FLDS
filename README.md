#### FLDS Tools

Two Perl scripts for a cleanup of FLDS raw reads and a determination of both termini of virus genome/segement.

=====

#### Script 1: Cleanup_FLDS_YT200318.pl

Raw Illumina reads are sequentially processed.

=====

## Manual

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

	 * -lst sample_list1.txt
	 * -seqDir Directory path containg fastq files
	 * -lib a kind of library, default FLDS
	 * -outdir Directory path of output, default cleanup

	# Skipping a step in script
     * -trimo  (1)/0 trim adaptor sequences and low-quality sequences, default 1
     * -bowtie2  (1)/0 remove the contamination of control library, default 1
     * -cut  (1)/0 trim cDNA synthesis adaptors, default 1
     * -low  (1)/0 exclude low-complexity reads, default 1
     * -dup  (T)/F exclude PCR duplicates, default T
     * -ribo (T)/F  remove rRNA-derived reads, default T

**Output:**
 This script generates the output files for each step in the cleanup processes.




