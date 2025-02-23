## Introduction

Minipileup is a simple pileup-based variant caller. It takes a reference FASTA
and one or multiple alignment BAM as input, and outputs a multi-sample VCF along with
allele counts:
```sh
samtools faidx ref.fa       # index FASTA; bgzip'd FASTA is not supported
minipileup -yf ref.fa -p.2 aln1.bam aln2.bam > var.vcf
```
You can adjust mapping quality, base quality, alignment length and allele count
thresholds, or specify regions on the command line.

Minipileup does not use realignment or local reassembly and cannot compete with
full-pledge variant callers such as GATK on typical input. Nonetheless, you
might find it useful for counting alleles in applications GATK is not designed
for. In principle, you can parse the samtools mpileup output to generate a VCF
but minipileup is faster and more convenient.

Minipileup is adapted from the [htsbox][htsbox] pileup command which was
initially implemented in 2012 and has been a tool I frequently use to
investigate alignment data.

## Methods

Minipileup briefly follows three steps:

1. Read filtering. This step controls reads used for pileup. You can filter
   reads by mapping quality (`-q`), alignment length (`-l` and `-S`) or by
   region (`-b`).

2. Allele counting. Initial pileup is generated with the htslib
   `bam_mplp_auto()` API. At each position, alleles across all input BAMs are
   grouped and the number of supporting reads for each allele is calculated.
   Bases of low quality (`-Q`) or close to reads ends (`-T`) may be ignored.
   Minipileup drops alleles supported by too few reads (`-s`, `-a` and `-p`).

3. Output. In the variant VCF mode (`-vc` and `-f`), a site is outputted if a
   non-reference allele remains. The VCF shows the number of supporting reads
   of each allele in each BAM. Strand information can be optionally obtained
   (`-C`). If the minimum allele fraction (`-p`) is specified, the minor allele
   fraction in a heterozygous genotype must be above the threshold.

[htsbox]: https://github.com/lh3/htsbox
