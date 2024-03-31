Minipileup is a simple pileup-based variant caller. It takes a reference FASTA
and one or multiple alignment BAM, and outputs a multi-sample VCF along with
allele counts:
```sh
minipileup -yf ref.fa aln1.bam aln2.bam > var.vcf
```
You can adjust mapping quality, base quality, alignment length and allele count
thresholds, or specify regions on the command line.

Minipileup does not use realignment or local reassembly and cannot compete with
full-pledge variant callers such as GATK on typical input. Nonetheless, you
might find it useful for counting alleles in applications GATK is not designed
for. In principle, you can parse the samtools mpileup output to generate a VCF
but minipileup is faster and more convenient.

I often see questions about how to get allele counts from alignments. The
[htsbox][htsbox] pileup command has implemented this functionality since 2012
but few are using it. Minipileup is derived from htsbox. I hope as a separate
repo, minipileup may become more useful.

[htsbox]: https://github.com/lh3/htsbox
