Minipileup is a simple pileup-based variant caller. It takes a reference FASTA
and one or multiple alignment BAM as input, and outputs a multi-sample VCF along with
allele counts:
```sh
samtools faidx ref.fa       # index FASTA; bgzip'd FASTA is not supported
minipileup -yf ref.fa aln1.bam aln2.bam > var.vcf
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

[htsbox]: https://github.com/lh3/htsbox
