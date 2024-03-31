// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include "sam.h"
#include "faidx.h"
#include "ksort.h"
#include "ketopt.h"

#define VERSION "0.1.0"

const char *hts_parse_reg(const char *s, int *beg, int *end);
void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

typedef struct {     // auxiliary data structure
	BGZF *fp;        // the file handler
	hts_itr_t *itr;  // NULL if a region not specified
	const bam_hdr_t *h;
	int min_mapQ, min_len; // mapQ filter; length filter
	int min_supp_len;
	void *bed;       // bedidx if not NULL
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->itr? bam_itr_next(aux->fp, aux->itr, b) : bam_read1(aux->fp, b);
	if (ret < 0) return ret;
	if (b->core.tid < 0) b->core.flag |= BAM_FUNMAP;
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) {
			b->core.flag |= BAM_FUNMAP;
		} else if (aux->min_len > 0 || aux->min_supp_len > 0 || aux->bed) {
			int k, qlen = 0, tlen = 0;
			const char *chr = aux->h->target_name[b->core.tid];
			const uint32_t *cigar = bam_get_cigar(b);
			for (k = 0; k < b->core.n_cigar; ++k) { // compute the query length in the alignment
				int op = bam_cigar_op(cigar[k]);
				int oplen = bam_cigar_oplen(cigar[k]);
				if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP)
					qlen += oplen;
				if (bam_cigar_type(op)&2)
					tlen += oplen;
			}
			if (qlen < aux->min_len) b->core.flag |= BAM_FUNMAP;
			if (qlen < aux->min_supp_len && (b->core.flag&BAM_FSUPP)) b->core.flag |= BAM_FUNMAP;
			if (aux->bed && !(b->core.flag&BAM_FUNMAP) && !bed_overlap(aux->bed, chr, b->core.pos, b->core.pos + tlen))
				b->core.flag |= BAM_FUNMAP;
		}
	}
	return ret;
}

typedef struct {
	uint32_t is_skip:1, is_rev:1, b:4, q:8, is_del:1, k:17; // b=base, q=quality, k=allele id
	int indel; // <0: deleteion; >0: insertion
	uint64_t hash;
	uint64_t pos; // i<<32|j: j-th read of the i-th sample
} allele_t;

#define allele_lt(a, b) ((a).hash < (b).hash || ((a).hash == (b).hash && (a).indel < (b).indel))
KSORT_INIT(allele, allele_t, allele_lt)

static inline allele_t pileup2allele(const bam_pileup1_t *p, int min_baseQ, uint64_t pos, int ref, int trim_len, int del_as_allele)
{ // collect allele information given a pileup1 record
	allele_t a;
	int i;
	const uint8_t *seq = bam_get_seq(p->b);
	a.k = (1<<17) - 1; // this will be set in count_alleles()
	a.q = bam_get_qual(p->b)[p->qpos];
	a.is_rev = bam_is_rev(p->b);
	if (del_as_allele) {
		a.is_skip = (p->is_refskip || a.q < min_baseQ);
		a.is_del = p->is_del;
	} else {
		a.is_skip = (p->is_del || p->is_refskip || a.q < min_baseQ);
		a.is_del = 0;
	}
	if (p->qpos < trim_len || p->b->core.l_qseq - p->qpos < trim_len) a.is_skip = 1;
	a.indel = p->indel;
	a.b = a.hash = bam_seqi(seq, p->qpos);
	a.pos = pos;
	if (a.is_del) a.hash = (uint64_t)-1;
	else if (p->indel > 0) // compute the hash for the insertion
		for (i = 0; i < p->indel; ++i)
			a.hash = (a.hash<<4) + a.hash + bam_seqi(seq, p->qpos + i + 1);
	a.hash = a.hash << 1 >> 1;
	if (p->indel != 0 || a.b != ref || ref == 15 || a.is_del) // the highest bit tells whether it is a reference allele or not
		a.hash |= 1ULL<<63;
	return a;
}

static inline void print_allele(const bam_pileup1_t *p, int l_ref, const char *ref, int pos, int max_del, int is_vcf, int del_as_allele)
{ // print the allele. The format depends on is_vcf.
	const uint8_t *seq = bam_get_seq(p->b);
	int i, rest = max_del;
	if (del_as_allele && p->is_del) {
		putchar('*');
		return;
	}
	putchar(seq_nt16_str[bam_seqi(seq, p->qpos)]);
	if (p->indel > 0) {
		if (!is_vcf) printf("+%d", p->indel);
		for (i = 1; i <= p->indel; ++i)
			putchar(seq_nt16_str[bam_seqi(seq, p->qpos + i)]);
	} else if (p->indel < 0) {
		if (!is_vcf) {
			printf("%d", p->indel);
			for (i = 1; i <= -p->indel; ++i)
				putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
		} else rest -= -p->indel, pos += -p->indel;
	}
	if (is_vcf)
		for (i = 1; i <= rest; ++i)
			putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
}

typedef struct {
	int n_a, n_alleles, max_del; // n_a: #reads used to compute quality sum; max_del: max deletion length
	int tot_dp, max_dp, n_cnt, max_cnt;
	allele_t *a; // allele of each read, of size $n_a
	int *cnt_strand, *cnt_supp; // cnt_strand: count of supporting reads on both strands; cnt_supp: sum of both strands
	int *support, *support_strand; // support across entire $a. It points to the last "row" of cnt_supp/cnt_strand
	int len, max_len;
	char *seq;
	int *depth;
} paux_t;

static void count_alleles(paux_t *pa, int n)
{
	allele_t *a = pa->a;
	int i;
	a[0].k = 0; // the first allele is given allele id 0
	pa->max_del = a[0].indel < 0? -a[0].indel : 0;
	for (i = pa->n_alleles = 1; i < pa->n_a; ++i) {
		if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) // change of allele
			++pa->n_alleles;
		a[i].k = pa->n_alleles - 1;
		pa->max_del = pa->max_del > -a[i].indel? pa->max_del : -a[i].indel; // max deletion
	}
	// collect per-BAM counts
	pa->n_cnt = pa->n_alleles * (n + 1);
	if (pa->n_cnt > pa->max_cnt) { // expand the arrays if necessary
		pa->max_cnt = pa->n_cnt;
		kroundup32(pa->max_cnt);
		pa->cnt_strand = (int*)realloc(pa->cnt_strand, pa->max_cnt * 2 * sizeof(int));
		pa->cnt_supp = (int*)realloc(pa->cnt_supp, pa->max_cnt * sizeof(int));
	}
	memset(pa->cnt_strand, 0, pa->n_cnt * 2 * sizeof(int));
	pa->support_strand = pa->cnt_strand + pa->n_alleles * n * 2;
	memset(pa->cnt_supp, 0, pa->n_cnt * sizeof(int));
	pa->support = pa->cnt_supp + pa->n_alleles * n; // points to the last row of cnt_supp
	for (i = 0; i < pa->n_a; ++i) { // compute counts and sums of qualities
		int j = (a[i].pos>>32)*pa->n_alleles + a[i].k;
		pa->cnt_strand[j<<1|a[i].is_rev]++;
		pa->cnt_supp[j]++;
		pa->support[a[i].k]++;
		pa->support_strand[a[i].k<<1|a[i].is_rev]++;
	}
}

int main(int argc, char *argv[])
{
	int i, j, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, l_ref = 0, min_support = 1, min_support_strand = 0, min_supp_len = 0;
	int is_vcf = 0, var_only = 0, show_2strand = 0, trim_len = 0, del_as_allele = 0;
	int last_tid;
	const bam_pileup1_t **plp;
	char *ref = 0, *reg = 0, *chr_end; // specified region
	char *fname = 0; // reference fasta
	faidx_t *fai = 0;
	bam_hdr_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	paux_t aux;
	bam_mplp_t mplp;
	void *bed = 0;
	ketopt_t o = KETOPT_INIT;

	// parse the command line
	while ((n = ketopt(&o, argc, argv, 1, "r:q:Q:l:f:vcCS:s:b:T:ea:y", 0)) >= 0) {
		if (n == 'f') { fname = o.arg; fai = fai_load(fname); }
		else if (n == 'b') bed = bed_read(o.arg);
		else if (n == 'l') min_len = atoi(o.arg); // minimum query length
		else if (n == 'r') reg = strdup(o.arg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(o.arg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(o.arg);    // mapping quality threshold
		else if (n == 's') min_support = atoi(o.arg);
		else if (n == 'a') min_support_strand = atoi(o.arg);
		else if (n == 'S') min_supp_len = atoi(o.arg);
		else if (n == 'v') var_only = 1;
		else if (n == 'c') is_vcf = var_only = 1;
		else if (n == 'C') show_2strand = 1;
		else if (n == 'T') trim_len = atoi(o.arg);
		else if (n == 'e') del_as_allele = 1;
		else if (n == 'y') mapQ = 30, baseQ = 20, min_support = 5, min_support_strand = 2, is_vcf = var_only = show_2strand = 1;
	}
	if (min_support < 1) min_support = 1;
	if (is_vcf && fai == 0) {
		fprintf(stderr, "[E::%s] with option -c, the reference genome must be provided.\n", __func__);
		return 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minipileup [options] in1.bam [in2.bam [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  General:\n");
		fprintf(stderr, "    -f FILE      reference genome [null]\n");
		fprintf(stderr, "    -v           show variants only\n");
		fprintf(stderr, "    -c           output in the VCF format (force -v)\n");
		fprintf(stderr, "    -C           show count of each allele on both strands\n");
		fprintf(stderr, "    -e           use '*' to mark deleted bases\n");
		fprintf(stderr, "  Filtering:\n");
		fprintf(stderr, "    -r STR       region in format of 'ctg:start-end' [null]\n");
		fprintf(stderr, "    -b FILE      BED or position list file to include [null]\n");
		fprintf(stderr, "    -q INT       minimum mapping quality [%d]\n", mapQ);
		fprintf(stderr, "    -Q INT       minimum base quality [%d]\n", baseQ);
		fprintf(stderr, "    -l INT       minimum query length [%d]\n", min_len);
		fprintf(stderr, "    -S INT       minimum supplementary alignment length [0]\n");
		fprintf(stderr, "    -V FLOAT     skip alignment with per-base divergence >FLOAT [1]\n");
		fprintf(stderr, "    -T INT       skip bases within INT-bp from either end of a read [0]\n");
		fprintf(stderr, "    -s INT       drop alleles with depth<INT [%d]\n", min_support);
		fprintf(stderr, "    -a INT       drop alleles with depth<INT on either strand [%d]\n", min_support_strand);
		fprintf(stderr, "    -y           variant calling mode (-vcC -a2 -s5 -q30 -Q20)\n");
		return 1;
	}

	// initialize the auxiliary data structures
	n = argc - o.ind; // the number of BAMs on the command line
	data = (aux_t**)calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	if (reg) {
		chr_end = (char*)hts_parse_reg(reg, &beg, &end);
		ref = fai? fai_fetch(fai, reg, &l_ref) : 0;
	} else chr_end = 0;

	// load the index or put the file position at the right place
	last_tid = -1;
	for (i = 0; i < n; ++i) {
		bam_hdr_t *htmp;
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = bgzf_open(argv[o.ind+i], "r"); // open BAM
		data[i]->min_mapQ = mapQ;                     // set the mapQ filter
		data[i]->min_len  = min_len;                  // set the qlen filter
		data[i]->min_supp_len = min_supp_len;
		data[i]->bed = bed;
		htmp = bam_hdr_read(data[i]->fp);             // read the BAM header
		if (i == 0 && chr_end) {
			char c = *chr_end;
			*chr_end = 0;
			last_tid = tid = bam_name2id(htmp, reg);
			*chr_end = c;
		}
		if (i) bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header
		else h = htmp; // keep the header of the 1st BAM
		if (tid >= 0) { // if a region is specified and parsed successfully
			hts_idx_t *idx = bam_index_load(argv[o.ind+i]); // load the index
			data[i]->itr = bam_itr_queryi(idx, tid, beg, end); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
		}
		data[i]->h = h;
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = (int*)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = (const bam_pileup1_t**)calloc(n, sizeof(const bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
	memset(&aux, 0, sizeof(paux_t));
	if (is_vcf) {
		puts("##fileformat=VCFv4.2");
		printf("##source=minipileup-%s\n", VERSION);
		if (fai) {
			printf("##reference=%s\n", fname);
			int i, n = faidx_fetch_nseq(fai);
			for (i=0; i<n; i++) {
				const char *seq = faidx_iseq(fai,i);
				int len = faidx_seq_len(fai, seq);
				printf("##contig=<ID=%s,length=%d>\n", seq, len);
			}
		}
		puts("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		if (show_2strand) {
			puts("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
			puts("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
		} else puts("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
		fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", stdout);
		for (i = 0; i < n; ++i) printf("\t%s", argv[o.ind+i]);
		putchar('\n');
	}
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && !bed_overlap(bed, h->target_name[tid], pos, pos + 1)) continue; // not overlapping BED
		for (i = aux.tot_dp = 0; i < n; ++i) aux.tot_dp += n_plp[i];
		if (last_tid != tid) {
			if (fai) { // switch of chromosomes
				free(ref);
				ref = fai_fetch(fai, h->target_name[tid], &l_ref);
			}
			last_tid = tid; aux.len = 0;
		}
		if (aux.tot_dp) {
			int k, r = 15, shift = 0, qual;
			allele_t *a;
			if (aux.tot_dp + 1 > aux.max_dp) { // expand array
				aux.max_dp = aux.tot_dp + 1;
				kroundup32(aux.max_dp);
				aux.a = (allele_t*)realloc(aux.a, aux.max_dp * sizeof(allele_t));
			}
			a = aux.a;
			// collect alleles
			r = (ref && pos - beg < l_ref)? seq_nt16_table[(int)ref[pos - beg]] : 15; // the reference allele
			for (i = aux.n_a = 0; i < n; ++i)
				for (j = 0; j < n_plp[i]; ++j) {
					a[aux.n_a] = pileup2allele(&plp[i][j], baseQ, (uint64_t)i<<32 | j, r, trim_len, del_as_allele);
					if (!a[aux.n_a].is_skip) ++aux.n_a;
				}
			if (aux.n_a == 0) continue; // no reads are good enough; zero effective coverage
			// count alleles
			ks_introsort(allele, aux.n_a, aux.a);
			count_alleles(&aux, n);
			// squeeze out weak alleles
			for (i = k = 0; i < aux.n_a; ++i)
				if (aux.support[a[i].k] >= min_support && aux.support_strand[a[i].k<<1] >= min_support_strand && aux.support_strand[a[i].k<<1|1] >= min_support_strand)
					a[k++] = a[i];
			if (k < aux.n_a) {
				if (k == 0) continue; // no alleles are good enough
				aux.n_a = k;
				count_alleles(&aux, n);
			}

			if (var_only && aux.n_alleles == 1 && a[0].hash>>63 == 0) continue; // var_only mode, but no ALT allele; skip
			if (var_only && aux.n_alleles <= 2 && del_as_allele) {
				int n_ref = 0, n_del = 0;
				for (i = 0; i < aux.n_a; ++i)
					if (a[i].is_del) ++n_del;
					else if (a[i].hash>>63 == 0) ++n_ref;
				if (n_ref + n_del == aux.n_a) continue;
			}
			// print VCF or allele summary
			fputs(h->target_name[tid], stdout); printf("\t%d", pos+1);
			if (is_vcf) {
				fputs("\t.\t", stdout);
				for (i = 0; i <= aux.max_del; ++i) // print the reference allele up to the longest deletion
					putchar(ref && pos + i < l_ref + beg? ref[pos + i - beg] : 'N');
				putchar('\t');
			} else printf("\t%c\t", ref && pos < l_ref + beg? ref[pos - beg] : 'N'); // print a single reference base
																					 // print alleles
			if (!is_vcf || a[0].hash>>63) { // print if there is no reference allele
				print_allele(&plp[a[0].pos>>32][(uint32_t)a[0].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf, del_as_allele);
				if (aux.n_alleles > 1) putchar(',');
			}
			for (i = k = 1; i < aux.n_a; ++i)
				if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) {
					print_allele(&plp[a[i].pos>>32][(uint32_t)a[i].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf, del_as_allele);
					if (++k != aux.n_alleles) putchar(',');
				}
			if (is_vcf && aux.n_alleles == 1 && a[0].hash>>63 == 0) putchar('.'); // print placeholder if there is only the reference allele
			// compute and print qual
			for (i = !(a[0].hash>>63), qual = 0; i < aux.n_alleles; ++i)
				qual = qual > aux.support[i]? qual : aux.support[i];
			if (is_vcf) printf("\t%d\t.\t.\tGT:%s", qual, show_2strand? "ADF:ADR" : "AD");
			// print counts
			shift = (is_vcf && a[0].hash>>63); // in VCF, if there is no ref allele, we need to shift the allele number
			for (i = k = 0; i < n; ++i, k += aux.n_alleles) {
				int max1 = 0, max2 = 0, a1 = -1, a2 = -1, *sum_q = &aux.cnt_supp[k];
				// estimate genotype
				for (j = 0; j < aux.n_alleles; ++j)
					if (sum_q[j] > max1) max2 = max1, a2 = a1, max1 = sum_q[j], a1 = j;
					else if (sum_q[j] > max2) max2 = sum_q[j], a2 = j;
				if (max1 == 0 || (min_support > 0 && max1 < min_support)) a1 = a2 = -1;
				else if (max2 == 0 || (min_support > 0 && max2 < min_support)) a2 = a1;
				// print genotypes
				if (a1 < 0) printf("\t./.:");
				else printf("\t%d/%d:", a1 + shift, a2 + shift);
				// print counts
				if (show_2strand) {
					if (shift) fputs("0,", stdout);
					for (j = 0; j < aux.n_alleles; ++j) {
						if (j) putchar(',');
						printf("%d", aux.cnt_strand[(k+j)<<1]);
					}
					putchar(':');
					if (shift) fputs("0,", stdout);
					for (j = 0; j < aux.n_alleles; ++j) {
						if (j) putchar(',');
						printf("%d", aux.cnt_strand[(k+j)<<1|1]);
					}
				} else {
					if (shift) fputs("0,", stdout);
					for (j = 0; j < aux.n_alleles; ++j) {
						if (j) putchar(',');
						printf("%d", aux.cnt_supp[k+j]);
					}
				}
			} // ~for(i)
			putchar('\n');
		} // ~if(aux.tot_dp)
	} // ~while()
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_hdr_destroy(h);
	for (i = 0; i < n; ++i) {
		bgzf_close(data[i]->fp);
		if (data[i]->itr) bam_itr_destroy(data[i]->itr);
		free(data[i]);
	}
	if (ref) free(ref);
	if (fai) fai_destroy(fai);
	free(aux.cnt_strand); free(aux.cnt_supp); free(aux.a);
	free(aux.seq); free(aux.depth);
	free(data); free(reg);
	if (bed) bed_destroy(bed);
	fprintf(stderr, "[M::%s] done\n", __func__);
	return 0;
}
