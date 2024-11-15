//https://github.com/lh3/kmer-cnt
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "khashl.h"
#include "kthread.h"
#include "ketopt.h"
#include "kmer.h"
#include "kseq.h"

///*** hash table kmer ***/
//#define MAX_KMER     31
#define ch_eq(a, b) ((a)>>SUF == (b)>>SUF) // lower 8 bits for counts; higher bits for k-mer
#define ch_hash(a) ((a)>>SUF)
KHASHL_MAP_INIT(, hashtable_kmer_t, ht, uint64_t, int32_t, ch_hash, ch_eq)

typedef struct {
	hashtable_kmer_t *h;
} hat_kmer_s;

typedef struct {
	int k, suf, n_hash, n_shift;
	uint64_t tot;
	hat_kmer_s *h;
} multi_hat_kmer_s;

/*** hash table ***/
multi_hat_kmer_s *ch_init(int k, int suf) {
	multi_hat_kmer_s *h;
	int i;
	//if (suf < TOTAL_BITS) return 0; //suf must be greater or equal (here is equal)
	mcalloc(h, 1);
	h->k = k, h->suf = suf;
	mcalloc(h->h, 1<<h->suf);
	for (i = 0; i < 1<<h->suf; ++i)
		h->h[i].h = ht_init();
	return h;
}

void ch_destroy(multi_hat_kmer_s *h) {
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->suf; ++i)
		ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

int hat_insert_kmers(multi_hat_kmer_s *h, int n, const uint64_t *a) {
	if (n == 0) return 0;
	int mask = (1<<h->suf) - 1, n_ins = 0;
	hat_kmer_s *g = &h->h[a[0]&mask];
	for (int j = 0; j < n; ++j) {
		int absent;
        //khint_t k = ht_put(g->h, a[j]&(~mask), &absent);
        khint_t k = ht_put(g->h, a[j], &absent);
        if (absent) ++n_ins;
        if (((uint32_t)(kh_val(g->h, k)&MASK_GENOME) >> BITS_GENOME) != ID_GENOME) {
            kh_val(g->h, k) = ((kh_val(g->h, k)+1) & (~MASK_GENOME)) | 
                              (ID_GENOME << BITS_GENOME);
        }
	}
	return n_ins;
}


typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} buf_kmer_s;

// insert a k-mer $y to a linear buffer
static inline void kmer_insert_buf(buf_kmer_s *buf, int suf, uint64_t y) {
    uint64_t y_suf = y & ((1<<suf) - 1);
	buf_kmer_s *b = &buf[y_suf];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		mrealloc(b->a, b->m);
	}
	b->a[b->n++] = y;
}

// insert canonical k-mers in $seq to linear buffer $buf
static void count_seq_buf_can(buf_kmer_s *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];      // Canonicals!
				kmer_insert_buf(buf, suf, hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// insert k-mers in $seq to linear buffer $buf
static void count_seq_buf(buf_kmer_s *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x, mask = (1ULL<<k*2) - 1;
	for (i = l = 0, x = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x = (x << 2 | c) & mask;                  // forward strand
			if (++l >= k) { // we find a k-mer
				kmer_insert_buf(buf, suf, hash64(x, mask));
			}
		} else l = 0, x = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const param_t *opt;
	kseq_t *ks;
	multi_hat_kmer_s *h;
} pldat_kmer_s;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_kmer_s *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	buf_kmer_s *buf;
} stepdat_kmer_s;

static void worker_for_kmer_insert_list(void *data, long i, int tid) { // callback for kt_for()
	stepdat_kmer_s *s = (stepdat_kmer_s*)data;
	buf_kmer_s *buf = &s->buf[i];
	multi_hat_kmer_s *h = s->p->h;
	buf->n_ins += hat_insert_kmers(h, buf->n, buf->a);
}

static void *worker_pipe_kmer(void *data, int step, void *in) { // callback for kt_pipeline()
	pldat_kmer_s *p = (pldat_kmer_s*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_kmer_s *s;
		mcalloc(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				mrealloc(s->len, s->m);
				mrealloc(s->seq, s->m);
			}
			mmalloc(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->ks->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->opt->k + 1;
			if (s->sum_len >= p->opt->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_kmer_s *s = (stepdat_kmer_s*)in;
		int i, n = 1<<p->opt->suf, m;
		mcalloc(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			mmalloc(s->buf[i].a, m);
		}
        if (p->opt->canonical) {
            for (i = 0; i < s->n; ++i) {
                count_seq_buf_can(s->buf, p->opt->k, p->opt->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        } else {
            for (i = 0; i < s->n; ++i) {
                count_seq_buf(s->buf, p->opt->k, p->opt->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        }
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_kmer_s *s = (stepdat_kmer_s*)in;
		int i, n = 1<<p->opt->suf;
		uint64_t n_ins = 0;
        kt_for(p->opt->n_thread, worker_for_kmer_insert_list, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		p->h->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M] processed %d sequences; %ld distinct k-mers in the hash table\n", s->n, (long)p->h->tot);
		free(s);
	}
	return 0;
}

multi_hat_kmer_s *count_kmer(const char *fn, const param_t *opt, multi_hat_kmer_s *h) {
	pldat_kmer_s pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.opt = opt;
    pl.h = (!h) ? ch_init(opt->k, opt->suf) : h;
	kt_pipeline(3, worker_pipe_kmer, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

multi_hat_kmer_s *count_kmer_file(const char *fn1, const param_t *opt, multi_hat_kmer_s *h2) {
	multi_hat_kmer_s *h;
	h = count_kmer(fn1, opt, h2); 
	return h;
}

/*** generate histogram ***/
typedef struct {
	uint64_t c[N_COUNTS];
} buf_hist_t;

typedef struct {
	const multi_hat_kmer_s *h;
	buf_hist_t *cnt;
} hist_kmer_s;

static void worker_hist_kmer(void *data, long i, int tid) // callback for kt_for()
{
	hist_kmer_s *a = (hist_kmer_s*)data;
	uint64_t *cnt = a->cnt[tid].c;
	hashtable_kmer_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
            ++cnt[kh_val(g, k)&MASK_COUNT];
}

void hist_kmer(const multi_hat_kmer_s *h, int64_t cnt[N_COUNTS], int n_thread) {
	hist_kmer_s a;
	int i, j;
	a.h = h;
	memset(cnt, 0, (NUM_GENOMES+1) * sizeof(uint64_t));
	mcalloc(a.cnt, n_thread); // count_kmer is divide for number of threads
	kt_for(n_thread, worker_hist_kmer, &a, 1<<h->suf);
	for (i = 0; i <= NUM_GENOMES; ++i) cnt[i] = 0; //This will be the total amount
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i <= NUM_GENOMES; ++i)
			cnt[i] += a.cnt[j].c[i]; //Combine all
	free(a.cnt);
}

void print_kmer_debug(multi_hat_kmer_s * ht) {
    int i;
    uint32_t j;
	for (i = 0; i < 1<<ht->suf; ++i) {
        hashtable_kmer_t *g = ht->h[i].h;
        for (j = 0; j < kh_end(g); ++j)
            if (kh_exist(g, j))
                printf("[%d] %s %d\n", i, bits2kmer(kh_key(g, j), ht->k), kh_val(g,j)&MASK_COUNT);
    }
}

void output_hist_kmer(int argc, char *argv[]){
	multi_hat_kmer_s *h = 0;
	int c;
	int i;
	param_t opt;
	ketopt_t o = KETOPT_INIT;
	param_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:K:t:i:b", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 's') opt.suf = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'b') opt.canonical = false;
		else if (c == 'i') opt.filelist = o.arg;
	}

	if (argc - o.ind < 1 && !opt.filelist) {
		fprintf(stderr, "Usage: pangrowth hist [options] <in.fa> [in.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -i PATH    file containing a list of fasta files on each line\n");
		fprintf(stderr, "  -b         turn off transformation into canonical [%d]\n", opt.canonical);
		fprintf(stderr, "  -s INT     suffix size for k-mer [%d]\n", opt.n_thread);
		//fprintf(stderr, "  -K INT     chunk size [100m]\n");
		return;
	}

    //if (opt.k == -1) {
	//	fprintf(stderr, "Calculating k\n");
    //}

    // Calculate necessary bits
    if (opt.filelist) {
        NUM_GENOMES = count_fasta(opt.filelist);
    } else {
        NUM_GENOMES = argc - o.ind;
    }
    ID_GENOME = 1;
    BITS_GENOME = ceil(log2(NUM_GENOMES+1));
    MASK_COUNT  = ((1 << BITS_GENOME) - 1);
    MASK_GENOME = ((1 << BITS_GENOME) - 1) << (BITS_GENOME);
    TOTAL_BITS  = BITS_GENOME*2;
    SUF = opt.suf;

    fprintf(stderr, "Number of genomes:\t%d\n", NUM_GENOMES);
    fprintf(stderr, "Number of threads:\t%d\n", opt.n_thread);
    fprintf(stderr, "Bits per genome:\t%d\n", BITS_GENOME);
    fprintf(stderr, "Bits per suffix:\t%d\n", SUF);
    fprintf(stderr, "Counting %s k-mers\n", opt.canonical ? "canonical" : "forward");
    if (opt.filelist) {
        FILE * fp;
        char * line = NULL;
        size_t len = 0;

        fp = fopen(opt.filelist, "r");
        if (fp == NULL) {
            fprintf(stderr, "Could not open file: %s",opt.filelist);
            exit(EXIT_FAILURE);
        }

        while ((getline(&line, &len, fp)) != -1) {
            chomp(line);
            if (h) {
                ID_GENOME++;
                h = count_kmer_file(line, &opt, h);
            } else {
                h = count_kmer_file(line, &opt, 0);
            }
        }
        fclose(fp);
        if (line) free(line);
    } else {
        h = count_kmer_file(argv[o.ind], &opt, 0);
        for (i = 1; i < NUM_GENOMES; i++){
            ID_GENOME++;
            h = count_kmer_file(argv[o.ind+i], &opt, h);
        }
    }
	fprintf(stderr, "[M::%s] %ld distinct k-mers\n", __func__, (long)h->tot);

	int64_t cnt[N_COUNTS];

    //print_kmer_debug(h);

	hist_kmer(h, cnt, opt.n_thread);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%d\t%lld\n", i,(long long)cnt[i]);
	for (i = 1; i <= NUM_GENOMES; ++i) printf("%lld\n", (long long)cnt[i]);
	ch_destroy(h);
}
