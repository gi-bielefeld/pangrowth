//https://github.com/lh3/kmer-cnt
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "khashl.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "kthread.h"
#include "ketopt.h"

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

//#define MAX_KMER     31
#define N_COUNTS     8192
#define MAX_COUNT    ((1<<TOTAL_BITS)-1)

int ID_GENOME = 1;
int NUM_GENOMES;
int TOTAL_BITS;
int BITS_GENOME;
int MASK_COUNT; 
int MASK_GENOME;
int SUF;

#define ch_eq(a, b) ((a)>>SUF == (b)>>SUF) // lower 8 bits for counts; higher bits for k-mer
#define ch_hash(a) ((a)>>SUF)
KHASHL_MAP_INIT(, ht_t, ht, uint64_t, int32_t, ch_hash, ch_eq)

typedef struct {
	int32_t k;
	int32_t suf;
	int32_t n_thread;
    bool canonical;
	int64_t chunk_size;
    char* filelist;
} copt_t;

typedef struct {
	ht_t *h;
} ch1_t;

typedef struct {
	int k, suf, n_hash, n_shift;
	uint64_t tot;
	ch1_t *h;
} ch_t;


static inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/*** hash table ***/
ch_t *ch_init(int k, int suf) {
	ch_t *h;
	int i;
	//if (suf < TOTAL_BITS) return 0; //suf must be greater or equal (here is equal)
	CALLOC(h, 1);
	h->k = k, h->suf = suf;
	CALLOC(h->h, 1<<h->suf);
	for (i = 0; i < 1<<h->suf; ++i)
		h->h[i].h = ht_init();
	return h;
}

void ch_destroy(ch_t *h) {
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->suf; ++i)
		ht_destroy(h->h[i].h);
	free(h->h); free(h);
}

int ch_insert_list(ch_t *h, int n, const uint64_t *a) {
	if (n == 0) return 0;
	int mask = (1<<h->suf) - 1, n_ins = 0;
	ch1_t *g = &h->h[a[0]&mask];
	for (int j = 0; j < n; ++j) {
		int absent;
        //khint_t k = ht_put(g->h, a[j]&(~mask), &absent);
        khint_t k = ht_put(g->h, a[j], &absent);
        if (absent) ++n_ins;
        if (((kh_val(g->h, k)&MASK_GENOME) >> BITS_GENOME) != ID_GENOME) {
            kh_val(g->h, k) = ((kh_val(g->h, k)+1) & (~MASK_GENOME)) | 
                              (ID_GENOME << BITS_GENOME);
        }
	}
	return n_ins;
}

/*** generate histogram ***/
typedef struct {
	uint64_t c[N_COUNTS];
} buf_cnt_t;

typedef struct {
	const ch_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	ht_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
            ++cnt[kh_val(g, k)&MASK_COUNT];
}

void ch_hist(const ch_t *h, int64_t cnt[N_COUNTS], int n_thread) {
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, (NUM_GENOMES+1) * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread); // count is divide for number of threads
	kt_for(n_thread, worker_hist, &a, 1<<h->suf);
	for (i = 0; i <= NUM_GENOMES; ++i) cnt[i] = 0; //This will be the total amount
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i <= NUM_GENOMES; ++i)
			cnt[i] += a.cnt[j].c[i]; //Combine all
	free(a.cnt);
}


void chomp(char *str) {
    size_t len = strlen(str);

    if (len == 0) return;

    size_t last_idx = len - 1;
    if( str[last_idx] == '\n' ) {
        str[last_idx] = '\0';
    }
}

/****************
 * From count.c *
 ****************/
#include <zlib.h>
#include <string.h>
#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void copt_init(copt_t *o) {
	memset(o, 0, sizeof(copt_t));
	o->k = 17;
	o->suf = 10;
	o->n_thread = 4;
	o->canonical = true;
	o->chunk_size = 10000000;
}

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} ch_buf_t;

// insert a k-mer $y to a linear buffer
static inline void ch_insert_buf(ch_buf_t *buf, int suf, uint64_t y) {
    uint64_t y_suf = y & ((1<<suf) - 1);
	ch_buf_t *b = &buf[y_suf];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

// insert canonical k-mers in $seq to linear buffer $buf
static void count_seq_buf_can(ch_buf_t *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];      // Canonicals!
				ch_insert_buf(buf, suf, hash64(y, mask));
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// insert k-mers in $seq to linear buffer $buf
static void count_seq_buf(ch_buf_t *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x, mask = (1ULL<<k*2) - 1;
	for (i = l = 0, x = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x = (x << 2 | c) & mask;                  // forward strand
			if (++l >= k) { // we find a k-mer
				ch_insert_buf(buf, suf, hash64(x, mask));
			}
		} else l = 0, x = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const copt_t *opt;
	kseq_t *ks;
	ch_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_t *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	ch_buf_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) { // callback for kt_for()
	stepdat_t *s = (stepdat_t*)data;
	ch_buf_t *buf = &s->buf[i];
	ch_t *h = s->p->h;
	buf->n_ins += ch_insert_list(h, buf->n, buf->a);
}

static void *worker_pipeline(void *data, int step, void *in) { // callback for kt_pipeline()
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		CALLOC(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (l < p->opt->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				REALLOC(s->len, s->m);
				REALLOC(s->seq, s->m);
			}
			MALLOC(s->seq[s->n], l);
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
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->suf, m;
		CALLOC(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			MALLOC(s->buf[i].a, m);
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
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->opt->suf;
		uint64_t n_ins = 0;
        kt_for(p->opt->n_thread, worker_for, s, n);
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

ch_t *count(const char *fn, const copt_t *opt, ch_t *h) {
	pldat_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.opt = opt;
    pl.h = (!h) ? ch_init(opt->k, opt->suf) : h;
	kt_pipeline(3, worker_pipeline, &pl, 3);
	kseq_destroy(pl.ks);
	gzclose(fp);
	return pl.h;
}

ch_t *count_file(const char *fn1, const copt_t *opt, ch_t *h2) {
	ch_t *h;
	h = count(fn1, opt, h2); 
	return h;
}

uint32_t count_fasta(char* filelist) {
    FILE * fp;
    char * line = NULL;
    uint32_t num_lines=0;
    size_t len = 0;

    fp = fopen(filelist, "r");
    if (fp == NULL) {
        fprintf(stderr, "Could not open file: %s",filelist);
        exit(EXIT_FAILURE);
    }

    while ((getline(&line, &len, fp)) != -1) num_lines++;
    fclose(fp);

    if (line) free(line);
    return num_lines;
}

void output_hist(int argc, char *argv[]){
	ch_t *h = 0;
	int c;
	int i;
	copt_t opt;
	ketopt_t o = KETOPT_INIT;
	copt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:s:K:t:b:i:", 0)) >= 0) {
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
    BITS_GENOME = ceil(log2(NUM_GENOMES+1));
    MASK_COUNT  = ((1 << BITS_GENOME) - 1);
    MASK_GENOME = ((1 << BITS_GENOME) - 1) << (BITS_GENOME);
    TOTAL_BITS  = BITS_GENOME*2;
    SUF = opt.suf;

    fprintf(stderr, "Number of genomes: %d\n", NUM_GENOMES);
    fprintf(stderr, "Number of bits per genome: %d\n", BITS_GENOME);
    fprintf(stderr, "Suffix size for k-mer: %d\n", SUF);
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
                h = count_file(line, &opt, h);
            } else {
                h = count_file(line, &opt, 0);
            }
        }
        fclose(fp);
        if (line) free(line);
    } else {
        h = count_file(argv[o.ind], &opt, 0);
        for (i = 1; i < NUM_GENOMES; i++){
            ID_GENOME++;
            h = count_file(argv[o.ind+i], &opt, h);
        }
    }
	fprintf(stderr, "[M::%s] %ld distinct k-mers\n", __func__, (long)h->tot);

	int64_t cnt[N_COUNTS];

	ch_hist(h, cnt, opt.n_thread);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%d\t%lld\n", i,(long long)cnt[i]);
	for (i = 1; i <= NUM_GENOMES; ++i) printf("%lld\n", (long long)cnt[i]);
	ch_destroy(h);
}
