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

#define PRINT_BITS(num) {int _ii; for(_ii = sizeof((__typeof__(num))num)*8-1; _ii>=0; _ii--) printf("%c", (num & (1ULL << _ii))?'1':'0'); printf("\n");}
#define mcalloc(pos, len) ((pos) = (__typeof__(pos))calloc((len), sizeof(*(pos))))
#define mmalloc(pos, len) ((pos) = (__typeof__(pos))malloc((len) * sizeof(*(pos))))
#define mrealloc(pos, len) ((pos) = (__typeof__(pos))realloc((pos), (len) * sizeof(*(pos))))

//#define MAX_KMER     31
#define N_COUNTS     8192
#define MAX_COUNT    ((1<<TOTAL_BITS)-1)

char* bits2kmer(uint64_t kmer_bits, int k) {
    char nucleotides[4] = {'A','C','G','T'};
    char *kmer_char;
    mcalloc(kmer_char, k+1);
    for (int i = 0; i < k; i++) {
        kmer_char[i] = nucleotides[(kmer_bits>>(2*k-(i+1)*2)) & 3];
    }
    kmer_char[k] ='\0';
    return kmer_char;
}

uint32_t ID_GENOME = 1;
int NUM_GENOMES;
int TOTAL_BITS;
int BITS_GENOME;
int MASK_COUNT; 
int MASK_GENOME;
int SUF;

typedef struct {
    uint32_t edges[16];
    uint8_t last_edge;
    uint32_t last_genome;
    uint32_t sigma; 
} infix_storage_t;

#define ch_eq(a, b) ((a) == (b))
#define ch_hash(a) ((a))
KHASHL_MAP_INIT(, hashtable_t, ht, uint64_t, infix_storage_t, ch_hash, ch_eq)

typedef struct {
	int32_t k;
	int32_t suf;
	int32_t n_thread;
    bool canonical;
	int64_t chunk_size;
    char* filelist;
} param_t;

typedef struct {
	hashtable_t *ht;
} ch1_t;

typedef struct {
	int k, suf, n_hash, n_shift;
	uint64_t tot;
	ch1_t *ht;
} multi_hashtable_t;

static inline uint64_t hash64(uint64_t key, uint64_t mask) { // invertible integer hash function
	//key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	//key = key ^ key >> 24;
	//key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	//key = key ^ key >> 14;
	//key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	//key = key ^ key >> 28;
	//key = (key + (key << 31)) & mask;
	return key;
}

static inline uint64_t get_infix (uint64_t kmer_bits, int k) {
    uint64_t mask = (1 << (2 * (k - 1))) - 1;
    return (kmer_bits >> 2) & mask;
}


/*** hash table ***/
multi_hashtable_t *ch_init(int k, int suf) {
	multi_hashtable_t *ht;
	int i;
	//if (suf < TOTAL_BITS) return 0; //suf must be greater or equal (here is equal)
	mcalloc(ht, 1);
	ht->k = k, ht->suf = suf;
	mcalloc(ht->ht, 1<<ht->suf);
	for (i = 0; i < 1<<ht->suf; ++i)
		ht->ht[i].ht = ht_init();
	return ht;
}

void ch_destroy(multi_hashtable_t *ht) {
	int i;
	if (ht == 0) return;
	for (i = 0; i < 1<<ht->suf; ++i)
		ht_destroy(ht->ht[i].ht);
	free(ht->ht); free(ht);
}

int ch_insert_list(multi_hashtable_t *ht, int n, const uint64_t *a) {
	if (n == 0) return 0;
    //uint64_t mask = ((1<<(ht->suf+2))-4);
    uint64_t mask = (1<<(ht->suf))-1;
    int n_ins = 0;
    int k = ht->k;
    //hashtable_t *h = ht->ht[(a[0]&mask)>>2].ht;
    hashtable_t *h = ht->ht[(get_infix(a[0],k)&mask)].ht;
	for (int j = 0; j < n; ++j) {
        uint64_t kp1mer = a[j];
        uint64_t infix = get_infix(kp1mer, k);
        uint8_t first_nt = (kp1mer >> (2 * k));
        uint8_t last_nt = kp1mer & 3ULL;
        uint8_t combined_nt = ((first_nt << 2) | last_nt);
        //printf("%s\n", bits2kmer(kp1mer,k+1)); 
        //printf("%s\n\n", bits2kmer(infix, k+1)); 
		int absent;
        khint_t z = ht_put(h, infix, &absent);
        //printf("%s %d\n", bits2kmer(infix,k-1), absent);
        if (absent) { 
            n_ins++;
            memset(kh_val(h,z).edges, 0, 16*sizeof(uint32_t));
            kh_val(h,z).last_edge = combined_nt;
            kh_val(h,z).edges[kh_val(h,z).last_edge]++;
            kh_val(h,z).last_genome = ID_GENOME;
            kh_val(h,z).sigma = 1;
        } else if (kh_val(h,z).last_genome == ID_GENOME && 
                   kh_val(h,z).last_edge != combined_nt && 
                   kh_val(h,z).last_edge != 255) {
            kh_val(h,z).edges[kh_val(h,z).last_edge]--;
            kh_val(h,z).last_edge = 255;
        } else if (kh_val(h,z).last_genome != ID_GENOME) {
            kh_val(h,z).last_edge = combined_nt;
            kh_val(h,z).edges[kh_val(h,z).last_edge]++;
            kh_val(h,z).last_genome = ID_GENOME;
            kh_val(h,z).sigma++;
        }
        //printf("%p %s(%d) %d %d\n", &kh_val(h,z), bits2kmer(kh_val(h,z).last_edge, 2), kh_val(h,z).last_edge, kh_val(h,z).last_genome, kh_val(h,z).sigma);
	}
	return n_ins;
}

/*** generate histogram ***/
//typedef struct {
//	uint64_t c[N_COUNTS];
//} buf_cnt_t;
//
//typedef struct {
//	const multi_hashtable_t *ht;
//	buf_cnt_t *cnt;
//} hist_aux_t;

//static void worker_hist(void *data, long i, int tid) { // callback for kt_for()
//	hist_aux_t *a = (hist_aux_t*)data;
//	uint64_t *cnt = a->cnt[tid].c;
//	hashtable_t *g = a->ht->ht[i].ht;
//	khint_t j;
//	//for (j = 0; j < kh_end(g); ++j)
//	//	if (kh_exist(g, j))
//    //        ++cnt[kh_val(g, j)&MASK_COUNT];
//}
//
//void hist(const multi_hashtable_t *ht, int64_t cnt[N_COUNTS], int n_thread) {
//	hist_aux_t a;
//	int i, j;
//	a.ht = ht;
//	memset(cnt, 0, (NUM_GENOMES+1) * sizeof(uint64_t));
//	mcalloc(a.cnt, n_thread); // count is divide for number of threads
//	kt_for(n_thread, worker_hist, &a, 1<<ht->suf);
//	for (i = 0; i <= NUM_GENOMES; ++i) cnt[i] = 0; //This will be the total amount
//	for (j = 0; j < n_thread; ++j)
//		for (i = 0; i <= NUM_GENOMES; ++i)
//			cnt[i] += a.cnt[j].c[i]; //Combine all
//	free(a.cnt);
//}


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

void copt_init(param_t *options) {
	memset(options, 0, sizeof(param_t));
	options->k = 17;
	options->suf = 10;
	options->n_thread = 4;
	options->canonical = true;
	options->chunk_size = 10000000;
}

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} ch_buf_t;

// insert a k-mer $y to a linear buffer
static inline void ch_insert_buf(ch_buf_t *buf, int suf, uint64_t y, int k) {
    //uint64_t mask = ((1<<(suf+2))-4);
    //uint64_t y_suf = (y & mask);
    //uint64_t y_suf = y & ((1<<suf) - 1);
    uint64_t infix = get_infix(y, k);
    uint64_t mask = (1<<(suf))-1;
    uint64_t y_suf = (infix & mask);
    //printf("%s\n",   bits2kmer(y, k+1));
    //printf("%s\n",   bits2kmer(infix, k+1));
    //printf("%s\n\n", bits2kmer(y_suf, k+1));
	ch_buf_t *b = &buf[y_suf];
	if (b->n == b->m) {
		b->m = b->m < 8 ? 8 : b->m + (b->m>>1);
		mrealloc(b->a, b->m);
	}
	b->a[b->n++] = y;
}

// insert canonical k-mers in $seq to linear buffer $buf
// k = k+1 here 
static void count_seq_buf_can(ch_buf_t *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x[2], mask = (1ULL<<(k+1)*2) - 1, shift = (k+1 - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k+1) { // we found a k+1-mer
                uint64_t infix[2] = {get_infix(x[0], k), get_infix(x[1], k)};
				uint64_t y = infix[0] < infix[1] ? x[0] : x[1];      // Canonicals!
				ch_insert_buf(buf, suf, hash64(y, mask), k);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// insert k-mers in $seq to linear buffer $buf
static void count_seq_buf(ch_buf_t *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x, mask = (1ULL<<(k+1)*2) - 1;
	for (i = l = 0, x = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x = (x << 2 | c) & mask;                  // forward strand
			if (++l >= k+1) { // we found a k+1-mer
				ch_insert_buf(buf, suf, hash64(x, mask), k);
			}
		} else l = 0, x = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const param_t *param;
	kseq_t *kseq;
	multi_hashtable_t *ht;
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
	multi_hashtable_t *ht = s->p->ht;
	buf->n_ins += ch_insert_list(ht, buf->n, buf->a);
}

static void *worker_pipeline(void *data, int step, void *in) { // callback for kt_pipeline()
	pldat_t *p = (pldat_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_t *s;
		mcalloc(s, 1);
		s->p = p;
		while ((ret = kseq_read(p->kseq)) >= 0) {
			int l = p->kseq->seq.l;
			if (l < p->param->k) continue;
			if (s->n == s->m) {
				s->m = s->m < 16? 16 : s->m + (s->n>>1);
				mrealloc(s->len, s->m);
				mrealloc(s->seq, s->m);
			}
			mmalloc(s->seq[s->n], l);
			memcpy(s->seq[s->n], p->kseq->seq.s, l);
			s->len[s->n++] = l;
			s->sum_len += l;
			s->nk += l - p->param->k + 1;
			if (s->sum_len >= p->param->chunk_size)
				break;
		}
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		stepdat_t *s = (stepdat_t*)in;
		int i, m;
        int n = 1<<p->param->suf;
		mcalloc(s->buf, n);
		m = (int)(s->nk * 1.2 / n) + 1;
		for (i = 0; i < n; ++i) {
			s->buf[i].m = m;
			mmalloc(s->buf[i].a, m);
		}
        if (p->param->canonical) {
            for (i = 0; i < s->n; ++i) {
                count_seq_buf_can(s->buf, p->param->k, p->param->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        } else {
            for (i = 0; i < s->n; ++i) {
                count_seq_buf(s->buf, p->param->k, p->param->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        }
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_t *s = (stepdat_t*)in;
		int i, n = 1<<p->param->suf;
		uint64_t n_ins = 0;
        kt_for(p->param->n_thread, worker_for, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		p->ht->tot += n_ins;
		free(s->buf);
		fprintf(stderr, "[M] processed %d sequences; %ld distinct k-mers in the hash table\n", s->n, (long)p->ht->tot);
		free(s);
	}
	return 0;
}

multi_hashtable_t *count(const char *fn, const param_t *param, multi_hashtable_t *ht) {
	pldat_t pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.kseq = kseq_init(fp);
	pl.param = param;
    pl.ht = (!ht) ? ch_init(param->k, param->suf) : ht;
	kt_pipeline(3, worker_pipeline, &pl, 3);
	kseq_destroy(pl.kseq);
	gzclose(fp);
	return pl.ht;
}

multi_hashtable_t *count_file(const char *fn1, const param_t *param, multi_hashtable_t *h2) {
	multi_hashtable_t *ht;
	ht = count(fn1, param, h2); 
	return ht;
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

void print_infix_debug(multi_hashtable_t* ht) {
    for (int i = 0; i < ht->k; i++) {
        printf(" ");
    }
    printf("\t");
    char nucleotides[4] = {'A','C','G','T'};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%c%c ",nucleotides[i], nucleotides[j]);
        }
    }
    printf(" sig\tlastG\tlastE\n");
	for (int i = 0; i < 1<<ht->suf; ++i) {
	    hashtable_t *h = ht->ht[i].ht;
        for (uint32_t z = 0; z < kh_end(h); ++z)
            if (kh_exist(h,z)) {
                printf("%s\t", bits2kmer(kh_key(h,z), ht->k-1));
                for (int j = 0; j < 16; j++) {
                    printf("%d  ", kh_val(h,z).edges[j]);
                }
                printf(" %d\t%d\t%s(%d)", kh_val(h,z).sigma, kh_val(h,z).last_genome, bits2kmer(kh_val(h,z).last_edge, 2), kh_val(h,z).last_edge);
                printf("\n");
            }
    }
}

void output_hist(int argc, char *argv[]){
	multi_hashtable_t *ht = 0;
	int c;
	int i;
	param_t param;
	ketopt_t options = KETOPT_INIT;
	copt_init(&param);
	while ((c = ketopt(&options, argc, argv, 1, "k:s:K:t:i:b", 0)) >= 0) {
		if (c == 'k') param.k = atoi(options.arg);
		else if (c == 's') param.suf = atoi(options.arg);
		else if (c == 'K') param.chunk_size = atoi(options.arg);
		else if (c == 't') param.n_thread = atoi(options.arg);
		else if (c == 'b') param.canonical = false;
		else if (c == 'i') param.filelist = options.arg;
	}

	if (argc - options.ind < 1 && !param.filelist) {
		fprintf(stderr, "Usage: pangrowth hist [options] <in.fa> [in.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", param.k);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", param.n_thread);
		fprintf(stderr, "  -i PATH    file containing a list of fasta files on each line\n");
		fprintf(stderr, "  -b         turn off transformation into canonical [%d]\n", param.canonical);
		fprintf(stderr, "  -s INT     suffix size for k-mer [%d]\n", param.suf);
		//fprintf(stderr, "  -K INT     chunk size [100m]\n");
		return;
	}

    //if (param.k == -1) {
	//	fprintf(stderr, "Calculating k\n");
    //}

    // Calculate necessary bits
    if (param.filelist) {
        NUM_GENOMES = count_fasta(param.filelist);
    } else {
        NUM_GENOMES = argc - options.ind;
    }
    BITS_GENOME = ceil(log2(NUM_GENOMES+1));
    MASK_COUNT  = ((1 << BITS_GENOME) - 1);
    MASK_GENOME = ((1 << BITS_GENOME) - 1) << (BITS_GENOME);
    TOTAL_BITS  = BITS_GENOME*2;
    SUF = param.suf;

    fprintf(stderr, "Number of genomes:\t%d\n", NUM_GENOMES);
    fprintf(stderr, "Number of threads:\t%d\n", param.n_thread);
    fprintf(stderr, "Bits per genome:\t%d\n", BITS_GENOME);
    fprintf(stderr, "Bits per suffix:\t%d\n", SUF);
    fprintf(stderr, "Counting %s k-mers\n", param.canonical ? "canonical" : "forward");

    if (param.filelist) {
        FILE * fp;
        char * line = NULL;
        size_t len = 0;

        fp = fopen(param.filelist, "r");
        if (fp == NULL) {
            fprintf(stderr, "Could not open file: %s",param.filelist);
            exit(EXIT_FAILURE);
        }

        while ((getline(&line, &len, fp)) != -1) {
            chomp(line);
            if (ht) {
                ID_GENOME++;
                ht = count_file(line, &param, ht);
            } else {
                ht = count_file(line, &param, 0);
            }
        }
        fclose(fp);
        if (line) free(line);
    } else {
        ht = count_file(argv[options.ind], &param, 0);
        for (i = 1; i < NUM_GENOMES; i++){
            ID_GENOME++;
            ht = count_file(argv[options.ind+i], &param, ht);
        }
    }
	fprintf(stderr, "[M::%s] %ld distinct k-mers\n", __func__, (long)ht->tot);

	print_infix_debug(ht);

	//int64_t cnt[N_COUNTS];
	//hist(ht, cnt, param.n_thread);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%d\t%lld\n", i,(long long)cnt[i]);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%lld\n", (long long)cnt[i]);
	ch_destroy(ht);
}
