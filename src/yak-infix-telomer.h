#pragma once
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

/*** hash table infix ***/
//#define MAX_KMER     31
typedef struct {
    uint16_t edges[24];
    uint8_t last_edge;
    uint32_t last_genome;
    uint16_t sigma; 
} infix_telomer_storage_t;
#define infix_eq(a, b) ((a) == (b))
#define infix_hash(a) ((a))
KHASHL_MAP_INIT(, hat_infix_tel_t, hti, uint64_t, infix_telomer_storage_t, infix_hash, infix_eq)

static inline uint64_t get_infix_tel (uint64_t kmer_bits, int k) {
    uint64_t mask = (1ULL << (2 * (k - 1))) - 1;
    return (kmer_bits >> 2) & mask;
}

typedef struct {
	hat_infix_tel_t *ht;
} hat_infix_tel_s;

typedef struct {
	int k, suf, n_hash, n_shift;
	uint64_t tot;
	hat_infix_tel_s *ht;
} multi_hat_infix_tel_s;

multi_hat_infix_tel_s *infix_tel_init(int k, int suf) {
	multi_hat_infix_tel_s *ht;
	int i;
	//if (suf < TOTAL_BITS) return 0; //suf must be greater or equal (here is equal)
	mcalloc(ht, 1);
	ht->k = k, ht->suf = suf;
	mcalloc(ht->ht, 1<<ht->suf);
	for (i = 0; i < 1<<ht->suf; ++i)
		ht->ht[i].ht = hti_init();
	return ht;
}

void infix_tel_destroy(multi_hat_infix_tel_s *ht) {
	int i;
	if (ht == 0) return;
	for (i = 0; i < 1<<ht->suf; ++i)
		hti_destroy(ht->ht[i].ht);
	free(ht->ht); free(ht);
}

int hat_insert_infix_tel(multi_hat_infix_tel_s *ht, int n, const uint64_t *a) {
	if (n == 0) return 0;
    //uint64_t mask = ((1<<(ht->suf+2))-4);
    uint64_t mask = (1<<(ht->suf))-1;
    int n_ins = 0;
    int k = ht->k;
    //hat_infix_tel_t *h = ht->ht[(a[0]&mask)>>2].ht;
    hat_infix_tel_t *h = ht->ht[(get_infix_tel(a[0],k)&mask)].ht;
	for (int j = 0; j < n; ++j) {
        uint64_t kp1mer = a[j];
        uint64_t infix = get_infix_tel(kp1mer, k);
        uint8_t first_nt = (kp1mer >> (2 * k));
        uint8_t last_nt = kp1mer & 3ULL;
        uint8_t combined_nt = ((first_nt << 2) | last_nt);
		int absent;
        khint_t z = hti_put(h, infix, &absent);
        //printf("%s\n", bits2kmer(kp1mer,k+1)); 
        //printf("%s %d\n", bits2kmer(infix,k+1), absent);
        if (absent) { 
            n_ins++;
            memset(kh_val(h,z).edges, 0, 24*sizeof(uint16_t));
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

/*** count infixes ***/
typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
} infix_tel_buf_s;

// insert a k-mer $y to a linear buffer
static inline void infix_tel_insert_buf(infix_tel_buf_s *buf, int suf, uint64_t y, int k) {
    //uint64_t mask = ((1<<(suf+2))-4);
    //uint64_t y_suf = (y & mask);
    //uint64_t y_suf = y & ((1<<suf) - 1);
    uint64_t infix = get_infix_tel(y, k);
    uint64_t mask = (1<<(suf))-1;
    uint64_t y_suf = (infix & mask);
    //printf("%s\n",   bits2kmer(y, k+1));
    //printf("%s\n",   bits2kmer(infix, k+1));
    //printf("%s\n\n", bits2kmer(y_suf, k+1));
	infix_tel_buf_s *b = &buf[y_suf];
	if (b->n == b->m) {
		b->m = b->m < 8 ? 8 : b->m + (b->m>>1);
		mrealloc(b->a, b->m);
	}
	b->a[b->n++] = y;
}

// insert canonical k-mers in $seq to linear buffer $buf
// k = k+1 here 
static void count_seq_buf_can_infix_tel(infix_tel_buf_s *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x[2], mask_kp1 = (1ULL<<(k+1)*2) - 1, shift = (k+1 - 1) * 2, mask_k = mask_kp1 >> 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask_kp1;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k+1 && (x[0] >> 2) != (x[0] & mask_k) && x[0] != x[1]) { 
                // we found a k+1-mer and it is not formed by the same two k-mers x != a^k
                // and x != yy^rc which means x != x^rc
                uint64_t infix[2] = {get_infix_tel(x[0], k), get_infix_tel(x[1], k)};
				uint64_t y = infix[0] < infix[1] ? x[0] : x[1];      // Canonicals!
				infix_tel_insert_buf(buf, suf, hash64(y, mask_kp1), k);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}

// insert k-mers in $seq to linear buffer $buf
static void count_seq_buf_infix_tel(infix_tel_buf_s *buf, int k, int suf, int len, const char *seq){ 
	int i, l;
	uint64_t x, mask_kp1 = (1ULL<<(k+1)*2) - 1, mask_k = mask_kp1 >> 2;
	for (i = l = 0, x = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x = (x << 2 | c) & mask_kp1;                  // forward strand
			if (++l >= k+1 && (x >> 2) != (x & mask_k)) { 
                // we found a k+1-mer and it is not formed by the same two k-mers x != a^k
				infix_tel_insert_buf(buf, suf, hash64(x, mask_kp1), k);
			}
		} else l = 0, x = 0; // if there is an "N", restart
	}
}

typedef struct { // global data structure for kt_pipeline()
	const param_t *param;
	kseq_t *kseq;
	multi_hat_infix_tel_s *ht;
} pldat_infix_tel_s;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_infix_tel_s *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
	infix_tel_buf_s *buf;
} stepdat_infix_tel_s;

static void worker_for_infix_tel_insert_list(void *data, long i, int tid) { // callback for kt_for()
	stepdat_infix_tel_s *s = (stepdat_infix_tel_s*)data;
	infix_tel_buf_s *buf = &s->buf[i];
	multi_hat_infix_tel_s *ht = s->p->ht;
	buf->n_ins += hat_insert_infix_tel(ht, buf->n, buf->a);
}

static void *worker_pipe_infix_tel(void *data, int step, void *in) { // callback for kt_pipeline()
	pldat_infix_tel_s *p = (pldat_infix_tel_s*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_infix_tel_s *s;
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
		stepdat_infix_tel_s *s = (stepdat_infix_tel_s*)in;
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
                count_seq_buf_can_infix_tel(s->buf, p->param->k, p->param->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        } else {
            for (i = 0; i < s->n; ++i) {
                count_seq_buf_infix_tel(s->buf, p->param->k, p->param->suf, s->len[i], s->seq[i]);
                free(s->seq[i]);
            }
        }
		free(s->seq); free(s->len);
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		stepdat_infix_tel_s *s = (stepdat_infix_tel_s*)in;
		int i, n = 1<<p->param->suf;
		uint64_t n_ins = 0;
        kt_for(p->param->n_thread, worker_for_infix_tel_insert_list, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			free(s->buf[i].a);
		}
		p->ht->tot += n_ins;
		free(s->buf);
        fprintf(stderr, "[M] processed %d sequences; %ld distinct %d-mers in the hash table\n", s->n, (long)p->ht->tot, p->param->k-1);
		free(s);
	}
	return 0;
}

multi_hat_infix_tel_s *count_infix_tel(const char *fn, const param_t *param, multi_hat_infix_tel_s *ht) {
	pldat_infix_tel_s pl;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.kseq = kseq_init(fp);
	pl.param = param;
    pl.ht = (!ht) ? infix_tel_init(param->k, param->suf) : ht;
	kt_pipeline(3, worker_pipe_infix_tel, &pl, 3);
	kseq_destroy(pl.kseq);
	gzclose(fp);
	return pl.ht;
}

multi_hat_infix_tel_s *count_infix_tel_file(const char *fn1, const param_t *param, multi_hat_infix_tel_s *h2) {
	multi_hat_infix_tel_s *ht;
	ht = count_infix_tel(fn1, param, h2); 
	return ht;
}

/*** generate histogram ***/
typedef struct {
	uint32_t *c;
} buf_hist_infix_tel_t;

typedef struct {
	const multi_hat_infix_tel_s *ht;
	buf_hist_infix_tel_t *cnt;
} hist_infix_tel_s;

static void worker_infix_tel_hist(void *data, long i, int tid) { // callback for kt_for()
	hist_infix_tel_s *a = (hist_infix_tel_s*)data;
	uint32_t *cnt = a->cnt[tid].c;
    //int k = a->ht->k;
	hat_infix_tel_t *g = a->ht->ht[i].ht;
	for (khint_t i = 0; i < kh_end(g); ++i) {
		if (kh_exist(g,i)) {
            //uint64_t infix = kh_key(g,i);
            for (int e = 0; e < 16; e++) {
                if(kh_val(g,i).edges[e]) {
                    uint64_t sigma = (uint64_t)(kh_val(g,i).sigma);
                    uint64_t edge = (uint64_t)(kh_val(g,i).edges[e]);
                    uint64_t idx = (sigma * (sigma-1))/2 + edge - 1;
                    cnt[idx]++;
                }
            }
        }
    }
}

void hist_infix_tel(const multi_hat_infix_tel_s *ht, uint32_t *cnt, int n_thread) {
	hist_infix_tel_s a;
	a.ht = ht;
	//memset(cnt, 0, (NUM_GENOMES+1) * sizeof(uint64_t));
	mcalloc(a.cnt, n_thread); // count is divide for number of threads
    int m = (NUM_GENOMES + 1) * NUM_GENOMES / 2;
    for (int i = 0; i < n_thread; i++) {
        mcalloc(a.cnt[i].c, m);
    }
	kt_for(n_thread, worker_infix_tel_hist, &a, 1<<ht->suf);
	//for (i = 0; i <= NUM_GENOMES; ++i) cnt[i] = 0; //This will be the total amount
	for (int j = 0; j < n_thread; ++j)
		for (int i = 0; i < m; ++i)
			cnt[i] += a.cnt[j].c[i]; //Combine all
	free(a.cnt);
}

void print_kmer_debug_infix_tel(multi_hat_infix_tel_s* ht) {
    int k = ht->k;
	for (int i = 0; i < 1<<ht->suf; ++i) {
        hat_infix_tel_t *g = ht->ht[i].ht;
        for (uint32_t j = 0; j < kh_end(g); ++j) {
            if (kh_exist(g,j)) {
                uint64_t infix = kh_key(g,j);
                for (int e = 0; e < 16; e++) {
                    if(kh_val(g,j).edges[e]) {
                        uint64_t left_kmer =  (((e & 12)) << (2*(k-1)-2)) | infix ; //12=0b1100
                        uint64_t right_kmer = (infix << 2) | (e & 3ULL);
                        printf("%s %d\n", bits2kmer(left_kmer, k), kh_val(g,j).edges[e]);
                        printf("%s %d\n", bits2kmer(right_kmer, k), kh_val(g,j).edges[e]);
                    }
                }
            }
        }
    }
}

void output_hist_infix_tel(int argc, char *argv[], param_t param, ketopt_t options){
    
	multi_hat_infix_tel_s *ht = 0;
	int i;

    // Calculate necessary bits
    if (param.filelist) {
        NUM_GENOMES = count_fasta(param.filelist);
    } else {
        NUM_GENOMES = argc - options.ind;
    }

    ID_GENOME = 1;
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
    fprintf(stderr, "Considering telomers as possible source of unitig break\n");

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
                ht = count_infix_tel_file(line, &param, ht);
            } else {
                ht = count_infix_tel_file(line, &param, 0);
            }
        }
        fclose(fp);
        if (line) free(line);
    } else {
        ht = count_infix_tel_file(argv[options.ind], &param, 0);
        for (i = 1; i < NUM_GENOMES; i++){
            ID_GENOME++;
            ht = count_infix_tel_file(argv[options.ind+i], &param, ht);
        }
    }
	fprintf(stderr, "[M::%s] %ld distinct %d-mers\n", __func__, (long)ht->tot, param.k-1);

	//print_infix_debug(ht);
	//print_kmer_debug(ht);

	uint32_t *cnt_infix;
    int m = (NUM_GENOMES + 1) * NUM_GENOMES / 2;
    mcalloc(cnt_infix, m);
	hist_infix_tel(ht, cnt_infix, param.n_thread);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%d\t%lld\n", i,(long long)cnt_infix[i]);
	for (i = 0; i < m; ++i) printf("%lld\n", (long long)cnt_infix[i]);
	infix_tel_destroy(ht);
}
