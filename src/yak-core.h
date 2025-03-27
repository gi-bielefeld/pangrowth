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
#include "yak-hist.h"

inline bool is_core(uint64_t kmer, uint64_t abs_quorum, multi_hat_kmer_s *h) {
    int mask = (1<<h->suf) - 1;
	hat_kmer_s *g = &h->h[kmer&mask];
    khint_t k = ht_get(g->h, kmer);
    //printf("C:%s %d %d \n",bits2kmer(kmer, 11), (kh_val(g->h, k)&MASK_COUNT), ((kh_val(g->h, k)&MASK_COUNT)) >= abs_quorum);
    return ((kh_val(g->h, k)&MASK_COUNT) >= abs_quorum);
}

static inline uint64_t stream_kmer_core(int k, int suf, int len, const char *seq, uint64_t abs_quorum, multi_hat_kmer_s *h){ 
    uint64_t core = 0;
	uint64_t i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    uint64_t first_core = true;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];      // Canonicals!
                if (is_core(y, abs_quorum, h)) {
                    if(first_core) core += k-1; //if it is the start add a plus k-1
                    core++;
                    first_core = false;
                } else {
                    first_core = true;
                }
			}
		} else l = 0, x[0] = x[1] = 0, first_core = true; // if there is an "N", restart
	}
    return core;
}

typedef struct { // global data structure for kt_pipeline()
	const param_t *opt;
	kseq_t *ks;
	multi_hat_kmer_s *h;
    uint64_t abs_quorum;
    uint64_t core_count;
} pldat_kmer_core_s;

typedef struct { // data structure for each step in kt_pipeline()
	pldat_kmer_core_s *p;
	int n, m, sum_len, nk;
	int *len;
	char **seq;
} stepdat_kmer_core_s;

static void *worker_pipe_kmer_core(void *data, int step, void *in) { // callback for kt_pipeline()
	pldat_kmer_core_s *p = (pldat_kmer_core_s*)data;
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		stepdat_kmer_core_s *s;
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
	} else if (step == 1) { // step 2: count core k-mers
		stepdat_kmer_core_s *s = (stepdat_kmer_core_s*)in;
		int i, n = 1<<p->opt->suf, m;
        if (p->opt->canonical) {
            for (i = 0; i < s->n; ++i) {
                uint64_t found_core = stream_kmer_core(p->opt->k, p->opt->suf, s->len[i], s->seq[i], s->p->abs_quorum, s->p->h);
                free(s->seq[i]);
		        p->core_count += found_core;
            }
        } 
        //else { //non-canonical
        //    for (i = 0; i < s->n; ++i) {
        //        count_seq_buf(s->buf, p->opt->k, p->opt->suf, s->len[i], s->seq[i]);
        //        free(s->seq[i]);
        //    }
        //}
        //fprintf(stderr, "[M] processed %d sequences; %ld core k-mers\n", s->n, p->core_count);
		free(s->seq); free(s->len);
		free(s);
	}
	return 0;
}

uint64_t count_kmer_core(const char *fn, const param_t *opt, multi_hat_kmer_s *h, double quorum) {
	pldat_kmer_core_s pl;
	gzFile fp;
    pl.core_count = 0;
    pl.abs_quorum = std::ceil(NUM_GENOMES*quorum);
	if ((fp = gzopen(fn, "r")) == 0) return 0;
	pl.ks = kseq_init(fp);
	pl.opt = opt;
    pl.h = h;
	kt_pipeline(2, worker_pipe_kmer_core, &pl, 2);
	kseq_destroy(pl.ks);
	gzclose(fp);
    return pl.core_count;
}

uint64_t count_kmer_core_file(const char *fn1, const param_t *opt, multi_hat_kmer_s *h, double quorum) {
	uint64_t core_kmer_count = count_kmer_core(fn1, opt, h, quorum); 
    printf("%s\t%ld\n", fn1, core_kmer_count);
    return core_kmer_count;
}

void output_kmer_core(int argc, char *argv[]){
	multi_hat_kmer_s *h = 0;
	int c;
	int i;
	param_t opt;
	ketopt_t o = KETOPT_INIT;
    double quorum = 0.9;
	param_init(&opt);

	while ((c = ketopt(&o, argc, argv, 1, "k:s:K:t:q:i:b", 0)) >= 0) {
		if (c == 'k') opt.k = atoi(o.arg);
		else if (c == 's') opt.suf = atoi(o.arg);
		else if (c == 'K') opt.chunk_size = atoi(o.arg);
		else if (c == 't') opt.n_thread = atoi(o.arg);
		else if (c == 'q') quorum = atof(o.arg);
		else if (c == 'i') opt.filelist = o.arg;
		else if (c == 'b') opt.canonical = false;
	}

	if (argc - o.ind < 1 && !opt.filelist) {
		fprintf(stderr, "Usage: pangrowth kmer_core [options] <in.fa> [in.fa]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT     k-mer size [%d]\n", opt.k);
		fprintf(stderr, "  -t INT     number of worker threads [%d]\n", opt.n_thread);
		fprintf(stderr, "  -i PATH    file containing a list of fasta files on each line\n");
		fprintf(stderr, "  -b         turn off transformation into canonical [%d]\n", opt.canonical);
		fprintf(stderr, "  -s INT     suffix size for k-mer [%d]\n", opt.n_thread);
		fprintf(stderr, "  -q FLOAT   quorum percentage between [0,1] for a k-mer to be core [%0.2f]\n", quorum);
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
    fprintf(stderr, "Quorum %0.2f\n", quorum);

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

    // Second phase: Core
    uint64_t core_count = 0;
    ID_GENOME = 1;
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
            core_count = count_kmer_core_file(line, &opt, h, quorum);
            ID_GENOME++;
        }
        fclose(fp);
        if (line) free(line);
    } else {
        for (i = 0; i < NUM_GENOMES; i++){
            core_count = count_kmer_core_file(argv[o.ind+i], &opt, h, quorum);
            ID_GENOME++;
        }
    }

    //print_kmer_debug(h);

	//hist_kmer(h, cnt, opt.n_thread);
	////for (i = 1; i <= NUM_GENOMES; ++i) printf("%d\t%lld\n", i,(long long)cnt[i]);
	//for (i = 1; i <= NUM_GENOMES; ++i) printf("%lld\n", (long long)cnt[i]);
	ch_destroy(h);
}
