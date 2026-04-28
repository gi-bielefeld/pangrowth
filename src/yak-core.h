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

#define kh_hash_generic(x) ((x))
KHASHL_SET_INIT(, kmer_core_t, ct, uint64_t, kh_hash_generic, kh_eq_generic);

kmer_core_t *extract_core_kmers(multi_hat_kmer_s *h, uint64_t abs_quorum) {
    kmer_core_t *core_ht = ct_init();

    for (uint64_t suf = 0; suf < (1ULL << h->suf); ++suf) {
        hat_kmer_s *g = &h->h[suf];
        for (khint_t i = 0; i != kh_end(g->h); ++i) {
            if (kh_exist(g->h, i)) {
                uint64_t kmer = kh_key(g->h, i);
                uint64_t count = kh_val(g->h, i) & MASK_COUNT;

                if (count >= abs_quorum) {
                    int absent;
                    ct_put(core_ht, kmer, &absent);
                }
            }
        }
    }

    return core_ht;
}

static inline uint64_t stream_kmer_core(int k, int suf, int len, const char *seq, uint64_t abs_quorum, kmer_core_t *core_ht){ 
    uint64_t core = 0;
    uint8_t overlap = 0;
	uint64_t i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				uint64_t y = x[0] < x[1]? x[0] : x[1];      // Canonicals!
                if (ct_get(core_ht, y) != kh_end(core_ht)) {
                    core += k - overlap;
                    overlap = k-1;
                } else {
                    overlap = max(0, overlap-1);
                }
			}
		} else l = 0, x[0] = x[1] = 0, overlap = 0;//, first_core = true; // if there is an "N", restart
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


uint64_t count_core_nucleotide(const char *fn, const param_t *opt, kmer_core_t *core_ht, uint64_t abs_quorum) {
    uint64_t core_nt = 0;
    gzFile fp = gzopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", fn);
        exit(1);
    }

    kseq_t* seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        const char* s = seq->seq.s;
        int len = seq->seq.l;
        core_nt += stream_kmer_core(opt->k, opt->suf, len, s, abs_quorum, core_ht);
    }

    kseq_destroy(seq);
    gzclose(fp);
    return core_nt;
}

uint64_t count_kmer_core_file(const char *fn1, const param_t *opt, kmer_core_t *core_ht, uint64_t abs_quorum) {
	uint64_t core_kmer_count = count_core_nucleotide(fn1, opt, core_ht, abs_quorum); 
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
    MASK_COUNT  = (1ULL << BITS_GENOME) - 1;
    MASK_GENOME = MASK_COUNT << BITS_GENOME;
    MASK_INTRA  = 7ULL << (2 * BITS_GENOME);
    TOTAL_BITS  = BITS_GENOME*2;
    MIN_COUNT   = opt.min_count;
    SUF = opt.suf;

    uint64_t abs_quorum = std::ceil(((double)NUM_GENOMES)*quorum);
    fprintf(stderr, "Number of genomes:\t%d\n", NUM_GENOMES);
    fprintf(stderr, "Number of threads:\t%d\n", opt.n_thread);
    fprintf(stderr, "Bits per genome:\t%d\n", BITS_GENOME);
    fprintf(stderr, "Bits per suffix:\t%d\n", SUF);
    fprintf(stderr, "Counting %s k-mers\n", opt.canonical ? "canonical" : "forward");
    fprintf(stderr, "Quorum %0.2f [%lu/%d]\n", quorum, (unsigned long)abs_quorum, NUM_GENOMES);//, NUM_GENOMES);

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
            if (h) ID_GENOME++;
            char *token = strtok(line, ",");
            while (token != NULL) {
                if (!h) h = count_kmer_file(token, &opt, 0);
                else    h = count_kmer_file(token, &opt, h);
                token = strtok(NULL, ",");
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
    //
	fprintf(stderr, "[M::%s] extracting core k-mers\n", __func__);
    kmer_core_t *core_ht = extract_core_kmers(h, abs_quorum);
    //print_kmer_debug(h);
    ch_destroy(h);
	fprintf(stderr, "[M::%s] counting core nucleotides\n", __func__);

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
            char *token = strtok(line, ",");
            while (token != NULL) {
                core_count = count_kmer_core_file(token, &opt, core_ht, abs_quorum);
                token = strtok(NULL, ",");
            }
            ID_GENOME++;
        }
        fclose(fp);
        if (line) free(line);
    } else {
        for (i = 0; i < NUM_GENOMES; i++){
            core_count = count_kmer_core_file(argv[o.ind+i], &opt, core_ht, abs_quorum);
            ID_GENOME++;
        }
    }
}
