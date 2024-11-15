#ifndef kmer_H
#define kmer_H

#include <zlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include "kseq.h" // FASTA/Q parser
                  
KSEQ_INIT(gzFile, gzread)

#define mcalloc(pos, len) ((pos) = (__typeof__(pos))calloc((len), sizeof(*(pos))))
#define mmalloc(pos, len) ((pos) = (__typeof__(pos))malloc((len) * sizeof(*(pos))))
#define mrealloc(pos, len) ((pos) = (__typeof__(pos))realloc((pos), (len) * sizeof(*(pos))))

#define N_COUNTS     8192
uint32_t ID_GENOME;
int NUM_GENOMES;
int TOTAL_BITS;
int BITS_GENOME;
int MASK_COUNT; 
int MASK_GENOME;
int SUF;

//*** utils ***/
void chomp(char *str) {
    size_t len = strlen(str);

    if (len == 0) return;

    size_t last_idx = len - 1;
    if( str[last_idx] == '\n' ) {
        str[last_idx] = '\0';
    }
}

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

static inline uint64_t hash64(uint64_t key, uint64_t mask) { // invertible integer hash function
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

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

//*** hist ***/
typedef struct {
	int32_t k;
	int32_t suf;
	int32_t n_thread;
    bool canonical;
	int64_t chunk_size;
    char* filelist;
} param_t;

void param_init(param_t *options) {
	memset(options, 0, sizeof(param_t));
	options->k = 17;
	options->suf = 10;
	options->n_thread = 4;
	options->canonical = true;
	options->chunk_size = 10000000;
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

#endif
