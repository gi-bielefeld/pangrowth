#include "size_pangenome.h"
#include "yak-hist.h"
#include "yak-infix.h"
#include "hill.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc == 1) {
        fprintf(stderr, "Usage: pangrowth <cmd> [options] arg...\n");
        fprintf(stderr, "Where <cmd> is one of: hist, growth or core\n");
        fprintf(stderr, "  hist:\n");
        fprintf(stderr, "       Outputs the histogram of k-mers from a list of fasta files\n");
        fprintf(stderr, "  hist_infix:\n");
        fprintf(stderr, "       Outputs the histogram of infix equivalents (k+1)-mers from a list of fasta files\n");
        fprintf(stderr, "  growth:\n");
        fprintf(stderr, "       Outputs the pangenome growth graph from the histogram (or panmatrix)\n");
        fprintf(stderr, "  hill:\n");
        fprintf(stderr, "       Outputs Hill's numbers for richness, exp. entropy and inv. gini-simpson index.\n"
                        "       It takes as input an histogram produced by 'hist'\n");
        fprintf(stderr, "  hill_cdbg:\n");
        fprintf(stderr, "       Outputs Hill's numbers for richness, exp. entropy and inv. gini-simpson index\n"
                        "       of the respective compacted de Bruijn graph. It takes as input two histograms\n"
                        "       representing the k-mer hist and the infix equivalents hist respectively\n");
        fprintf(stderr, "  core:\n");
        fprintf(stderr, "       Outputs the pangenome core graph from the histogram (or panmatrix)\n");
        return 0;
    }

    if (strcmp(argv[1],"hist") == 0) {  
        output_hist_kmer(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "hist_infix") == 0) { 
        output_hist_infix(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "growth") == 0) { 
        output_pangenome(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "hill") == 0) { 
        output_hill(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "hill_cdbg") == 0) { 
        output_hill_cdbg(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "core") == 0) { 
        output_core(argc-1, argv+sizeof(char));
    }

    return 0;
}
