#include "size_pangenome.h"
#include "yak-hist.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc == 1) {
        fprintf(stderr, "Usage: pangrowth <cmd> [options] arg...\n");
        fprintf(stderr, "Where <cmd> is one of: hist, growth or core\n");
        fprintf(stderr, "  hist:\n");
        fprintf(stderr, "       Outputs the histogram of k-mers from a list of fasta files\n");
        fprintf(stderr, "  growth:\n");
        fprintf(stderr, "       Outputs the pangenome growth graph from the histogram (or panmatrix)\n");
        fprintf(stderr, "  core:\n");
        fprintf(stderr, "       Outputs the pangenome core graph from the histogram (or panmatrix)\n");
        return 0;
    }

    if (strcmp(argv[1],"hist") == 0) {  
        output_hist(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "growth") == 0) { 
        output_pangenome(argc-1, argv+sizeof(char));
    } else if(strcmp(argv[1], "core") == 0) { 
        output_core(argc-1, argv+sizeof(char));
    }

    return 0;
}
