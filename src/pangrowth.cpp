#include "size_pangenome.h"
#include "yak-hist.h"
#include "yak-infix.h"
#include "yak-infix-telomer.h"
#include "yak-core.h"
#include "hill.h"
#include "ggcat-infix.h"

#include <string>
#include <vector>
#include <unistd.h>

using namespace std;

typedef void (*pangrowth_cmd_fn)(int, char**);

static bool streq(const char *a, const char *b) {
    return strcmp(a, b) == 0;
}

static void print_hist_usage() {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  pangrowth hist [options] <in.fa> [in.fa] > hist.txt\n");
    fprintf(stderr, "  pangrowth hist --ggcat [options] <ggcat_k_graph.fa> > hist.txt\n");
    fprintf(stderr, "  pangrowth hist --cdbg [-o <prefix>] [options] <in.fa> [in.fa]\n");
    fprintf(stderr, "  pangrowth hist --cdbg --ggcat [-o <prefix>] [options] <ggcat_k_graph.fa>\n");
    fprintf(stderr, "Modes:\n");
    fprintf(stderr, "  --cdbg     output both k-mer and infix-equivalent histograms\n");
    fprintf(stderr, "  --ggcat    read a colored ggcat k-mer graph instead of FASTA files\n");
    fprintf(stderr, "  -o STR     cdbg output prefix [out]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -k INT     k-mer size [17]\n");
    fprintf(stderr, "  -t INT     number of worker threads [4]\n");
    fprintf(stderr, "FASTA options:\n");
    fprintf(stderr, "  -i PATH    file containing one FASTA path per line\n");
    fprintf(stderr, "  -b         turn off canonical k-mers\n");
    fprintf(stderr, "  -s INT     suffix size for k-mer [4]\n");
    fprintf(stderr, "  -c INT     minimum k-mer count within a file [1]\n");
    fprintf(stderr, "  -T         account for telomeres in cdbg mode\n");
    fprintf(stderr, "cdbg output: <prefix>_hist.txt and <prefix>_hist_infix.txt\n");
}

static vector<char*> filtered_hist_argv(int argc, char **argv, bool for_kmer) {
    vector<char*> out;
    out.push_back(argv[0]);
    for (int i = 1; i < argc; ++i) {
        if (streq(argv[i], "--cdbg") || streq(argv[i], "--ggcat")) continue;
        if (streq(argv[i], "-o") || streq(argv[i], "--out-prefix")) {
            ++i;
            continue;
        }
        if (strncmp(argv[i], "-o", 2) == 0 && argv[i][2] != '\0') continue;
        if (strncmp(argv[i], "--out-prefix=", 13) == 0) continue;
        if (for_kmer && (streq(argv[i], "-T") || streq(argv[i], "-D"))) continue;
        out.push_back(argv[i]);
    }
    return out;
}

static void run_command_to_file(pangrowth_cmd_fn fn,
                                vector<char*>& args,
                                const string& path) {
    fflush(stdout);
    int saved_stdout = dup(fileno(stdout));
    if (saved_stdout < 0) {
        fprintf(stderr, "Could not save stdout\n");
        exit(EXIT_FAILURE);
    }

    FILE *out = fopen(path.c_str(), "w");
    if (out == 0) {
        fprintf(stderr, "Could not open output file: %s\n", path.c_str());
        close(saved_stdout);
        exit(EXIT_FAILURE);
    }
    if (dup2(fileno(out), fileno(stdout)) < 0) {
        fprintf(stderr, "Could not redirect stdout to: %s\n", path.c_str());
        fclose(out);
        close(saved_stdout);
        exit(EXIT_FAILURE);
    }

    fn((int)args.size(), args.data());
    fflush(stdout);

    if (dup2(saved_stdout, fileno(stdout)) < 0) {
        fprintf(stderr, "Could not restore stdout\n");
        fclose(out);
        close(saved_stdout);
        exit(EXIT_FAILURE);
    }
    close(saved_stdout);
    fclose(out);
}

static void output_hist(int argc, char **argv) {
    if (argc == 1 || (argc == 2 && (streq(argv[1], "-h") || streq(argv[1], "--help")))) {
        print_hist_usage();
        return;
    }

    bool cdbg = false;
    bool ggcat = false;
    bool telomeres = false;
    const char *out_prefix = "out";
    bool has_out_prefix = false;

    for (int i = 1; i < argc; ++i) {
        if (streq(argv[i], "--cdbg")) {
            cdbg = true;
        } else if (streq(argv[i], "--ggcat")) {
            ggcat = true;
        } else if (streq(argv[i], "-T")) {
            telomeres = true;
        } else if (streq(argv[i], "-o") || streq(argv[i], "--out-prefix")) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -o requires an output prefix.\n");
                print_hist_usage();
                exit(EXIT_FAILURE);
            }
            out_prefix = argv[++i];
            has_out_prefix = true;
        } else if (strncmp(argv[i], "-o", 2) == 0 && argv[i][2] != '\0') {
            out_prefix = argv[i] + 2;
            has_out_prefix = true;
        } else if (strncmp(argv[i], "--out-prefix=", 13) == 0) {
            out_prefix = argv[i] + 13;
            has_out_prefix = true;
        }
    }

    if (has_out_prefix && !cdbg) {
        fprintf(stderr, "Error: -o is only used with hist --cdbg.\n");
        print_hist_usage();
        exit(EXIT_FAILURE);
    }
    if (ggcat && telomeres) {
        fprintf(stderr, "Error: -T is not yet supported with --ggcat.\n");
        exit(EXIT_FAILURE);
    }
    if (!cdbg && !ggcat) {
        output_hist_kmer(argc, argv);
        return;
    }

    vector<char*> args = filtered_hist_argv(argc, argv, false);
    if (ggcat && !cdbg) {
        output_hist_ggcat((int)args.size(), args.data());
        return;
    }

    if (ggcat) {
        output_hist_cdbg_ggcat((int)args.size(), args.data(), out_prefix);
        return;
    }

    vector<char*> kmer_args = filtered_hist_argv(argc, argv, true);
    vector<char*> infix_args = filtered_hist_argv(argc, argv, false);
    run_command_to_file(output_hist_kmer, kmer_args, string(out_prefix) + "_hist.txt");
    run_command_to_file(output_hist_infix, infix_args, string(out_prefix) + "_hist_infix.txt");
}

int main(int argc, char** argv) {
    if (argc == 1) {
        fprintf(stderr, "Usage: pangrowth <cmd> [options] arg...\n");
        fprintf(stderr, "Where <cmd> is one of: hist, growth, hill or core\n");
        fprintf(stderr, "  hist:\n");
        fprintf(stderr, "       Outputs k-mer histograms; use --cdbg for k-mer and infix histograms\n");
        fprintf(stderr, "  kmer_core:\n");
        fprintf(stderr, "       Outputs the pangenome average core from the fasta files\n");
        fprintf(stderr, "  growth:\n");
        fprintf(stderr, "       Outputs the pangenome growth graph from the histogram (or panmatrix)\n");
        fprintf(stderr, "  hill:\n");
        fprintf(stderr, "       Outputs Hill's numbers for richness, exp. entropy and inv. gini-simpson index.\n"
                        "       It takes one k-mer histogram, or k-mer and infix histograms for cdbg.\n"
                        "       Options: -p <int>  Number of points to sample (default: 30)\n");
        fprintf(stderr, "  core:\n");
        fprintf(stderr, "       Outputs the pangenome core graph from the histogram (or panmatrix)\n");
        fprintf(stderr, "Legacy aliases: hist_infix, hist_infix_ggcat, hill_cdbg\n");
        return 0;
    }

    if (strcmp(argv[1],"hist") == 0) {  
        output_hist(argc-1, argv+1);
    } else if(strcmp(argv[1], "kmer_core") == 0) { 
        output_kmer_core(argc-1, argv+1);
    } else if(strcmp(argv[1], "hist_infix") == 0) { 
        output_hist_infix(argc-1, argv+1);
    } else if(strcmp(argv[1], "hist_infix_ggcat") == 0) { 
        output_hist_infix_ggcat(argc-1, argv+1);
    } else if(strcmp(argv[1], "growth") == 0) { 
        output_pangenome(argc-1, argv+1);
    } else if(strcmp(argv[1], "hill") == 0) { 
        output_hill(argc-1, argv+1);
    } else if(strcmp(argv[1], "hill_cdbg") == 0) { 
        output_hill_cdbg(argc-1, argv+1);
    } else if(strcmp(argv[1], "core") == 0) { 
        output_core(argc-1, argv+1);
    } else {
        fprintf(stderr, "Unknown command: %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    return 0;
}
