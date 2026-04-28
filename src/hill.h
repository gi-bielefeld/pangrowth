#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <fstream> // std::ifstream
#include <cstdlib> // For atoi, exit
#include "Rmath.h"
#include "ketopt.h"

#define EPSILON 0.000001
//#define EPSILON numeric_limits<double>::epsilon()
using namespace std;

struct param_hill_t{ 
    bool err = false;
    int o_index;
    //points
    int num_points = 30;
    vector<int> points;
    char *points_file;
    bool use_points_file = false;
    //genomes
    uint32_t num_genomes = -1;
    char *kmer_hist_filename;
};

struct param_hill_cdbg_t : param_hill_t { 
    char *infix_hist_filename;
};

vector<double> read_histogram_1based_index(const char* file){
    vector<double> H {0}; //1-based index

    cerr << "Reading file " << file << ".."<< flush;
    ifstream fin(file);
    if (!fin.is_open()) {
        cerr << "Error opening file: " << file << '\n';
        exit(EXIT_FAILURE);
    }

    string line;
    stringstream ss; 
    while (getline(fin, line, '\n')) {
        double tmp;
        ss.str(line);
        ss >> tmp;
        H.push_back(tmp);
        ss.clear();
    }
    cerr << "complete\n" << flush;
    return H;
}

void update_params(param_hill_t &params, vector<double> &h_kmer) {
    uint32_t num_genomes = h_kmer.size()-1;
    params.num_genomes = num_genomes;
    //if 0 return all the points
    if (params.num_points == 0) { params.num_points = 3*num_genomes; }
    params.num_points = min(int(3*num_genomes), params.num_points);
}

void generate_sample_points(param_hill_t &params) {
    uint32_t n = params.num_genomes;
    int num_points = params.num_points;

    if (n==1) {
        cerr << "Warn: number of genomes needs to be greater than 1\n" << flush;
        return;
    }

    double end = 3.0 * n;

    if (num_points <= 3) {
        cerr << "Warn: number of points needs to be at least 3. Set to 3.\n" << flush;
        params.points.push_back(1);
        params.points.push_back(n);
        params.points.push_back(int(end));
        return;
    }

    if (num_points > end) {
        params.points.reserve(end);
        for (int i = 1; i <= end; i++) {
            params.points.push_back(i);
        }
        return;
    }

    params.points.reserve(num_points);
    double step = end / num_points;

    for (double i = 1.0; i <= end; i+=step) {
        params.points.push_back(int(floor(i)));
    }


    for (int j = 1; j < num_points; j++) {
        if (params.points[j-1] < n && params.points[j] > n) {
            params.points[j-1] = n;
            break;
        }
    }

    if (!params.points.empty()) {
        params.points.back() = int(floor(end));
    }

    if (params.points.size() != num_points) {
        cout << "Warn: outputting " << params.points.size() << " instead of " << num_points << '\n' << flush;
    }
}

void sample_points(param_hill_t &params) {
    if (params.use_points_file) {
        // Ignore -p if a points file is provided
        std::ifstream fin(params.points_file);
        if (!fin) {
            cerr << "Error: Cannot open points file: " << params.points_file << '\n';
            return;
        }
        int pt;
        while (fin >> pt) {
            if (pt > 0) params.points.push_back(pt);
        }
        if (params.points.empty()) {
            cerr << "Error: No valid points found in file.\n";
            return;
        }
        sort(params.points.begin(), params.points.end());
    } else {
        generate_sample_points(params);
    }
}

double static inline 
richness(int m, vector<double> &h) {
    double richness = 0;
    for (int i = 1; i <= m; i++) {
        richness += h[i];
    }
    return richness;
}

double static inline
exp_entropy(int m, vector<double> &h) {
    double entropy = 0;
    double U_hat_t = 0;
    for (int i = 1; i <= m; i++) {
        U_hat_t += i*h[i];
    }
    for (int i = 1; i <= m; i++) {
        double x = (double)i/(double)U_hat_t;
        entropy -= x*log(x)*h[i];
    }
    return exp(entropy);
}

double static inline
inv_gini_simpson(int m, vector<double> &h){
    double gini_simpson = 0;
    double U_hat_t = 0;
    for (int i = 1; i <= m; i++) {
        U_hat_t += i*h[i];
    }
    for (int i = 1; i <= m; i++) {
        double x = (double)i/(double)U_hat_t;
        gini_simpson += x*x*h[i];
    }
    return 1/gini_simpson;
}

//double static inline 
//lchoose(int n, int k) {
//    //if (k < 0 || k > n) return R_NegInf; // log(C(n,k)) where C(n,k)=0
//    if (k == 0 || k == n) return 0.0;    // log(C(n,k)) where C(n,k)=1
//    if (k > n / 2) k = n - k;           // Symmetry: C(n, k) = C(n, n-k)
//    return lgammaf(n + 1.0) - lgammaf(k + 1.0) - lgammaf(n - k + 1.0);
//}

double static inline 
est_h(vector<double>& h, int i, int n, int m, double lchoose_nm) {
    double tot = 0;

    for (int j = i; j <= n - m + i; j++) {
        double log_hj = log(h[j]);
        double term = exp(log_hj + lchoose(j, i) + lchoose(n - j, m - i) - lchoose_nm);
        tot += term;
    }

    return tot;
}

void static inline 
est_h_hill(vector<double>& h, int n, int m, double lchoose_nm, vector<double>& h_hat_kmer) {
    vector<double> log_h(n+1);
    for (int i = 1; i <= n; ++i) log_h[i] = log(h[i]);

    double lfact_i = 0;
    for (int i = 1; i <= m; i++) {
        //fprintf(stderr,"%d\n",i);
        h_hat_kmer[i] = 0;
        //double lfact_i = lgamma(i + 1);
        double lfact_m_i = lgamma(m - i + 1);
        //lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
        double lfact_j = lfact_i;
        double lfact_j_i = 0;
        for (int j = i; j <= n - m + i; j++) {
            double lchoose_j_i = lfact_j - lfact_i - lfact_j_i;

            double lfact_n_j = lgamma(n - j + 1);
            double lfact_n_j_m_i = lgamma(n - j - m + i + 1);
            double lchoose_n_j_m_i = lfact_n_j - lfact_m_i - lfact_n_j_m_i;

            h_hat_kmer[i] += exp(log_h[j] + lchoose_j_i + lchoose_n_j_m_i - lchoose_nm);
            lfact_j += log(j+1);
            lfact_j_i += log(j-i+1);
        }
        lfact_i += log(i+1);
    }
}

double static inline
est_h_unimer(vector<double> &hbar, int i, int n, int m) {
    double tot = 0;
    for (int sigma = i; sigma <= n; sigma++) {
        for (int j = i; j <= sigma; j++) {
            tot += exp(log(hbar[(sigma*(sigma-1)/2)+j]) + lchoose(j,i) + lchoose(n-sigma, m-i) - lchoose(n,m));
        }
    }
    return tot;
}

//double static inline 
//est_h_unimer(vector<double> &hbar, int i, int n, int m) {
//    double tot = 0;
//    for (int sigma = i; sigma < i-m+n; sigma++) {
//    //for (int sigma = i; sigma < n; sigma++) {
//        for (int j = i; j <= sigma; j++) {
//            tot += exp(log(hbar[(sigma*(sigma-1)/2)+j]) + lchoose(j,i) + lchoose(n-sigma, m-i) - lchoose(n,m));
//        }
//    }
//    return tot;
//}

double static inline 
est_h_extra(vector<double> &h, int i, int n, int x) {
    double tot = 0;
    double p_j;
    for (int j = max(1,i-x); j <= min(i,n); j++) {
        if (h[j] == 0) continue;
        //if (j == 0) p_j = 1.0/(double)(n+x);
        p_j = (double)j/(double)n;
        if (p_j < EPSILON) break;
        tot += ((double)h[j])*exp(lchoose(x,i-j) + (i-j)*log(p_j) + (x-i+j)*log(1-p_j));
    }
    return tot;
}

double static inline
est_h_unimer_extra(vector<double> &h, double pi_broken, int i, int n, int x) {
    double tot = 0;
    double p_j;
    for (int j = max(0,i-x); j <= min(i,n-1); j++) {
        if (h[j] == 0) continue;
        if (j == 0) p_j = 1.0/(double)(n+x);
        else p_j = (double)j/(double)n;
        if (p_j < EPSILON) break;
        tot += ((double)h[j])*exp(lchoose(x,i-j) + (i-j)*log(p_j) + (x-i+j)*log(1-p_j)) * (pow(1 - pi_broken, x));
    }
    return tot;
}

void static inline
print_hill(int m, vector<double> &h) {
    printf("%.2f\t", richness(m, h));
    fflush(stdout);
    printf("%.2f\t", exp_entropy(m, h));
    fflush(stdout);
    printf("%.2f\n", inv_gini_simpson(m, h));
    fflush(stdout);
}

void 
hill_cdbg(vector<double> &h_kmer, vector<double> &h_infix_eq, vector<int> &points) {
    double n = h_kmer.size()-1;
    int m = points[points.size()-1];
    size_t p = 0;
    double sum_h = 0;
    for (int i = 1; i <= n; ++i) {
        sum_h += h_kmer[i];
    }
    //1-based index
    vector<double> h_hat_unitig(m+1);
    vector<double> h_unimer(n+1);
    for (int i = 1; i <= n; i++) h_unimer[i] = h_infix_eq[((i+1)*i/2)];

    printf("fit\tm\trichness\texp_entropy\tinv_gini_simp\n");
    //** Interpolation **//
    while (points[p] < n) {
        int m = points[p++];
        
        double lchoose_nm = lchoose(n, m);
        for (int i = 1; i <= m; i++) {
            h_hat_unitig[i] = est_h(h_kmer, i, n, m, lchoose_nm) - est_h_unimer(h_infix_eq, i, n, m);
        }
        printf("int\t%d\t", m);
        print_hill(m, h_hat_unitig);
    }

    //** Observed **//
    if (points[p] == n) {
        p++;
        printf("obs\t%d\t", (int)n);
        for (int i = 1; i <= n; i++) {
            h_hat_unitig[i] = h_kmer[i] - h_unimer[i];
        }
        print_hill(m, h_hat_unitig);
    }

    //** Extrapolation **//
    // Unseen k-mers
    double Q1 = h_kmer[1];
    double Q2 = h_kmer[2];
    double Q0hat;
    if (Q2 == 0) {
        Q0hat = ((n - 1)/n) * Q1*(Q1 - 1)/2;
    } else {
        Q0hat = ((n - 1)/n) * Q1*Q1/(2*Q2);
    }
    // Unseen uni-mers
    double Q1_unimer = h_infix_eq[1];
    double Q2_unimer = h_infix_eq[3];
    double Q0hat_unimer;
    if (Q2_unimer == 0) {
        Q0hat_unimer = ((n - 1)/n) * Q1_unimer*(Q1_unimer - 1)/2;
    } else {
        Q0hat_unimer = ((n - 1)/n) * Q1_unimer*Q1_unimer/(2*Q2_unimer);
    }
    // Broken uni-mers
    double Q_unimer = 0;
    double Q_broken_unimer = 0;
    for (int i = 1; i < n; i++) {
        Q_unimer += h_infix_eq[((i+1)*i)/2];
        Q_broken_unimer += h_infix_eq[((i+1)*i)/2+1];
    }
    Q_unimer += h_infix_eq[h_infix_eq.size()-1];
    double pi_broken = Q_broken_unimer/(Q_broken_unimer+n*Q_unimer);

    while (p < points.size()) {
        int m = points[p++];
        int x = m-n;
        //printf("%d]",m);
        for (int i = 1; i < m; i++) {
            h_hat_unitig[i] = est_h_extra(h_kmer, i, n, x) - est_h_unimer_extra(h_unimer, pi_broken, i, n, x);
        }
        h_hat_unitig[1] += Q0hat * (1 - pow(1 - Q1/(Q1+n*Q0hat), x));
        h_hat_unitig[1] -= Q0hat_unimer * (1 - pow(1 - Q1_unimer/(Q1_unimer+n*Q0hat_unimer), x));
        h_hat_unitig[m] = h_kmer[n] - h_unimer[n];

        printf("ext\t%d\t", m);
        print_hill(m, h_hat_unitig);
    }
}

void 
hill(vector<double> &h_kmer, vector<int> &points) {
    double n = h_kmer.size()-1;
    int m = points[points.size()-1];
    size_t p = 0;
    double sum_h = 0;
    for (int i = 1; i <= n; ++i) {
        sum_h += h_kmer[i];
    }
    //1-based index
    vector<double> h_hat_kmer(m+1);

    printf("fit\tm\trichness\texp_entropy\tinv_gini_simp\n");
    //** Interpolation **//
    while (points[p] < n) {
        int m = points[p++];
        
        double lchoose_nm = lchoose(n, m);
        est_h_hill(h_kmer, n, m, lchoose_nm, h_hat_kmer);
        printf("int\t%d\t", m);
        print_hill(m, h_hat_kmer);
    }

    //** Observed **//
    if (points[p] == n) {
        p++;
        printf("obs\t%d\t", (int)n);
        for (int i = 1; i <= n; i++) {
            h_hat_kmer[i] = h_kmer[i];
        }
        print_hill(m, h_hat_kmer);
    }

    //** Extrapolation **//
    // Unseen k-mers
    double Q1 = h_kmer[1];
    double Q2 = h_kmer[2];
    double Q0hat;
    if (Q2 == 0) {
        Q0hat = ((n - 1)/n) * Q1*(Q1 - 1)/2;
    } else {
        Q0hat = ((n - 1)/n) * Q1*Q1/(2*Q2);
    }

    while (p < points.size()) {
        int m = points[p++];
        int x = m-n;
        //printf("%d]",m);
        for (int i = 1; i < m; i++) {
            h_hat_kmer[i] = est_h_extra(h_kmer, i, n, x);
        }
        h_hat_kmer[1] += Q0hat * (1 - pow(1 - Q1/(Q1+n*Q0hat), x));
        h_hat_kmer[m] = h_kmer[n];

        printf("ext\t%d\t", m);
        print_hill(m, h_hat_kmer);
    }
}

param_hill_t parse_common_cli(int argc, char *argv[]){
    param_hill_cdbg_t params;

    int c;
    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "p:f:", 0)) >= 0) {
        if (c == 'p') params.num_points = atoi(o.arg); 
        if (c == 'f') {
            params.points_file = o.arg;
            params.use_points_file = true;
        }
    }

    //num_points
    if (params.num_points < 0) {
        fprintf(stderr, "Error: Number of points for -p must be positive.\n");
        params.err = true;
    }

    if (params.num_points == 0) {
        fprintf(stderr, "Info: running on all points.\n");
    }

    params.o_index = o.ind;

    return params;
}

param_hill_t cli_hill(int argc, char *argv[]) {
    param_hill_t params = parse_common_cli(argc, argv);

    //input_hist_kmer
    if (argc - params.o_index != 1) {
        fprintf(stderr, "Error: Expected 1 file argument (k-mer histogram).\n");
        params.err = true;
    }

    //helper
    if (params.err) {
        fprintf(stderr, "Usage: pangrowth %s [-p <num_points>, -f <file_points>] <hist_kmer>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    params.kmer_hist_filename = argv[params.o_index];
    return params;
}

param_hill_cdbg_t cli_hill_cdbg(int argc, char *argv[]) {
    param_hill_cdbg_t params;
    static_cast<param_hill_t&>(params) =  parse_common_cli(argc, argv);

    //input_hist_kmer
    if (argc - params.o_index != 2) {
        fprintf(stderr, "Error: Expected 2 files argument (k-mer and infix).\n");
        params.err = true;
    }

    //helper
    if (params.err) {
        fprintf(stderr, "Usage: pangrowth %s [-p <num_points>, -f <file_points>] <hist_kmer> <hist_infix>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    params.kmer_hist_filename = argv[params.o_index];
    params.infix_hist_filename = argv[params.o_index+1];
    return params;
}

void output_hill(int argc, char *argv[]) {
    param_hill_t params = cli_hill(argc, argv);
    vector<double> h_kmer = read_histogram_1based_index(params.kmer_hist_filename);
    update_params(params, h_kmer);
    sample_points(params);
    hill(h_kmer, params.points);
}

void output_hill_cdbg(int argc, char *argv[]) {
    param_hill_cdbg_t params = cli_hill_cdbg(argc, argv);
    vector<double> h_kmer = read_histogram_1based_index(params.kmer_hist_filename);
    vector<double> h_infix_eq = read_histogram_1based_index(params.infix_hist_filename);
    update_params(params, h_kmer);
    sample_points(params);
    hill_cdbg(h_kmer, h_infix_eq, params.points);
}
