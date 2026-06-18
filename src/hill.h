#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <fstream> // std::ifstream
#include <cstdlib> // For atoi, exit
#include <limits>
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
    // cdbg interpolation only
    bool use_bernoulli_unimer = false;
    int bernoulli_exact_tail = 0;
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

double static inline
log_choose_fast(const vector<double>& log_fact, int n, int k) {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    return log_fact[n] - log_fact[k] - log_fact[n-k];
}

vector<double> static inline
make_log_factorials(int n) {
    vector<double> log_fact(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) log_fact[i] = log_fact[i - 1] + log((double)i);
    return log_fact;
}

void static inline
log_add_exp(double& acc, double x) {
    if (std::isinf(x) && x < 0) return;
    if (std::isinf(acc) && acc < 0) {
        acc = x;
    } else if (acc < x) {
        acc = x + log1p(exp(acc - x));
    } else {
        acc = acc + log1p(exp(x - acc));
    }
}

double static inline
est_h_unimer_exact_i(vector<double> &hbar,
                    int i,
                    int n,
                    int m,
                    const vector<double>& log_fact,
                    double lchoose_nm) {
    double log_total = -std::numeric_limits<double>::infinity();
    int max_sigma = min(n, n - m + i);

    for (int sigma = i; sigma <= max_sigma; ++sigma) {
        double lchoose_n_sigma_m_i = log_choose_fast(log_fact, n - sigma, m - i);
        for (int j = i; j <= sigma; ++j) {
            double h = hbar[(sigma * (sigma - 1) / 2) + j];
            if (h <= 0.0) continue;
            double log_term = log(h)
                + log_choose_fast(log_fact, j, i)
                + lchoose_n_sigma_m_i
                - lchoose_nm;
            log_add_exp(log_total, log_term);
        }
    }

    if (std::isinf(log_total) && log_total < 0) return 0.0;
    return exp(log_total);
}

void static inline
est_h_unimer_bernoulli_hybrid(vector<double> &hbar,
                              int n,
                              int m,
                              int exact_tail,
                              const vector<double>& log_fact,
                              vector<double>& h_hat_unimer) {
    fill(h_hat_unimer.begin(), h_hat_unimer.end(), 0.0);

    double p = (double)m / (double)n;
    double log_p = log(p);
    double log_q = log1p(-p);
    vector<double> log_A(n + 1, -std::numeric_limits<double>::infinity());

    for (int j = 1; j <= n; ++j) {
        double log_sum = -std::numeric_limits<double>::infinity();
        for (int sigma = j; sigma <= n; ++sigma) {
            double h = hbar[(sigma * (sigma - 1) / 2) + j];
            if (h <= 0.0) continue;
            log_add_exp(log_sum, log(h) + (double)(sigma - j) * log_q);
        }
        log_A[j] = log_sum;
    }

    for (int i = 1; i <= m; ++i) {
        double log_total = -std::numeric_limits<double>::infinity();
        for (int j = i; j <= n; ++j) {
            if (std::isinf(log_A[j]) && log_A[j] < 0) continue;
            double log_term = log_A[j]
                + log_choose_fast(log_fact, j, i)
                + (double)i * log_p
                + (double)(j - i) * log_q;
            log_add_exp(log_total, log_term);
        }
        if (!(std::isinf(log_total) && log_total < 0)) h_hat_unimer[i] = exp(log_total);
    }

    if (exact_tail > 0) {
        int start_i = max(1, m - exact_tail + 1);
        double lchoose_nm = log_choose_fast(log_fact, n, m);
        for (int i = start_i; i <= m; ++i) {
            h_hat_unimer[i] = est_h_unimer_exact_i(hbar, i, n, m, log_fact, lchoose_nm);
        }
    }
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
hill_cdbg(vector<double> &h_kmer,
          vector<double> &h_infix_eq,
          vector<int> &points,
          int bernoulli_exact_tail = -1) {
    int n = (int)h_kmer.size()-1;
    int m = points[points.size()-1];
    size_t p = 0;
    bool use_bernoulli_unimer = bernoulli_exact_tail >= 0;
    double sum_h = 0;
    for (int i = 1; i <= n; ++i) {
        sum_h += h_kmer[i];
    }
    //1-based index
    vector<double> h_hat_unitig(m+1);
    vector<double> h_hat_unimer(m+1);
    vector<double> h_unimer(n+1);
    vector<double> log_fact;
    if (use_bernoulli_unimer) log_fact = make_log_factorials(n);
    for (int i = 1; i <= n; i++) h_unimer[i] = h_infix_eq[((i+1)*i/2)];

    printf("fit\tm\trichness\texp_entropy\tinv_gini_simp\n");
    //** Interpolation **//
    while (p < points.size() && points[p] < n) {
        int m = points[p++];
        
        double lchoose_nm = use_bernoulli_unimer ? log_choose_fast(log_fact, n, m) : lchoose(n, m);
        if (use_bernoulli_unimer) {
            est_h_unimer_bernoulli_hybrid(h_infix_eq, n, m, bernoulli_exact_tail, log_fact, h_hat_unimer);
        }

        for (int i = 1; i <= m; i++) {
            double h_unimer_i = use_bernoulli_unimer ? h_hat_unimer[i] : est_h_unimer(h_infix_eq, i, n, m);
            h_hat_unitig[i] = est_h(h_kmer, i, n, m, lchoose_nm) - h_unimer_i;
        }
        printf("int\t%d\t", m);
        print_hill(m, h_hat_unitig);
    }

    //** Observed **//
    if (p < points.size() && points[p] == n) {
        p++;
        printf("obs\t%d\t", (int)n);
        for (int i = 1; i <= n; i++) {
            h_hat_unitig[i] = h_kmer[i] - h_unimer[i];
        }
        print_hill(n, h_hat_unitig);
    }

    //** Extrapolation **//
    // Unseen k-mers
    double Q1 = h_kmer[1];
    double Q2 = h_kmer[2];
    double Q0hat;
    if (Q2 == 0) {
        Q0hat = ((n - 1)/(double)n) * Q1*(Q1 - 1)/2;
    } else {
        Q0hat = ((n - 1)/(double)n) * Q1*Q1/(2*Q2);
    }
    // Unseen uni-mers
    double Q1_unimer = h_infix_eq[1];
    double Q2_unimer = h_infix_eq[3];
    double Q0hat_unimer;
    if (Q2_unimer == 0) {
        Q0hat_unimer = ((n - 1)/(double)n) * Q1_unimer*(Q1_unimer - 1)/2;
    } else {
        Q0hat_unimer = ((n - 1)/(double)n) * Q1_unimer*Q1_unimer/(2*Q2_unimer);
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
    while ((c = ketopt(&o, argc, argv, 1, "p:f:B:", 0)) >= 0) {
        if (c == 'p') params.num_points = atoi(o.arg); 
        if (c == 'f') {
            params.points_file = o.arg;
            params.use_points_file = true;
        }
        if (c == 'B') {
            params.use_bernoulli_unimer = true;
            params.bernoulli_exact_tail = atoi(o.arg);
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

    if (params.use_bernoulli_unimer) {
        fprintf(stderr, "Error: -B is only supported by hill_cdbg.\n");
        params.err = true;
    }

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

    if (params.bernoulli_exact_tail < 0) {
        fprintf(stderr, "Error: -B must be >= 0.\n");
        params.err = true;
    }

    //input_hist_kmer
    if (argc - params.o_index != 2) {
        fprintf(stderr, "Error: Expected 2 files argument (k-mer and infix).\n");
        params.err = true;
    }

    //helper
    if (params.err) {
        fprintf(stderr, "Usage: pangrowth %s [-p <num_points>, -f <file_points>, -B <exact_tail>] <hist_kmer> <hist_infix>\n", argv[0]);
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
    if (params.use_bernoulli_unimer) {
        hill_cdbg(h_kmer, h_infix_eq, params.points, params.bernoulli_exact_tail);
    } else {
        hill_cdbg(h_kmer, h_infix_eq, params.points);
    }
}
