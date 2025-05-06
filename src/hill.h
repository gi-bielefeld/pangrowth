#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <cstdlib> // For atoi, exit
#include <cstring> // For strcmp (though parsing part is removed, ketopt might use it)
#include "Rmath.h"
#include "ketopt.h"

#define EPSILON 0.00001
//#define EPSILON numeric_limits<double>::epsilon()
using namespace std;

// Helper function to parse command line options for hill commands
// Returns o.ind (index of the first non-option argument) on success.
// Returns a negative value on error:
//  -1 for unknown option or missing argument.
//  -2 for invalid argument to -p.
// Error messages related to option parsing are printed by this helper.
static int parse_hill_cli_options(int argc, char *argv[], int *num_points_ptr, const char* cmd_name) {
    ketopt_t o = KETOPT_INIT;
    int c;
    *num_points_ptr = 30; // Default value for number of points

    // argv[0] is the command name (e.g., "hill"), options start from argv[1]
    while ((c = ketopt(&o, argc, argv, 1, "p:", 0)) >= 0) {
        if (c == 'p') {
            *num_points_ptr = atoi(o.arg);
            if (*num_points_ptr < 0) { // Allow 0 for "all points"
                fprintf(stderr, "Error: Number of points for -p must be non-negative.\n");
                fprintf(stderr, "Usage: pangrowth %s [-p <num_points>] ...\n", cmd_name);
                return -2; // Invalid argument error
            }
        } else if (c == '?') { // Unknown option
            fprintf(stderr, "Error: unknown option `-%c'.\n", o.optopt);
            fprintf(stderr, "Usage: pangrowth %s [-p <num_points>] ...\n", cmd_name);
            return -1; // Unknown option error
        } else if (c == ':') { // Missing argument for an option
            fprintf(stderr, "Error: option `-%c' requires an argument.\n", o.optopt);
            fprintf(stderr, "Usage: pangrowth %s [-p <num_points>] ...\n", cmd_name);
            return -1; // Missing argument error
        }
    }
    return o.ind; // Index of the first non-option argument
}

void read_file(const char* file, vector<double>& R){
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
        R.push_back(tmp);
        ss.clear();
    }
    cerr << "complete\n" << flush;
}

void get_points(int n, int num_points, vector<int> &points) {
    points.clear(); // Ensure points vector is empty before filling

    double end_val = 3.0 * n; // Use a distinct name for the double value
    int int_end_val = static_cast<int>(floor(end_val));

    if (n==1) { // Original check for n=1
        cerr << "WARN: number of genomes needs to be greater than 1\n" << flush;
        return; // Original behavior: return if n=1, points remains empty.
    }
    
    if (int_end_val <= 0) { // If max point (3*n) is not positive, no points to generate.
        return;
    }

    // If num_points effectively requests all points up to 3*n.
    // This handles the case where num_points (from CLI via output_hill*) is set to int_end_val for "-p 0"
    if (num_points >= int_end_val) {
        points.reserve(int_end_val);
        for (int i = 1; i <= int_end_val; i++) {
            points.push_back(i);
        }
        // The original warning about size mismatch at the end of this function will reflect this.
        return; // All points generated
    }

    // Original logic for num_points <= 3 (this block is now after the "all points" check)
    if (num_points <= 3) {
        points.push_back(1);
        points.push_back(n); 
        points.push_back(int_end_val); // Use the calculated int_end_val
    } else { // Original sampling logic for num_points > 3 and < int_end_val
      points.reserve(num_points);
      double step = end_val / num_points; // Use end_val (double) for step calculation

      for (double i = 1.0; i <= end_val; i+=step) {
        points.push_back(int(floor(i)));
    }


    for (int j = 1; j < num_points; j++) {
        if (points[j-1] < n && points[j] > n) {
            points[j-1] = n;
            break;
        }
    }

    points[num_points-1] = int(floor(end));

    if (points.size() != num_points) {
        cout << "WARN: outputting " << points.size() << " instead of " << num_points << '\n' << flush;
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
//    if (k < 0 || k > n) return R_NegInf; // log(C(n,k)) where C(n,k)=0
//    if (k == 0 || k == n) return 0.0;    // log(C(n,k)) where C(n,k)=1
//    if (k > n / 2) k = n - k;           // Symmetry: C(n, k) = C(n, n-k)
//    return lgammafn(n + 1.0) - lgammafn(k + 1.0) - lgammafn(n - k + 1.0);
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
    for (int sigma = 1; sigma <= n; sigma++) {
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
    p++;
    printf("obs\t%d\t", (int)n);
    for (int i = 1; i <= n; i++) {
        h_hat_unitig[i] = h_kmer[i] - h_unimer[i];
    }
    print_hill(m, h_hat_unitig);

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
    p++;
    printf("obs\t%d\t", (int)n);
    for (int i = 1; i <= n; i++) {
        h_hat_kmer[i] = h_kmer[i];
    }
    print_hill(m, h_hat_kmer);

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

void output_hill_cdbg(int argc, char *argv[]) {
    vector<double> h_kmer{0}, h_infix_eq{0}; //1-based index
    vector<int> points;
    int num_points = 30; // Default
    char *kmer_hist_filename = nullptr;
    char *infix_hist_filename = nullptr;
    // num_points will be set to default (30) by helper if -p not given, or to parsed value.

    int opt_idx = parse_hill_cli_options(argc, argv, &num_points, argv[0]);

    if (opt_idx < 0) { // Error in parsing options, message already printed by helper.
        return;       // Specific usage for this command is part of helper's message.
    }

    if (argc - opt_idx != 2) {
        fprintf(stderr, "Error: Expected 2 file arguments (k-mer histogram and infix histogram).\n");
        fprintf(stderr, "Usage: pangrowth %s [-p <num_points>] <hist_kmer> <hist_infix>\n", argv[0]);
        return;
    }
    kmer_hist_filename = argv[o.ind];
    infix_hist_filename = argv[o.ind + 1];

    read_file(kmer_hist_filename, h_kmer);
    read_file(infix_hist_filename, h_infix_eq);
    int n = h_kmer.size()-1;

    int final_num_points;
    if (num_points == 0) { // User specified -p 0, meaning all points up to 3*n
        final_num_points = (n > 0) ? static_cast<int>(floor(3.0 * n)) : 1; // If n=0, 3*n=0, so 1 point.
        // Ensure at least 1 point if 3*n was < 1 (e.g. n is very small positive, like 0.1 if n could be double)
        // Since n is int from h_kmer.size()-1, if n>0, 3*n >= 3. So floor(3.0*n) >=3.
        // If n=0, final_num_points = 1.
        // If n is from hist, n >=0. If hist empty, n=-1. If hist has 1 el, n=0.
        // Let's assume n is number of genomes, so n >= 1 for meaningful analysis.
        // If n (from hist size) is 0 or less, default to 1 point.
        if (n <= 0) final_num_points = 1;
        else final_num_points = static_cast<int>(floor(3.0 * n));

        if (final_num_points <= 0) final_num_points = 1; // Ensure at least 1 point.
    } else {
        final_num_points = num_points;
        if (n > 0) { // Cap user-provided num_points if it's too large, only if n is positive
            int cap = static_cast<int>(floor(3.0 * n));
            if (cap <= 0) cap = 1; // Cap must be at least 1
            final_num_points = min(cap, final_num_points);
        } else { // if n <= 0
             final_num_points = 1; // If no genomes, just 1 point.
        }
    }
    if (final_num_points <= 0) final_num_points = 1; // Overall safety: ensure at least 1 point

    get_points(n, final_num_points, points);
    hill_cdbg(h_kmer, h_infix_eq, points);
}

void output_hill(int argc, char *argv[]) {
    vector<double> h_kmer {0}; //1-based index
    vector<int> points;
    int num_points = 30; // Will be updated by helper, or remains default.
    char *kmer_hist_filename = nullptr;

    int opt_idx = parse_hill_cli_options(argc, argv, &num_points, argv[0]);

    if (opt_idx < 0) { // Error in parsing options, message already printed by helper.
        return;       // Specific usage for this command is part of helper's message.
    }

    if (argc - opt_idx != 1) {
        fprintf(stderr, "Error: Expected 1 file argument (k-mer histogram).\n");
        fprintf(stderr, "Usage: pangrowth %s [-p <num_points>] <hist_kmer>\n", argv[0]);
        return;
    }
    kmer_hist_filename = argv[o.ind];

    read_file(kmer_hist_filename, h_kmer);
    int n = h_kmer.size()-1;

    int final_num_points;
    if (num_points == 0) { // User specified -p 0, meaning all points up to 3*n
        if (n <= 0) final_num_points = 1; // If no genomes, default to 1 point.
        else final_num_points = static_cast<int>(floor(3.0 * n));

        if (final_num_points <= 0) final_num_points = 1; // Ensure at least 1 point.
    } else {
        final_num_points = num_points;
        if (n > 0) { // Cap user-provided num_points if it's too large, only if n is positive
            int cap = static_cast<int>(floor(3.0 * n));
            if (cap <= 0) cap = 1; // Cap must be at least 1
            final_num_points = min(cap, final_num_points);
        } else { // if n <= 0
            final_num_points = 1; // If no genomes, just 1 point.
        }
    }
    if (final_num_points <= 0) final_num_points = 1; // Overall safety: ensure at least 1 point

    get_points(n, final_num_points, points);
    hill(h_kmer, points);
}
