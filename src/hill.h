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
    bool force_exact_cdbg = false;
    int bernoulli_exact_tail = 0;
    bool bernoulli_tail_set = false;
    bool use_adaptive_bernoulli_tail = false;
    double bernoulli_tail_tolerance = 0.001;
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

double static inline
est_h_fast(vector<double>& h,
           int i,
           int n,
           int m,
           double lchoose_nm,
           const vector<double>& log_fact) {
    double total = 0.0;
    for (int j = i; j <= n - m + i; ++j) {
        if (h[j] <= 0.0) continue;
        total += h[j] * exp(log_choose_fast(log_fact, j, i)
            + log_choose_fast(log_fact, n - j, m - i)
            - lchoose_nm);
    }
    return total;
}

vector<double> static inline
make_log_factorials(int n) {
    vector<double> log_fact(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) log_fact[i] = log_fact[i - 1] + log((double)i);
    return log_fact;
}

void static inline
add_binomial_mixture(vector<double>& out,
                     double weight,
                     int j,
                     int m,
                     double p,
                     double log_p,
                     double log_q,
                     double odds,
                     double inverse_odds,
                     const vector<double>& log_fact) {
    if (weight <= 0.0) return;
    int upper = min(j, m);
    int start = (int)floor((j + 1) * p);
    if (start < 1) start = 1;
    if (start > upper) start = upper;

    double prob = exp(log_choose_fast(log_fact, j, start)
        + (double)start * log_p
        + (double)(j - start) * log_q);
    out[start] += weight * prob;

    double p_i = prob;
    for (int i = start - 1; i >= 1; --i) {
        p_i *= ((double)(i + 1) / (double)(j - i)) * inverse_odds;
        out[i] += weight * p_i;
    }

    p_i = prob;
    for (int i = start + 1; i <= upper; ++i) {
        p_i *= ((double)(j - i + 1) / (double)i) * odds;
        out[i] += weight * p_i;
    }
}

void static inline
est_h_kmer_bernoulli_hybrid(vector<double>& h,
                            int n,
                            int m,
                            int exact_tail,
                            const vector<double>& log_fact,
                            vector<double>& h_hat) {
    fill(h_hat.begin(), h_hat.begin() + m + 1, 0.0);

    double p = (double)m / (double)n;
    double log_p = log(p);
    double log_q = log1p(-p);
    double q = 1.0 - p;
    double odds = p / q;
    double inverse_odds = q / p;

    for (int j = 1; j <= n; ++j) {
        add_binomial_mixture(h_hat, h[j], j, m, p, log_p, log_q, odds, inverse_odds, log_fact);
    }

    if (exact_tail > 0) {
        int start_i = max(1, m - exact_tail + 1);
        double lchoose_nm = log_choose_fast(log_fact, n, m);
        for (int i = start_i; i <= m; ++i) {
            h_hat[i] = est_h(h, i, n, m, lchoose_nm);
        }
    }
}

void static inline
est_h_unimer_bernoulli(vector<double> &hbar,
                       int n,
                       int m,
                       const vector<double>& log_fact,
                       vector<double>& h_hat_unimer) {
    fill(h_hat_unimer.begin(), h_hat_unimer.begin() + m + 1, 0.0);

    double p = (double)m / (double)n;
    double log_p = log(p);
    double log_q = log1p(-p);
    double q = 1.0 - p;
    double odds = p / q;
    double inverse_odds = q / p;
    vector<double> A(n + 1, 0.0);

    // A[j] = sum_sigma h(j,sigma) * (1-p)^(sigma-j)
    for (int j = 1; j <= n; ++j) {
        double q_pow = 1.0;
        double sum = 0.0;
        for (int sigma = j; sigma <= n; ++sigma) {
            double h = hbar[(sigma * (sigma - 1) / 2) + j];
            if (h > 0.0) sum += h * q_pow;
            q_pow *= q;
        }
        A[j] = sum;
    }

    // Add A[j] * Binomial(j,p)[i] for all needed i using recurrences
    // around the mode. This avoids O(N^2) log-sum-exp work.
    for (int j = 1; j <= n; ++j) {
        add_binomial_mixture(h_hat_unimer, A[j], j, m, p, log_p, log_q, odds, inverse_odds, log_fact);
    }
}

void static inline
exact_unimer_aggregate_at_depth(vector<double>& hbar,
                                int n,
                                int d,
                                vector<double>& aggregate) {
    fill(aggregate.begin(), aggregate.end(), 0.0);

    // Normalize by C(n-j,d). The remaining factor is a hypergeometric
    // probability, so neither intermediate can overflow.
    for (int j = 1; j <= n - d; ++j) {
        double weight = 1.0;
        double total = 0.0;
        for (int sigma = j; sigma <= n - d; ++sigma) {
            double h = hbar[(size_t)sigma * (sigma - 1) / 2 + j];
            if (h > 0.0) total += h * weight;
            if (sigma < n - d) {
                weight *= (double)(n - sigma - d) / (double)(n - sigma);
            }
        }
        aggregate[j] = total;
    }
}

double static inline
exact_unimer_at_depth(const vector<double>& aggregate,
                      int n,
                      int m,
                      int d,
                      double lchoose_nm,
                      const vector<double>& log_fact) {
    int i = m - d;
    double total = 0.0;
    for (int j = i; j <= n - d; ++j) {
        if (aggregate[j] <= 0.0) continue;
        total += aggregate[j] * exp(log_choose_fast(log_fact, j, i)
            + log_choose_fast(log_fact, n - j, d)
            - lchoose_nm);
    }
    return total;
}

vector<double> static inline
precompute_exact_unimer_tail(vector<double>& hbar,
                             int n,
                             const vector<int>& points,
                             size_t interpolation_points,
                             int tail_width,
                             const vector<double>& log_fact) {
    vector<double> exact(interpolation_points * (size_t)tail_width, 0.0);
    vector<double> aggregate(n + 1, 0.0);

    for (int d = 0; d < tail_width; ++d) {
        exact_unimer_aggregate_at_depth(hbar, n, d, aggregate);
        for (size_t q = 0; q < interpolation_points; ++q) {
            int m = points[q];
            if (d >= m) continue;
            double lchoose_nm = log_choose_fast(log_fact, n, m);
            exact[q * (size_t)tail_width + d] = exact_unimer_at_depth(
                aggregate, n, m, d, lchoose_nm, log_fact);
        }
    }

    return exact;
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

struct hill_sums_t {
    double richness = 0.0;
    double incidence = 0.0;
    double incidence_log_frequency = 0.0;
    double squared_frequency = 0.0;
};

void static inline
add_hill_bin(hill_sums_t& sums, int i, double value) {
    sums.richness += value;
    sums.incidence += (double)i * value;
    if (i > 1) {
        sums.incidence_log_frequency += (double)i * log((double)i) * value;
    }
    sums.squared_frequency += (double)i * (double)i * value;
}

void static inline
print_hill_sums(const hill_sums_t& sums) {
    double entropy = log(sums.incidence)
        - sums.incidence_log_frequency / sums.incidence;
    double inverse_gini_simpson = sums.incidence * sums.incidence
        / sums.squared_frequency;
    printf("%.2f\t%.2f\t%.2f\n",
           sums.richness, exp(entropy), inverse_gini_simpson);
    fflush(stdout);
}

struct adaptive_tail_point_t {
    int m = 0;
    double lchoose_nm = 0.0;
    int tail_width = 0;
    hill_sums_t sums;
    int similar_bins = 0;
    int exact_bins = 0;
    int negative_bins = 0;
    double negative_mass = 0.0;
    bool active = true;
    bool converged = false;
};

void static inline
replace_adaptive_bin(adaptive_tail_point_t& point,
                     int i,
                     double bernoulli,
                     double exact) {
    add_hill_bin(point.sums, i, exact - bernoulli);
    if (bernoulli < -EPSILON) {
        --point.negative_bins;
        point.negative_mass += bernoulli;
    }
    if (exact < -EPSILON) {
        ++point.negative_bins;
        point.negative_mass -= exact;
    }
}

bool static inline
adaptive_hill_values_are_close(const hill_sums_t& before,
                               const hill_sums_t& after,
                               double tolerance) {
    double before_entropy = exp(log(before.incidence)
        - before.incidence_log_frequency / before.incidence);
    double after_entropy = exp(log(after.incidence)
        - after.incidence_log_frequency / after.incidence);
    double before_gini = before.incidence * before.incidence
        / before.squared_frequency;
    double after_gini = after.incidence * after.incidence
        / after.squared_frequency;

    double richness_scale = max(1.0, max(abs(before.richness), abs(after.richness)));
    double entropy_scale = max(1.0, max(abs(before_entropy), abs(after_entropy)));
    double gini_scale = max(1.0, max(abs(before_gini), abs(after_gini)));
    return abs(after.richness - before.richness) <= tolerance * richness_scale
        && abs(after_entropy - before_entropy) <= tolerance * entropy_scale
        && abs(after_gini - before_gini) <= tolerance * gini_scale;
}

vector<adaptive_tail_point_t> static inline
adaptive_bernoulli_interpolation(vector<double>& h_kmer,
                                 vector<double>& h_infix_eq,
                                 int n,
                                 const vector<int>& points,
                                 size_t interpolation_points,
                                 int max_tail,
                                 double tolerance,
                                 const vector<double>& log_fact) {
    const int required_similar_bins = 3;
    vector<adaptive_tail_point_t> states(interpolation_points);
    vector<double> bernoulli_tails(
        interpolation_points * (size_t)max_tail, 0.0);
    vector<double> h_hat_kmer(n + 1, 0.0);
    vector<double> h_hat_unimer(n + 1, 0.0);

    for (size_t q = 0; q < interpolation_points; ++q) {
        adaptive_tail_point_t& state = states[q];
        state.m = points[q];
        state.lchoose_nm = log_choose_fast(log_fact, n, state.m);
        state.tail_width = min(max_tail, state.m);
        state.active = state.tail_width > 0;

        est_h_kmer_bernoulli_hybrid(
            h_kmer, n, state.m, 0, log_fact, h_hat_kmer);
        est_h_unimer_bernoulli(
            h_infix_eq, n, state.m, log_fact, h_hat_unimer);

        for (int i = 1; i <= state.m; ++i) {
            double unitig = h_hat_kmer[i] - h_hat_unimer[i];
            add_hill_bin(state.sums, i, unitig);
            if (unitig < -EPSILON) {
                ++state.negative_bins;
                state.negative_mass -= unitig;
            }
            int d = state.m - i;
            if (d < state.tail_width) {
                bernoulli_tails[q * (size_t)max_tail + d] = unitig;
            }
        }
    }

    size_t active_points = 0;
    for (size_t q = 0; q < states.size(); ++q) {
        if (states[q].active) ++active_points;
    }

    vector<double> aggregate(n + 1, 0.0);
    for (int d = 0; d < max_tail && active_points > 0; ++d) {
        bool any_at_depth = false;
        for (size_t q = 0; q < states.size(); ++q) {
            if (states[q].active && d < states[q].tail_width) {
                any_at_depth = true;
                break;
            }
        }
        if (!any_at_depth) break;

        exact_unimer_aggregate_at_depth(h_infix_eq, n, d, aggregate);
        for (size_t q = 0; q < states.size(); ++q) {
            adaptive_tail_point_t& state = states[q];
            if (!state.active || d >= state.tail_width) continue;

            int i = state.m - d;
            double exact_kmer = est_h_fast(
                h_kmer, i, n, state.m, state.lchoose_nm, log_fact);
            double exact_unimer = exact_unimer_at_depth(
                aggregate, n, state.m, d, state.lchoose_nm, log_fact);
            double exact_unitig = exact_kmer - exact_unimer;
            double bernoulli_unitig = bernoulli_tails[
                q * (size_t)max_tail + d];
            hill_sums_t before = state.sums;
            replace_adaptive_bin(state, i, bernoulli_unitig, exact_unitig);
            ++state.exact_bins;

            if (adaptive_hill_values_are_close(before, state.sums, tolerance)) {
                ++state.similar_bins;
            } else {
                state.similar_bins = 0;
            }

            if (state.similar_bins >= required_similar_bins) {
                state.active = false;
                state.converged = true;
                --active_points;
            } else if (state.exact_bins == state.tail_width) {
                state.active = false;
                state.converged = state.exact_bins == state.m;
                --active_points;
            }
        }
    }

    if (states.empty()) return states;

    int min_used = max_tail;
    int max_used = 0;
    int reached_cap = 0;
    double total_used = 0.0;
    for (size_t q = 0; q < states.size(); ++q) {
        min_used = min(min_used, states[q].exact_bins);
        max_used = max(max_used, states[q].exact_bins);
        total_used += states[q].exact_bins;
        if (!states[q].converged) ++reached_cap;
    }
    double mean_used = total_used / (double)states.size();
    fprintf(stderr,
            "Info: adaptive Bernoulli tail used %d-%d exact bins (mean %.2f); %d/%zu points reached their cap.\n",
            min_used, max_used, mean_used, reached_cap, states.size());

    return states;
}

void 
hill_cdbg(vector<double> &h_kmer,
          vector<double> &h_infix_eq,
          vector<int> &points,
          int bernoulli_exact_tail = -1,
          double adaptive_tolerance = -1.0) {
    int n = (int)h_kmer.size()-1;
    int m = points[points.size()-1];
    size_t p = 0;
    bool use_bernoulli_unimer = bernoulli_exact_tail >= 0;
    bool use_adaptive_tail = use_bernoulli_unimer && adaptive_tolerance > 0.0;
    double sum_h = 0;
    for (int i = 1; i <= n; ++i) {
        sum_h += h_kmer[i];
    }
    //1-based index
    vector<double> h_hat_unitig(m+1);
    vector<double> h_hat_kmer(m+1);
    vector<double> h_hat_unimer(m+1);
    vector<double> h_unimer(n+1);
    vector<double> log_fact;
    if (use_bernoulli_unimer) log_fact = make_log_factorials(n);
    for (int i = 1; i <= n; i++) h_unimer[i] = h_infix_eq[((i+1)*i/2)];

    size_t interpolation_points = 0;
    while (interpolation_points < points.size() && points[interpolation_points] < n) {
        ++interpolation_points;
    }
    int exact_tail_width = use_bernoulli_unimer ? min(bernoulli_exact_tail, n) : 0;
    vector<double> exact_unimer_tail;
    if (exact_tail_width > 0 && !use_adaptive_tail) {
        exact_unimer_tail = precompute_exact_unimer_tail(
            h_infix_eq, n, points, interpolation_points, exact_tail_width, log_fact);
    }

    printf("fit\tm\trichness\texp_entropy\tinv_gini_simp\n");
    //** Interpolation **//
    if (use_adaptive_tail) {
        vector<adaptive_tail_point_t> states = adaptive_bernoulli_interpolation(
            h_kmer, h_infix_eq, n, points, interpolation_points,
            exact_tail_width, adaptive_tolerance, log_fact);
        for (size_t q = 0; q < states.size(); ++q) {
            adaptive_tail_point_t& state = states[q];
            if (state.negative_bins > 0) {
                fprintf(stderr,
                        "Warn: Bernoulli cdbg interpolation has %d negative bins at m=%d (negative mass %.6g); Hill numbers may be invalid.\n",
                        state.negative_bins, state.m, state.negative_mass);
            }
            printf("int\t%d\t", state.m);
            print_hill_sums(state.sums);
        }
        p = interpolation_points;
    }

    while (!use_adaptive_tail && p < points.size() && points[p] < n) {
        size_t point_index = p;
        int m = points[p++];
        
        double lchoose_nm = use_bernoulli_unimer ? log_choose_fast(log_fact, n, m) : lchoose(n, m);
        if (use_bernoulli_unimer) {
            est_h_kmer_bernoulli_hybrid(h_kmer, n, m, bernoulli_exact_tail, log_fact, h_hat_kmer);
            est_h_unimer_bernoulli(h_infix_eq, n, m, log_fact, h_hat_unimer);
            for (int d = 0; d < exact_tail_width && d < m; ++d) {
                h_hat_unimer[m - d] = exact_unimer_tail[
                    point_index * (size_t)exact_tail_width + d];
            }
        }

        for (int i = 1; i <= m; i++) {
            if (use_bernoulli_unimer) {
                h_hat_unitig[i] = h_hat_kmer[i] - h_hat_unimer[i];
            } else {
                h_hat_unitig[i] = est_h(h_kmer, i, n, m, lchoose_nm) - est_h_unimer(h_infix_eq, i, n, m);
            }
        }
        if (use_bernoulli_unimer) {
            int negative_bins = 0;
            double negative_mass = 0.0;
            for (int i = 1; i <= m; ++i) {
                if (h_hat_unitig[i] < -EPSILON) {
                    ++negative_bins;
                    negative_mass += -h_hat_unitig[i];
                }
            }
            if (negative_bins > 0) {
                fprintf(stderr,
                        "Warn: Bernoulli cdbg interpolation has %d negative bins at m=%d (negative mass %.6g); Hill numbers may be invalid.\n",
                        negative_bins, m, negative_mass);
            }
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
    while ((c = ketopt(&o, argc, argv, 1, "p:f:B:A:E", 0)) >= 0) {
        if (c == 'p') params.num_points = atoi(o.arg); 
        if (c == 'f') {
            params.points_file = o.arg;
            params.use_points_file = true;
        }
        if (c == 'B') {
            params.use_bernoulli_unimer = true;
            params.bernoulli_tail_set = true;
            params.bernoulli_exact_tail = atoi(o.arg);
        }
        if (c == 'A') {
            char *end = NULL;
            params.use_bernoulli_unimer = true;
            params.use_adaptive_bernoulli_tail = true;
            params.bernoulli_tail_tolerance = strtod(o.arg, &end);
            if (end == o.arg || *end != '\0') {
                params.err = true;
            }
        }
        if (c == 'E') params.force_exact_cdbg = true;
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

void print_hill_usage() {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  pangrowth hill [options] <hist_kmer>\n");
    fprintf(stderr, "  pangrowth hill [options] <hist_kmer> <hist_infix>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -p INT     number of sample points; 0 outputs all points [30]\n");
    fprintf(stderr, "  -f PATH    file containing one sample point per line\n");
    fprintf(stderr, "cdbg interpolation options:\n");
    fprintf(stderr, "  -B INT     maximum exact right-tail bins [5; 20 with -A]\n");
    fprintf(stderr, "  -A FLOAT   stop after 3 tail corrections each change all Hill numbers by <= FLOAT\n");
    fprintf(stderr, "  -E         use exact interpolation instead of the Bernoulli hybrid\n");
}

void validate_hill_params(param_hill_t &params, int argc, char *argv[]) {
    if (params.use_bernoulli_unimer || params.force_exact_cdbg) {
        fprintf(stderr, "Error: -B and -E require two input histograms for cdbg mode.\n");
        params.err = true;
    }

    if (argc - params.o_index != 1) {
        fprintf(stderr, "Error: Expected 1 file argument (k-mer histogram).\n");
        params.err = true;
    }

    //helper
    if (params.err) {
        print_hill_usage();
        exit(EXIT_FAILURE);
    }

    params.kmer_hist_filename = argv[params.o_index];
}

param_hill_t cli_hill(int argc, char *argv[]) {
    param_hill_t params = parse_common_cli(argc, argv);
    validate_hill_params(params, argc, argv);
    return params;
}

void validate_hill_cdbg_params(param_hill_cdbg_t &params, int argc, char *argv[]) {
    if (params.force_exact_cdbg && params.use_bernoulli_unimer) {
        fprintf(stderr, "Error: -E cannot be used with -B or -A.\n");
        params.err = true;
    }
    if (!params.force_exact_cdbg && !params.use_bernoulli_unimer) {
        params.use_bernoulli_unimer = true;
        params.bernoulli_exact_tail = 5;
    }
    if (params.use_adaptive_bernoulli_tail && !params.bernoulli_tail_set) {
        params.bernoulli_exact_tail = 20;
    }
    if (params.bernoulli_exact_tail < 0) {
        fprintf(stderr, "Error: -B must be >= 0.\n");
        params.err = true;
    }
    if (params.use_adaptive_bernoulli_tail &&
        (!(params.bernoulli_tail_tolerance > 0.0) ||
         !(params.bernoulli_tail_tolerance < 1.0))) {
        fprintf(stderr, "Error: -A must be > 0 and < 1.\n");
        params.err = true;
    }
    if (params.use_adaptive_bernoulli_tail && params.bernoulli_exact_tail == 0) {
        fprintf(stderr, "Error: adaptive interpolation requires -B > 0.\n");
        params.err = true;
    }

    //input_hist_kmer
    if (argc - params.o_index != 2) {
        fprintf(stderr, "Error: Expected 2 files argument (k-mer and infix).\n");
        params.err = true;
    }

    //helper
    if (params.err) {
        print_hill_usage();
        exit(EXIT_FAILURE);
    }

    params.kmer_hist_filename = argv[params.o_index];
    params.infix_hist_filename = argv[params.o_index+1];
}

param_hill_cdbg_t cli_hill_cdbg(int argc, char *argv[]) {
    param_hill_cdbg_t params;
    static_cast<param_hill_t&>(params) =  parse_common_cli(argc, argv);
    validate_hill_cdbg_params(params, argc, argv);
    return params;
}

void run_hill_cdbg(param_hill_cdbg_t &params) {
    vector<double> h_kmer = read_histogram_1based_index(params.kmer_hist_filename);
    vector<double> h_infix_eq = read_histogram_1based_index(params.infix_hist_filename);
    update_params(params, h_kmer);
    sample_points(params);
    if (params.use_bernoulli_unimer) {
        double adaptive_tolerance = params.use_adaptive_bernoulli_tail
            ? params.bernoulli_tail_tolerance : -1.0;
        hill_cdbg(h_kmer, h_infix_eq, params.points,
                  params.bernoulli_exact_tail, adaptive_tolerance);
    } else {
        hill_cdbg(h_kmer, h_infix_eq, params.points);
    }
}

void output_hill(int argc, char *argv[]) {
    if (argc == 1 || (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))) {
        print_hill_usage();
        return;
    }

    param_hill_t params = parse_common_cli(argc, argv);
    int n_files = argc - params.o_index;
    if (n_files == 2) {
        param_hill_cdbg_t cdbg_params;
        static_cast<param_hill_t&>(cdbg_params) = params;
        validate_hill_cdbg_params(cdbg_params, argc, argv);
        run_hill_cdbg(cdbg_params);
        return;
    }

    validate_hill_params(params, argc, argv);
    vector<double> h_kmer = read_histogram_1based_index(params.kmer_hist_filename);
    update_params(params, h_kmer);
    sample_points(params);
    hill(h_kmer, params.points);
}

void output_hill_cdbg(int argc, char *argv[]) {
    param_hill_cdbg_t params = cli_hill_cdbg(argc, argv);
    run_hill_cdbg(params);
}
