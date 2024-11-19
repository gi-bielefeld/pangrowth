#include <iostream>
#include <vector>
#include <string>
#include "Rmath.h"
//cout << lchoose(10,0) << '\n' << flush; // 0
//cout << lchoose(0,0) << '\n' << flush;  // 0
//cout << lchoose(0,10) << '\n' << flush; // -inf
//#define EPSILON 0.0000000001
#define EPSILON numeric_limits<double>::epsilon()
using namespace std;

void read_file(const char* file, std::vector<double>& R){
    std::cerr << "Reading file " << file << ".."<< std::flush;
    std::ifstream fin(file);
    if (!fin.is_open()) {
        std::cerr << "Error opening file: " << file << '\n';
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    std::stringstream ss; 
    while (std::getline(fin, line, '\n')) {
        double tmp;
        ss.str(line);
        ss >> tmp;
        R.push_back(tmp);
        ss.clear();
    }
    std::cerr << "complete\n" << std::flush;
}

void get_points(int n, int m, vector<int> &points) {
    double start = 1.0;
    double end = 2.0 * n;
    double step = (end - start) / (m - 1);

    for (int i = 0; i < m; ++i) {
        points.push_back(start + i * step);
        if (points[i-1] < n && points[i] > n) points[i] = n;
    }
}

double static inline 
richness(int m, vector<double> &h_hat_kmer, vector<double> &h_hat_unimer) {
    double richness = 0;
    for (int i = 1; i <= m; i++) {
        richness += h_hat_kmer[i];
        richness -= h_hat_unimer[i];
    }
    return richness;
}

double static inline
exp_entropy(int m, vector<double> &h_hat_kmer, vector<double> &h_hat_unimer) {
    double entropy = 0;
    double U_hat_t = 0;
    for (int i = 1; i <= m; i++) {
        U_hat_t += i*(h_hat_kmer[i]-h_hat_unimer[i]);
    }
    for (int i = 1; i <= m; i++) {
        double x = (double)i/(double)U_hat_t;
        entropy -= x*log(x)*(h_hat_kmer[i]-h_hat_unimer[i]);
    }
    return exp(entropy);
}

double static inline
inv_gini_simpson(int m, vector<double> &h_hat_kmer, vector<double> &h_hat_unimer){
    double gini_simpson = 0;
    double U_hat_t = 0;
    for (int i = 1; i <= m; i++) {
        U_hat_t += i*(h_hat_kmer[i]-h_hat_unimer[i]);
    }
    for (int i = 1; i <= m; i++) {
        double x = (double)i/(double)U_hat_t;
        gini_simpson += x*x*(h_hat_kmer[i]-h_hat_unimer[i]);
    }
    return 1/gini_simpson;
}

double static inline 
est_h(vector<double> &h,int i,int n, int m) {
    double tot = 0;
    for (int j = i; j <= n-m+i; j++) {
        tot += exp(log((double)h[j]) + lchoose(j,i) + lchoose(n-j,m-i) - lchoose(n,m));
    }
    return tot;
}

double static inline 
est_h_unimer(vector<double> &hbar, int i, int n, int m) {
    double tot = 0;
    //for (int sigma = 1; sigma <= n; sigma++) {
    //untested, in case of error, go back to sigma=1
    for (int sigma = i; sigma <= n; sigma++) {
        for (int j = i; j <= sigma; j++) {
            tot += exp(log(hbar[(sigma*(sigma-1)/2)+j]) + lchoose(j,i) + lchoose(n-sigma, m-i) - lchoose(n,m));
        }
    }
    return tot;
}

double static inline 
est_h_extra(vector<double> &h, int i, int n, int x) {
    double tot = 0;
    double p_j;
    for (int j = max(0,i-x); j <= min(i,n); j++) {
        if (h[j] == 0) continue;
        if (j == 0) p_j = 1.0/(double)(n+x);
        else p_j = (double)j/(double)n;
        if (p_j < EPSILON) break;
        tot += ((double)h[j])*exp(lchoose(x,i-j) + (i-j)*log(p_j) + (x-i+j)*log(1-p_j));
    }
    return tot;
}

double static inline 
est_h_unimer_extra(vector<double> &h, vector<double> &h_infix_eq, double pi_broken, int i, int n, int x) {
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
print_hill(int m, vector<double> &h_hat_kmer, vector<double> &h_hat_unimer) {
    printf("%.2f\t", richness(m, h_hat_kmer, h_hat_unimer));
    fflush(stdout);
    printf("%.2f\t", exp_entropy(m, h_hat_kmer, h_hat_unimer));
    fflush(stdout);
    printf("%.2f\n", inv_gini_simpson(m, h_hat_kmer, h_hat_unimer));
}

void 
hill_cdbg(vector<double> &h_kmer, vector<double> &h_infix_eq, vector<int> &points) {
    double n = h_kmer.size()-1;
    int m = points[points.size()-1];
    int p = 0;
    //1-based index
    vector<double> h_hat_kmer(m+1);
    vector<double> h_hat_unimer(m+1);
    vector<double> h_unimer(n+1);
    for (int i = 1; i <= n; i++) h_unimer[i] = h_infix_eq[((i+1)*i/2)];

    printf("fit\tm\trichness\texp_entropy\tinv_gini_simp\n");
    //** Interpolation **//
    while (points[p] < n) {
        int m = points[p++];
        
        for (int i = 1; i <= m; i++) {
            h_hat_kmer[i] = est_h(h_kmer, i, n, m);
            h_hat_unimer[i] = est_h_unimer(h_infix_eq, i, n, m);
        }
        printf("int\t%d\t", m);
        print_hill(m, h_hat_kmer, h_hat_unimer);
    }

    //** Observed **//
    p++;
    printf("obs\t%d\t", (int)n);
    print_hill(n, h_kmer, h_unimer);

    //** Extrapolation **//
    // Unseen k-mers
    double Q1 = h_kmer[1];
    double Q2 = h_kmer[2];
    if (Q2 == 0) {
        h_kmer[0] = ((n - 1)/n) * Q1*(Q1 - 1)/2;
    } else {
        h_kmer[0] = ((n - 1)/n) * Q1*Q1/(2*Q2);
    }
    // Unseen uni-mers
    double Q1_infix = h_infix_eq[1];
    double Q2_infix = h_infix_eq[3];
    if (Q2_infix == 0) {
        h_unimer[0] = ((n - 1)/n) * Q1_infix*(Q1_infix - 1)/2;
    } else {
        h_unimer[0] = ((n - 1)/n) * Q1_infix*Q1_infix/(2*Q2_infix);
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
            h_hat_kmer[i] = est_h_extra(h_kmer, i, n, x);
            h_hat_unimer[i] = est_h_unimer_extra(h_unimer, h_infix_eq, pi_broken, i, n, x);
            //h_hat_unimer[i] *= n/(n+x);
            //printf("%d %.2f %.2f\n", i, h_hat_kmer[i], h_hat_unimer[i]);
            //printf("%.2f ", h_hat_unimer[i]);
        }
        h_hat_kmer[m] = h_kmer[n];
        h_hat_unimer[m] = h_unimer[n];//*n/(n+x);
        //printf("%.2f\n", h_hat_unimer[m]);
        //printf("%d %.2f %.2f\n", m, h_hat_kmer[m], h_hat_unimer[m]);

        printf("ext\t%d\t", m);
        print_hill(m, h_hat_kmer, h_hat_unimer);
    }
}

void output_hill_cdbg(int argc, char *argv[]) {
    vector<double> h_kmer {0}, h_infix_eq{0}; //1-based index
    vector<int> points;
    //points = {1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 100};
    int n_points = 40;

    if (argc < 3) {
        fprintf(stderr, "Usage: pangrowth hill_cdbg <hist_kmer> <hist_infix> \n");
        return;
    }

    read_file(argv[1], h_kmer);
    read_file(argv[2], h_infix_eq);

    int n = h_kmer.size()-1;
    //vector<double> h_unitig(n+1);
    //for (int i = 1; i <= n; i++) {
    //    h_unitig[i] = h_kmer[i] - h_infix_eq[((i+1)*i/2)];
    //}

    get_points(n, n_points, points);

    hill_cdbg(h_kmer, h_infix_eq, points);
}
