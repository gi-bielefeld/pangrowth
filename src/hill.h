#include <iostream>
#include <vector>
#include <string>
#include "Rmath.h"
//cout << lchoose(10,0) << '\n' << flush; // 0
//cout << lchoose(0,0) << '\n' << flush;  // 0
//cout << lchoose(0,10) << '\n' << flush; // -inf

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
        tot += exp(log((double)h[j-1]) + lchoose(j,i) + lchoose(n-j,m-i) - lchoose(n,m));
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
            tot += exp(log(hbar[(sigma*(sigma-1)/2)+j-1]) + lchoose(j,i) + lchoose(n-sigma, m-i) - lchoose(n,m));
        }
    }
    return tot;
}

double static inline 
est_h_extra(vector<double> &h,int i,int n, int m) {
    double tot = 0;
    for (int j = i; j <= n-m+i; j++) {
        tot += exp(log((double)h[j-1]) + lchoose(j,i) + lchoose(n-j,m-i) - lchoose(n,m));
    }
    return tot;
}

double static inline 
est_h_unimer_extrj(vector<double> &hbar, int i, int n, int m) {
    double tot = 0;
    //for (int sigma = 1; sigma <= n; sigma++) {
    //untested, in case of error, go back to sigma=1
    for (int sigma = i; sigma <= n; sigma++) {
        for (int j = i; j <= sigma; j++) {
            tot += exp(log(hbar[(sigma*(sigma-1)/2)+j-1]) + lchoose(j,i) + lchoose(n-sigma, m-i) - lchoose(n,m));
        }
    }
    return tot;
}

void 
hill(vector<double> &h_unitig, vector<double> &h_kmer, vector<double> &h_infix_eq, vector<int> &points) {
    double n = h_unitig.size();

    //double U = 0.0;
    //for (int i = 0; i < n; i++) {
    //    U += (i+1)*h_unitig[i];
    //}
    
    printf("m\trichness\texp_entropy\tinv_gini_simp\n");
    int p = 0;
    //** Interpolation **//
    vector<double> h_hat_kmer(n);
    vector<double> h_hat_unimer(n);
    while (points[p] <= n) {
        int m = points[p++];
        
        for (int i = 1; i <= m; i++) {
            h_hat_kmer[i] = est_h(h_kmer, i, n, m);
            h_hat_unimer[i] = est_h_unimer(h_infix_eq, i, n, m);
            //printf("%d %.2f %.2f\n", i, h_hat_kmer[i], h_hat_unimer[i]);
        }
        printf("%d\t", m);
        printf("%.2f\t", richness(m, h_hat_kmer, h_hat_unimer));
        fflush(stdout);
        printf("%.2f\t", exp_entropy(m, h_hat_kmer, h_hat_unimer));
        fflush(stdout);
        printf("%.2f\n", inv_gini_simpson(m, h_hat_kmer, h_hat_unimer));
    }
    //** Extrapolation **//
    //while (p < points.size()) {
    //    int m = points[p++];
    //    for (int i = 1; i <= m; i++) {
    //        h_hat_kmer[i] = est_h(h_kmer, i, n, m);
    //        h_hat_unimer[i] = est_h_unimer(h_infix_eq, i, n, m);
    //        //printf("%d %.2f %.2f\n", i, h_hat_kmer[i], h_hat_unimer[i]);
    //    }
    //    printf("%d\t", m);
    //    printf("%.2f\t", richness(m, h_hat_kmer, h_hat_unimer));
    //    fflush(stdout);
    //    printf("%.2f\t", exp_entropy(m, U, h_hat_kmer, h_hat_unimer));
    //    fflush(stdout);
    //    printf("%.2f\n", inv_gini_simpson(m, U, h_hat_kmer, h_hat_unimer));
    //}
}

void output_hill_cdbg(int argc, char *argv[]) {
    //double* h = 0;
    vector<double> h_kmer, h_infix_eq;
    vector<int> points;
    //points = {1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 100};
    int n_points = 40;

    if (argc < 3) {
        fprintf(stderr, "Usage: pangrowth hill_cdbg <hist_kmer> <hist_infix> \n");
        return;
    }

    read_file(argv[1], h_kmer);
    read_file(argv[2], h_infix_eq);

    int n = h_kmer.size();
    vector<double> h_unitig(n);
    for (int i = 1; i <= n; i++) {
        h_unitig[i-1] = h_kmer[i-1] - h_infix_eq[((i+1)*i/2)-1];
    }

    get_points(n, n_points, points);

    hill(h_unitig, h_kmer, h_infix_eq, points);
}
