#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h> 
#include <stdint.h> 
#include <limits> 
#include <string.h> 

#define EPS numeric_limits<double>::epsilon()

using namespace std;

#define ll unsigned long long

double choose_log(uint32_t n, uint32_t k) {
    double res = 0.0;
    if (k > n) {
        return 0.0;
    }

    k = ((k > n - k) ? n - k : k);
    for (uint32_t i = 0; i < k; i++) {
        res += log(n-i);
        res -= log(i+1);
    }
    return res;
}

void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

void read_panmatrix(const char* file, double* h, int n){
    fprintf(stderr,"Reading panmatrix...");
    ifstream fin(file);
    if (!fin.is_open()) {
        cout << "Error opening file: " << file << '\n';
        exit(EXIT_FAILURE);
    }

    stringstream ss; 
    string line;
    while (getline(fin, line, '\n')) {
        char separator = ((line.find('\t') !=-1ULL) ? '\t': ' ');
        vector<string> line_split;
        split(line, separator, line_split);

        int sum = 0;
        for (int i = 0; i < n; i++) {
            sum += stoi(line_split[i]);
        }
        h[sum]++;
    }
    fprintf(stderr,"ok\n");
}

void read_hist(const char* file, double* h, int n){
    fprintf(stderr,"Reading histogram...");
    ifstream fin(file);
    if (!fin.is_open()) {
        cout << "Error opening file: " << file << '\n';
        exit(EXIT_FAILURE);
    }

    string line;
    stringstream ss; 

    int i = 1;
    while (getline(fin, line, '\n')) {
        int val;
        ss.str(line);
        ss >> val;
        h[i++] = val;
        ss.clear();
    }
    fprintf(stderr,"ok\n");
}

void get_pangenome_size(double* h, uint32_t n){
    double tot = 0;
    double n_fall_m = 0;
    for (uint32_t i = 0; i < n+1; i++) tot += h[i];

    double* F = new double[n+1];
    for (uint32_t i = 0; i < n+1; i++) F[i] = 0;

    for (uint32_t m = 1; m <= n; m++) {
        double y = 0;
        n_fall_m += log(n-m+1);
        for (uint32_t i = 1; i <= n-m; i++) {
            F[i] += log((double)n-(double)m-(double)i+1);
            y += exp(log(h[i]) + F[i] - n_fall_m);
        }
        if (m > 1) printf(" ");
        printf("%.2f",tot - y);
    }
    cout <<  '\n' << flush;
}

void get_pangenome_core(double* h, uint32_t n){
    double n_fall_m = 0;

    double* F = new double[n+1];
    for (uint32_t i = 0; i < n+1; i++) F[i] = 0;

    for (uint32_t m = 1; m <= n; m++) {
        double y = 0;
        n_fall_m += log(n-m+1);
        for (uint32_t i = m; i <= n; i++) {
            F[i] += (log((double)i-(double)m+1));
            y += exp(log(h[i]) + F[i] - n_fall_m);
        }
        if (m > 1) printf(" ");
        printf("%.2f",y);
    }
    cout <<  '\n' << flush;
}

void get_pangenome_corequ(double* h, uint32_t n, double qu){
    double n_fall_m = 0.0;
    double m_fact = 0.0;

    double* F = new double[n+1]();
    double** q = new double*[n+1];
    for (uint32_t i = 0; i < n+1; i++) q[i] = new double[n+1]();

    for (uint32_t m = 1; m <= n; m++) {
        m_fact += log(m);
        //items present in all genomes
        double yl = 0.0;
        n_fall_m += log(n-m+1);
        for (uint32_t i = ceil(m*qu); i <= n; i++) {
            F[i] += log(i-m+1);
            yl += exp(log(h[i]) + F[i] - n_fall_m);
        }

        //items present in exactly qu*n genomes 
        double yr = 0.0;
        //std::cout << "m:" << m << '\n' << std::flush;
        for (uint32_t i = ceil(m*qu); i <= n; i++) {
            double sum_q = 0.0;
            bool add = false;
            for (uint32_t j = ceil(m*qu); j <= m-1; j++) {
                if (n + j + 1 > i + m && j <= i) {
                    if (q[i][j] == 0.0) {
                        q[i][j] = choose_log(i,j);
                    }
                    q[i][j] += log(n-i-m+1+j) - log(m-j);
                    sum_q += exp(q[i][j] + m_fact - n_fall_m);
                    add = true;
                } 
            }
            if (add) {
                yr += exp(log(h[i]) + log(sum_q));
            }
        }
        if (m > 1) printf(" ");
        printf("%.2f", yl+yr);
    }
    printf("\n");
}

int get_n_hist(const char* file){
    ifstream fin(file);
    string line;             // line == ""
    uint32_t n = 0;
    while (getline(fin, line, '\n')) {
        n++;  
    }
    return n;
}

int get_n_panmatrix(const char* file){
    ifstream fin(file);
    if (!fin.is_open()) {
        cout << "Error opening file: " << file << '\n';
        exit(EXIT_FAILURE);
    }

    stringstream ss; 
    string line;
    getline(fin, line, '\n'); 
    vector<string> line_split;
    char separator = ((line.find('\t') !=-1ULL) ? '\t': ' ');
    split(line, separator, line_split);
    return line_split.size();
}

void print_info() {
    fprintf(stderr, "Usage: pangrowth core [-q] <-h|-p> <input> \n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h         if <input> is an histogram\n");
    fprintf(stderr, "  -p         if <input> is a panmatrix\n");
    fprintf(stderr, "  -q         [0-1], default: 1.0. Quorum for an item to be in the core. \n\
             When 1, the item is counted as core if present at least \n\
             once in all genomes. If set to 0.9 the item needs to be \n\
             in 90%% of the genomes.\n");
}

void output_pangenome(int argc, char *argv[]){
    uint32_t n=0;
    double* h = 0;
    string input;
    bool panmatrix = false;
    bool hist = false;

    if (argc < 3) {
        fprintf(stderr, "Usage: pangrowth growth <-h|-p> <input> \n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -h         if <input> is an histogram\n");
        fprintf(stderr, "  -p         if <input> is a panmatrix\n");
        return;
    }

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-')
            continue;
        if (strncmp(argv[i],"-h",3) == 0) {  // Use histogram h
            hist = true;
            input = string(argv[i+1]);
        } else if(strncmp(argv[i], "-p",3) == 0) { // Use panmatrix
            panmatrix = true;
            input = string(argv[i+1]);
        } 
    }

    if (hist) {
        n = get_n_hist(input.c_str());
        h = new double[n+1]();
        read_hist(input.c_str(), h, n);
    }
    if (panmatrix) {
        n = get_n_panmatrix(input.c_str());
        h = new double[n+1]();
        read_panmatrix(input.c_str(), h, n);
    }

    get_pangenome_size(h, n); 

}

void output_core(int argc, char *argv[]){
    uint32_t n=0;
    double* h = 0;
    double quorum = 1.0;
    string input;
    bool panmatrix = false;
    bool hist = false;

    if (argc < 3) {
        print_info();
        return; 
    }

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-')
            continue;
        if (strncmp(argv[i],"-h",3) == 0) {  // Use histogram h
            hist = true;
            input = string(argv[i+1]);
        } else if(strncmp(argv[i], "-p",3) == 0) { // Use panmatrix
            panmatrix = true;
            input = string(argv[i+1]);
        } else if(strncmp(argv[i], "-q",3) == 0) {
            quorum = stod(argv[i+1]);
            if (fabs(quorum - 1.0) > EPS && quorum > 1.0) {
                quorum = 1.0;
            }
        } 

    }
    if (!panmatrix && ! hist) {
        print_info();
        return;
    }

    if (hist) {
        n = get_n_hist(input.c_str());
        h = new double[n+1]();
        read_hist(input.c_str(), h, n);
    }
    if (panmatrix) {
        n = get_n_panmatrix(input.c_str());
        h = new double[n+1]();
        read_panmatrix(input.c_str(), h, n);
    }

    if (fabs(quorum - 1.0) > EPS) {
        fprintf(stderr, "Quorum set to: %.2f\n", quorum);
        get_pangenome_corequ(h, n, quorum); 
    } else 
        get_pangenome_core(h, n); 
}
