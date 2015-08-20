#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <omp.h>

#include <armadillo>

using namespace std;

void usage(int rc) {
    cerr << 
        "Output Pearson correlations between input matrix columns to stdout\n"
        "USAGE: pcor [-n num]\n\n"
        "Arguments:\n"
        "-n <N> :\n"
        "   If -n is specified, output the top N most correlated columns labels.\n"
        "   If not, output all pairwise correlations as a large matrix\n"
        "-m <N> :\n"
        "   Pairwise columns must have at least N samples that are mutually\n"
        "   non-NaN, otherwise the correlation will not be computed.\n"
        "-p <N> :\n"
        "   Use N OpenMP threads (defaults to the OMP_NUM_THREADS environment variable\n"
        "   if set, otherwise implementation-dependent (usually number of hyperthreads)\n";
    exit(rc);            
}

vector<string> 
split(const string &text, char sep) {
    vector<string> tokens;
    int start = 0, end = 0;
    while ((end = text.find(sep, start)) != string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

struct Matrix {
    vector<string> index, columns;
    arma::mat data;
};

Matrix
read_matrix(std::istream& in) {
    Matrix o;
    vector<vector<double> > data;
    string line;

    getline(in, line);
    vector<string> columns = split(line, '\t');
    for (int i=1; i<columns.size(); i++) {
        o.columns.push_back(columns[i]);
    }

    while (getline(in, line)) {
        vector<double> v;
        vector<string> tokens = split(line, '\t');
        o.index.push_back(tokens[0]);
        for (int i=1; i<tokens.size(); i++) {
            v.push_back(atof(tokens[i].c_str()));
        }
        data.push_back(v);
    }

    o.data = arma::mat(data.size(), data[0].size());
    for (int i=0; i<data.size(); i++) {
        for (int j=0; j<data[0].size(); j++) {
            o.data(i,j) = data[i][j];
        }
    }
    return o;
}

arma::vec
pairwise_correlate(const arma::mat& X, const arma::vec& v1, size_t min_samples) {
    size_t nc = X.n_cols;
    size_t nr = X.n_rows;
    arma::vec o(nc);
    for (int j=0; j<nc; j++)
        o(j) = nan("");

    for (int j=0; j<nc; j++) {
        arma::vec v2 = X.col(j);
        vector<double> vv1, vv2;
        size_t n = 0;
        for (int i=0; i<nr; i++) {
            double x1 = v1(i);
            double x2 = v2(i);
            if ((x1 == x1) && (x2 == x2)) {
                n++;
                vv1.push_back(x1);
                vv2.push_back(x2);
            }
        }
        if (n >= min_samples) {
            arma::mat av1 = arma::mat(vv1);
            arma::mat av2 = arma::mat(vv2);
            arma::mat rs = arma::cor(av1, av2);
            double r = rs(0,0);
            o(j) = r;
        }
    }
    return o;
}

int main(int argc, char* argv[]) {
    bool show_usage = false;
    int top_n = 0;
    size_t n_cpu = 0;
    size_t min_samples = 3;
    int c;
    while ((c = getopt(argc, argv, "hn:m:")) != -1) {
        switch (c) {
            case 'h':
                show_usage = true;
                break;
            case 'n':
                top_n = atoi(optarg);
                break;
            case 'm':
                min_samples = atol(optarg);
                break;
            case 'p':
                n_cpu = atoi(optarg);
                break;
        }
    }

    if (n_cpu > 0)
        omp_set_num_threads(n_cpu);

    Matrix X = read_matrix(cin);
    if ((top_n > 0) && (top_n > X.data.n_cols)) {
        cerr << "ERROR: -n=" << top_n 
            << " is greater than # matrix columns (" << X.data.n_cols << ")\n";
        exit(1);
    }

    if (top_n == 0) {
        for (auto& c : X.columns) {
            cout << "\t" << c;
        }
        cout << endl;
    }

    size_t nc = X.data.n_cols;

    #pragma omp parallel for shared(X)
    for (int i=0; i<nc; i++) {
        auto x = X.data.col(i);
        auto r = pairwise_correlate(X.data, x, min_samples);
        auto ix = arma::sort_index(r);
        #pragma omp critical
        {
            cout << X.columns[i];
            if (top_n > 0) {
                size_t no = 0, j = nc;
                while ((no < top_n) && (j > 0)) {
                    j--;
                    double rr = r(j);
                    if (rr == rr) {
                        cout << "\t" << X.columns[ix(j)];
                        no++;
                    }
                }
            } else {
                for (int j=0; j<nc; j++) {
                    cout << "\t" << r(j);
                }
            }
            cout << endl;
        }
    }
}
