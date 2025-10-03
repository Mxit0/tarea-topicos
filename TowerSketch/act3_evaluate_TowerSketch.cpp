#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <cmath>

using namespace std;

unordered_map<string, int> load_kmers(const string& filename) {
    unordered_map<string, int> kmers;
    ifstream in(filename);
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t tab = line.find('\t');
        if (tab == string::npos) continue;
        string kmer = line.substr(0, tab);
        int freq = stoi(line.substr(tab + 1));
        kmers[kmer] = freq;
    }
    return kmers;
}

int main() {
    auto real = load_kmers("HH/HH_k31.txt");
    auto approx = load_kmers("output/towersketch/k31/phi6/TS_k31.txt");

    int tp = 0, fp = 0, fn = 0;
    for (const auto& [kmer, _] : approx) {
        if (real.count(kmer)) tp++;
        else fp++;
    }
    for (const auto& [kmer, _] : real) {
        if (!approx.count(kmer)) fn++;
    }
    double precision = tp ? (double)tp / (tp + fp) : 0.0;
    double recall = tp ? (double)tp / (tp + fn) : 0.0;
    double f1 = (precision + recall) ? 2 * precision * recall / (precision + recall) : 0.0;

    cout << "\n--- Métricas de candidatos ---\n";
    cout << "True Positives:  " << tp << endl;
    cout << "False Positives: " << fp << endl;
    cout << "False Negatives: " << fn << endl;
    cout << "Precision:       " << precision << endl;
    cout << "Recall:          " << recall << endl;
    cout << "F1 Score:        " << f1 << endl;

//------------------------------------------------------//

    double total_abs_err = 0, total_rel_err = 0;
    int count = 0;

    for (const auto& [kmer, real_freq] : real) {
        auto it = approx.find(kmer);
        if (it != approx.end()) {
            int est = it->second;
            double abs_err = abs(est - real_freq);
            double rel_err = abs_err / real_freq;
            total_abs_err += abs_err;
            total_rel_err += rel_err;
            count++;
            cout << kmer << "\tReal: " << real_freq << "\tAprox: " << est
                 << "\tAbsErr: " << abs_err << "\tRelErr: " << rel_err << endl;
        } else {
            total_abs_err += real_freq;
            total_rel_err += 1.0;
            count++;
            cout << kmer << "\tReal: " << real_freq << "\tAprox: 0"
                 << "\tAbsErr: " << real_freq << "\tRelErr: 1" << endl;
        }
    }

    cout << "\nError absoluto promedio: " << (total_abs_err / count) << endl;
    cout << "Error relativo promedio: " << (total_rel_err / count) << endl;
    return 0;
}