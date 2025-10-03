#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <cmath>

using namespace std;

unordered_map<string, int> load_kmers(const string& filename) {
    unordered_map<string, int> kmers;
    ifstream in(filename);
    if (!in) {
        cerr << "No se pudo abrir archivo: " << filename << endl;
        return kmers;
    }
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t tab = line.find('\t');
        if (tab == string::npos) {
            size_t pos = line.find_last_of(" ");
            if (pos == string::npos) continue;
            string kmer = line.substr(0, pos);
            string num = line.substr(pos + 1);
            try { kmers[kmer] = stoi(num); }
            catch(...) { continue; }
        } else {
            string kmer = line.substr(0, tab);
            int freq = 0;
            try { freq = stoi(line.substr(tab + 1)); }
            catch(...) { freq = 0; }
            kmers[kmer] = freq;
        }
    }
    return kmers;
}

int main(int argc, char* argv[]) {
    string real_path = "HH/HH_k31.txt";
    string approx_path = "output/countsketch/k31/phi6/CS_k31.txt";
    string out_path = "evaluateCountSketchApprox_results31.txt";

    if (argc >= 2) real_path = argv[1];
    if (argc >= 3) approx_path = argv[2];
    if (argc >= 4) out_path = argv[3];

    auto real = load_kmers(real_path);
    auto approx = load_kmers(approx_path);

    ofstream fout(out_path);
    if (!fout) {
        cerr << "No se pudo crear archivo de salida: " << out_path << endl;
        return 1;
    }

    int tp = 0, fp = 0, fn = 0;
    for (const auto& [kmer, _] : approx) {
        if (real.count(kmer)) tp++;
        else fp++;
    }
    for (const auto& [kmer, _] : real) {
        if (!approx.count(kmer)) fn++;
    }
    double precision = (tp + fp) ? (double)tp / (tp + fp) : 0.0;
    double recall = (tp + fn) ? (double)tp / (tp + fn) : 0.0;
    double f1 = (precision + recall) ? 2 * precision * recall / (precision + recall) : 0.0;

    fout << "--- Métricas de candidatos ---\n";
    fout << "True Positives:  " << tp << "\n";
    fout << "False Positives: " << fp << "\n";
    fout << "False Negatives: " << fn << "\n";
    fout << "Precision:       " << precision << "\n";
    fout << "Recall:          " << recall << "\n";
    fout << "F1 Score:        " << f1 << "\n\n";

    cout << "\n--- Métricas de candidatos (en consola) ---\n";
    cout << "True Positives:  " << tp << endl;
    cout << "False Positives: " << fp << endl;
    cout << "False Negatives: " << fn << endl;
    cout << "Precision:       " << precision << endl;
    cout << "Recall:          " << recall << endl;
    cout << "F1 Score:        " << f1 << endl;

    double total_abs_err = 0, total_rel_err = 0;
    int count = 0;

    fout << "kmer\tReal\tAprox\tAbsErr\tRelErr\n";

    for (const auto& [kmer, real_freq] : real) {
        auto it = approx.find(kmer);
        if (it != approx.end()) {
            int est = it->second;
            double abs_err = abs(est - real_freq);
            double rel_err = real_freq ? (abs_err / real_freq) : 0.0;
            total_abs_err += abs_err;
            total_rel_err += rel_err;
            count++;
            fout << kmer << "\t" << real_freq << "\t" << est << "\t" << abs_err << "\t" << rel_err << "\n";
        } else {
            total_abs_err += real_freq;
            total_rel_err += 1.0;
            count++;
            fout << kmer << "\t" << real_freq << "\t0\t" << real_freq << "\t1\n";
        }
    }

    if (count > 0) {
        fout << "\nError absoluto promedio: " << (total_abs_err / count) << "\n";
        fout << "Error relativo promedio: " << (total_rel_err / count) << "\n";
    } else {
        fout << "\nNo hay k-mers para calcular errores.\n";
    }

    fout.close();
    cout << "\nResultados guardados en: " << out_path << endl;
    return 0;
}

