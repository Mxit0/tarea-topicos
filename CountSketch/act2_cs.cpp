#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cmath>
#include <iomanip>
#include <sys/resource.h>

namespace fs = std::filesystem;
using namespace std;

class CountSketch {
private:
    int d; 
    int w; 
    vector<vector<int64_t>> C;
    vector<uint64_t> seeds_h;
    vector<uint64_t> seeds_g;

    uint64_t hash64(uint64_t key) {
        key = (~key) + (key << 21);
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8);
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }

    uint64_t h(int j, uint64_t item) { return hash64(item ^ seeds_h[j]) % w; }
    int g(int j, uint64_t item) { return (hash64(item ^ seeds_g[j]) % 2) * 2 - 1; }

public:
    CountSketch(int depth, int width) : d(depth), w(width) {
        C.assign(d, vector<int64_t>(w, 0));
        mt19937_64 rng(1337);
        for (int i = 0; i < d; ++i) {
            seeds_h.push_back(rng());
            seeds_g.push_back(rng());
        }
    }

    void update(uint64_t item, int64_t count = 1) {
        for (int j = 0; j < d; ++j) {
            C[j][h(j, item)] += g(j, item) * count;
        }
    }

    int64_t estimate(uint64_t item) {
        vector<int64_t> estimates;
        for (int j = 0; j < d; ++j) {
            estimates.push_back(g(j, item) * C[j][h(j, item)]);
        }
        sort(estimates.begin(), estimates.end());
        return estimates[d / 2];
    }
    
    size_t size_in_bytes() const { return d * w * sizeof(int64_t); }
};


vector<string> list_fasta_files(const string &folder);
uint64_t canonical_kmer(uint64_t kmer, uint k);
inline uint8_t base_to_bits(char c);
string bits_to_sequence(uint64_t kmer, uint k);
long report_memory();
uint64_t extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts);


int main() {
    string folder = "Genomas/";
    vector<string> fasta_files = list_fasta_files(folder);
    cout << "Archivos encontrados: " << fasta_files.size() << endl;

    size_t subset_size = min((size_t)25, fasta_files.size()); 
    vector<string> subset_files(fasta_files.begin(), fasta_files.begin() + subset_size);
    cout << "Procesando subconjunto de " << subset_files.size() << " archivos para ground truth y calibración." << endl;

    unordered_map<uint64_t, uint64_t> ground_truth_counts;
    uint64_t N_total = 0;
    const size_t k = 21;

    cout << "\n--- Generando Ground Truth (k=" << k << ") ---" << endl;
    for (const auto &fp : subset_files) {
        cout << "Procesando: " << fp << endl;
        N_total += extract_kmers_bits(fp, k, ground_truth_counts);
    }
    cout << "Ground Truth generado. Total de k-mers: " << N_total << ". K-mers unicos: " << ground_truth_counts.size() << endl;

    double phi = 1e-6;
    uint64_t threshold = static_cast<uint64_t>(phi * N_total);
    if (threshold < 2) threshold = 2;
    
    vector<pair<uint64_t, uint64_t>> heavy_hitters;
    for (const auto& [kmer, count] : ground_truth_counts) {
        if (count >= threshold) {
            heavy_hitters.push_back({kmer, count});
        }
    }
    cout << "Se usarán " << heavy_hitters.size() << " heavy hitters reales (frec >= " << threshold << ") para la calibración." << endl;

    ofstream outfile("calibration.txt");
    if (!outfile.is_open()) {
        cerr << "Error: No se pudo abrir el archivo calibration.txt para escribir." << endl;
        return 1; 
    }

    stringstream header_stream;
    header_stream << "--------------------------------------------------------------------------------" << endl;
    header_stream << setw(5) << "d" << setw(10) << "w" << setw(15) << "Tamaño (KB)"
                  << setw(20) << "Error Relativo Prom" << setw(20) << "Error Absoluto Prom" << endl;
    header_stream << "--------------------------------------------------------------------------------" << endl;

    cout << "\n--- Calibrando CountSketch (error vs tamaño) ---" << endl;
    cout << header_stream.str();
    outfile << "Tabla de Calibración de CountSketch (error vs tamaño)\n";
    outfile << header_stream.str();

    vector<int> d_values = {3, 5, 7, 9};
    vector<int> w_values = {500000, 1000000, 2500000, 5000000, 10000000};

    for (int d : d_values) {
        for (int w : w_values) {
            CountSketch cs(d, w);
            
            for (const auto& [kmer, count] : ground_truth_counts) {
                cs.update(kmer, count);
            }

            if (heavy_hitters.empty()) {
                if (w == w_values[0] && d == d_values[0])
                    cerr << "Advertencia: No se encontraron heavy hitters para calibrar. Prueba bajando el valor de phi." << endl;
                continue;
            }

            double total_relative_error = 0.0;
            double total_absolute_error = 0.0;
            for (const auto& [kmer, true_count] : heavy_hitters) {
                int64_t estimated_count = cs.estimate(kmer);
                double absolute_error = std::abs(estimated_count - (int64_t)true_count);
                total_absolute_error += absolute_error;
                total_relative_error += absolute_error / true_count;
            }

            double avg_relative_error = total_relative_error / heavy_hitters.size();
            double avg_absolute_error = total_absolute_error / heavy_hitters.size();
            
            stringstream data_line_stream;
            data_line_stream << fixed << setprecision(4);
            data_line_stream << setw(5) << d << setw(10) << w 
                             << setw(15) << cs.size_in_bytes() / 1024.0
                             << setw(20) << avg_relative_error
                             << setw(20) << avg_absolute_error << endl;

            cout << data_line_stream.str();
            outfile << data_line_stream.str();
        }
    }
    
    cout << "--------------------------------------------------------------------------------" << endl;
    outfile << "--------------------------------------------------------------------------------" << endl;

    outfile.close();
    cout << "\nResultados de la calibración guardados en 'calibration.txt'" << endl;
    
    long mem_kb_final = report_memory();
    cout << "Memoria máxima usada (proceso): " << mem_kb_final/1024.0 << " MB" << endl;
    cout << "Calibración completa." << endl;

    return 0;
}


vector<string> list_fasta_files(const string &folder) {
    vector<string> files;
    for (const auto &entry : fs::directory_iterator(folder)) {
        string filename = entry.path().filename().string();
        if (filename.rfind("GCA_", 0) == 0 && filename.size() > 4 && filename.substr(filename.size()-4) == ".fna") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

uint64_t canonical_kmer(uint64_t kmer, uint k) {
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;
    kmer = ((kmer >> 2) & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8) & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32 ) | ( kmer << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
    return (b_kmer < reverse) ? b_kmer : reverse;
}

inline uint8_t base_to_bits(char c) {
    switch(toupper(c)) {
        case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3;
        default: return 4;
    }
}

string bits_to_sequence(uint64_t kmer, uint k) {
    string seq(k, 'N');
    for (int i = k-1; i >= 0; --i) {
        uint8_t b = kmer & 0x3ULL;
        switch(b) {
            case 0: seq[i] = 'A'; break; case 1: seq[i] = 'C'; break;
            case 2: seq[i] = 'G'; break; case 3: seq[i] = 'T'; break;
        }
        kmer >>= 2;
    }
    return seq;
}

long report_memory() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_maxrss;
}

uint64_t extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts) {
    ifstream infile(filepath);
    if (!infile) {
        cerr << "No se pudo abrir: " << filepath << endl;
        return 0;
    }
    string line, seq;
    uint64_t N_total_file = 0;
    
    while (getline(infile, line)) {
        if (line.empty() || line[0] == '>') continue;
        seq.append(line);
    }

    if (!seq.empty()) {
        uint64_t kmer = 0;
        uint64_t mask = (1ULL << (k*2)) - 1;
        size_t valid_bases = 0;
        for (char c : seq) {
            uint8_t bits = base_to_bits(c);
            if (bits > 3) { valid_bases = 0; continue; }
            kmer = ((kmer << 2) | bits) & mask;
            if (++valid_bases >= k) {
                counts[canonical_kmer(kmer, k)]++;
                N_total_file++;
            }
        }
    }
    return N_total_file;
}