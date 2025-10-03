#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <random>
#include <iomanip>
#include "towersketch.cpp"

namespace fs = std::filesystem;
using namespace std;

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

uint64_t canonical_kmer(uint64_t kmer, uint k = 31) {
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;
    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32 ) | ( kmer << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));
    return (b_kmer < reverse) ? b_kmer : reverse;
}

inline uint8_t base_to_bits(char c) {
    switch(toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4; 
    }
}

void extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts) {
    ifstream infile(filepath);
    if (!infile) {
        cerr << "No se pudo abrir: " << filepath << endl;
        return;
    }
    string line, seq;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty()) {
                uint64_t kmer = 0;
                size_t bases = 0;
                for (char c : seq) {
                    uint8_t bits = base_to_bits(c);
                    if (bits > 3) { kmer = 0; bases = 0; continue; }
                    kmer = ((kmer << 2) | bits) & ((1ULL << (k*2)) - 1);
                    if (++bases >= k) {
                        uint64_t cano = canonical_kmer(kmer, k);
                        counts[cano]++;
                        bases--;
                    }
                }
                seq.clear();
            }
        } else {
            for (char &c : line) seq.push_back(toupper(c));
        }
    }
    if (!seq.empty()) {
        uint64_t kmer = 0;
        size_t bases = 0;
        for (char c : seq) {
            uint8_t bits = base_to_bits(c);
            if (bits > 3) { kmer = 0; bases = 0; continue; }
            kmer = ((kmer << 2) | bits) & ((1ULL << (k*2)) - 1);
            if (++bases >= k) {
                uint64_t cano = canonical_kmer(kmer, k);
                counts[cano]++;
                bases--;
            }
        }
    }
}

size_t sketch_size_bytes(uint32_t d, uint32_t w, uint32_t levels) {

    return size_t(d) * size_t(w) * size_t(levels) * sizeof(uint32_t);
}

int main(int argc, char* argv[]) {
    string folder = "Genomas/";
    size_t k = 31;

    vector<string> fasta_files = list_fasta_files(folder);
    if (fasta_files.size() < 25) {
        cerr << "Se requieren al menos 25 archivos fasta en la carpeta.\n";
        return 1;
    }

    random_device rd;
    mt19937 g(rd());
    shuffle(fasta_files.begin(), fasta_files.end(), g);
    vector<string> sample_files(fasta_files.begin(), fasta_files.begin() + 25);

    unordered_map<uint64_t, uint64_t> counts;
    for (const auto& file : sample_files) {
        extract_kmers_bits(file, k, counts);
    }

    vector<uint32_t> d_values = {5, 7, 9};
    vector<uint32_t> w_values = {2500000, 5000000, 10000000};
    vector<uint32_t> levels_values = {5, 6, 7};

    cout << "---------------------------------------------------------------------------------------------\n";
    cout << "    d         w   levels   Tamaño (KB) Error Relativo Prom Error Absoluto Prom\n";
    cout << "---------------------------------------------------------------------------------------------\n";
    cout << fixed << setprecision(4);

    for (auto d : d_values) {
        for (auto w : w_values) {
            for (auto levels : levels_values) {
                TowerSketchCMCU sketch(d, w, levels);
                for (const auto& [kmer, freq] : counts) {
                    for (uint64_t i = 0; i < freq; ++i) {
                        sketch.insert(kmer, k);
                    }
                }
                double total_abs_err = 0, total_rel_err = 0;
                int n = 0;
                for (const auto& [kmer, real_freq] : counts) {
                    uint32_t est = sketch.estimate(kmer, k);
                    double abs_err = abs((double)est - (double)real_freq);
                    double rel_err = real_freq > 0 ? abs_err / real_freq : 0.0;
                    total_abs_err += abs_err;
                    total_rel_err += rel_err;
                    n++;
                }
                double mean_abs_err = n ? total_abs_err / n : 0.0;
                double mean_rel_err = n ? total_rel_err / n : 0.0;
                double size_kb = sketch_size_bytes(d, w, levels) / 1024.0;

                cout << setw(5) << d << setw(11) << w << setw(9) << levels
                     << setw(14) << size_kb
                     << setw(21) << mean_rel_err
                     << setw(21) << mean_abs_err << endl;
            }
        }
    }
    cout << "---------------------------------------------------------------------------------------------\n";
    return 0;
}