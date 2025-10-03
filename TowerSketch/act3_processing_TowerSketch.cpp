#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <random>
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

string bits_to_sequence(uint64_t kmer, uint k) {
    string seq(k, 'N');
    for (int i = k-1; i >= 0; --i) {
        uint8_t b = kmer & 0x3ULL;
        switch(b) {
            case 0: seq[i] = 'A'; break;
            case 1: seq[i] = 'C'; break;
            case 2: seq[i] = 'G'; break;
            case 3: seq[i] = 'T'; break;
        }
        kmer >>= 2;
    }
    return seq;
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

void save_approx_kmers(
    const unordered_map<uint64_t,uint64_t> &counts, 
    TowerSketchCMCU &sketch, uint k, double phi, const string &out_folder)
{
    fs::create_directories(out_folder);
    string filename = out_folder + "/TS_k" + to_string(k) + ".txt";
    ofstream out(filename);

    uint64_t N = 0;
    for (auto &[kmer, c] : counts) N += c;
    uint64_t threshold = static_cast<uint64_t>(phi * N);
    size_t hh_count = 0;

    out << "#Heavy hitters TowerSketch k=" << k << " (threshold=" << threshold << ")\n";
    for (auto &[kmer, c] : counts) {
        uint32_t est = sketch.estimate(kmer, k);
        if (est >= threshold) {
            out << bits_to_sequence(kmer, k) << "\t" << est << "\n";
            hh_count++;
        }
    }
    cout << "  [TowerSketch] k=" << k
         << " | Total k-mers=" << N
         << " | HH encontrados=" << hh_count
         << " | Threshold=" << threshold
         << " | Guardado en: " << filename << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0] << " <carpeta_genomas> <k>\n";
        return 1;
    }
    string folder = argv[1];
    size_t k = std::stoi(argv[2]);
    vector<string> fasta_files = list_fasta_files(folder);
  
    random_device rd;
    mt19937 g(rd());
    shuffle(fasta_files.begin(), fasta_files.end(), g);
    vector<string> gen_files(fasta_files.begin(), fasta_files.begin() + 50);

    cout << "Genomas seleccionados para procesamiento:\n";
    for (const auto& file : gen_files) cout << "  " << file << "\n";

    unordered_map<uint64_t, uint64_t> counts;
    for (const auto& file : gen_files) {
        cout << "Procesando: " << file << endl;
        extract_kmers_bits(file, k, counts);
    }
    cout << "Total de k-mers distintos: " << counts.size() << endl;
    uint32_t d = 9, w = 10000000, levels = 6;
    TowerSketchCMCU sketch(d, w, levels);

    for (const auto& [kmer, freq] : counts) {
        for (uint64_t i = 0; i < freq; ++i) {
            sketch.insert(kmer, k);
        }
    }
    double phi = 1e-6;
    save_approx_kmers(counts, sketch, k, phi, "output/towersketch/k" + to_string(k) + "/phi" + to_string(int(-log10(phi))));
    return 0;
}