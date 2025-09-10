#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>

namespace fs = std::filesystem;
using namespace std;

// --- Listar archivos FASTA ---
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

// --- Función para obtener k-mer canónico en bits ---
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

// --- Convertir base a 2 bits ---
inline uint8_t base_to_bits(char c) {
    switch(toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4; // base inválida (N u otro)
    }
}

// --- Convertir k-mer de bits a secuencia ADN ---
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

// --- Extraer k-mers canónicos y contarlos ---
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

// --- Guardar heavy hitters como secuencias legibles ---
void save_heavy_hitters(const unordered_map<uint64_t,uint64_t> &counts, double phi, const string &filename, uint k) {
    ofstream out(filename);
    uint64_t N = 0;
    for (auto &[kmer, c] : counts) N += c;
    uint64_t threshold = static_cast<uint64_t>(phi * N);

    for (auto &[kmer, c] : counts) {
        if (c >= threshold) {
            out << bits_to_sequence(kmer, k) << "\n";
        }
    }
}

int main() {
    string folder = "Genomas/";
    vector<string> fasta_files = list_fasta_files(folder);
    cout << "Archivos encontrados: " << fasta_files.size() << endl;

    size_t subset_size = min((size_t)10, fasta_files.size());
    vector<string> subset_files(fasta_files.begin(), fasta_files.begin() + subset_size);
    cout << "Procesando subconjunto representativo de " << subset_files.size() << " archivos." << endl;

    unordered_map<uint64_t,uint64_t> counts21, counts31;
    for (const auto &fp : subset_files) {
        cout << "Procesando: " << fp << endl;
        extract_kmers_bits(fp, 21, counts21);
        extract_kmers_bits(fp, 31, counts31);
    }

    double phi = 0; // umbral de heavy hitters
    save_heavy_hitters(counts21, phi, "heavy_hitters_k21.txt", 21);
    save_heavy_hitters(counts31, phi, "heavy_hitters_k31.txt", 31);

    cout << "Heavy hitters guardados" << endl;
    return 0;
}
