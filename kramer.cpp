#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>

namespace fs = std::filesystem;
using namespace std;

// --- Listar archivos ---
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
        default: return 4; // base inválida
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

// --- Extraer k-mers canónicos y contarlos, devolviendo N ---
uint64_t extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts) {
    ifstream infile(filepath);
    if (!infile) {
        cerr << "No se pudo abrir: " << filepath << endl;
        return 0;
    }

    string line, seq;
    uint64_t N_total = 0;

    while (getline(infile, line)) {
        if (line[0] == '>') { // línea de encabezado FASTA
            if (!seq.empty()) {
                uint64_t kmer = 0;
                size_t bases = 0;
                for (char c : seq) {
                    uint8_t bits = base_to_bits(c);
                    if (bits > 3) { kmer = 0; bases = 0; continue; } // reiniciar si base inválida
                    kmer = ((kmer << 2) | bits) & ((1ULL << (k*2)) - 1); // desplazar bits
                    if (++bases >= k) {
                        counts[canonical_kmer(kmer, k)]++; // contar k-mer canónico
                        N_total++;
                        bases--;
                    }
                }
                seq.clear();
            }
        } else {
            for (char &c : line) seq.push_back(toupper(c));
        }
    }

    // procesar la última secuencia
    if (!seq.empty()) {
        uint64_t kmer = 0;
        size_t bases = 0;
        for (char c : seq) {
            uint8_t bits = base_to_bits(c);
            if (bits > 3) { kmer = 0; bases = 0; continue; }
            kmer = ((kmer << 2) | bits) & ((1ULL << (k*2)) - 1);
            if (++bases >= k) {
                counts[canonical_kmer(kmer, k)]++;
                N_total++;
                bases--;
            }
        }
    }

    return N_total;
}

// --- Guardar heavy hitters como secuencia\tfrecuencia---
void save_heavy_hitters(const unordered_map<uint64_t,uint64_t> &counts, double phi, const string &filename, uint k, uint64_t N_total) {
    ofstream out(filename);
    uint64_t nivel_minimo = static_cast<uint64_t>(phi * N_total);

    out << "# phi=" << phi << " nivel_minimo=" << nivel_minimo << " N_total=" << N_total << "\n";
    for (auto &[kmer, c] : counts) {
        if (c >= nivel_minimo) { // solo guardar heavy hitters
            string seq = bits_to_sequence(kmer, k);
            //cout << seq << "\t" << c << endl;  
            //out << seq << "\t" << c << "\n";  
        }
    }
}

int main() {
    string folder = "Genomas/";
    vector<string> fasta_files = list_fasta_files(folder);
    cout << "Archivos encontrados: " << fasta_files.size() << endl;

    size_t subset_size = min((size_t)20, fasta_files.size());
    vector<string> subset_files(fasta_files.begin(), fasta_files.begin() + subset_size);
    cout << "Procesando subconjunto de " << subset_files.size() << " archivos." << endl;

    unordered_map<uint64_t,uint64_t> counts21, counts31;
    uint64_t N21 = 0, N31 = 0;

    for (const auto &fp : subset_files) {
        cout << "Procesando: " << fp << endl;
        N21 += extract_kmers_bits(fp, 21, counts21);
        N31 += extract_kmers_bits(fp, 31, counts31);
    }

    double phi = 1e-7;

    save_heavy_hitters(counts21, phi, "heavy_hitters_k21.txt", 21, N21);
    save_heavy_hitters(counts31, phi, "heavy_hitters_k31.txt", 31, N31);

    cout << "Heavy hitters guardados con phi=" << phi << endl;
    return 0;
}
