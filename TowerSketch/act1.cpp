#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <sys/resource.h> 
#include <cmath>

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

void save_heavy_hitters_phi(
    const unordered_map<uint64_t,uint64_t> &counts, uint k,
    double phi, const string &out_folder)
{
    fs::create_directories(out_folder);
    string filename = out_folder + "/HH_k" + to_string(k) + ".txt";
    ofstream out(filename);

    uint64_t N = 0;
    for (auto &[kmer, c] : counts) N += c;
    uint64_t threshold = static_cast<uint64_t>(phi * N);
    size_t hh_count = 0;

    out << "#Heavy hitters k=" << k << " (threshold=" << threshold << ")\n";
    for (auto &[kmer, c] : counts) {
        if (c >= threshold) {
            out << bits_to_sequence(kmer, k) << "\t" << c << "\n";
            hh_count++;
        }
    }

    cout << "  [Phi=" << phi << "] k=" << k
         << " | Total k-mers=" << N
         << " | HH encontrados=" << hh_count
         << " | Threshold=" << threshold
         << " | Guardado en: " << filename << endl;
}

long report_memory() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_maxrss;
}


int main() {
    string folder = "Genomas/";
    vector<string> fasta_files = list_fasta_files(folder);
    cout << "Archivos encontrados: " << fasta_files.size() << endl;

    unordered_map<uint64_t,uint64_t> counts21, counts31;
    size_t file_counter = 0;
    for (const auto &fp : fasta_files) {
        file_counter++;
        cout << "------------------------------------------------------------" << endl;
        cout << "Procesando archivo " << file_counter << "/" << fasta_files.size() << ": " << fp << endl;

        size_t before21 = counts21.size();
        size_t before31 = counts31.size();

        extract_kmers_bits(fp, 21, counts21);
        extract_kmers_bits(fp, 31, counts31);

        cout << "  -> K=21: Total k-mers únicos hasta ahora: " << counts21.size() << " (+" << counts21.size()-before21 << ")" << endl;
        cout << "  -> K=31: Total k-mers únicos hasta ahora: " << counts31.size() << " (+" << counts31.size()-before31 << ")" << endl;

        long mem_kb = report_memory();
        cout << "  -> Memoria actual (proceso): " << mem_kb/1024.0 << " MB" << endl;
    }

    // Calcular memoria estimada de las estructuras
    size_t memory_counts21 = counts21.size() * (sizeof(uint64_t) + sizeof(uint64_t));
    size_t memory_counts31 = counts31.size() * (sizeof(uint64_t) + sizeof(uint64_t));
    cout << "------------------------------------------------------------" << endl;
    cout << "Memoria estimada de counts21: " << memory_counts21 / 1024.0 / 1024.0 << " MB" << endl;
    cout << "Memoria estimada de counts31: " << memory_counts31 / 1024.0 / 1024.0 << " MB" << endl;

    // --- Múltiples phi ---
    vector<double> phi_values = {1e-6};
    for (double phi : phi_values) {
        cout << "============================================================" << endl;
        cout << "Guardando HH para phi=" << phi << endl;

        save_heavy_hitters_phi(counts21, 21, phi, "output/k21/phi" + to_string(int(-log10(phi))));
        save_heavy_hitters_phi(counts31, 31, phi, "output/k31/phi" + to_string(int(-log10(phi))));
    }

    long mem_kb_final = report_memory();
    cout << "------------------------------------------------------------" << endl;
    cout << "Memoria máxima usada (proceso): " << mem_kb_final/1024.0 << " MB" << endl;
    cout << "Procesamiento completo." << endl;

    return 0;
}
