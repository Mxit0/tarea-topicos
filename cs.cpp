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

namespace fs = std::filesystem;
using namespace std;

class CountSketch {
private:
    int d; // Número de filas (hash functions)
    int w; // Número de columnas (contadores)
    vector<vector<int64_t>> C; // Matriz del sketch, renombrada a C
    vector<uint64_t> seeds_h;   // Semillas para la familia de hash h
    vector<uint64_t> seeds_g;   // Semillas para la familia de hash g

    // Función de hash base (privada)
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

    // Familia de funciones hash h_j(item) -> [0, w-1]
    uint64_t h(int j, uint64_t item) {
        return hash64(item ^ seeds_h[j]) % w;
    }

    // Familia de funciones hash g_j(item) -> {-1, 1}
    int g(int j, uint64_t item) {
        return (hash64(item ^ seeds_g[j]) % 2) * 2 - 1;
    }

public:
    // --- Inicialización ---
    CountSketch(int depth, int width) : d(depth), w(width) {
        // C[1..d][1..w] = 0
        C.assign(d, vector<int64_t>(w, 0));

        // Inicializar las "semillas" para simular d funciones hash h_j y g_j
        mt19937_64 rng(1337); // Semilla fija para reproducibilidad
        for (int i = 0; i < d; ++i) {
            seeds_h.push_back(rng());
            seeds_g.push_back(rng());
        }
    }

    // --- Actualización (corresponde a insertCS) ---
    // El parámetro 'count' corresponde a la constante 'c' en tu pseudocódigo
    void update(uint64_t item, int64_t count = 1) {
        for (int j = 0; j < d; ++j) {
            // C[j][h_j[item]] += g_j(item) * count
            C[j][h(j, item)] += g(j, item) * count;
        }
    }

    // --- Estimación (corresponde a estimarFreq) ---
    int64_t estimate(uint64_t item) {
        vector<int64_t> estimates;
        for (int j = 0; j < d; ++j) {
            // Coleccionar g_j(item) * C[j][h_j(item)]
            estimates.push_back(g(j, item) * C[j][h(j, item)]);
        }

        // Devolver la mediana de las estimaciones
        sort(estimates.begin(), estimates.end());
        return estimates[d / 2];
    }
    
    // Devuelve el tamaño del sketch en bytes
    size_t size_in_bytes() const {
        return d * w * sizeof(int64_t);
    }
};

vector<string> list_fasta_files(const string &folder);
uint64_t canonical_kmer(uint64_t kmer, uint k);
inline uint8_t base_to_bits(char c);
string bits_to_sequence(uint64_t kmer, uint k);
uint64_t extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts);


int main() {
    string folder = "Genomas/";
    vector<string> fasta_files = list_fasta_files(folder);
    cout << "Archivos encontrados: " << fasta_files.size() << endl;

    size_t subset_size = min((size_t)20, fasta_files.size());
    vector<string> subset_files(fasta_files.begin(), fasta_files.begin() + subset_size);
    cout << "Procesando subconjunto de " << subset_files.size() << " archivos para ground truth y calibración." << endl;

    // --- PASO 1: Generar Ground Truth (conteo exacto) ---
    unordered_map<uint64_t, uint64_t> ground_truth_counts;
    uint64_t N_total = 0;
    const size_t k = 31;

    cout << "\n--- Generando Ground Truth (k=" << k << ") ---" << endl;
    for (const auto &fp : subset_files) {
        cout << "Procesando: " << fp << endl;
        N_total += extract_kmers_bits(fp, k, ground_truth_counts);
    }
    cout << "Ground Truth generado. Total de k-mers: " << N_total << ". K-mers unicos: " << ground_truth_counts.size() << endl;

    // Identificar los Heavy Hitters reales para la evaluación
    double phi = 1e-7; // ** Umbral ajustado para encontrar HH **
    uint64_t threshold = static_cast<uint64_t>(phi * N_total);
    if (threshold < 2) threshold = 2;
    
    vector<pair<uint64_t, uint64_t>> heavy_hitters;
    for (const auto& [kmer, count] : ground_truth_counts) {
        if (count >= threshold) {
            heavy_hitters.push_back({kmer, count});
        }
    }
    cout << "Se usarán " << heavy_hitters.size() << " heavy hitters reales (frec >= " << threshold << ") para la calibración." << endl;

    // --- PASO 2: Calibración de CountSketch ---
    cout << "\n--- Calibrando CountSketch (error vs tamaño) ---" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << setw(5) << "d" << setw(10) << "w" << setw(15) << "Tamaño (KB)"
         << setw(20) << "Error Relativo Prom" << setw(20) << "Error Absoluto Prom" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;

    vector<int> d_values = {3, 5, 7};
    vector<int> w_values = {50000, 100000, 250000, 500000, 1000000};

    for (int d : d_values) {
        for (int w : w_values) {
            CountSketch cs(d, w);
            
            for (const auto& [kmer, count] : ground_truth_counts) {
                cs.update(kmer, count);
            }

            if (heavy_hitters.empty()) {
                if (w == w_values[0] && d == d_values[0]) // Imprimir advertencia solo una vez
                    cerr << "Advertencia: No se encontraron heavy hitters para calibrar." << endl;
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
            
            cout << fixed << setprecision(4);
            cout << setw(5) << d << setw(10) << w 
                 << setw(15) << cs.size_in_bytes() / 1024.0
                 << setw(20) << avg_relative_error
                 << setw(20) << avg_absolute_error << endl;
        }
    }
    cout << "--------------------------------------------------------------------------------" << endl;

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
    kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8) & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
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
uint64_t extract_kmers_bits(const string &filepath, size_t k, unordered_map<uint64_t,uint64_t> &counts) {
    ifstream infile(filepath);
    if (!infile) {
        cerr << "No se pudo abrir: " << filepath << endl;
        return 0;
    }
    string line, seq;
    uint64_t N_total = 0;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
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
                        N_total++;
                    }
                }
                seq.clear();
            }
        } else {
            seq.append(line);
        }
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
                N_total++;
            }
        }
    }
    return N_total;
}