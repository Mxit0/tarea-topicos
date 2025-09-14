#pragma once
#include <vector>
#include <random>
#include <algorithm>
#include <cstdint>

class TowerSketchCMCU {
public:
    TowerSketchCMCU(uint32_t d, uint32_t w, uint32_t levels)
        : d_(d), w_(w), levels_(levels) {
        sketches_.resize(levels_);
        for (uint32_t l = 0; l < levels_; ++l) {
            sketches_[l] = std::vector<std::vector<uint32_t>>(d_, std::vector<uint32_t>(w_, 0));
        }
        std::mt19937_64 gen(42);
        for (uint32_t i = 0; i < d_; ++i) hash_seeds_.push_back(gen());
    }

    void insert(uint64_t kmer) {
        for (uint32_t l = 0; l < levels_; ++l) {
            uint32_t min_val = UINT32_MAX;
            std::vector<uint32_t> idxs(d_);
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = hash(kmer, hash_seeds_[i]) % w_;
                idxs[i] = idx;
                min_val = std::min(min_val, sketches_[l][i][idx]);
            }
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = idxs[i];
                if (sketches_[l][i][idx] == min_val) sketches_[l][i][idx]++;
            }
        }
    }

    uint32_t estimate(uint64_t kmer) const {
        uint32_t est = 0;
        for (uint32_t l = 0; l < levels_; ++l) {
            uint32_t min_val = UINT32_MAX;
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = hash(kmer, hash_seeds_[i]) % w_;
                min_val = std::min(min_val, sketches_[l][i][idx]);
            }
            est = std::max(est, min_val); 
        }
        return est;
    }

private:
    uint32_t d_, w_, levels_;
    std::vector<std::vector<std::vector<uint32_t>>> sketches_;
    std::vector<uint64_t> hash_seeds_;

    static uint32_t hash(uint64_t key, uint64_t seed) {
        key ^= seed;
        key ^= (key >> 33);
        key *= 0xff51afd7ed558ccdULL;
        key ^= (key >> 33);
        key *= 0xc4ceb9fe1a85ec53ULL;
        key ^= (key >> 33);
        return static_cast<uint32_t>(key);
    }
};