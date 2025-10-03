#pragma once
#include <vector>
#include <random>
#include <algorithm>
#include <cstdint>
#include <climits>
#include <iostream>

class TowerSketchCMCU {
public:
    TowerSketchCMCU(uint32_t d, uint32_t w, uint32_t levels)
        : d_(d), w_(w), levels_(levels) 
    {
        sketches_.resize(levels_);
        for (uint32_t l = 0; l < levels_; ++l) {
            sketches_[l] = std::vector<std::vector<uint32_t>>(d_, std::vector<uint32_t>(w_, 0));
        }
        std::mt19937_64 gen(42);
        for (uint32_t i = 0; i < d_; ++i) hash_seeds_.push_back(gen());
    }

    void insert(uint64_t kmer, uint k = 31) {
        uint64_t canon = canonical_kmer(kmer, k);
        for (uint32_t l = 0; l < levels_; ++l) {
            uint32_t min_val = UINT32_MAX;
            std::vector<uint32_t> idxs(d_);
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = hash(canon, hash_seeds_[i]) % w_;
                idxs[i] = idx;
                min_val = std::min(min_val, sketches_[l][i][idx]);
            }
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = idxs[i];
                if (sketches_[l][i][idx] == min_val) sketches_[l][i][idx]++;
            }
        }
    }

    uint32_t estimate(uint64_t kmer, uint k = 31) const {
        uint64_t canon = canonical_kmer(kmer, k);
        uint32_t est = 0;
        for (uint32_t l = 0; l < levels_; ++l) {
            uint32_t min_val = UINT32_MAX;
            for (uint32_t i = 0; i < d_; ++i) {
                uint32_t idx = hash(canon, hash_seeds_[i]) % w_;
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

    static uint64_t canonical_kmer(uint64_t kmer, uint k) {
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
};
