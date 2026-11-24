// Verify behaviour of db_sobol_int for sampleID beyond the 2^16 construction
// limit: no memory out-of-bounds, but does the sequence stay valid (distinct)?
#include "mask.h"
#include <cstdint>
#include <cstdio>
#include <unordered_set>

int main() {
    // 1) Build the set of distinct values for the valid domain [0, 2^16).
    std::unordered_set<uint32_t> valid;
    for (uint32_t m = 0; m < (1u << 16); ++m) valid.insert(db_sobol_int(m, 0));
    std::printf("distinct values for m in [0, 2^16): %zu / 65536\n", valid.size());

    // 2) For m just above 2^16, do we collide with an already-used point?
    int coll = 0, novel = 0;
    for (uint32_t m = (1u << 16); m < (1u << 16) + 64u; ++m) {
        uint32_t v = db_sobol_int(m, 0);
        if (valid.count(v)) ++coll; else ++novel;
        if (coll + novel <= 6)
            std::printf("  m=%u -> %08x  %s\n", m, v,
                        valid.count(v) ? "DUPLICATE" : "novel");
    }
    std::printf("m in [2^16, 2^16+64): collisions=%d, novel=%d\n", coll, novel);

    // 3) Show that the whole image of db_sobol_int is at most 2^16 values
    //    (a bijection only on the low 16 bits).
    std::unordered_set<uint32_t> all;
    for (uint32_t m = 0; m < (1u << 16); ++m) all.insert(db_sobol_int(m, 0));
    for (uint32_t m = (1u << 16); m < (1u << 17); ++m) all.insert(db_sobol_int(m, 0));
    std::printf("distinct over [0, 2^17): %zu\n", all.size());

    // 4) Concretely: what does sample(i,j,sampleID,d) return for a large index?
    for (int s = 65536; s <= 65536 + 2; ++s)
        std::printf("sample(0,0,%d,0) = %.6f\n", s, (double)sample(0, 0, s, 0));
    return 0;
}
