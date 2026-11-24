// Verify the dynamic mask: its embedded generator db_sobol_int(m, d) must equal
// the static table sobol_sequence[m][d] (the integer underlying
// sobol_burley::sample(m, d, seed=0)) for every m < 4096, d < 256, and must
// also produce valid, distinct points for m in (4096, 2^16).

#include "mask_dynamic.h"            // dynamic mask (self-contained)
#include "sobol_4096spp_256d.h"      // static table for cross-check (only for this test)

#include <cstdint>
#include <cstdio>

int main() {
    // 1) Dynamic == static over the whole table range.
    uint64_t mismatches = 0;
    for (uint32_t m = 0; m < 4096; ++m) {
        for (uint32_t d = 0; d < 256; ++d) {
            if (db_sobol_int(m, d) != sobol_sequence[m][d]) {
                if (mismatches < 8)
                    std::printf("MISMATCH m=%u d=%u dyn=%08x table=%08x\n",
                                m, d, db_sobol_int(m, d), sobol_sequence[m][d]);
                ++mismatches;
            }
        }
    }
    std::printf("dynamic vs static: mismatches=%llu / %llu\n",
                (unsigned long long)mismatches,
                (unsigned long long)4096ULL * 256ULL);

    // 2) Indices beyond 4096: valid [0,2^32) values, distinct within a prefix,
    //    and the first row of a 8192-sample prefix stays in range.
    const uint32_t probe[] = {4096, 5000, 8191, 10000, 32768, 65535};
    for (uint32_t m : probe) {
        std::printf("m=%u d=0 -> %08x\n", m, db_sobol_int(m, 0));
    }
    const uint32_t dims[] = {0u, 7u, 15u, 255u};
    for (uint32_t k = 0; k < 4; ++k) {
        std::printf("m=65535 d=%u -> %08x\n", dims[k], db_sobol_int(65535, dims[k]));
    }

    // A prefix of 8192 samples, dim 0: all distinct and in range (a Sobol net
    // prefix must be a permutation of 2^13 points in the first 13 bits).
    bool dup[8192] = {};
    bool ok = true;
    for (uint32_t m = 0; m < 8192; ++m) {
        uint32_t v = db_sobol_int(m, 0);
        uint32_t hi = v >> 19;  // top 13 bits = the net stratum index
        if (hi >= 8192 || dup[hi]) { ok = false; std::printf("dup/bad at m=%u\n", m); break; }
        dup[hi] = true;
    }
    std::printf("8192-sample prefix is a (1,13)-net permutation of top-13 bits: %s\n",
                ok ? "YES" : "NO");

    return (mismatches == 0 && ok) ? 0 : 1;
}
