// Verifies the generated sobol_4096spp_256d.h against the live library.
//
// Includes the generated header and checks that every entry equals
// sample_int(i, d) computed directly from sobol_burley's building blocks.
//
// Usage: verify_table <path-to-sobol_4096spp_256d.h>

#include "sobol_burley.hpp"

#include <cstdio>
#include <cstdint>
#include <cstdlib>

using namespace sobol_burley;

namespace {

inline uint32_t sample_int(uint32_t sample_index, uint32_t dimension) {
    constexpr uint32_t kSeed = 0u;
    const uint32_t shuffled_rev_index =
        owen_scramble_rev(reverse_bits(sample_index), hash(kSeed ^ 0x79c68e4au));
    const uint32_t sobol = sobol_rev(shuffled_rev_index, dimension);
    static constexpr std::array<uint32_t, 4> kScrambleConsts = {
        0x912f69bau, 0x174f18abu, 0x691e72cau, 0xb40cc1b8u};
    const uint32_t seed_ = kSeed * 0x9c8f2d3bu;
    const uint32_t ds = dimension >> 2;
    const uint32_t scramble = ds ^ seed_ ^ kScrambleConsts[dimension & 0b11u];
    const uint32_t sobol_owen_rev = owen_scramble_rev(sobol, hash(scramble));
    return reverse_bits(sobol_owen_rev);
}

}  // namespace

#include "sobol_4096spp_256d.h"

int main(int argc, char** argv) {
    constexpr uint32_t kSamples = 4096;
    constexpr uint32_t kDims = 256;

    // Sanity: correct array dimensions.
    static_assert(sizeof(sobol_sequence) / sizeof(sobol_sequence[0]) == kSamples,
                  "unexpected sample count");
    static_assert(sizeof(sobol_sequence[0]) / sizeof(sobol_sequence[0][0]) == kDims,
                  "unexpected dimension count");

    // Check every entry against the live library.
    uint64_t mismatches = 0;
    for (uint32_t i = 0; i < kSamples; ++i) {
        for (uint32_t d = 0; d < kDims; ++d) {
            const uint32_t expected = sample_int(i, d);
            if (sobol_sequence[i][d] != expected) {
                if (mismatches < 10) {
                    std::printf("mismatch at [%u][%u]: table=0x%08x live=0x%08x\n", i, d,
                                sobol_sequence[i][d], expected);
                }
                ++mismatches;
            }
        }
    }

    if (mismatches == 0) {
        std::printf("verify_table: ALL %u x %u ENTRIES MATCH LIVE LIBRARY\n", kSamples, kDims);
        return 0;
    }
    std::printf("verify_table: %llu MISMATCHES\n", (unsigned long long)mismatches);
    return 1;
}
