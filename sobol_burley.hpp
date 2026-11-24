// A seedable Owen-scrambled Sobol sequence.
//
// A C++17 port of the Rust crate `sobol_burley` (v0.5.0) by Nathan Vegdahl,
// based on the paper "Practical Hash-based Owen Scrambling" by Brent Burley,
// with an improved hash from "Building a Better LK Hash" and more dimensions
// due to Kuo et al.
//
// Design constraints (kept identical to the original):
//   * The maximum sequence length is 2^16.
//   * The maximum number of dimensions is 256 (workaround: reseed to reuse
//     dimensions, i.e. "padding").
//   * Only `f32` output is supported.
//
// This is a scalar (non-SIMD) implementation. The original crate's `sample_4d`
// routine used SIMD on x86-64; here it is computed with plain arithmetic, but
// produces bit-identical results to the Rust `sample_4d` (and to calling
// `sample` four times).
//
// License: MIT OR Apache-2.0.  The Sobol direction numbers embedded in
// `direction_numbers_embedded.hpp` are adapted from work by Stephen Joe and
// Frances Y. Kuo and are under the 3-clause BSD license (see licenses/JOE_KUO.txt).

#pragma once

#include <array>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <system_error>
#include <vector>

#include "direction_numbers_embedded.hpp"

namespace sobol_burley {

//---------------------------------------------------------------------------
// Constants
//---------------------------------------------------------------------------

/// The number of available dimensions.
constexpr uint32_t NUM_DIMENSIONS = 256;

/// The number of available 4d dimension sets.  Just `NUM_DIMENSIONS / 4`.
constexpr uint32_t NUM_DIMENSION_SETS_4D = NUM_DIMENSIONS / 4;

/// Bits per direction vector element.
constexpr int SOBOL_BITS = 32;

/// Number of direction vector elements per dimension.
constexpr int SOBOL_DEPTH = 16;

/// Storage layout mirrors the Rust crate: REV_VECTORS[dimension_set][depth][sub_dimension].
/// Each stored element already has its bits reversed.
using RevVectors =
    std::array<std::array<std::array<uint32_t, 4>, SOBOL_DEPTH>, NUM_DIMENSION_SETS_4D>;

//---------------------------------------------------------------------------
// Internal helpers
//---------------------------------------------------------------------------

/// Reverse the order of the 32 bits of `n`.
inline uint32_t reverse_bits(uint32_t n) {
    n = ((n & 0x55555555u) << 1) | ((n >> 1) & 0x55555555u);
    n = ((n & 0x33333333u) << 2) | ((n >> 2) & 0x33333333u);
    n = ((n & 0x0f0f0f0fu) << 4) | ((n >> 4) & 0x0f0f0f0fu);
    n = ((n & 0x00ff00ffu) << 8) | ((n >> 8) & 0x00ff00ffu);
    n = (n << 16) | (n >> 16);
    return n;
}

namespace detail {

/// Parse the direction numbers from a single line of the Joe & Kuo text.
/// Returns false (and leaves `a`/`m` untouched) when the line is not a valid
/// data line (header, comments, blank, or unparseable).
inline bool parse_direction_numbers(std::string_view line, uint32_t& a,
                                    std::vector<uint32_t>& m) {
    // Mirror the Rust parser's tokenization.
    std::vector<std::string_view> tokens;
    {
        size_t pos = 0;
        while (pos <= line.size()) {
            size_t start = pos;
            while (pos < line.size() &&
                   !(line[pos] == ' ' || line[pos] == '\t' || line[pos] == '\r')) {
                ++pos;
            }
            if (pos > start) tokens.push_back(line.substr(start, pos - start));
            ++pos;  // skip the delimiter (may run past the end harmlessly)
        }
    }

    if (tokens.size() < 4 || (!line.empty() && line.front() == '#')) return false;

    // tokens[0] = dimension, tokens[1] = count  (both discarded here)
    // a = tokens[2], m = tokens[3..]
    a = 0;
    {
        size_t pos = 0;
        auto r = std::from_chars(tokens[2].data(), tokens[2].data() + tokens[2].size(), a);
        if (r.ec != std::errc{} || r.ptr != tokens[2].data() + tokens[2].size()) return false;
    }
    m.clear();
    m.reserve(tokens.size() - 3);
    for (size_t i = 3; i < tokens.size(); ++i) {
        uint32_t v = 0;
        auto r = std::from_chars(tokens[i].data(), tokens[i].data() + tokens[i].size(), v);
        if (r.ec != std::errc{} || r.ptr != tokens[i].data() + tokens[i].size()) return false;
        m.push_back(v);
    }
    return true;
}

/// Generate the full set of Sobol direction vectors (mirrors `build.rs`),
/// then interleave them into the reversed, SIMD-friendly storage layout.
inline const RevVectors& rev_vectors() {
    static const RevVectors table = []() {
        RevVectors table{};

        // Generate the raw direction vectors, one `SOBOL_DEPTH`-vector per
        // dimension, exactly as the Rust `generate_direction_vectors`.
        std::array<std::array<uint32_t, SOBOL_DEPTH>, NUM_DIMENSIONS> vectors{};

        // Dimension 0: the van der Corput sequence.
        for (int i = 0; i < SOBOL_DEPTH; ++i) {
            vectors[0][i] = 1u << (SOBOL_BITS - 1 - i);
        }

        // Remaining dimensions from the direction numbers text.
        {
            size_t start = 0;
            size_t dim = 1;
            while (dim < NUM_DIMENSIONS) {
                size_t nl = DIRECTION_NUMBERS_TEXT.find('\n', start);
                std::string_view line = DIRECTION_NUMBERS_TEXT.substr(
                    start, nl == std::string_view::npos ? std::string_view::npos : nl - start);

                uint32_t a = 0;
                std::vector<uint32_t> m;
                if (parse_direction_numbers(line, a, m)) {
                    const int s = static_cast<int>(m.size());

                    // Fill in the leading direction numbers.
                    for (int i = 0; i < std::min(s, SOBOL_DEPTH); ++i) {
                        vectors[dim][i] = m[static_cast<size_t>(i)] << (SOBOL_BITS - (i + 1));
                    }

                    // Extend with the recurrence if fewer than SOBOL_DEPTH numbers were given.
                    if (s < SOBOL_DEPTH) {
                        for (int i = s; i < SOBOL_DEPTH; ++i) {
                            vectors[dim][i] = vectors[dim][i - s] ^ (vectors[dim][i - s] >> s);
                            for (int k = 1; k < s; ++k) {
                                vectors[dim][i] ^=
                                    ((a >> (s - 1 - k)) & 1u) * vectors[dim][i - k];
                            }
                        }
                    }

                    ++dim;
                }

                if (nl == std::string_view::npos) break;
                start = nl + 1;
            }
        }

        // Interleave each set of four dimensions, writing elements with
        // reversed bits and reversed depth order, to match the Rust `REV_VECTORS`.
        for (uint32_t set = 0; set < NUM_DIMENSION_SETS_4D; ++set) {
            for (int row = 0; row < SOBOL_DEPTH; ++row) {
                const int depth = SOBOL_DEPTH - 1 - row;
                for (int sub = 0; sub < 4; ++sub) {
                    table[set][row][sub] =
                        reverse_bits(vectors[set * 4 + sub][depth]);
                }
            }
        }

        return table;
    }();
    return table;
}

}  // namespace detail

//---------------------------------------------------------------------------
// Building blocks
//---------------------------------------------------------------------------

/// Compute one dimension of a single sample in the Sobol sequence
/// (result has reversed bits).
inline uint32_t sobol_rev(uint32_t sample_index_rev, uint32_t dimension) {
    assert(dimension < NUM_DIMENSIONS);

    const uint32_t dimension_set = dimension >> 2;
    const uint32_t sub_dimension = dimension & 0b11u;

    const auto& vecs = detail::rev_vectors()[dimension_set];
    uint32_t sobol = 0u;
    for (int i = 16; i < 32; ++i) {
        const uint32_t mask = 0u - ((sample_index_rev >> i) & 1u);
        sobol ^= mask & vecs[static_cast<size_t>(i - (32 - SOBOL_DEPTH))][sub_dimension];
    }
    return sobol;
}

/// Same as `sobol_rev` but returns four dimensions at once.
///
/// `dimension_set` indexes into sets of four dimensions:
///   * `0` -> [dim0, dim1, dim2, dim3]
///   * `1` -> [dim4, dim5, dim6, dim7]
///   * etc.
inline std::array<uint32_t, 4> sobol_int4_rev(uint32_t sample_index_rev,
                                              uint32_t dimension_set) {
    assert(dimension_set < NUM_DIMENSION_SETS_4D);

    const auto& vecs = detail::rev_vectors()[dimension_set];
    std::array<uint32_t, 4> sobol = {0u, 0u, 0u, 0u};
    for (int i = 16; i < 32; ++i) {
        const uint32_t mask = 0u - ((sample_index_rev >> i) & 1u);
        for (int sub = 0; sub < 4; ++sub) {
            sobol[sub] ^= mask & vecs[static_cast<size_t>(i - (32 - SOBOL_DEPTH))][sub];
        }
    }
    return sobol;
}

/// Scramble `n` using a hash function that closely approximates a reverse-bit
/// Owen scramble.  `scramble` must already be well randomized (e.g. pass it
/// through `hash()` first, or take it from a random source).
inline uint32_t owen_scramble_rev(uint32_t n_rev, uint32_t scramble) {
    n_rev ^= n_rev * 0x3d20adeau;
    n_rev += scramble;
    n_rev *= (scramble >> 16) | 1u;
    n_rev ^= n_rev * 0x05526c56u;
    n_rev ^= n_rev * 0x53a22864u;
    return n_rev;
}

/// Same as `owen_scramble_rev`, but on four integers at a time.
inline std::array<uint32_t, 4> owen_scramble_int4_rev(
    std::array<uint32_t, 4> n_rev, std::array<uint32_t, 4> scramble) {
    for (int sub = 0; sub < 4; ++sub) {
        n_rev[sub] ^= n_rev[sub] * 0x3d20adeau;
        n_rev[sub] += scramble[sub];
        n_rev[sub] *= (scramble[sub] >> 16) | 1u;
        n_rev[sub] ^= n_rev[sub] * 0x05526c56u;
        n_rev[sub] ^= n_rev[sub] * 0x53a22864u;
    }
    return n_rev;
}

/// A fast 32-bit hash function (from <https://github.com/skeeto/hash-prospector>).
inline uint32_t hash(uint32_t n) {
    n ^= 0xe6fe3bebu;  // So zero doesn't map to zero.
    n ^= n >> 16;
    n *= 0x7feb352du;
    n ^= n >> 15;
    n *= 0x846ca68bu;
    n ^= n >> 16;
    return n;
}

/// Same as `hash`, but on four numbers at once.
inline std::array<uint32_t, 4> hash_int4(std::array<uint32_t, 4> n) {
    for (int sub = 0; sub < 4; ++sub) {
        n[sub] ^= 0xe6fe3bebu;  // So zero doesn't map to zero.
        n[sub] ^= n[sub] >> 16;
        n[sub] *= 0x7feb352du;
        n[sub] ^= n[sub] >> 15;
        n[sub] *= 0x846ca68bu;
        n[sub] ^= n[sub] >> 16;
    }
    return n;
}

/// Convert a `u32` to a float in [0.0, 1.0).  Maps the full range of `u32`
/// onto the [0, 1) range.
inline float u32_to_f32_norm(uint32_t n) {
    const uint32_t bits = (n >> 9) | 0x3f800000u;
    float f;
    std::memcpy(&f, &bits, sizeof(f));
    return f - 1.0f;
}

//---------------------------------------------------------------------------
// Public API
//---------------------------------------------------------------------------

/// Compute one dimension of a single sample in the Sobol sequence.
///
/// `sample_index` specifies which sample in the Sobol sequence to compute
/// (a maximum of 2^16 samples is supported).  `dimension` specifies which
/// dimension to compute.  `seed` produces statistically independent Sobol
/// sequences.  Returns a number in the interval [0, 1).
inline float sample(uint32_t sample_index, uint32_t dimension, uint32_t seed) {
    assert(sample_index < (1u << 16));

    // Shuffle the index using the given seed to produce a unique statistically
    // independent Sobol sequence.
    const uint32_t shuffled_rev_index =
        owen_scramble_rev(reverse_bits(sample_index), hash(seed ^ 0x79c68e4au));

    const uint32_t sobol = sobol_rev(shuffled_rev_index, dimension);

    // Compute the scramble value for doing Owen scrambling.  The multiply on
    // `seed` is to avoid accidental cancellation with `dimension` on an
    // incrementing or otherwise structured seed.
    uint32_t scramble;
    {
        const uint32_t seed_ = seed * 0x9c8f2d3bu;
        const uint32_t ds = dimension >> 2;
        static constexpr std::array<uint32_t, 4> kScrambleConsts = {
            0x912f69bau, 0x174f18abu, 0x691e72cau, 0xb40cc1b8u};
        scramble = ds ^ seed_ ^ kScrambleConsts[dimension & 0b11u];
    }

    const uint32_t sobol_owen_rev = owen_scramble_rev(sobol, hash(scramble));

    return u32_to_f32_norm(reverse_bits(sobol_owen_rev));
}

/// Compute four dimensions of a single sample in the Sobol sequence.
///
/// Identical to calling `sample()` four times, but computes all four
/// dimensions at once.  `dimension_set` specifies which four dimensions to
/// compute: `0` yields the first four, `1` the second four, and so on.
inline std::array<float, 4> sample_4d(uint32_t sample_index, uint32_t dimension_set,
                                      uint32_t seed) {
    assert(sample_index < (1u << 16));

    // Shuffle the index using the given seed to produce a unique statistically
    // independent Sobol sequence.
    const uint32_t shuffled_rev_index =
        owen_scramble_rev(reverse_bits(sample_index), hash(seed ^ 0x79c68e4au));

    const std::array<uint32_t, 4> sobol = sobol_int4_rev(shuffled_rev_index, dimension_set);

    // Compute the scramble values for doing Owen scrambling.
    static constexpr std::array<uint32_t, 4> kScrambleConsts = {
        0x912f69bau, 0x174f18abu, 0x691e72cau, 0xb40cc1b8u};
    std::array<uint32_t, 4> scramble;
    const uint32_t seed_ = seed * 0x9c8f2d3bu;
    for (int sub = 0; sub < 4; ++sub) {
        scramble[sub] = dimension_set ^ seed_ ^ kScrambleConsts[sub];
    }

    const std::array<uint32_t, 4> sobol_owen_rev =
        owen_scramble_int4_rev(sobol, hash_int4(scramble));

    std::array<float, 4> out;
    for (int sub = 0; sub < 4; ++sub) {
        out[sub] = u32_to_f32_norm(reverse_bits(sobol_owen_rev[sub]));
    }
    return out;
}

}  // namespace sobol_burley
