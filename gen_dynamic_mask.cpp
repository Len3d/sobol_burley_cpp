// Generator for a self-contained dynamic `mask.h`.
//
// Reads the optimized keys (scramblingKeys, rankingKeys) from an existing
// mask.h, and emits a new mask.h that:
//   * embeds the reversed-bit Sobol direction numbers (Joe & Kuo) and a compact
//     dynamic generator for the seedable Owen-scrambled Sobol sequence
//     (sobol_burley, seed = 0), so the static sobol_4096spp_256d.h table is no
//     longer needed;
//   * supports any sample index in [0, 2^16) (the construction's bit depth),
//     i.e. far beyond the 4096 rows of the static table;
//   * keeps the exact same scramblingKeys / rankingKeys values.
//
// The dynamic sample function reproduces, for m < 4096, bit-identical integers
// to sobol_sequence[m][d] in the static table (both are the integer underlying
// sobol_burley::sample(m, d, seed = 0)).

#include "sobol_burley.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace sobol_burley;

static const char* kHelper = R"HLP(
static inline uint32_t db_reverse_bits(uint32_t n) {
    n = ((n & 0x55555555u) << 1) | ((n >> 1) & 0x55555555u);
    n = ((n & 0x33333333u) << 2) | ((n >> 2) & 0x33333333u);
    n = ((n & 0x0f0f0f0fu) << 4) | ((n >> 4) & 0x0f0f0f0fu);
    n = ((n & 0x00ff00ffu) << 8) | ((n >> 8) & 0x00ff00ffu);
    n = (n << 16) | (n >> 16);
    return n;
}

static inline uint32_t db_hash(uint32_t n) {
    n ^= 0xe6fe3bebu;  // So zero doesn't map to zero.
    n ^= n >> 16;
    n *= 0x7feb352du;
    n ^= n >> 15;
    n *= 0x846ca68bu;
    n ^= n >> 16;
    return n;
}

static inline uint32_t db_owen_scramble_rev(uint32_t n_rev, uint32_t scramble) {
    n_rev ^= n_rev * 0x3d20adeau;
    n_rev += scramble;
    n_rev *= (scramble >> 16) | 1u;
    n_rev ^= n_rev * 0x05526c56u;
    n_rev ^= n_rev * 0x53a22864u;
    return n_rev;
}

static inline uint32_t db_sobol_rev(uint32_t sample_index_rev, uint32_t dimension) {
    const uint32_t set = dimension >> 2;
    const uint32_t sub = dimension & 0b11u;
    uint32_t sobol = 0u;
    for (int i = 16; i < 32; ++i) {
        const uint32_t mask = 0u - ((sample_index_rev >> i) & 1u);
        sobol ^= mask & db_rev_vectors[set][i - 16][sub];
    }
    return sobol;
}

// The 32-bit integer underlying sobol_burley::sample(m, d, seed = 0).
// Valid for any m in [0, 2^16) and d in [0, 256).
static inline uint32_t db_sobol_int(uint32_t m, uint32_t d) {
    const uint32_t shuffled_rev_index =
        db_owen_scramble_rev(db_reverse_bits(m), db_hash(0x79c68e4au));
    const uint32_t sobol = db_sobol_rev(shuffled_rev_index, d);
    static const uint32_t kScrambleConsts[4] = {
        0x912f69bau, 0x174f18abu, 0x691e72cau, 0xb40cc1b8u};
    const uint32_t ds = d >> 2;
    const uint32_t scramble = ds ^ 0u ^ kScrambleConsts[d & 0b11u];
    const uint32_t sobol_owen_rev = db_owen_scramble_rev(sobol, db_hash(scramble));
    return db_reverse_bits(sobol_owen_rev);
}
)HLP";

// Copy a `static const ... = { ... };` block from the source text, starting at
// the declaration substring `decl`.
static std::string extract_block(const std::string& src, const std::string& decl) {
    size_t p = src.find(decl);
    if (p == std::string::npos) {
        std::fprintf(stderr, "error: declaration not found: %s\n", decl.c_str());
        std::exit(1);
    }
    size_t open = src.find('{', p);
    if (open == std::string::npos) {
        std::fprintf(stderr, "error: '{' not found after %s\n", decl.c_str());
        std::exit(1);
    }
    size_t close = src.find("};", open);
    if (close == std::string::npos) {
        std::fprintf(stderr, "error: '};' not found after %s\n", decl.c_str());
        std::exit(1);
    }
    return src.substr(p, close - p + 2);  // 'static const ... { ... };'
}

int main(int argc, char** argv) {
    const char* src_path = (argc > 1) ? argv[1] : "mask.h";
    const char* out_path = (argc > 2) ? argv[2] : "mask_dynamic.h";

    std::ifstream fin(src_path, std::ios::binary);
    if (!fin) {
        std::fprintf(stderr, "error: cannot open %s\n", src_path);
        return 1;
    }
    std::stringstream ss;
    ss << fin.rdbuf();
    std::string src = ss.str();
    fin.close();

    std::string scram = extract_block(src, "static const uint32_t scramblingKeys[128][128][16]");
    std::string rank  = extract_block(src, "static const uint32_t rankingKeys[128][128]");

    // Emit the reversed direction-number table.
    const auto& rv = detail::rev_vectors();
    std::stringstream tab;
    tab << "static const uint32_t db_rev_vectors[64][16][4] = {\n";
    for (int set = 0; set < 64; ++set) {
        tab << "  {";
        for (int row = 0; row < 16; ++row) {
            tab << "{";
            for (int sub = 0; sub < 4; ++sub) {
                tab << "0x" << std::hex;
                char buf[16];
                std::snprintf(buf, sizeof(buf), "%08x", rv[set][row][sub]);
                tab << buf << std::dec << (sub == 3 ? "" : ",");
            }
            tab << "}" << (row == 15 ? "" : ",");
        }
        tab << "}" << (set == 63 ? "" : ",") << "\n";
    }
    tab << "};\n";

    std::ofstream fout(out_path, std::ios::binary);
    if (!fout) {
        std::fprintf(stderr, "error: cannot open %s for writing\n", out_path);
        return 1;
    }

    fout << "#pragma once\n\n";
    fout << "// Auto-generated by gen_dynamic_mask. Do not edit.\n";
    fout << "//\n";
    fout << "// Self-contained dynamic mask: the per-pixel sample sets are the seedable\n";
    fout << "// Owen-scrambled Sobol sequence (sobol_burley, seed = 0), generated on the\n";
    fout << "// fly from the embedded reversed-bit direction numbers (Joe & Kuo). It needs\n";
    fout << "// no external sobol_4096spp_256d.h table and supports any sample index in\n";
    fout << "// [0, 2^16) (the construction's 16-bit depth).\n";
    fout << "\n#include <cstdint>\n\n";
    fout << tab.str() << "\n";
    fout << kHelper << "\n";
    fout << scram << "\n\n";
    fout << rank << "\n\n";
    fout << "// Per-pixel sample evaluation (Algorithm 1 of Heitz et al. 2019):\n";
    fout << "//   m = xor(sampleID, rankingKey)   -- permutes the sample order within the\n";
    fout << "//       power-of-two block (rankingKey < 16), preserving prefix low discrepancy;\n";
    fout << "//   p = xor(sobol_int(m), scramble) -- Owen scrambling with the per-pixel key.\n";
    fout << "// sampleID may be any index in [0, 2^16); d is used modulo 16 (the mask stores\n";
    fout << "// 16 scrambling-key dimensions).\n";
    fout << "static inline float sample(int i, int j, int sampleID, int d) {\n";
    fout << "    i = i & 127;\n";
    fout << "    j = j & 127;\n";
    fout << "    d = d & 15;\n";
    fout << "    uint32_t m = (uint32_t)sampleID ^ rankingKeys[i][j];\n";
    fout << "    uint32_t val = db_sobol_int(m, (uint32_t)d) ^ scramblingKeys[i][j][d];\n";
    fout << "    return (val + 0.5f) / 4294967296ULL;\n";
    fout << "}\n";
    fout.close();

    std::printf("Wrote dynamic mask to %s\n", out_path);
    return 0;
}
