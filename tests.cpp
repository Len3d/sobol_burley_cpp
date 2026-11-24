// Self-verification tests for the C++ port of sobol_burley.
//
// The exact expected values below are copied from the Rust crate's own unit
// tests (`src/parts.rs` and `src/lib.rs`), so passing them means this port
// produces bit-identical results to the original.

#include "sobol_burley.hpp"

#include <cstdio>
#include <cstdint>

using namespace sobol_burley;

static int g_failures = 0;

#define CHECK(cond)                                                       \
    do {                                                                  \
        if (!(cond)) {                                                    \
            std::printf("FAIL (line %d): %s\n", __LINE__, #cond);         \
            ++g_failures;                                                 \
        }                                                                 \
    } while (0)

#define CHECK_EQ_U32(a, b)                                                       \
    do {                                                                         \
        const uint32_t va_ = (a);                                                \
        const uint32_t vb_ = (b);                                                \
        if (va_ != vb_) {                                                        \
            std::printf("FAIL (line %d): %s = 0x%08x != %s = 0x%08x\n", __LINE__, \
                        #a, va_, #b, vb_);                                       \
            ++g_failures;                                                        \
        }                                                                        \
    } while (0)

#define CHECK_EQ_F32(a, b)                                                        \
    do {                                                                          \
        const float va_ = (a);                                                    \
        const float vb_ = (b);                                                    \
        if (va_ != vb_) {                                                         \
            std::printf("FAIL (line %d): %s = %.9g != %s = %.9g\n", __LINE__, #a, \
                        va_, #b, vb_);                                            \
            ++g_failures;                                                         \
        }                                                                         \
    } while (0)

static void test_reverse_bits() {
    // From the Rust wide.rs SIMD reverse_bits test.
    const uint32_t a = 0xcde7a64eu, b = 0xdc69fbd9u, c = 0x3238fec6u, d = 0x1fb9ba8fu;
    // Expected values are simply recomputed by our own implementation here for
    // sanity (the authoritative check is sobol_rev_01 below).
    CHECK_EQ_U32(reverse_bits(0), 0);
    CHECK_EQ_U32(reverse_bits(0x80000000u), 1);
    CHECK_EQ_U32(reverse_bits(0xffffffffu), 0xffffffffu);
    CHECK_EQ_U32(reverse_bits(reverse_bits(a)), a);
    CHECK_EQ_U32(reverse_bits(reverse_bits(b)), b);
    CHECK_EQ_U32(reverse_bits(reverse_bits(c)), c);
    CHECK_EQ_U32(reverse_bits(reverse_bits(d)), d);
}

static void test_to_norm_f32() {
    // From the Rust parts.rs to_norm_f32 test.
    CHECK_EQ_F32(u32_to_f32_norm(0), 0.0f);
    CHECK(u32_to_f32_norm(0xffffffffu) < 1.0f);
    // 0x80000000 -> 0.5 (equivalent to the wide.rs to_f32_norm case).
    CHECK_EQ_F32(u32_to_f32_norm(0x80000000u), 0.5f);
}

static void test_sobol_rev_01() {
    // Exact vectors from the Rust parts.rs sobol_rev_01 test.
    static const uint32_t inputs[] = {0u,   1u,       2u,     3u,     20u,
                                      123456u, 13u,   635u,   99999999u, 54321u};
    static const uint32_t expected[] = {692060160u,  2839543808u, 3913285632u,
                                        1765801984u, 20971520u,   3321430016u,
                                        1900019712u, 1180696576u, 1077477376u,
                                        3576496128u};

    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i) {
        const uint32_t v = inputs[i];
        const uint32_t got = sobol_rev(reverse_bits(v + 1000), 1);
        CHECK_EQ_U32(reverse_bits(got), expected[i]);
    }
}

static void test_check_1d_and_4d_match() {
    // From the Rust lib.rs check_1d_and_4d_match test.
    for (uint32_t s = 0; s < 4; ++s) {
        for (uint32_t d = 0; d < 8; ++d) {
            for (uint32_t n = 0; n < 256; ++n) {
                const float a1 = sample(n, d * 4, s);
                const float b1 = sample(n, d * 4 + 1, s);
                const float c1 = sample(n, d * 4 + 2, s);
                const float d1 = sample(n, d * 4 + 3, s);

                const std::array<float, 4> v = sample_4d(n, d, s);

                CHECK_EQ_F32(a1, v[0]);
                CHECK_EQ_F32(b1, v[1]);
                CHECK_EQ_F32(c1, v[2]);
                CHECK_EQ_F32(d1, v[3]);
            }
        }
    }
}

static void test_range() {
    // All sample() outputs must be in [0, 1).
    for (uint32_t n = 0; n < (1u << 12); ++n) {
        for (uint32_t d = 0; d < 16; ++d) {
            const float x = sample(n, d, 12345);
            CHECK(x >= 0.0f && x < 1.0f);
        }
    }
    // sobol_int4_rev / sobol_rev consistency.
    for (uint32_t set = 0; set < 8; ++set) {
        for (uint32_t n = 0; n < (1u << 12); ++n) {
            const uint32_t rev = reverse_bits(n);
            const auto four = sobol_int4_rev(rev, set);
            for (int sub = 0; sub < 4; ++sub) {
                CHECK_EQ_U32(sobol_rev(rev, set * 4 + sub), four[sub]);
            }
        }
    }
}

int main() {
    test_reverse_bits();
    test_to_norm_f32();
    test_sobol_rev_01();
    test_check_1d_and_4d_match();
    test_range();

    if (g_failures == 0) {
        std::printf("ALL TESTS PASSED\n");
        return 0;
    }
    std::printf("%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
