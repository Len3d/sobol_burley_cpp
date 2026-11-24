// Example program that reproduces the Rust crate's `examples/plots.rs`.
//
// It writes plain-Sobol scatter plots (`sobol.pbm`, `sobol_0.pbm` ... `sobol_19.pbm`)
// and Owen-scrambled, seed-shuffled Sobol scatter plots (`sobol_burley_0.pbm` ...
// `sobol_burley_19.pbm`) as binary PBM images (P1, 512x512).
//
// Usage: plots [output_directory]   (defaults to the current directory)

#include "sobol_burley.hpp"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

using namespace sobol_burley;

namespace {

constexpr size_t X_RES = 512;
constexpr size_t Y_RES = 512;
constexpr uint32_t POINTS = 1u << 10;

// Plain (unscrambled) Sobol sample in [0, 1).
inline float sobol(uint32_t i, uint32_t dimension) {
    const uint32_t sobol_int = sobol_rev(reverse_bits(i), dimension);
    return u32_to_f32_norm(reverse_bits(sobol_int));
}

template <typename XFn, typename YFn>
void plot(XFn x_fn, YFn y_fn, const std::string& path) {
    std::vector<uint8_t> image(X_RES * Y_RES, 1u);
    for (uint32_t i = 0; i < POINTS; ++i) {
        const size_t x = static_cast<size_t>(x_fn(i) * static_cast<float>(X_RES - 1));
        const size_t y = static_cast<size_t>(y_fn(i) * static_cast<float>(Y_RES - 1));
        image[y * X_RES + x] = 0u;
    }

    FILE* f = fopen(path.c_str(), "wb");
    if (!f) {
        std::fprintf(stderr, "error: cannot open %s for writing\n", path.c_str());
        std::exit(1);
    }
    std::fprintf(f, "P1\n%zu %zu\n\n", X_RES, Y_RES);

    // Write the image as lines of up to 80 bits, one bit per character,
    // followed by a newline after each line (matching the Rust example).
    for (size_t i = 0; i < image.size(); i += 80) {
        const size_t end = std::min(i + 80, image.size());
        for (size_t j = i; j < end; ++j) {
            std::fputc(image[j] == 0u ? '0' : '1', f);
        }
        std::fputc('\n', f);
    }

    std::fclose(f);
}

}  // namespace

int main(int argc, char** argv) {
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    const auto path = [&outdir](const std::string& name) { return outdir + "/" + name; };

    // Plain sobol.
    plot([&](uint32_t i) { return sobol(i, 0); },
         [&](uint32_t i) { return sobol(i, 1); }, path("sobol.pbm"));
    for (int d = 0; d < 20; ++d) {
        plot([&](uint32_t i) { return sobol(i, static_cast<uint32_t>(d)); },
             [&](uint32_t i) { return sobol(i, static_cast<uint32_t>(d + 1)); },
             path("sobol_" + std::to_string(d) + ".pbm"));
    }

    // Scrambled shuffled sobol (Burley).
    for (int d = 0; d < 20; ++d) {
        plot([&](uint32_t i) { return sample(i, static_cast<uint32_t>(d), 0); },
             [&](uint32_t i) { return sample(i, static_cast<uint32_t>(d + 1), 0); },
             path("sobol_burley_" + std::to_string(d) + ".pbm"));
    }

    std::printf("Wrote 41 PBM images to %s\n", outdir.c_str());
    return 0;
}
