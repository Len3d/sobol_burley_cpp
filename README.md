# sobol_burley (C++17)

A C++17 port of the Rust crate
[`sobol_burley`](https://github.com/cessen/sobol_burley) (v0.5.0) by Nathan
Vegdahl — a seedable Owen-scrambled Sobol sequence, based on Brent Burley's
paper *[Practical Hash-based Owen
Scrambling](http://www.jcgt.org/published/0009/04/01/)*, the improved hash from
*[Building a Better LK Hash](https://psychopath.io/post/2021_01_30_building_a_better_lk_hash)*,
and more dimensions due to [Kuo et al.](http://web.maths.unsw.edu.au/~fkuo/sobol/).

This port is **scalar only** (the original crate's SIMD branch is not ported).
`sample_4d()` still works and produces bit-identical results to the original,
just without the 4x SIMD speed-up.

## Design constraints (identical to the original)

* Maximum sequence length is 2^16.
* Maximum number of dimensions is 256 (workaround: reseed to "pad" dimensions).
* Only `f32` output is supported.

## Files

| File | Purpose |
|------|---------|
| `sobol_burley.hpp` | Header-only library (constants, direction-vector generation, `sample`, `sample_4d`, and the `parts`-style building blocks). |
| `direction_numbers_embedded.hpp` | Joe & Kuo direction numbers, embedded (generated from `direction_numbers/new-joe-kuo-6.256.txt`). |
| `plots.cpp` | Example reproducing the original `examples/plots.rs`; writes PBM scatter plots. |
| `tests.cpp` | Self-verification tests (exact expected values taken from the Rust crate's own unit tests). |
| `direction_numbers/` | Original Joe & Kuo data file (reference). |
| `licenses/` | MIT, Apache-2.0, and Joe & Kuo (3-clause BSD) licenses. |

## Usage

```cpp
#include "sobol_burley.hpp"

// Print 1024 3-dimensional points.
for (uint32_t i = 0; i < 1024; ++i) {
    float x = sobol_burley::sample(i, 0, 0);
    float y = sobol_burley::sample(i, 1, 0);
    float z = sobol_burley::sample(i, 2, 0);
    // ...
}
```

The first argument is the sample index, the second is the dimension index
(both zero-indexed), and outputs are in the interval [0, 1). The third argument
is a seed that produces statistically independent Sobol sequences.

## Building and running

```sh
cmake -S . -B build
cmake --build build --config Release
./build/plots            # writes sobol*.pbm into the current directory
ctest --test-dir build   # runs the self-verification tests
```

`plots` takes an optional output directory: `plots out/`.

## Verification

The unit tests in `tests.cpp` assert the exact values from the Rust crate's
own tests (`sobol_rev_01`, `to_norm_f32`, `check_1d_and_4d_match`), so passing
them confirms bit-identical output. The generated `sobol_*.pbm` plots can be
compared pixel-for-pixel against the reference images committed in the
original Rust repository.

## License

MIT OR Apache-2.0.  The Sobol direction numbers are adapted from work by
Stephen Joe and Frances Y. Kuo, licensed under the 3-clause BSD license (see
`licenses/JOE_KUO.txt`).
