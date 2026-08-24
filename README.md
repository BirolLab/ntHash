[![Release](https://img.shields.io/github/release/BirolLab/ntHash.svg)](https://github.com/BirolLab/ntHash/releases)
[![Downloads](https://img.shields.io/github/downloads/BirolLab/ntHash/total?logo=github)](https://github.com/BirolLab/ntHash/archive/master.zip)
[![Issues](https://img.shields.io/github/issues/BirolLab/ntHash.svg)](https://github.com/BirolLab/ntHash/issues)

![Logo](nthash-logo.png)

ntHash is an efficient rolling hash function for k-mers and spaced seeds.

# Installation

You can simply download `nthash.hpp` and drop it directly into your project's include path. No compiling or linking is required.

If you prefer to install it system-wide using [Meson](https://mesonbuild.com/), run the following in the project root:

```shell
meson setup --buildtype=release --prefix=<PREFIX> build
meson install -C build
```

This will copy `nthash.hpp` to your specified `<PREFIX>/include` directory.

# Usage

To use ntHash in a C++ project:
- Include the library: `#include <nthash.hpp>`
- Compile your code with `-std=c++17` (and preferably `-O3`)

Refer to the [docs](https://birollab.github.io/ntHash/) for more information.

# Examples

## Object-Oriented API

Generally, the `nthash::NtHash` and `nthash::SeedNtHash` classes are the easiest way to hash sequences:

```C++
nthash::NtHash nth("TGACTGATCGAGTCGTACTAG", 1, 5);  // 1 hash per 5-mer
while (nth.roll()) {
    // use nth.hashes() for canonical hashes
    // nth.get_forward_hash() returns forward strand hashes
    // nth.get_reverse_hash() returns reverse strand hashes
}
```

```C++
std::vector<std::string> seeds = {"10101", "11011"};
nthash::SeedNtHash nth("TGACTGATCGAGTCGTACTAG", seeds, 3, 5);
while (nth.roll()) {
    // nth.hashes()[0] = "T#A#T"'s first hash
    // nth.hashes()[1] = "T#A#T"'s second hash
    // nth.hashes()[2] = "T#A#T"'s third hash
    // nth.hashes()[3] = "TG#CT"'s first hash
}
```

## Stateless Functional API

You can also use the underlying internal functions directly.

```C++
#include <nthash.hpp>

unsigned k = 5;
uint64_t hash = nthash::kmer::base_forward_hash("TGACT", k);

// Roll forward by dropping 'T' and adding 'G'
hash = nthash::kmer::next_forward_hash(hash, k, 'T', 'G');
```

## Functional API with Precomputated Tables

For maximum performance, you can generate a thread-local cache table to achieve O(1) rolling operations.

```C++
#include <nthash.hpp>

unsigned k = 5;
uint64_t hash = nthash::kmer::base_forward_hash("TGACT", k);

// Generate the table for k=5 (cached for k per-thread)
const auto& table = nthash::kmer::generate_rollk_table(k);

// Fast roll using the precomputed masks
hash = nthash::kmer::next_forward_hash(hash, 'T', 'G', table);
```

# For developers

If you would like to contribute to the development of ntHash, after forking/cloning the repo, create the `build` directory:

```
meson setup build
```

Compile the tests and benchmarking scripts using:

```
meson compile -C build
```

**Note:** ntHash is distributed as a single header file. If you make changes to any of the files in the src/ directory, you must regenerate the single header before committing. You can do this by running `ninja nthash.hpp -C build` in the project root.

Before sending a PR, please make sure that:

- the single header is up-to-date by running `ninja nthash.hpp` in build
- tests pass by running `meson test -v -C build` in the project directory
- code is formatted properly by running `ninja clang-format` in the `build` folder (requires `clang-format`)
- coding standards have been met by making sure running `ninja clang-tidy-check` in `build` returns no errors (requires `clang-tools`)
- documentation is up-to-date by running `ninja docs` in `build` (requires [doxygen](https://www.doxygen.nl/))

# Publications

Parham Kazemi, Johnathan Wong, Vladimir Nikolić, Hamid Mohamadi, René L Warren, Inanç Birol, ntHash2: recursive spaced seed hashing for nucleotide sequences, Bioinformatics, 2022;, btac564, [https://doi.org/10.1093/bioinformatics/btac564](https://doi.org/10.1093/bioinformatics/btac564)

Hamid Mohamadi, Justin Chu, Benjamin P Vandervalk, and Inanc Birol.
ntHash: recursive nucleotide hashing.
*Bioinformatics* (2016) 32 (22): 3492-3494.
[doi:10.1093/bioinformatics/btw397](http://dx.doi.org/10.1093/bioinformatics/btw397)
