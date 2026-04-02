[![Release](https://img.shields.io/github/release/BirolLab/ntHash.svg)](https://github.com/BirolLab/ntHash/releases)
[![Downloads](https://img.shields.io/github/downloads/BirolLab/ntHash/total?logo=github)](https://github.com/BirolLab/ntHash/archive/master.zip)
[![Issues](https://img.shields.io/github/issues/BirolLab/ntHash.svg)](https://github.com/BirolLab/ntHash/issues)

![Logo](nthash-logo.png)

ntHash is an efficient rolling hash function for k-mers and spaced seeds.

# Installation

Make sure [Meson](https://mesonbuild.com/) is installed on the system.

Download the repo (either from the releases section or close using `git clone https://github.com/BirolLab/ntHash`). Setup meson in an arbitrary directory (e.g. `build`), by running the following command in the project's root (include `--prefix=PREFIX` set the installation prefix to `PREFIX`):

```shell
meson setup --buildtype=release --prefix=PREFIX build
```

Then, install the project and its dependencies using:

```shell
meson install -C build 
```

This will install `include/nthash` and `lib/libnthash.a` to the installation prefix.

# Usage

To use ntHash in a C++ project:
- Import ntHash in the code using `#include <nthash/nthash.hpp>`
- Access ntHash classes from the `nthash` namespace
- Add the `include` directory (pass `-IPREFIX/include` to the compiler)
- Link the code with `libnthash.a` (i.e. pass `-LPREFIX/lib -lnthash` to the compiler, where `PREFIX` is the installation prefix)
- Compile your code with `-std=c++17` (and preferably `-O3`) enabled

Refer to [docs](https://birollab.github.io/ntHash/) for more information.

# Examples

Generally, the `nthash::NtHash` and `nthash::SeedNtHash` classes are used for hashing sequences:

```C++
nthash::NtHash nth("TGACTGATCGAGTCGTACTAG", 1, 5);  // 1 hash per 5-mer
while (nth.roll()) {
    // use nth.hashes() for canonical hashes
    //     nth.get_forward_hash() for forward strand hashes
    //     nth.get_reverse_hash() for reverse strand hashes
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

# For developers

If you would like to contribute to the development of ntHash, after forking/cloning the repo, create the `build` directory without the release flag:

```
meson setup build
```

Compile the code, tests, and benchmarking script using:

```
meson compile -C build
```

If compilation is successful, `libnthash.a` will be available in the `build` folder. The benchmarking script is also compiled as the `bench` binary file in `build`.

Before sending a PR, make sure that:

- tests pass by running `meson test` in the project directory
- code is formatted properly by running `ninja clang-format` in the `build` folder (requires `clang-format` to be available)
- coding standards have been met by making sure running `ninja clang-tidy-check` in `build` returns no errors (requires `clang-tools` to be installed)
- documentation is up-to-date by running `ninja docs` in `build` (requires [doxygen](https://www.doxygen.nl/))

# Publications

Parham Kazemi, Johnathan Wong, Vladimir Nikolić, Hamid Mohamadi, René L Warren, Inanç Birol, ntHash2: recursive spaced seed hashing for nucleotide sequences, Bioinformatics, 2022;, btac564, [https://doi.org/10.1093/bioinformatics/btac564](https://doi.org/10.1093/bioinformatics/btac564)

Hamid Mohamadi, Justin Chu, Benjamin P Vandervalk, and Inanc Birol.
ntHash: recursive nucleotide hashing.
*Bioinformatics* (2016) 32 (22): 3492-3494.
[doi:10.1093/bioinformatics/btw397](http://dx.doi.org/10.1093/bioinformatics/btw397)
