#pragma once

#include <cstdint>
#include <limits>
#include <string_view>

#include "kmer/blind.hpp"
#include "kmer/nthash.hpp"
#include "seed/blind.hpp"
#include "seed/seed.hpp"

namespace nthash {

/**
 * String representing the hash function's name. Only change if hash outputs
 * are different from the previous version. Useful for tracking differences in
 * saved hashes, e.g., in Bloom filters.
 */
inline constexpr std::string_view NTHASH_FN_NAME = "ntHash_v3.0";

static_assert(std::numeric_limits<uint64_t>::max() + 1 == 0,
              "Integers don't overflow on this platform which is necessary for "
              "ntHash canonical hash computation.");

// Expose core classes
using kmer::BlindNtHash;
using kmer::NtHash;
using seed::BlindSeedNtHash;
using seed::SeedNtHash;

} // namespace nthash
