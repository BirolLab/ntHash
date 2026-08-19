#pragma once

#include <cstdint>
#include <limits>

namespace nthash {

static_assert(std::numeric_limits<uint64_t>::max() + 1 == 0,
              "Integers don't overflow on this platform which is necessary for "
              "ntHash canonical hash computation.");

/**
 * String representing the hash function's name. Only change if hash outputs
 * are different from the previous version. Useful for tracking differences in
 * saved hashes, e.g., in Bloom filters.
 */
static const char* const NTHASH_FN_NAME = "ntHash_v2.5";

} // namespace nthash
