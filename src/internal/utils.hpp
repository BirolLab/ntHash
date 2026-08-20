#pragma once

#include <climits>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace nthash::internal {

using HASH_TYPE = uint64_t;
constexpr unsigned int HASH_BITS = sizeof(HASH_TYPE) * CHAR_BIT;

/**
 * This lets us minimize object size. Good for performance if it's
 * copied in, e.g., DBG traversal.
 */
using K_TYPE = uint16_t;

// number of rotations per roll
constexpr unsigned int ROT_R = 7;

// number of shifts per roll
constexpr unsigned int SHIFT_C = 27; // C has to be >= 22

// shift for generating multiple hash values
constexpr unsigned int MULTISHIFT = 27;

// seed for generating multiple hash values
constexpr HASH_TYPE MULTISEED = 0x90b45d39fb6da1fa;

[[nodiscard]] inline constexpr HASH_TYPE
canonical(const HASH_TYPE fwd, const HASH_TYPE rev) noexcept
{
  return fwd + rev;
}

[[nodiscard]] inline constexpr HASH_TYPE
rotl(HASH_TYPE x, unsigned int r) noexcept
{
  return (x << (r & (HASH_BITS - 1))) | (x >> ((-r) & (HASH_BITS - 1)));
}

[[nodiscard]] inline constexpr HASH_TYPE
rotr(HASH_TYPE x, unsigned int r) noexcept
{
  return (x >> (r & (HASH_BITS - 1))) | (x << ((-r) & (HASH_BITS - 1)));
}

/**
 * Steps the hash forward (to the right) by d positions.
 */
[[nodiscard]] inline constexpr HASH_TYPE
roll_next(HASH_TYPE hash_value) noexcept
{
  hash_value = rotl(hash_value, ROT_R);
  return hash_value ^ (hash_value << SHIFT_C);
}

[[nodiscard]] inline constexpr HASH_TYPE
roll_next(HASH_TYPE hash_value, unsigned d) noexcept
{
  for (unsigned i = 0; i < d; i++) {
    hash_value = roll_next(hash_value);
  }
  return hash_value;
}

/**
 * Steps the hash backward (to the left) by d positions.
 */
[[nodiscard]] inline constexpr HASH_TYPE
roll_back(HASH_TYPE hash_value) noexcept
{
  HASH_TYPE y = hash_value;
  for (unsigned s = SHIFT_C; s < HASH_BITS; s *= 2) {
    y ^= (y << s);
  }
  return rotr(y, ROT_R);
}

[[nodiscard]] inline constexpr HASH_TYPE
roll_back(HASH_TYPE hash_value, unsigned d) noexcept
{
  for (unsigned i = 0; i < d; i++) {
    hash_value = roll_back(hash_value);
  }
  return hash_value;
}

inline void
extend_hashes(HASH_TYPE fwd_hash,
              HASH_TYPE rev_hash,
              HASH_TYPE k_mult,
              std::vector<HASH_TYPE>& hash_array) noexcept
{
  HASH_TYPE t_val;
  const auto base_hash = canonical(fwd_hash, rev_hash);
  hash_array[0] = base_hash;
  for (unsigned i = 1; i < hash_array.size(); i++) {
    t_val = base_hash * (i ^ k_mult);
    t_val ^= t_val >> MULTISHIFT;
    hash_array[i] = t_val;
  }
}

/**
 * Extend hash array using a base hash value.
 * @param fwd_hash Forward hash value
 * @param rev_hash Reverse hash value
 * @param k k-mer size
 * @param h Size of the resulting hash array (number of extra hashes minus one)
 * @param h_val Array of size h for storing the output hashes
 */
inline void
extend_hashes(HASH_TYPE fwd_hash,
              HASH_TYPE rev_hash,
              K_TYPE k,
              std::vector<HASH_TYPE>& hash_array) noexcept
{
  const auto k_mult = static_cast<HASH_TYPE>(k) * MULTISEED;
  return extend_hashes(fwd_hash, rev_hash, k_mult, hash_array);
}

} // namespace nthash::internal