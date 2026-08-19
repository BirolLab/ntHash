#pragma once

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "tables.hpp"

namespace nthash::utils {

using HASH_TYPE = uint64_t;
constexpr unsigned int HASH_BITS = sizeof(HASH_TYPE) * CHAR_BIT;

/**
 * This lets us minimize object size. Good for performance if it's
 * copied in, e.g., DBG traversal.
 */
using K_TYPE = uint16_t;
using NUM_HASHES_TYPE = uint8_t;

// number of rotations per roll
constexpr unsigned int ROT_R = 7;

// number of shifts per roll
constexpr unsigned int SHIFT_C = 27; // C has to be >= 22

// shift for generating multiple hash values
constexpr unsigned int MULTISHIFT = 27;

// seed for generating multiple hash values
constexpr HASH_TYPE MULTISEED = 0x90b45d39fb6da1fa;

inline HASH_TYPE
canonical(const HASH_TYPE fwd, const HASH_TYPE rev)
{
  return fwd + rev;
}

inline HASH_TYPE
rotl(HASH_TYPE x, unsigned int r)
{
  return (x << (r & (HASH_BITS - 1))) | (x >> ((-r) & (HASH_BITS - 1)));
}

inline HASH_TYPE
rotr(HASH_TYPE x, unsigned int r)
{
  return (x >> (r & (HASH_BITS - 1))) | (x << ((-r) & (HASH_BITS - 1)));
}

/**
 * Steps the hash forward (to the right) by d positions.
 */
inline HASH_TYPE
roll_next(HASH_TYPE hash_value)
{
  hash_value = rotl(hash_value, ROT_R);
  hash_value ^= (hash_value << SHIFT_C);
  return hash_value;
}

inline HASH_TYPE
roll_next(HASH_TYPE hash_value, unsigned d)
{
  for (unsigned i = 0; i < d; i++) {
    hash_value = roll_next(hash_value);
  }
  return hash_value;
}

/**
 * Steps the hash backward (to the left) by d positions.
 */
inline HASH_TYPE
roll_back(HASH_TYPE hash_value)
{
  HASH_TYPE y = hash_value;
  for (unsigned s = SHIFT_C; s < HASH_BITS; s *= 2) {
    y ^= (y << s);
  }
  return rotr(y, ROT_R);
}

inline HASH_TYPE
roll_back(HASH_TYPE hash_value, unsigned d)
{
  for (unsigned i = 0; i < d; ++i) {
    hash_value = roll_back(hash_value);
  }
  return hash_value;
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
              NUM_HASHES_TYPE h,
              HASH_TYPE* hash_array)
{
  HASH_TYPE t_val;
  hash_array[0] = canonical(fwd_hash, rev_hash);
  const auto k_mult = static_cast<HASH_TYPE>(k) * MULTISEED;
  for (unsigned i = 1; i < h; i++) {
    t_val = hash_array[0] * (i ^ k_mult);
    t_val ^= t_val >> MULTISHIFT;
    hash_array[i] = t_val;
  }
}

/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
inline bool
is_invalid_kmer(const char* seq, unsigned k, size_t& pos_n)
{
  for (size_t i = k; i-- > 0;) {
    if (tables::SEED_TAB[(unsigned char)seq[i]] == tables::SEED_N) {
      pos_n = i;
      return true;
    }
  }
  return false;
}

} // namespace nthash