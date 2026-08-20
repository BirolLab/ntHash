#pragma once

#include "../internal/internal.hpp"

#include <cstring>
#include <memory>
#include <stdexcept>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;

/**
 * Container for the two strand-specific forward-out and reverse-in masks.
 * Forward-out mask: k-th rotation of each base's seed.
 * Reverse-in mask: (k-1)-th rotation of each base's complement's seed.
 * Masks are cached for k, since k doesn't normally change in a single run.
 */
struct StrandMasks
{
  std::array<HASH_TYPE, 4> fwd_out{};
  std::array<HASH_TYPE, 4> rev_in{};

  constexpr StrandMasks() noexcept = default;

  StrandMasks(K_TYPE k)
  {
    static K_TYPE cached_k = 0;
    static std::array<HASH_TYPE, 4> cached_fwd;
    static std::array<HASH_TYPE, 4> cached_rev;
    if (k != cached_k) {
      constexpr std::array<HASH_TYPE, 4> seeds{
        internal::SEED_A, internal::SEED_C, internal::SEED_G, internal::SEED_T
      };
      for (size_t i = 0; i < seeds.size(); i++) {
        fwd_out[i] = internal::roll_next(seeds[i], k);
        rev_in[i] = internal::roll_next(seeds[seeds.size() - 1 - i], k - 1);
      }
      cached_fwd = fwd_out;
      cached_rev = rev_in;
      cached_k = k;
    } else {
      fwd_out = cached_fwd;
      rev_in = cached_rev;
    }
  }
};

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
[[nodiscard]] inline constexpr HASH_TYPE
base_forward_hash(const char* seq, K_TYPE k) noexcept
{
  HASH_TYPE h_val = 0;
  for (K_TYPE i = 0; i < k; i++) {
    const auto index = static_cast<unsigned char>(seq[i]);
    h_val = internal::roll_next(h_val) ^ internal::SEED_TAB[index];
  }
  return h_val;
}

/**
 * Perform a roll operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param out_mask Character mask to be removed (k-times roll of char_out)
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_forward_hash(HASH_TYPE fh_val,
                  HASH_TYPE out_mask,
                  unsigned char char_in) noexcept
{
  const auto in_mask = internal::SEED_TAB[char_in];
  return nthash::internal::roll_next(fh_val) ^ in_mask ^ out_mask;
}

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param char_out Character to be removed
 * @param in_mask Character mask to be included (k-th rotation of char_in)
 * @return Forward hash value rolled back
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_forward_hash(HASH_TYPE fh_val,
                  unsigned char char_out,
                  HASH_TYPE in_mask) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out];
  return nthash::internal::roll_back(fh_val ^ out_mask) ^ in_mask;
}

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
[[nodiscard]] inline constexpr HASH_TYPE
base_reverse_hash(const char* seq, K_TYPE k) noexcept
{
  HASH_TYPE h_val = 0;
  for (K_TYPE i = 0; i < k; i++) {
    const auto index =
      static_cast<unsigned char>(seq[k - 1 - i]) & internal::CP_OFF;
    h_val = nthash::internal::roll_next(h_val) ^ internal::SEED_TAB[index];
  }
  return h_val;
}

/**
 * Perform a roll operation on the reverse-complement.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param char_out Character to be removed
 * @param in_mask Character mask to be included (k-th rotation of char_in)
 * @return Rolled hash value for the reverse-complement
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_reverse_hash(HASH_TYPE rh_val,
                  unsigned char char_out,
                  HASH_TYPE in_mask) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out & internal::CP_OFF];
  return nthash::internal::roll_back(rh_val ^ out_mask) ^ in_mask;
}

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param out_mask Character mask to be removed (k-times rotation of char_out)
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_reverse_hash(HASH_TYPE rh_val,
                  HASH_TYPE out_mask,
                  unsigned char char_in) noexcept
{
  const auto in_mask = internal::SEED_TAB[char_in & internal::CP_OFF];
  return nthash::internal::roll_next(rh_val ^ out_mask) ^ in_mask;
}

inline HASH_TYPE
blind_fwd_out_mask(unsigned char c,
                   K_TYPE k,
                   const std::array<HASH_TYPE, 4>& fwd_out_mask)
{
  const auto loc = internal::CONVERT_TAB[c];
  return loc < 4 ? fwd_out_mask[loc]
                 : internal::roll_next(internal::SEED_TAB[c],
                                       static_cast<unsigned>(k));
}

inline HASH_TYPE
blind_rev_in_mask(unsigned char c,
                  K_TYPE k,
                  const std::array<HASH_TYPE, 4>& rev_in_mask)
{
  const auto loc = internal::CONVERT_TAB[c];
  return loc < 4 ? rev_in_mask[loc]
                 : internal::roll_next(internal::SEED_TAB[c & internal::CP_OFF],
                                       static_cast<unsigned>(k) - 1);
}

inline HASH_TYPE
blind_fwd_in_mask(unsigned char c,
                  K_TYPE k,
                  const std::array<HASH_TYPE, 4>& rev_in_mask)
{
  const auto loc = internal::CONVERT_TAB[c];
  return loc < 4 ? rev_in_mask[3 - loc]
                 : internal::roll_next(internal::SEED_TAB[c],
                                       static_cast<unsigned>(k) - 1);
}

inline HASH_TYPE
blind_rev_out_mask(unsigned char c,
                   K_TYPE k,
                   const std::array<HASH_TYPE, 4>& fwd_out_mask)
{
  const auto loc = internal::CONVERT_TAB[c];
  return loc < 4 ? fwd_out_mask[3 - loc]
                 : internal::roll_next(internal::SEED_TAB[c & internal::CP_OFF],
                                       static_cast<unsigned>(k));
}

} // namespace nthash::kmer
