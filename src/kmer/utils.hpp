#pragma once

#include "../internal/internal.hpp"

#include <array>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;
using RollKTable = std::array<internal::HASH_TYPE, internal::ASCII_SIZE>;

/**
 * Generates a thread-safe, lock-free lookup table for roll^k values.
 *
 * This function utilizes thread-local storage to cache precomputed rolls.
 * The table is only recalculated if the k-mer size differs from the
 * previously cached size for the calling thread, leading to zero allocation
 * overhead and O(1) lookups during hot loops.
 * @param k The k-mer size used to compute the roll^k masks
 * @return A const reference to the thread-local precomputed roll^k table
 */
[[nodiscard]] inline const RollKTable&
generate_rollk_table(K_TYPE k) noexcept
{
  thread_local K_TYPE cached_k = 0;
  thread_local RollKTable table{};
  if (k != cached_k) {
    const auto mask_A = internal::roll_next(internal::SEED_A, k);
    const auto mask_C = internal::roll_next(internal::SEED_C, k);
    const auto mask_G = internal::roll_next(internal::SEED_G, k);
    const auto mask_T = internal::roll_next(internal::SEED_T, k);
    table['A'] = table['a'] = mask_A;
    table['C'] = table['c'] = mask_C;
    table['G'] = table['g'] = mask_G;
    table['T'] = table['t'] = table['U'] = table['u'] = mask_T;
    table['A' & internal::CP_OFF] = mask_T;
    table['C' & internal::CP_OFF] = mask_G;
    table['T' & internal::CP_OFF] = mask_A;
    table['U' & internal::CP_OFF] = mask_A;
    table['G' & internal::CP_OFF] = mask_C;
    cached_k = k;
  }
  return table;
}

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
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous forward hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character leaving the sliding window
 * @param char_in Character entering the sliding window
 * @return Rolled forward hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_forward_hash(HASH_TYPE fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto out_mask = internal::roll_next(internal::SEED_TAB[char_out], k);
  const auto in_mask = internal::SEED_TAB[char_in];
  return nthash::internal::roll_next(fh_val) ^ in_mask ^ out_mask;
}

/**
 * Perform a roll operation on the forward strand using a precomputed roll^k
 * table.
 * @param fh_val Previous forward hash value computed for the sequence
 * @param char_out Character leaving the sliding window
 * @param char_in Character entering the sliding window
 * @param rollk_table roll^k table generated using generate_rollk_table(k)
 * @return Rolled forward hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_forward_hash(HASH_TYPE fh_val,
                  unsigned char char_out,
                  unsigned char char_in,
                  const RollKTable& rollk_table) noexcept
{
  const auto out_mask = rollk_table[char_out];
  const auto in_mask = internal::SEED_TAB[char_in];
  return nthash::internal::roll_next(fh_val) ^ in_mask ^ out_mask;
}

/**
 * Perform a roll operation on the reverse-complement strand by removing
 * char_out and including char_in from the forward sequence.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character leaving the forward sliding window
 * @param char_in Character entering the forward sliding window
 * @return Rolled reverse-complement hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_reverse_hash(HASH_TYPE rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out & internal::CP_OFF];
  const auto in_mask =
    internal::roll_next(internal::SEED_TAB[char_in & internal::CP_OFF], k);
  return nthash::internal::roll_back(rh_val ^ out_mask ^ in_mask);
}

/**
 * Perform a roll operation on the reverse-complement strand using a precomputed
 * roll^k table.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param char_out Character leaving the forward sliding window
 * @param char_in Character entering the forward sliding window
 * @param rollk_table roll^k table generated using generate_rollk_table(k)
 * @return Rolled reverse-complement hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
next_reverse_hash(HASH_TYPE rh_val,
                  unsigned char char_out,
                  unsigned char char_in,
                  const RollKTable& rollk_table) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out & internal::CP_OFF];
  const auto in_mask = rollk_table[char_in & internal::CP_OFF];
  return nthash::internal::roll_back(rh_val ^ out_mask ^ in_mask);
}

/**
 * Perform a backward roll operation on the forward strand by removing char_out
 * and including char_in.
 * @param fh_val Previous forward hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character leaving the sliding window
 * @param char_in Character entering the sliding window
 * @return Rolled forward hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_forward_hash(HASH_TYPE fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out];
  const auto in_mask = internal::roll_next(internal::SEED_TAB[char_in], k);
  return nthash::internal::roll_back(fh_val ^ out_mask ^ in_mask);
}

/**
 * Perform a backward roll operation on the forward strand using a precomputed
 * roll^k table.
 * @param fh_val Previous forward hash value computed for the sequence
 * @param char_out Character leaving the sliding window
 * @param char_in Character entering the sliding window
 * @param rollk_table roll^k table generated using generate_rollk_table(k)
 * @return Rolled forward hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_forward_hash(HASH_TYPE fh_val,
                  unsigned char char_out,
                  unsigned char char_in,
                  const RollKTable& rollk_table) noexcept
{
  const auto out_mask = internal::SEED_TAB[char_out];
  const auto in_mask = rollk_table[char_in];
  return nthash::internal::roll_back(fh_val ^ out_mask ^ in_mask);
}

/**
 * Perform a backward roll operation on the reverse-complement strand by
 * removing char_out and including char_in from the forward sequence.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character leaving the forward sliding window
 * @param char_in Character entering the forward sliding window
 * @return Rolled reverse-complement hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_reverse_hash(HASH_TYPE rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto out_mask =
    internal::roll_next(internal::SEED_TAB[char_out & internal::CP_OFF], k);
  const auto in_mask = internal::SEED_TAB[char_in & internal::CP_OFF];
  return nthash::internal::roll_next(rh_val) ^ in_mask ^ out_mask;
}

/**
 * Perform a backward roll operation on the reverse-complement strand using a
 * precomputed roll^k table.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param char_out Character leaving the forward sliding window
 * @param char_in Character entering the forward sliding window
 * @param rollk_table roll^k table generated using generate_rollk_table(k)
 * @return Rolled reverse-complement hash value
 */
[[nodiscard]] inline constexpr HASH_TYPE
prev_reverse_hash(HASH_TYPE rh_val,
                  unsigned char char_out,
                  unsigned char char_in,
                  const RollKTable& rollk_table) noexcept
{
  const auto out_mask = rollk_table[char_out & internal::CP_OFF];
  const auto in_mask = internal::SEED_TAB[char_in & internal::CP_OFF];
  return nthash::internal::roll_next(rh_val) ^ in_mask ^ out_mask;
}

} // namespace nthash::kmer
