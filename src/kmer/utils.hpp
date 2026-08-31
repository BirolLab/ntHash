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
    const auto rollk_a = internal::roll_next(internal::SEED_A, k);
    const auto rollk_c = internal::roll_next(internal::SEED_C, k);
    const auto rollk_g = internal::roll_next(internal::SEED_G, k);
    const auto rollk_t = internal::roll_next(internal::SEED_T, k);
    table['A'] = table['a'] = rollk_a;
    table['C'] = table['c'] = rollk_c;
    table['G'] = table['g'] = rollk_g;
    table['T'] = table['t'] = table['U'] = table['u'] = rollk_t;
    table['A' & internal::CP_OFF] = rollk_t;
    table['C' & internal::CP_OFF] = rollk_g;
    table['T' & internal::CP_OFF] = rollk_a;
    table['U' & internal::CP_OFF] = rollk_a;
    table['G' & internal::CP_OFF] = rollk_c;
    cached_k = k;
  }
  return table;
}

alignas(64) inline constexpr auto DIMER_TAB =
  internal::generate_kmer_table<2, internal::roll_next>();
alignas(64) inline constexpr auto TRIMER_TAB =
  internal::generate_kmer_table<3, internal::roll_next>();
alignas(64) inline constexpr auto TETRAMER_TAB =
  internal::generate_kmer_table<4, internal::roll_next>();

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
[[nodiscard]] inline constexpr HASH_TYPE
base_forward_hash(const char* seq, K_TYPE k) noexcept
{
  HASH_TYPE hash = 0;
  size_t i = 0;
  for (; i + 4 <= k; i += 4) {
    const uint8_t idx =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 6) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])] << 4) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 2])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 3])];
    hash = internal::roll_next(hash, 4) ^ TETRAMER_TAB[idx];
  }
  const std::size_t rem = k - i;
  if (rem == 3) {
    const uint8_t idx =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 4) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 2])];
    hash = internal::roll_next(hash, rem) ^ TRIMER_TAB[idx];
  } else if (rem == 2) {
    const uint8_t idx =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])];
    hash = internal::roll_next(hash, rem) ^ DIMER_TAB[idx];
  } else if (rem == 1) {
    const uint8_t idx = static_cast<unsigned char>(seq[i]);
    hash = internal::roll_next(hash) ^ internal::SEED_TAB[idx];
  }
  return hash;
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
  HASH_TYPE hash = 0;
  const size_t remainder = k % 4;
  if (remainder == 3) {
    const uint8_t idx =
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[k - 1])] << 4) |
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[k - 2])] << 2) |
      internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[k - 3])];
    hash = TRIMER_TAB[idx];
  } else if (remainder == 2) {
    const uint8_t idx =
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[k - 1])] << 2) |
      internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[k - 2])];
    hash = DIMER_TAB[idx];
  } else if (remainder == 1) {
    const auto b0 = static_cast<unsigned char>(seq[k - 1]);
    hash = internal::SEED_TAB[b0 & internal::CP_OFF];
  }
  for (int i = static_cast<int>(k - remainder) - 1; i >= 3; i -= 4) {
    const uint8_t index =
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 6) |
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[i - 1])] << 4) |
      (internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[i - 2])] << 2) |
      internal::RC_CONVERT_TAB[static_cast<unsigned char>(seq[i - 3])];
    hash = internal::roll_next(hash, 4) ^ TETRAMER_TAB[index];
  }
  return hash;
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
