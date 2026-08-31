#pragma once

#include "../internal/internal.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nthash::legacy {

using internal::HASH_TYPE;
using internal::K_TYPE;

alignas(64) inline constexpr auto DIMER_TAB =
  internal::generate_kmer_table<2, internal::rotl>();
alignas(64) inline constexpr auto TRIMER_TAB =
  internal::generate_kmer_table<3, internal::rotl>();
alignas(64) inline constexpr auto TETRAMER_TAB =
  internal::generate_kmer_table<4, internal::rotl>();

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
  std::size_t i = 0;
  for (; i + 4 <= k; i += 4) {
    const uint8_t index =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 6) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])] << 4) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 2])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 3])];
    hash = internal::rotl(hash, 4) ^ TETRAMER_TAB[index];
  }
  const auto rem = k - i;
  if (rem == 3) {
    const uint8_t index =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 4) |
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 2])];
    hash = internal::rotl(hash, 3) ^ TRIMER_TAB[index];
  } else if (rem == 2) {
    const uint8_t idx =
      (internal::CONVERT_TAB[static_cast<unsigned char>(seq[i])] << 2) |
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[i + 1])];
    hash = internal::rotl(hash, 2) ^ DIMER_TAB[idx];
  } else if (rem == 1) {
    hash = internal::rotl(hash) ^
           internal::SEED_TAB[static_cast<unsigned char>(seq[i])];
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
  const auto remainder = k % 4;

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
    hash = internal::rotl(hash, 4) ^ TETRAMER_TAB[index];
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
                  K_TYPE k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto h_out = internal::rotl(internal::SEED_TAB[char_out], k);
  return internal::rotl(fh_val) ^ h_out ^ internal::SEED_TAB[char_in];
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
                  K_TYPE k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto h_out = internal::SEED_TAB[char_out & internal::CP_OFF];
  const auto h_in = internal::SEED_TAB[char_in & internal::CP_OFF];
  return internal::rotr(rh_val ^ h_out) ^ internal::rotl(h_in, k - 1);
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
                  K_TYPE k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto h_out = internal::SEED_TAB[char_out];
  const auto h_in = internal::rotl(internal::SEED_TAB[char_in], k - 1);
  return internal::rotr(fh_val ^ h_out) ^ h_in;
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
                  K_TYPE k,
                  unsigned char char_out,
                  unsigned char char_in) noexcept
{
  const auto h_out = internal::SEED_TAB[char_out & internal::CP_OFF];
  const auto h_in = internal::SEED_TAB[char_in & internal::CP_OFF];
  return internal::rotl(rh_val ^ internal::rotl(h_out, k - 1)) ^ h_in;
}

} // namespace nthash::legacy