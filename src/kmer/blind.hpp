#pragma once

#include "utils.hpp"

#include <stdexcept>
#include <vector>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;

class BlindNtHash
{

public:
  /**
   * Construct an ntHash object for hashing k-mers on-the-fly.
   * @param seq C-string of the data. Only \p k characters will be used,
   * starting from \p pos.
   * @param hash_num Number of hashes to generate per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  BlindNtHash(const char* seq, unsigned num_hashes, K_TYPE k, ssize_t pos = 0)
    : buffer(seq + pos, seq + pos + k)
    , buffer_idx(0)
    , pos(pos)
    , hash_arr(num_hashes)
    , rollk_tab(generate_rollk_table(k))
    , k_mult(static_cast<HASH_TYPE>(k) * internal::MULTISEED)
  {
    if (k == 0) {
      throw std::invalid_argument("BlindNtHash: k must be greater than 0");
    }
    fwd_hash = kmer::base_forward_hash(seq + pos, k);
    rev_hash = kmer::base_reverse_hash(seq + pos, k);
    internal::extend_hashes(fwd_hash, rev_hash, k_mult, hash_arr);
  }

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindNtHash object was constructed on, the provided character
   * \p char_in is used as the next base. Useful if you want to query for
   * possible paths in an implicit de Bruijn graph graph.
   */
  void roll(char char_in)
  {
    const auto char_out = static_cast<unsigned char>(buffer[buffer_idx]);
    const auto uchar_in = static_cast<unsigned char>(char_in);
    fwd_hash = kmer::next_forward_hash(fwd_hash, char_out, uchar_in, rollk_tab);
    rev_hash = kmer::next_reverse_hash(rev_hash, char_out, uchar_in, rollk_tab);
    internal::extend_hashes(fwd_hash, rev_hash, k_mult, hash_arr);
    buffer[buffer_idx] = char_in;
    if (++buffer_idx == get_k()) {
      buffer_idx = 0;
    }
    ++pos;
  }

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in)
  {
    const auto back_idx = (buffer_idx == 0) ? get_k() - 1 : buffer_idx - 1;
    const auto char_out = static_cast<unsigned char>(buffer[back_idx]);
    const auto uchar_in = static_cast<unsigned char>(char_in);
    fwd_hash = kmer::prev_forward_hash(fwd_hash, char_out, uchar_in, rollk_tab);
    rev_hash = kmer::prev_reverse_hash(rev_hash, char_out, uchar_in, rollk_tab);
    internal::extend_hashes(fwd_hash, rev_hash, k_mult, hash_arr);
    buffer[back_idx] = char_in;
    buffer_idx = back_idx;
    --pos;
  }

  /**
   * Like NtHash::peek(), but as if roll(char char_in) was called.
   */
  void peek(char char_in)
  {
    const auto char_out = static_cast<unsigned char>(buffer[buffer_idx]);
    const auto uchar_in = static_cast<unsigned char>(char_in);
    internal::extend_hashes(
      kmer::next_forward_hash(fwd_hash, char_out, uchar_in, rollk_tab),
      kmer::next_reverse_hash(rev_hash, char_out, uchar_in, rollk_tab),
      k_mult,
      hash_arr);
  }

  /**
   * Like peek(char char_in), but as if roll_back(char char_in) was called.
   */
  void peek_back(char char_in)
  {
    const auto back_idx = (buffer_idx == 0) ? get_k() - 1 : buffer_idx - 1;
    const auto char_out = static_cast<unsigned char>(buffer[back_idx]);
    const auto uchar_in = static_cast<unsigned char>(char_in);
    internal::extend_hashes(
      kmer::prev_forward_hash(fwd_hash, char_out, uchar_in, rollk_tab),
      kmer::prev_reverse_hash(rev_hash, char_out, uchar_in, rollk_tab),
      k_mult,
      hash_arr);
  }

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  [[nodiscard]] const HASH_TYPE* hashes() const noexcept
  {
    return hash_arr.data();
  }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  [[nodiscard]] ssize_t get_pos() const noexcept { return pos; }

  /**
   * Get the number of hashes generated per k-mer.
   * @return Number of hashes per k-mer
   */
  [[nodiscard]] unsigned get_hash_num() const noexcept
  {
    return hash_arr.size();
  }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  [[nodiscard]] K_TYPE get_k() const noexcept
  {
    return static_cast<K_TYPE>(buffer.size());
  }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  [[nodiscard]] HASH_TYPE get_forward_hash() const noexcept { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  [[nodiscard]] HASH_TYPE get_reverse_hash() const noexcept { return rev_hash; }

private:
  std::string buffer;
  size_t buffer_idx;
  ssize_t pos;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  std::vector<HASH_TYPE> hash_arr;
  const RollKTable& rollk_tab;
  const HASH_TYPE k_mult;
};

} // namespace nthash::kmer