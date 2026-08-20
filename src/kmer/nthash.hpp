#pragma once

#include "utils.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;

class NtHash
{

public:
  /**
   * Construct an ntHash object for k-mers.
   * @param seq C-string containing sequence data
   * @param seq_len Length of the sequence
   * @param num_hashes Number of hashes to generate per k-mer
   * @param k K-mer size
   * @param pos Position in the sequence to start hashing from
   */
  NtHash(const char* seq,
         size_t seq_len,
         unsigned num_hashes,
         K_TYPE k,
         size_t pos = 0)
    : seq(seq, seq_len)
    , k(k)
    , pos(pos)
    , initialized(false)
    , masks(k)
    , hash_arr(num_hashes)
  {
    if (k == 0) {
      throw std::invalid_argument("NtHash: k must be greater than 0");
    }
    if (this->seq.size() < k) {
      std::ostringstream err;
      err << "NtHash: sequence length (" << this->seq.size() << ") ";
      err << "is smaller than k (" << k << ")";
      throw std::invalid_argument(err.str());
    }
    if (pos > this->seq.size() - k) {
      std::ostringstream err;
      err << "NtHash: passed position (" << pos << ") ";
      err << "is larger than sequence length (" << this->seq.size() << ")";
      throw std::invalid_argument(err.str());
    }
  }

  /**
   * Construct an ntHash object for k-mers.
   * @param seq Sequence string
   * @param num_hashes Number of hashes to produce per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  NtHash(std::string_view seq, unsigned num_hashes, K_TYPE k, size_t pos = 0)
    : NtHash(seq.data(), seq.size(), num_hashes, k, pos)
  {
  }

  /**
   * Calculate the hash values of current k-mer and advance to the next k-mer.
   * NtHash advances one nucleotide at a time until it finds a k-mer with valid
   * characters (ACGTU) and skips over those with invalid characters (non-ACGTU,
   * including N). This method must be called before hashes() is accessed, for
   * the first and every subsequent hashed kmer. get_pos() may be called at any
   * time to obtain the position of last hashed k-mer or the k-mer to be hashed
   * if roll() has never been called on this NtHash object. It is important to
   * note that the number of roll() calls is NOT necessarily equal to get_pos(),
   * if there are N's or invalid characters in the hashed sequence.
   * @return \p true on success and \p false otherwise
   */
  bool roll()
  {
    if (!initialized) {
      return init();
    }
    if (pos >= seq.size() - k) {
      return false;
    }
    const auto char_in = static_cast<unsigned char>(seq[pos + k]);
    if (internal::SEED_TAB[char_in] == internal::SEED_N) {
      pos += k;
      return init();
    }
    const auto char_out = static_cast<unsigned char>(seq[pos]);
    const auto mask_out = masks.fwd_out[internal::CONVERT_TAB[char_out]];
    fwd_hash = kmer::next_forward_hash(fwd_hash, mask_out, char_in);
    const auto mask_in = masks.rev_in[internal::CONVERT_TAB[char_in]];
    rev_hash = kmer::next_reverse_hash(rev_hash, char_out, mask_in);
    internal::extend_hashes(fwd_hash, rev_hash, k, hash_arr);
    ++pos;
    return true;
  }

  /**
   * Like the roll() function, but advance backwards.
   * @return \p true on success and \p false otherwise
   */
  bool roll_back()
  {
    if (!initialized) {
      return init();
    }
    if (pos == 0) {
      return false;
    }
    const auto char_in = static_cast<unsigned char>(seq[pos - 1]);
    if (internal::SEED_TAB[char_in] == internal::SEED_N) {
      if (pos >= k) {
        pos -= k;
        return init();
      }
      return false;
    }
    const auto char_out = static_cast<unsigned char>(seq[pos + k - 1]);
    const auto mask_in = masks.rev_in[3 - internal::CONVERT_TAB[char_in]];
    fwd_hash = kmer::prev_forward_hash(fwd_hash, char_out, mask_in);
    const auto mask_out = masks.rev_in[internal::CONVERT_TAB[char_out]];
    rev_hash = kmer::prev_reverse_hash(rev_hash, mask_out, char_in);
    internal::extend_hashes(fwd_hash, rev_hash, k, hash_arr);
    --pos;
    return true;
  }

  /**
   * Peeks the hash values as if roll() was called (without advancing the
   * NtHash object. The peeked hash values can be obtained through the
   * hashes() method.
   * @return \p true on success and \p false otherwise
   */
  bool peek()
  {
    if (pos >= seq.size() - k) {
      return false;
    }
    return peek(seq[pos + k]);
  }

  /**
   * Like peek(), but as if roll_back() was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back()
  {
    if (pos == 0) {
      return false;
    }
    return peek_back(seq[pos - 1]);
  }

  /**
   * Peeks the hash values as if roll() was called for char_in (without
   * advancing the NtHash object. The peeked hash values can be obtained through
   * the hashes() method.
   * @return \p true on success and \p false otherwise
   */
  bool peek(char char_in)
  {
    if (!initialized && !init()) {
      return false;
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (internal::SEED_TAB[char_in_u] == internal::SEED_N) {
      return false;
    }
    const auto char_out = static_cast<unsigned char>(seq[pos]);
    const auto mask_out = masks.fwd_out[internal::CONVERT_TAB[char_out]];
    const auto fwd = kmer::next_forward_hash(fwd_hash, mask_out, char_in_u);
    const auto mask_in = masks.rev_in[internal::CONVERT_TAB[char_in_u]];
    const auto rev = kmer::next_reverse_hash(rev_hash, char_out, mask_in);
    internal::extend_hashes(fwd, rev, k, hash_arr);
    return true;
  }

  /**
   * Like peek(), but as if roll_back on char_in was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back(char char_in)
  {
    if (!initialized && !init()) {
      return false;
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (internal::SEED_TAB[char_in_u] == internal::SEED_N) {
      return false;
    }
    const auto char_out = static_cast<unsigned char>(seq[pos + k - 1]);
    const auto mask_in = masks.rev_in[3 - internal::CONVERT_TAB[char_in_u]];
    const auto fwd = kmer::prev_forward_hash(fwd_hash, char_out, mask_in);
    const auto mask_out = masks.rev_in[internal::CONVERT_TAB[char_out]];
    const auto rev = kmer::prev_reverse_hash(rev_hash, mask_out, char_in_u);
    internal::extend_hashes(fwd, rev, k, hash_arr);
    return true;
  }

  /**
   * Get the array of current canonical hash values (length = \p get_hash_num())
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
  [[nodiscard]] size_t get_pos() const noexcept { return pos; }

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
  [[nodiscard]] K_TYPE get_k() const noexcept { return k; }

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
  std::string_view seq;
  K_TYPE k;
  size_t pos;
  bool initialized;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  kmer::StrandMasks masks;
  std::vector<HASH_TYPE> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init()
  {
    while (pos <= seq.size() - k) {
      bool valid = true;
      for (size_t i = k; i > 0 && valid; --i) {
        const auto c = static_cast<unsigned char>(seq[pos + i - 1]);
        if (internal::SEED_TAB[c] == internal::SEED_N) {
          pos += i;
          valid = false;
        }
      }
      if (valid) {
        fwd_hash = kmer::base_forward_hash(seq.data() + pos, k);
        rev_hash = kmer::base_reverse_hash(seq.data() + pos, k);
        internal::extend_hashes(fwd_hash, rev_hash, k, hash_arr);
        initialized = true;
        return true;
      }
    }
    return false;
  }
};

} // namespace nthash::kmer