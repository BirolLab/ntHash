#pragma once

#include "utils.hpp"

#include <cstring>
#include <memory>
#include <stdexcept>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;
using internal::NUM_HASHES_TYPE;

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
         NUM_HASHES_TYPE num_hashes,
         K_TYPE k,
         size_t pos = 0)
    : seq(seq, seq_len)
    , num_hashes(num_hashes)
    , k(k)
    , pos(pos)
    , initialized(false)
    , masks(k)
    , hash_arr(new HASH_TYPE[num_hashes])
  {
    if (k == 0) {
      throw std::invalid_argument("NtHash: k must be greater than 0");
    }
    if (this->seq.size() < k) {
      throw std::invalid_argument(
        "NtHash: sequence length (" + std::to_string(this->seq.size()) +
        ") is smaller than k (" + std::to_string(k) + ")");
    }
    if (pos > this->seq.size() - k) {
      throw std::invalid_argument("NtHash: passed position (" +
                                  std::to_string(pos) +
                                  ") is larger than sequence length (" +
                                  std::to_string(this->seq.size()) + ")");
    }
  }

  /**
   * Construct an ntHash object for k-mers.
   * @param seq Sequence string
   * @param num_hashes Number of hashes to produce per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  NtHash(std::string_view seq,
         NUM_HASHES_TYPE num_hashes,
         K_TYPE k,
         size_t pos = 0)
    : NtHash(seq.data(), seq.size(), num_hashes, k, pos)
  {
  }

  NtHash(const NtHash& obj)
    : seq(obj.seq)
    , num_hashes(obj.num_hashes)
    , k(obj.k)
    , pos(obj.pos)
    , initialized(obj.initialized)
    , fwd_hash(obj.fwd_hash)
    , rev_hash(obj.rev_hash)
    , masks(obj.masks)
  {
    hash_arr = std::make_unique<HASH_TYPE[]>(num_hashes);
    const auto copy_size = num_hashes * sizeof(HASH_TYPE);
    std::memcpy(hash_arr.get(), obj.hash_arr.get(), copy_size);
  }

  NtHash(NtHash&&) = default;

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
    if (internal::SEED_TAB[(unsigned char)seq[pos + k]] == internal::SEED_N) {
      pos += k;
      return init();
    }
    const auto out_loc =
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[pos])];
    fwd_hash =
      kmer::next_forward_hash(fwd_hash,
                              masks.fwd_out[out_loc],
                              static_cast<unsigned char>(seq[pos + k]));
    const auto in_loc =
      internal::CONVERT_TAB[static_cast<unsigned char>(seq[pos + k])];
    rev_hash = kmer::next_reverse_hash(
      rev_hash, static_cast<unsigned char>(seq[pos]), masks.rev_in[in_loc]);
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
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
    const auto in_loc = internal::CONVERT_TAB[char_in];
    const auto out_loc = internal::CONVERT_TAB[char_out];

    fwd_hash =
      kmer::prev_forward_hash(fwd_hash, char_out, masks.rev_in[3 - in_loc]);
    rev_hash =
      kmer::prev_reverse_hash(rev_hash, masks.rev_in[out_loc], char_in);
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
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
      return init();
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (internal::SEED_TAB[char_in_u] == internal::SEED_N) {
      return false;
    }

    const auto char_out = static_cast<unsigned char>(seq[pos]);
    const auto out_loc = internal::CONVERT_TAB[char_out];
    const auto in_loc = internal::CONVERT_TAB[char_in_u];

    const auto fwd =
      kmer::next_forward_hash(fwd_hash, masks.fwd_out[out_loc], char_in_u);
    const auto rev =
      kmer::next_reverse_hash(rev_hash, char_out, masks.rev_in[in_loc]);
    internal::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
    return true;
  }

  /**
   * Like peek(), but as if roll_back on char_in was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back(char char_in)
  {
    if (!initialized && !init()) {
      return init();
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (internal::SEED_TAB[char_in_u] == internal::SEED_N) {
      return false;
    }

    const auto char_out = static_cast<unsigned char>(seq[pos + k - 1]);
    const auto in_loc = internal::CONVERT_TAB[char_in_u];
    const auto out_loc = internal::CONVERT_TAB[char_out];

    const auto fwd =
      kmer::prev_forward_hash(fwd_hash, char_out, masks.rev_in[3 - in_loc]);
    const auto rev =
      kmer::prev_reverse_hash(rev_hash, masks.rev_in[out_loc], char_in_u);
    internal::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
    return true;
  }

  /**
   * Get the array of current canonical hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const HASH_TYPE* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  size_t get_pos() const { return pos; }

  /**
   * Get the number of hashes generated per k-mer.
   * @return Number of hashes per k-mer
   */
  NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  K_TYPE get_k() const { return k; }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  HASH_TYPE get_forward_hash() const { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  HASH_TYPE get_reverse_hash() const { return rev_hash; }

private:
  std::string_view seq;
  NUM_HASHES_TYPE num_hashes;
  K_TYPE k;
  size_t pos;
  bool initialized;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  kmer::StrandMasks masks;
  std::unique_ptr<HASH_TYPE[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init()
  {
    size_t pos_n = 0;
    while (pos <= seq.size() - k &&
           internal::is_invalid_kmer(seq.data() + pos, k, pos_n)) {
      pos += pos_n + 1;
    }
    if (pos > seq.size() - k) {
      return false;
    }
    fwd_hash = kmer::base_forward_hash(seq.data() + pos, k);
    rev_hash = kmer::base_reverse_hash(seq.data() + pos, k);
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    initialized = true;
    return true;
  }
};

} // namespace nthash::kmer