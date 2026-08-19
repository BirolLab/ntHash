#pragma once

#include "tables.hpp"
#include "utils.hpp"

#include <cstring>
#include <deque>
#include <memory>
#include <stdexcept>

namespace nthash {

using utils::HASH_TYPE;
using utils::K_TYPE;
using utils::NUM_HASHES_TYPE;

namespace kmer {

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
inline HASH_TYPE
base_forward_hash(const char* seq, K_TYPE k)
{
  HASH_TYPE h_val = 0;
  for (size_t i = 0; i < k; i++) {
    h_val = nthash::utils::roll_next(h_val);
    h_val ^= nthash::tables::SEED_TAB[static_cast<unsigned char>(seq[i])];
  }
  return h_val;
}

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous hash value computed for the sequence
 * @param out_mask Character mask to be removed (k-times rotation of char_out)
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
inline HASH_TYPE
next_forward_hash(HASH_TYPE fh_val, HASH_TYPE out_mask, unsigned char char_in)
{
  auto h_val = nthash::utils::roll_next(fh_val);
  h_val ^= out_mask;
  h_val ^= nthash::tables::SEED_TAB[char_in];
  return h_val;
}

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param char_out Character to be removed
 * @param in_mask Character mask to be included (k-th rotation of char_in)
 * @return Forward hash value rolled back
 */
inline HASH_TYPE
prev_forward_hash(HASH_TYPE fh_val, unsigned char char_out, HASH_TYPE in_mask)
{
  auto h_val = fh_val ^ nthash::tables::SEED_TAB[char_out];
  h_val = nthash::utils::roll_back(h_val);
  h_val ^= in_mask;
  return h_val;
}

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
inline HASH_TYPE
base_reverse_hash(const char* seq, K_TYPE k)
{
  HASH_TYPE h_val = 0;
  for (K_TYPE i = 0; i < k; i++) {
    h_val = nthash::utils::roll_next(h_val);
    const auto index = seq[k - 1 - i] & nthash::tables::CP_OFF;
    h_val ^= nthash::tables::SEED_TAB[index];
  }
  return h_val;
}

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param char_out Character to be removed
 * @param in_mask Character mask to be included (k-th rotation of char_in)
 * @return Rolled hash value for the reverse-complement
 */
inline HASH_TYPE
next_reverse_hash(HASH_TYPE rh_val, unsigned char char_out, HASH_TYPE in_mask)
{
  const auto index = char_out & nthash::tables::CP_OFF;
  auto h_val = rh_val ^ nthash::tables::SEED_TAB[index];
  h_val = nthash::utils::roll_back(h_val);
  h_val ^= in_mask;
  return h_val;
}

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param out_mask Character mask to be removed (k-times rotation of char_out)
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
inline HASH_TYPE
prev_reverse_hash(HASH_TYPE rh_val, HASH_TYPE out_mask, unsigned char char_in)
{
  auto h_val = rh_val ^ out_mask;
  h_val = nthash::utils::roll_next(h_val);
  h_val ^= nthash::tables::SEED_TAB[char_in & nthash::tables::CP_OFF];
  return h_val;
}

inline void
init_strand_masks(K_TYPE k,
                  HASH_TYPE (&fwd_out_mask)[4],
                  HASH_TYPE (&rev_in_mask)[4])
{
  constexpr HASH_TYPE BASE_SEEDS[4] = {
    tables::SEED_A, tables::SEED_C, tables::SEED_G, tables::SEED_T
  };
  for (unsigned i = 0; i < 4; i++) {
    fwd_out_mask[i] = utils::roll_next(BASE_SEEDS[i], k);
    rev_in_mask[i] = utils::roll_next(BASE_SEEDS[3 - i], k - 1);
  }
}

inline HASH_TYPE
blind_fwd_out_mask(unsigned char c,
                   K_TYPE k,
                   const HASH_TYPE (&fwd_out_mask)[4])
{
  const auto loc = tables::CONVERT_TAB[c];
  return loc < 4
           ? fwd_out_mask[loc]
           : utils::roll_next(tables::SEED_TAB[c], static_cast<unsigned>(k));
}

inline HASH_TYPE
blind_rev_in_mask(unsigned char c, K_TYPE k, const HASH_TYPE (&rev_in_mask)[4])
{
  const auto loc = tables::CONVERT_TAB[c];
  return loc < 4 ? rev_in_mask[loc]
                 : utils::roll_next(tables::SEED_TAB[c & tables::CP_OFF],
                                    static_cast<unsigned>(k) - 1);
}

inline HASH_TYPE
blind_fwd_in_mask(unsigned char c, K_TYPE k, const HASH_TYPE (&rev_in_mask)[4])
{
  const auto loc = tables::CONVERT_TAB[c];
  return loc < 4 ? rev_in_mask[3 - loc]
                 : utils::roll_next(tables::SEED_TAB[c],
                                    static_cast<unsigned>(k) - 1);
}

inline HASH_TYPE
blind_rev_out_mask(unsigned char c,
                   K_TYPE k,
                   const HASH_TYPE (&fwd_out_mask)[4])
{
  const auto loc = tables::CONVERT_TAB[c];
  return loc < 4 ? fwd_out_mask[3 - loc]
                 : utils::roll_next(tables::SEED_TAB[c & tables::CP_OFF],
                                    static_cast<unsigned>(k));
}

} // namespace kmer

class NtHash
{

private:
  std::string_view seq;
  NUM_HASHES_TYPE num_hashes;
  K_TYPE k;
  size_t pos;
  bool initialized;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  HASH_TYPE fwd_out_mask[4];
  HASH_TYPE rev_in_mask[4];
  std::unique_ptr<HASH_TYPE[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init()
  {
    size_t pos_n = 0;
    while (pos <= seq.size() - k &&
           utils::is_invalid_kmer(seq.data() + pos, k, pos_n)) {
      pos += pos_n + 1;
    }
    if (pos > seq.size() - k) {
      return false;
    }
    fwd_hash = kmer::base_forward_hash(seq.data() + pos, k);
    rev_hash = kmer::base_reverse_hash(seq.data() + pos, k);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    initialized = true;
    return true;
  }

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
    kmer::init_strand_masks(k, fwd_out_mask, rev_in_mask);
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
    , hash_arr(new HASH_TYPE[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(HASH_TYPE));
    std::memcpy(fwd_out_mask, obj.fwd_out_mask, sizeof(fwd_out_mask));
    std::memcpy(rev_in_mask, obj.rev_in_mask, sizeof(rev_in_mask));
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
    if (tables::SEED_TAB[(unsigned char)seq[pos + k]] == tables::SEED_N) {
      pos += k;
      return init();
    }
    const auto out_loc =
      tables::CONVERT_TAB[static_cast<unsigned char>(seq[pos])];
    fwd_hash =
      kmer::next_forward_hash(fwd_hash,
                              fwd_out_mask[out_loc],
                              static_cast<unsigned char>(seq[pos + k]));
    const auto in_loc =
      tables::CONVERT_TAB[static_cast<unsigned char>(seq[pos + k])];
    rev_hash = kmer::next_reverse_hash(
      rev_hash, static_cast<unsigned char>(seq[pos]), rev_in_mask[in_loc]);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
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
    if (tables::SEED_TAB[char_in] == tables::SEED_N) {
      if (pos >= k) {
        pos -= k;
        return init();
      }
      return false;
    }

    const auto char_out = static_cast<unsigned char>(seq[pos + k - 1]);
    const auto in_loc = tables::CONVERT_TAB[char_in];
    const auto out_loc = tables::CONVERT_TAB[char_out];

    fwd_hash =
      kmer::prev_forward_hash(fwd_hash, char_out, rev_in_mask[3 - in_loc]);
    rev_hash = kmer::prev_reverse_hash(rev_hash, rev_in_mask[out_loc], char_in);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
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
    if (!initialized) {
      return init();
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (tables::SEED_TAB[char_in_u] == tables::SEED_N) {
      return false;
    }

    const auto char_out = static_cast<unsigned char>(seq[pos]);
    const auto out_loc = tables::CONVERT_TAB[char_out];
    const auto in_loc = tables::CONVERT_TAB[char_in_u];

    const auto fwd =
      kmer::next_forward_hash(fwd_hash, fwd_out_mask[out_loc], char_in_u);
    const auto rev =
      kmer::next_reverse_hash(rev_hash, char_out, rev_in_mask[in_loc]);
    utils::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
    return true;
  }

  /**
   * Like peek(), but as if roll_back on char_in was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back(char char_in)
  {
    if (!initialized) {
      return init();
    }
    const auto char_in_u = static_cast<unsigned char>(char_in);
    if (tables::SEED_TAB[char_in_u] == tables::SEED_N) {
      return false;
    }

    const auto char_out = static_cast<unsigned char>(seq[pos + k - 1]);
    const auto in_loc = tables::CONVERT_TAB[char_in_u];
    const auto out_loc = tables::CONVERT_TAB[char_out];

    const auto fwd =
      kmer::prev_forward_hash(fwd_hash, char_out, rev_in_mask[3 - in_loc]);
    const auto rev =
      kmer::prev_reverse_hash(rev_hash, rev_in_mask[out_loc], char_in_u);
    utils::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
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
};

class BlindNtHash
{

public:
  /**
   * Construct an ntHash object for hashing k-mers on-the-fly.
   * @param seq C-string of the data. Only the first \p k characters will be
   * used, starting from \p pos.
   * @param hash_num Number of hashes to generate per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  BlindNtHash(const char* seq,
              NUM_HASHES_TYPE num_hashes,
              K_TYPE k,
              ssize_t pos = 0)
    : seq(seq + pos, seq + pos + k)
    , num_hashes(num_hashes)
    , pos(pos)
    , hash_arr(new HASH_TYPE[num_hashes])
  {
    if (k == 0) {
      throw std::invalid_argument("BlindNtHash: k must be greater than 0");
    }
    kmer::init_strand_masks(k, fwd_out_mask, rev_in_mask);
    fwd_hash = kmer::base_forward_hash(seq + pos, k);
    rev_hash = kmer::base_reverse_hash(seq + pos, k);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  }

  BlindNtHash(const BlindNtHash& obj)
    : seq(obj.seq)
    , num_hashes(obj.num_hashes)
    , pos(obj.pos)
    , fwd_hash(obj.fwd_hash)
    , rev_hash(obj.rev_hash)
    , hash_arr(new HASH_TYPE[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(HASH_TYPE));
    std::memcpy(fwd_out_mask, obj.fwd_out_mask, sizeof(fwd_out_mask));
    std::memcpy(rev_in_mask, obj.rev_in_mask, sizeof(rev_in_mask));
  }

  BlindNtHash(BlindNtHash&&) = default;

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindNtHash object was constructed on, the provided character
   * \p char_in is used as the next base. Useful if you want to query for
   * possible paths in an implicit de Bruijn graph graph.
   */
  void roll(char char_in)
  {
    const auto k = static_cast<K_TYPE>(seq.size());
    const auto char_out_u = static_cast<unsigned char>(seq.front());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto out_mask = kmer::blind_fwd_out_mask(char_out_u, k, fwd_out_mask);
    const auto in_mask = kmer::blind_rev_in_mask(char_in_u, k, rev_in_mask);

    fwd_hash = kmer::next_forward_hash(fwd_hash, out_mask, char_in_u);
    rev_hash = kmer::next_reverse_hash(rev_hash, char_out_u, in_mask);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    seq.pop_front();
    seq.push_back(char_in);
    ++pos;
  }

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in)
  {
    const auto k = static_cast<K_TYPE>(seq.size());
    const auto char_out_u = static_cast<unsigned char>(seq.back());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto in_mask = kmer::blind_fwd_in_mask(char_in_u, k, rev_in_mask);
    const auto out_mask = kmer::blind_rev_in_mask(char_out_u, k, rev_in_mask);

    fwd_hash = kmer::prev_forward_hash(fwd_hash, char_out_u, in_mask);
    rev_hash = kmer::prev_reverse_hash(rev_hash, out_mask, char_in_u);
    utils::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    seq.pop_back();
    seq.push_front(char_in);
    --pos;
  }

  /**
   * Like NtHash::peek(), but as if roll(char char_in) was called.
   */
  void peek(char char_in)
  {
    const auto k = static_cast<K_TYPE>(seq.size());
    const auto char_out_u = static_cast<unsigned char>(seq.front());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto out_mask = kmer::blind_fwd_out_mask(char_out_u, k, fwd_out_mask);
    const auto in_mask = kmer::blind_rev_in_mask(char_in_u, k, rev_in_mask);

    const auto fwd = kmer::next_forward_hash(fwd_hash, out_mask, char_in_u);
    const auto rev = kmer::next_reverse_hash(rev_hash, char_out_u, in_mask);
    utils::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  }

  /**
   * Like peek(char char_in), but as if roll_back(char char_in) was called.
   */
  void peek_back(char char_in)
  {
    const auto k = static_cast<K_TYPE>(seq.size());
    const auto char_out_u = static_cast<unsigned char>(seq.back());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto in_mask = kmer::blind_fwd_in_mask(char_in_u, k, rev_in_mask);
    const auto out_mask = kmer::blind_rev_in_mask(char_out_u, k, rev_in_mask);

    const auto fwd = kmer::prev_forward_hash(fwd_hash, char_out_u, in_mask);
    const auto rev = kmer::prev_reverse_hash(rev_hash, out_mask, char_in_u);
    utils::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  }

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const HASH_TYPE* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  ssize_t get_pos() const { return pos; }

  /**
   * Get the number of hashes generated per k-mer.
   * @return Number of hashes per k-mer
   */
  NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  K_TYPE get_k() const { return static_cast<K_TYPE>(seq.size()); }

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
  std::deque<char> seq;
  NUM_HASHES_TYPE num_hashes;
  ssize_t pos;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  HASH_TYPE fwd_out_mask[4];
  HASH_TYPE rev_in_mask[4];
  std::unique_ptr<HASH_TYPE[]> hash_arr;
};

} // namespace nthash