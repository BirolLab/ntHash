#pragma once

#include "utils.hpp"

#include <cstring>
#include <memory>
#include <stdexcept>

namespace nthash::kmer {

using internal::HASH_TYPE;
using internal::K_TYPE;
using internal::NUM_HASHES_TYPE;

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
    : kmer(seq + pos, seq + pos + k)
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
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
  }

  BlindNtHash(const BlindNtHash& obj)
    : kmer(obj.kmer)
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
    const auto k = static_cast<K_TYPE>(kmer.size());
    const auto char_out_u = static_cast<unsigned char>(kmer.front());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto out_mask = kmer::blind_fwd_out_mask(char_out_u, k, fwd_out_mask);
    const auto in_mask = kmer::blind_rev_in_mask(char_in_u, k, rev_in_mask);

    fwd_hash = kmer::next_forward_hash(fwd_hash, out_mask, char_in_u);
    rev_hash = kmer::next_reverse_hash(rev_hash, char_out_u, in_mask);
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    kmer.erase(0, 1);
    kmer.push_back(char_in);
    ++pos;
  }

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in)
  {
    const auto k = static_cast<K_TYPE>(kmer.size());
    const auto char_out_u = static_cast<unsigned char>(kmer.back());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto in_mask = kmer::blind_fwd_in_mask(char_in_u, k, rev_in_mask);
    const auto out_mask = kmer::blind_rev_in_mask(char_out_u, k, rev_in_mask);

    fwd_hash = kmer::prev_forward_hash(fwd_hash, char_out_u, in_mask);
    rev_hash = kmer::prev_reverse_hash(rev_hash, out_mask, char_in_u);
    internal::extend_hashes(fwd_hash, rev_hash, k, num_hashes, hash_arr.get());
    kmer.pop_back();
    kmer.insert(kmer.begin(), char_in);
    --pos;
  }

  /**
   * Like NtHash::peek(), but as if roll(char char_in) was called.
   */
  void peek(char char_in)
  {
    const auto k = static_cast<K_TYPE>(kmer.size());
    const auto char_out_u = static_cast<unsigned char>(kmer.front());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto out_mask = kmer::blind_fwd_out_mask(char_out_u, k, fwd_out_mask);
    const auto in_mask = kmer::blind_rev_in_mask(char_in_u, k, rev_in_mask);

    const auto fwd = kmer::next_forward_hash(fwd_hash, out_mask, char_in_u);
    const auto rev = kmer::next_reverse_hash(rev_hash, char_out_u, in_mask);
    internal::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
  }

  /**
   * Like peek(char char_in), but as if roll_back(char char_in) was called.
   */
  void peek_back(char char_in)
  {
    const auto k = static_cast<K_TYPE>(kmer.size());
    const auto char_out_u = static_cast<unsigned char>(kmer.back());
    const auto char_in_u = static_cast<unsigned char>(char_in);
    const auto in_mask = kmer::blind_fwd_in_mask(char_in_u, k, rev_in_mask);
    const auto out_mask = kmer::blind_rev_in_mask(char_out_u, k, rev_in_mask);

    const auto fwd = kmer::prev_forward_hash(fwd_hash, char_out_u, in_mask);
    const auto rev = kmer::prev_reverse_hash(rev_hash, out_mask, char_in_u);
    internal::extend_hashes(fwd, rev, k, num_hashes, hash_arr.get());
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
  K_TYPE get_k() const { return static_cast<K_TYPE>(kmer.size()); }

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
  std::string kmer;
  NUM_HASHES_TYPE num_hashes;
  ssize_t pos;
  HASH_TYPE fwd_hash = 0;
  HASH_TYPE rev_hash = 0;
  HASH_TYPE fwd_out_mask[4];
  HASH_TYPE rev_in_mask[4];
  std::unique_ptr<HASH_TYPE[]> hash_arr;
};

} // namespace nthash::kmer