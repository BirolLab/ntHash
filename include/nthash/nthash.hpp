#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace nthash {

/**
 * String representing the hash function's name. Only change if hash outputs
 * are different from the previous version. Useful for tracking differences in
 * saved hashes, e.g., in Bloom filters.
 */
static const char* const NTHASH_FN_NAME = "ntHash_v2.5";

/**
 * This lets us minimize NtHash object size. Good for performance if it's
 * copied in, e.g., DBG traversal.
 */
namespace typedefs {
using NUM_HASHES_TYPE = uint8_t;
using K_TYPE = uint16_t;
using SpacedSeedBlocks = std::vector<std::array<unsigned, 2>>;
using SpacedSeedMonomers = std::vector<unsigned>;
} // namespace typedefs

/**
 * Normal k-mer hashing.
 */
class NtHash;

/**
 * Similar to the NtHash class, but instead of rolling on a predefined sequence,
 * BlindNtHash needs to be fed the new character on each roll. This is useful
 * when traversing an implicit de Bruijn Graph, as we need to query all bases
 * to know the possible extensions.
 */
class BlindNtHash;

/**
 * Spaced seed hashing.
 */
class SeedNtHash;

/**
 * Similar to the SeedNtHash class, but instead of rolling on a predefined
 * sequence, BlindSeedNtHash needs to be fed the new character on each roll.
 */
class BlindSeedNtHash;

/**
 * Parse each spaced seed pattern into lists of don't care positions. Legacy
 * function used in btllib Bloom filters.
 */
std::vector<std::vector<unsigned>>
parse_seeds(const std::vector<std::string>& seed_strings);

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
         typedefs::NUM_HASHES_TYPE num_hashes,
         typedefs::K_TYPE k,
         size_t pos = 0);

  /**
   * Construct an ntHash object for k-mers.
   * @param seq Sequence string
   * @param num_hashes Number of hashes to produce per k-mer
   * @param k K-mer size
   * @param pos Position in sequence to start hashing from
   */
  NtHash(const std::string& seq,
         typedefs::NUM_HASHES_TYPE num_hashes,
         typedefs::K_TYPE k,
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
    , hash_arr(new uint64_t[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(uint64_t));
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
  bool roll();

  /**
   * Like the roll() function, but advance backwards.
   * @return \p true on success and \p false otherwise
   */
  bool roll_back();

  /**
   * Peeks the hash values as if roll() was called (without advancing the
   * NtHash object. The peeked hash values can be obtained through the
   * hashes() method.
   * @return \p true on success and \p false otherwise
   */
  bool peek();

  /**
   * Like peek(), but as if roll_back() was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back();

  /**
   * Peeks the hash values as if roll() was called for char_in (without
   * advancing the NtHash object. The peeked hash values can be obtained through
   * the hashes() method.
   * @return \p true on success and \p false otherwise
   */
  bool peek(char char_in);

  /**
   * Like peek(), but as if roll_back on char_in was called.
   * @return \p true on success and \p false otherwise
   */
  bool peek_back(char char_in);

  /**
   * Get the array of current canonical hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

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
  typedefs::NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  typedefs::K_TYPE get_k() const { return k; }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  uint64_t get_forward_hash() const { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  uint64_t get_reverse_hash() const { return rev_hash; }

private:
  std::string_view seq;
  typedefs::NUM_HASHES_TYPE num_hashes;
  typedefs::K_TYPE k;
  size_t pos;
  bool initialized;
  uint64_t fwd_hash = 0;
  uint64_t rev_hash = 0;
  std::unique_ptr<uint64_t[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init();
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
              typedefs::NUM_HASHES_TYPE num_hashes,
              typedefs::K_TYPE k,
              ssize_t pos = 0);

  BlindNtHash(const BlindNtHash& obj)
    : seq(obj.seq)
    , num_hashes(obj.num_hashes)
    , pos(obj.pos)
    , fwd_hash(obj.fwd_hash)
    , rev_hash(obj.rev_hash)
    , hash_arr(new uint64_t[obj.num_hashes])
  {
    std::memcpy(
      hash_arr.get(), obj.hash_arr.get(), num_hashes * sizeof(uint64_t));
  }

  BlindNtHash(BlindNtHash&&) = default;

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindNtHash object was constructed on, the provided character
   * \p char_in is used as the next base. Useful if you want to query for
   * possible paths in an implicit de Bruijn graph graph.
   */
  void roll(char char_in);

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in);

  /**
   * Like NtHash::peek(), but as if roll(char char_in) was called.
   */
  void peek(char char_in);

  /**
   * Like peek(char char_in), but as if roll_back(char char_in) was called.
   */
  void peek_back(char char_in);

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

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
  typedefs::NUM_HASHES_TYPE get_hash_num() const { return num_hashes; }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  typedefs::K_TYPE get_k() const { return seq.size(); }

  /**
   * Get the hash value of the forward strand.
   * @return Forward hash value
   */
  uint64_t get_forward_hash() const { return fwd_hash; }

  /**
   * Get the hash value of the reverse strand.
   * @return Reverse-complement hash value
   */
  uint64_t get_reverse_hash() const { return rev_hash; }

private:
  std::deque<char> seq;
  typedefs::NUM_HASHES_TYPE num_hashes;
  ssize_t pos;
  uint64_t fwd_hash = 0;
  uint64_t rev_hash = 0;
  std::unique_ptr<uint64_t[]> hash_arr;
};

class SeedNtHash
{

public:
  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq C-string of the sequence to be hashed
   * @param seq_len Length of the sequence
   * @param seeds Vector of spaced seed patterns as strings (1s as cares, 0s
   * as don't cares, must be of size \p k)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::string>& seeds,
             typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
             typedefs::K_TYPE k,
             size_t pos = 0);

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq String of the sequence to be hashed
   * @param seeds Vector of spaced seed patterns as strings (1s as cares, 0s
   * as don't cares, must be of size \p k)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const std::string& seq,
             const std::vector<std::string>& seeds,
             typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
             typedefs::K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq C-string of the sequence to be hashed
   * @param seq_len Length of the sequence
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::vector<unsigned>>& seeds,
             typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
             typedefs::K_TYPE k,
             size_t pos = 0);

  /**
   * Construct an ntHash object for spaced seeds.
   * @param seq String of the sequence to be hashed
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  SeedNtHash(const std::string& seq,
             const std::vector<std::vector<unsigned>>& seeds,
             typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
             typedefs::K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  SeedNtHash(const SeedNtHash& obj)
    : seq(obj.seq)
    , num_hashes_per_seed(obj.num_hashes_per_seed)
    , k(obj.k)
    , pos(obj.pos)
    , initialized(obj.initialized)
    , blocks(obj.blocks)
    , monomers(obj.monomers)
    , fwd_hash_nomonos(new uint64_t[obj.blocks.size()])
    , rev_hash_nomonos(new uint64_t[obj.blocks.size()])
    , fwd_hash(new uint64_t[obj.blocks.size()])
    , rev_hash(new uint64_t[obj.blocks.size()])
    , hash_arr(new uint64_t[obj.num_hashes_per_seed * obj.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                obj.fwd_hash_nomonos.get(),
                obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos.get(),
                obj.rev_hash_nomonos.get(),
                obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(
      fwd_hash.get(), obj.fwd_hash.get(), obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(
      rev_hash.get(), obj.rev_hash.get(), obj.blocks.size() * sizeof(uint64_t));
    std::memcpy(hash_arr.get(),
                obj.hash_arr.get(),
                obj.num_hashes_per_seed * obj.blocks.size() * sizeof(uint64_t));
  }

  SeedNtHash(SeedNtHash&&) = default;

  /**
   * Calculate the next hash value. Refer to \ref NtHash::roll() for more
   * information.
   * @return \p true on success and \p false otherwise.
   */
  bool roll();

  /**
   * Like the roll() function, but advance backwards.
   * @return \p true on success and \p false otherwise.
   */
  bool roll_back();

  /**
   * Peeks the hash values as if roll() was called. Refer to
   * \ref NtHash::peek() for more information.
   * @return \p true on success and \p false otherwise.
   */
  bool peek();

  /**
   * Like peek(), but as if roll_back() was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek_back();

  /**
   * Like peek(), but as if roll(char char_in) was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek(char char_in);

  /**
   * Like peek(), but as if roll_back(char char_in) was called.
   * @return \p true on success and \p false otherwise.
   */
  bool peek_back(char char_in);

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  size_t get_pos() const { return pos; }

  /**
   * Get the length of the hash array.
   * @return Number of seeds times \p get_hash_num_per_seed()
   */
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }

  /**
   * Get the number of hashes generated per seed.
   * @return Number of hashes per seed
   */
  typedefs::NUM_HASHES_TYPE get_hash_num_per_seed() const
  {
    return num_hashes_per_seed;
  }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  typedefs::K_TYPE get_k() const { return k; }

  /**
   * Get the hash values of the forward strand.
   * @return Array of forward hash value arrays for each seed
   */
  uint64_t* get_forward_hash() const { return fwd_hash.get(); }

  /**
   * Get the hash values of the reverse strand.
   * @return Array of reverse-complement hash value arrays for each seed
   */
  uint64_t* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::string_view seq;
  typedefs::NUM_HASHES_TYPE num_hashes_per_seed;
  typedefs::K_TYPE k;
  size_t pos;
  bool initialized;
  std::vector<typedefs::SpacedSeedBlocks> blocks;
  std::vector<typedefs::SpacedSeedMonomers> monomers;
  std::unique_ptr<uint64_t[]> fwd_hash_nomonos;
  std::unique_ptr<uint64_t[]> rev_hash_nomonos;
  std::unique_ptr<uint64_t[]> fwd_hash;
  std::unique_ptr<uint64_t[]> rev_hash;
  std::unique_ptr<uint64_t[]> hash_arr;

  /**
   * Initialize the internal state of the iterator
   * @return \p true if successful, \p false otherwise
   */
  bool init();
};

class BlindSeedNtHash
{

public:
  /**
   * Construct an ntHash object for hashing spaced seeds on-the-fly.
   * @param seq C-string of the data. Only the first \p k characters will be
   * used, starting from \p pos.
   * @param seeds Vector of parsed spaced seed patterns (vectors of don't care
   * positions)
   * @param num_hashes_per_seed Number of hashes to generate per seed
   * @param k K-mer size
   * @param pos Position in seq to start hashing from
   */
  BlindSeedNtHash(const char* seq,
                  const std::vector<std::string>& seeds,
                  typedefs::NUM_HASHES_TYPE num_hashes_per_seed,
                  typedefs::K_TYPE k,
                  ssize_t pos = 0);

  BlindSeedNtHash(const BlindSeedNtHash& seed_nthash)
    : seq(seed_nthash.seq)
    , num_hashes_per_seed(seed_nthash.num_hashes_per_seed)
    , k(seed_nthash.k)
    , pos(seed_nthash.pos)
    , blocks(seed_nthash.blocks)
    , monomers(seed_nthash.monomers)
    , fwd_hash_nomonos(new uint64_t[seed_nthash.blocks.size()])
    , rev_hash_nomonos(new uint64_t[seed_nthash.blocks.size()])
    , fwd_hash(new uint64_t[seed_nthash.blocks.size()])
    , rev_hash(new uint64_t[seed_nthash.blocks.size()])
    , hash_arr(new uint64_t[num_hashes_per_seed * seed_nthash.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                seed_nthash.fwd_hash_nomonos.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash_nomonos.get(),
                seed_nthash.rev_hash_nomonos.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(fwd_hash.get(),
                seed_nthash.fwd_hash.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(rev_hash.get(),
                seed_nthash.rev_hash.get(),
                seed_nthash.blocks.size() * sizeof(uint64_t));
    std::memcpy(hash_arr.get(),
                seed_nthash.hash_arr.get(),
                num_hashes_per_seed * seed_nthash.blocks.size() *
                  sizeof(uint64_t));
  }

  BlindSeedNtHash(BlindSeedNtHash&&) = default;

  /**
   * Like the NtHash::roll() function, but instead of advancing in the
   * sequence BlindSeedNtHash object was constructed on, the provided character
   * \p char_in is used as the next base.
   */
  void roll(char char_in);

  /**
   * Like the roll(char char_in) function, but advance backwards.
   */
  void roll_back(char char_in);

  /**
   * Get the array of current hash values (length = \p get_hash_num())
   * @return Pointer to the hash array
   */
  const uint64_t* hashes() const { return hash_arr.get(); }

  /**
   * Get the position of last hashed k-mer or the k-mer to be hashed if roll()
   * has never been called on this NtHash object.
   * @return Position of the most recently hashed k-mer's first base-pair
   */
  ssize_t get_pos() const { return pos; }

  /**
   * Get the length of the hash array.
   * @return Number of seeds times \p get_hash_num_per_seed()
   */
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }

  /**
   * Get the number of hashes generated per seed.
   * @return Number of hashes per seed
   */
  typedefs::NUM_HASHES_TYPE get_hash_num_per_seed() const
  {
    return num_hashes_per_seed;
  }

  /**
   * Get the length of the k-mers.
   * @return \p k
   */
  typedefs::K_TYPE get_k() const { return k; }

  /**
   * Get the hash values of the forward strand.
   * @return Array of forward hash value arrays for each seed
   */
  uint64_t* get_forward_hash() const { return fwd_hash.get(); }

  /**
   * Get the hash values of the reverse strand.
   * @return Array of reverse-complement hash value arrays for each seed
   */
  uint64_t* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::deque<char> seq;
  typedefs::NUM_HASHES_TYPE num_hashes_per_seed;
  typedefs::K_TYPE k;
  ssize_t pos;
  std::vector<typedefs::SpacedSeedBlocks> blocks;
  std::vector<typedefs::SpacedSeedMonomers> monomers;
  std::unique_ptr<uint64_t[]> fwd_hash_nomonos;
  std::unique_ptr<uint64_t[]> rev_hash_nomonos;
  std::unique_ptr<uint64_t[]> fwd_hash;
  std::unique_ptr<uint64_t[]> rev_hash;
  std::unique_ptr<uint64_t[]> hash_arr;
};

} // namespace nthash
