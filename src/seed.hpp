#pragma once

#include "tables.hpp"
#include "utils.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nthash {

using utils::HASH_TYPE;
using utils::K_TYPE;
using utils::NUM_HASHES_TYPE;
using SpacedSeedBlocks = std::vector<std::array<unsigned, 2>>;
using SpacedSeedMonomers = std::vector<unsigned>;

namespace seed {

constexpr HASH_TYPE BASE_SEEDS[4] = { tables::SEED_A,
                                      tables::SEED_C,
                                      tables::SEED_G,
                                      tables::SEED_T };

inline void
get_blocks(const std::vector<std::string>& seed_strings,
           std::vector<SpacedSeedBlocks>& blocks,
           std::vector<SpacedSeedMonomers>& monomers)
{
  blocks.clear();
  monomers.clear();
  blocks.reserve(seed_strings.size());
  monomers.reserve(seed_strings.size());

  for (const auto& seed_string : seed_strings) {
    const char pad = seed_string.back() == '1' ? '0' : '1';
    const std::string padded_string = seed_string + pad;
    SpacedSeedBlocks care_blocks;
    SpacedSeedBlocks ignore_blocks;
    std::vector<unsigned> care_monos;
    std::vector<unsigned> ignore_monos;

    unsigned i_start = 0;
    bool is_care_block = padded_string[0] == '1';
    for (unsigned pos = 0; pos < padded_string.length(); pos++) {
      if (is_care_block && padded_string[pos] == '0') {
        if (pos - i_start == 1) {
          care_monos.push_back(i_start);
        } else {
          care_blocks.push_back({ i_start, pos });
        }
        i_start = pos;
        is_care_block = false;
      } else if (!is_care_block && padded_string[pos] == '1') {
        if (pos - i_start == 1) {
          ignore_monos.push_back(i_start);
        } else {
          ignore_blocks.push_back({ i_start, pos });
        }
        i_start = pos;
        is_care_block = true;
      }
    }

    const unsigned num_cares = care_blocks.size() * 2 + care_monos.size();
    const unsigned num_ignores =
      ignore_blocks.size() * 2 + ignore_monos.size() + 2;

    if (num_ignores < num_cares) {
      ignore_blocks.push_back(
        { 0, static_cast<unsigned>(seed_string.length()) });
      blocks.push_back(ignore_blocks);
      monomers.push_back(ignore_monos);
    } else {
      blocks.push_back(care_blocks);
      monomers.push_back(care_monos);
    }
  }
}

inline void
parsed_seeds_to_blocks(const std::vector<std::vector<unsigned>>& seeds,
                       K_TYPE k,
                       std::vector<SpacedSeedBlocks>& blocks,
                       std::vector<SpacedSeedMonomers>& monomers)
{
  std::vector<std::string> seed_strings;
  seed_strings.reserve(seeds.size());

  for (const auto& seed : seeds) {
    std::string seed_string(k, '1');
    for (const auto i : seed) {
      if (i >= k) {
        throw std::invalid_argument(
          "SeedNtHash: parsed seed index out of range");
      }
      seed_string[i] = '0';
    }
    seed_strings.push_back(std::move(seed_string));
  }

  get_blocks(seed_strings, blocks, monomers);
}

inline void
check_seeds(const std::vector<std::string>& seeds, K_TYPE k)
{
  for (const auto& seed : seeds) {
    if (seed.length() != k) {
      throw std::invalid_argument("SeedNtHash: spaced seed string length (" +
                                  std::to_string(seed.length()) +
                                  ") not equal to k=" + std::to_string(k) +
                                  " in " + seed);
    }
  }
}

inline bool
ntmsm64_init(const char* kmer_seq,
             const std::vector<SpacedSeedBlocks>& seeds_blocks,
             const std::vector<SpacedSeedMonomers>& seeds_monomers,
             K_TYPE k,
             unsigned m,
             NUM_HASHES_TYPE m2,
             const std::vector<std::array<HASH_TYPE, 4>>& fwd_shift_table,
             const std::vector<std::array<HASH_TYPE, 4>>& rev_shift_table,
             HASH_TYPE* fh_nomonos,
             HASH_TYPE* rh_nomonos,
             HASH_TYPE* fh_val,
             HASH_TYPE* rh_val,
             HASH_TYPE* h_val)
{
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    HASH_TYPE fh_seed = 0;
    HASH_TYPE rh_seed = 0;

    for (const auto& block : seeds_blocks[i_seed]) {
      for (unsigned pos = block[0]; pos < block[1]; pos++) {
        const auto c = static_cast<unsigned char>(kmer_seq[pos]);
        const auto idx = tables::CONVERT_TAB[c];
        fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
        rh_seed ^= rev_shift_table[pos][idx];
      }
    }

    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = static_cast<unsigned char>(kmer_seq[pos]);
      const auto idx = tables::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = utils::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t =
        h_val[base] * (i_hash ^ static_cast<HASH_TYPE>(k) * utils::MULTISEED);
      t ^= t >> utils::MULTISHIFT;
      h_val[base + i_hash] = t;
    }
  }

  return true;
}

template<typename CharAt>
inline void
ntmsm64_forward_core(
  CharAt char_at,
  const std::vector<SpacedSeedBlocks>& seeds_blocks,
  const std::vector<SpacedSeedMonomers>& seeds_monomers,
  K_TYPE k,
  unsigned m,
  NUM_HASHES_TYPE m2,
  const std::vector<std::array<HASH_TYPE, 4>>& fwd_shift_table,
  const std::vector<std::array<HASH_TYPE, 4>>& rev_shift_table,
  HASH_TYPE* fh_nomonos,
  HASH_TYPE* rh_nomonos,
  HASH_TYPE* fh_val,
  HASH_TYPE* rh_val,
  HASH_TYPE* h_val)
{
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    HASH_TYPE fh_seed = utils::roll_next(fh_nomonos[i_seed]);
    HASH_TYPE rh_seed = rh_nomonos[i_seed];

    for (const auto& block : seeds_blocks[i_seed]) {
      const auto i_in = block[1];
      const auto i_out = block[0];
      const auto char_in = char_at(i_in);
      const auto char_out = char_at(i_out);

      const auto idx_in = tables::CONVERT_TAB[char_in];
      const auto idx_out = tables::CONVERT_TAB[char_out];

      fh_seed ^= fwd_shift_table[k - i_out][idx_out];
      fh_seed ^= fwd_shift_table[k - i_in][idx_in];
      rh_seed ^= rev_shift_table[i_out][idx_out];
      rh_seed ^= rev_shift_table[i_in][idx_in];
    }

    rh_seed = utils::roll_back(rh_seed);
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = char_at(pos + 1);
      const auto idx = tables::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = utils::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t =
        h_val[base] * (i_hash ^ static_cast<HASH_TYPE>(k) * utils::MULTISEED);
      t ^= t >> utils::MULTISHIFT;
      h_val[base + i_hash] = t;
    }
  }
}

template<typename CharAt>
inline void
ntmsm64_backward_core(
  CharAt char_at,
  const std::vector<SpacedSeedBlocks>& seeds_blocks,
  const std::vector<SpacedSeedMonomers>& seeds_monomers,
  K_TYPE k,
  unsigned m,
  NUM_HASHES_TYPE m2,
  const std::vector<std::array<HASH_TYPE, 4>>& fwd_shift_table,
  const std::vector<std::array<HASH_TYPE, 4>>& rev_shift_table,
  HASH_TYPE* fh_nomonos,
  HASH_TYPE* rh_nomonos,
  HASH_TYPE* fh_val,
  HASH_TYPE* rh_val,
  HASH_TYPE* h_val)
{
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    HASH_TYPE fh_seed = fh_nomonos[i_seed];
    HASH_TYPE rh_seed = utils::roll_next(rh_nomonos[i_seed]);

    for (const auto& block : seeds_blocks[i_seed]) {
      const auto i_in = block[0];
      const auto i_out = block[1];
      const auto char_in = char_at(i_in);
      const auto char_out = char_at(i_out);

      const auto idx_in = tables::CONVERT_TAB[char_in];
      const auto idx_out = tables::CONVERT_TAB[char_out];

      fh_seed ^= fwd_shift_table[k - i_out][idx_out];
      fh_seed ^= fwd_shift_table[k - i_in][idx_in];
      rh_seed ^= rev_shift_table[i_out][idx_out];
      rh_seed ^= rev_shift_table[i_in][idx_in];
    }

    fh_seed = utils::roll_back(fh_seed);
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = char_at(pos); // Fixed backward monomer index offset
      const auto idx = tables::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = utils::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t =
        h_val[base] * (i_hash ^ static_cast<HASH_TYPE>(k) * utils::MULTISEED);
      t ^= t >> utils::MULTISHIFT;
      h_val[base + i_hash] = t;
    }
  }
}

inline std::vector<std::vector<unsigned>>
parse_seeds(const std::vector<std::string>& seed_strings)
{
  std::vector<std::vector<unsigned>> seed_set;
  seed_set.reserve(seed_strings.size());

  for (const auto& seed_string : seed_strings) {
    std::vector<unsigned> seed;
    for (size_t pos = 0; pos < seed_string.size(); ++pos) {
      if (seed_string[pos] != '1') {
        seed.push_back(static_cast<unsigned>(pos));
      }
    }
    seed_set.push_back(std::move(seed));
  }

  return seed_set;
}

} // namespace seed

class SeedNtHash
{
public:
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::string>& seeds,
             NUM_HASHES_TYPE num_hashes_per_seed,
             K_TYPE k,
             size_t pos = 0)
    : seq(seq, seq_len)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , initialized(false)
    , fwd_hash_nomonos(new HASH_TYPE[seeds.size()])
    , rev_hash_nomonos(new HASH_TYPE[seeds.size()])
    , fwd_hash(new HASH_TYPE[seeds.size()])
    , rev_hash(new HASH_TYPE[seeds.size()])
    , hash_arr(new HASH_TYPE[num_hashes_per_seed * seeds.size()])
  {
    if (k == 0)
      throw std::invalid_argument("SeedNtHash: k must be greater than 0");
    if (this->seq.size() < k)
      throw std::invalid_argument("SeedNtHash: sequence smaller than k");
    if (pos > this->seq.size() - k)
      throw std::invalid_argument("SeedNtHash: invalid pos");
    if (seeds.empty())
      throw std::invalid_argument("SeedNtHash: empty seeds");

    seed::check_seeds(seeds, k);
    seed::get_blocks(seeds, blocks, monomers);
    init_shift_tables();
  }

  SeedNtHash(std::string_view seq,
             const std::vector<std::string>& seeds,
             NUM_HASHES_TYPE num_hashes_per_seed,
             K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::vector<unsigned>>& seeds,
             NUM_HASHES_TYPE num_hashes_per_seed,
             K_TYPE k,
             size_t pos = 0)
    : seq(seq, seq_len)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , initialized(false)
    , fwd_hash_nomonos(new HASH_TYPE[seeds.size()])
    , rev_hash_nomonos(new HASH_TYPE[seeds.size()])
    , fwd_hash(new HASH_TYPE[seeds.size()])
    , rev_hash(new HASH_TYPE[seeds.size()])
    , hash_arr(new HASH_TYPE[num_hashes_per_seed * seeds.size()])
  {
    if (k == 0)
      throw std::invalid_argument("SeedNtHash: k must be greater than 0");
    if (this->seq.size() < k)
      throw std::invalid_argument("SeedNtHash: sequence smaller than k");
    if (pos > this->seq.size() - k)
      throw std::invalid_argument("SeedNtHash: invalid pos");
    if (seeds.empty())
      throw std::invalid_argument("SeedNtHash: empty seeds");

    seed::parsed_seeds_to_blocks(seeds, k, blocks, monomers);
    init_shift_tables();
  }

  SeedNtHash(std::string_view seq,
             const std::vector<std::vector<unsigned>>& seeds,
             NUM_HASHES_TYPE num_hashes_per_seed,
             K_TYPE k,
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
    , fwd_shift_table(obj.fwd_shift_table)
    , rev_shift_table(obj.rev_shift_table)
    , fwd_hash_nomonos(new HASH_TYPE[obj.blocks.size()])
    , rev_hash_nomonos(new HASH_TYPE[obj.blocks.size()])
    , fwd_hash(new HASH_TYPE[obj.blocks.size()])
    , rev_hash(new HASH_TYPE[obj.blocks.size()])
    , hash_arr(new HASH_TYPE[obj.num_hashes_per_seed * obj.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                obj.fwd_hash_nomonos.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos.get(),
                obj.rev_hash_nomonos.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(fwd_hash.get(),
                obj.fwd_hash.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(rev_hash.get(),
                obj.rev_hash.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(hash_arr.get(),
                obj.hash_arr.get(),
                obj.num_hashes_per_seed * obj.blocks.size() *
                  sizeof(HASH_TYPE));
  }

  SeedNtHash(SeedNtHash&&) = default;

  bool roll()
  {
    if (!initialized)
      return init();
    if (pos >= seq.size() - k)
      return false;

    size_t pos_n = 0;
    if (utils::is_invalid_kmer(seq.data() + pos + 1, k, pos_n)) {
      pos += pos_n + 1;
      return init();
    }

    seed::ntmsm64_forward_core(
      [this](unsigned idx) {
        return static_cast<unsigned char>(seq[pos + idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(blocks.size()),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos.get(),
      rev_hash_nomonos.get(),
      fwd_hash.get(),
      rev_hash.get(),
      hash_arr.get());
    ++pos;
    return true;
  }

  bool roll_back()
  {
    if (!initialized)
      return init();
    if (pos == 0)
      return false;

    size_t pos_n = 0;
    if (utils::is_invalid_kmer(seq.data() + pos - 1, k, pos_n)) {
      return false; // Safe exit instead of infinite loop
    }

    seed::ntmsm64_backward_core(
      [this](unsigned idx) {
        return static_cast<unsigned char>(seq[pos - 1 + idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(blocks.size()),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos.get(),
      rev_hash_nomonos.get(),
      fwd_hash.get(),
      rev_hash.get(),
      hash_arr.get());
    --pos;
    return true;
  }

  bool peek()
  {
    if (pos >= seq.size() - k)
      return false;
    return peek(seq[pos + k]);
  }

  bool peek(char char_in)
  {
    if (!initialized)
      return init();

    const auto n = blocks.size();
    auto fwd_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto fwd_hash_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_cpy = std::make_unique<HASH_TYPE[]>(n);

    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));

    seed::ntmsm64_forward_core(
      [this, char_in](unsigned idx) {
        if (idx == k)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[pos + idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(n),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos_cpy.get(),
      rev_hash_nomonos_cpy.get(),
      fwd_hash_cpy.get(),
      rev_hash_cpy.get(),
      hash_arr.get());
    return true;
  }

  bool peek_back()
  {
    if (pos == 0)
      return false;
    return peek_back(seq[pos - 1]);
  }

  bool peek_back(char char_in)
  {
    if (!initialized)
      return init();

    const auto n = blocks.size();
    auto fwd_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto fwd_hash_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_cpy = std::make_unique<HASH_TYPE[]>(n);

    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));

    seed::ntmsm64_backward_core(
      [this, char_in](unsigned idx) {
        if (idx == 0)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[pos - 1 + idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(n),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos_cpy.get(),
      rev_hash_nomonos_cpy.get(),
      fwd_hash_cpy.get(),
      rev_hash_cpy.get(),
      hash_arr.get());
    return true;
  }

  const HASH_TYPE* hashes() const { return hash_arr.get(); }
  size_t get_pos() const { return pos; }
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }
  NUM_HASHES_TYPE get_hash_num_per_seed() const { return num_hashes_per_seed; }
  K_TYPE get_k() const { return k; }
  const HASH_TYPE* get_forward_hash() const { return fwd_hash.get(); }
  const HASH_TYPE* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::string_view seq;
  NUM_HASHES_TYPE num_hashes_per_seed;
  K_TYPE k;
  size_t pos;
  bool initialized;
  std::vector<SpacedSeedBlocks> blocks;
  std::vector<SpacedSeedMonomers> monomers;

  // O(1) Precomputed shift tables: [rotation_distance][base_index (0-3)]
  std::vector<std::array<HASH_TYPE, 4>> fwd_shift_table;
  std::vector<std::array<HASH_TYPE, 4>> rev_shift_table;

  std::unique_ptr<HASH_TYPE[]> fwd_hash_nomonos;
  std::unique_ptr<HASH_TYPE[]> rev_hash_nomonos;
  std::unique_ptr<HASH_TYPE[]> fwd_hash;
  std::unique_ptr<HASH_TYPE[]> rev_hash;
  std::unique_ptr<HASH_TYPE[]> hash_arr;

  void init_shift_tables()
  {
    fwd_shift_table.resize(k + 1);
    rev_shift_table.resize(k + 1);
    for (size_t d = 0; d <= k; ++d) {
      for (unsigned i = 0; i < 4; ++i) {
        HASH_TYPE seed_f = seed::BASE_SEEDS[i];
        HASH_TYPE seed_r = seed::BASE_SEEDS[3 - i]; // RC mapping
        fwd_shift_table[d][i] = utils::roll_next(seed_f, d);
        rev_shift_table[d][i] = utils::roll_next(seed_r, d);
      }
    }
  }

  bool init()
  {
    size_t pos_n = 0;
    while (pos <= seq.size() - k) {
      if (utils::is_invalid_kmer(seq.data() + pos, k, pos_n)) {
        pos += pos_n + 1;
        continue;
      }
      seed::ntmsm64_init(seq.data() + pos,
                         blocks,
                         monomers,
                         k,
                         static_cast<unsigned>(blocks.size()),
                         num_hashes_per_seed,
                         fwd_shift_table,
                         rev_shift_table,
                         fwd_hash_nomonos.get(),
                         rev_hash_nomonos.get(),
                         fwd_hash.get(),
                         rev_hash.get(),
                         hash_arr.get());
      initialized = true;
      return true;
    }
    return false;
  }
};

class BlindSeedNtHash
{
public:
  /**
   * Construct a BlindSeedNtHash object for spaced seeds on-the-fly.
   */
  BlindSeedNtHash(const char* seq_ptr,
                  const std::vector<std::string>& seeds,
                  NUM_HASHES_TYPE num_hashes_per_seed,
                  K_TYPE k,
                  ssize_t pos = 0)
    : seq(seq_ptr + pos, seq_ptr + pos + k)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , fwd_hash_nomonos(new HASH_TYPE[seeds.size()])
    , rev_hash_nomonos(new HASH_TYPE[seeds.size()])
    , fwd_hash(new HASH_TYPE[seeds.size()])
    , rev_hash(new HASH_TYPE[seeds.size()])
    , hash_arr(new HASH_TYPE[num_hashes_per_seed * seeds.size()])
  {
    if (k == 0)
      throw std::invalid_argument("BlindSeedNtHash: k must be greater than 0");
    if (seeds.empty())
      throw std::invalid_argument("BlindSeedNtHash: empty seeds");

    size_t pos_n = 0;
    if (utils::is_invalid_kmer(seq_ptr + pos, k, pos_n)) {
      throw std::invalid_argument(
        "BlindSeedNtHash: initial k-mer contains invalid characters");
    }

    seed::check_seeds(seeds, k);
    seed::get_blocks(seeds, blocks, monomers);
    init_shift_tables();

    seed::ntmsm64_init(seq_ptr + pos,
                       blocks,
                       monomers,
                       k,
                       static_cast<unsigned>(blocks.size()),
                       num_hashes_per_seed,
                       fwd_shift_table,
                       rev_shift_table,
                       fwd_hash_nomonos.get(),
                       rev_hash_nomonos.get(),
                       fwd_hash.get(),
                       rev_hash.get(),
                       hash_arr.get());
  }

  BlindSeedNtHash(const char* seq_ptr,
                  const std::vector<std::vector<unsigned>>& seeds,
                  NUM_HASHES_TYPE num_hashes_per_seed,
                  K_TYPE k,
                  ssize_t pos = 0)
    : seq(seq_ptr + pos, seq_ptr + pos + k)
    , num_hashes_per_seed(num_hashes_per_seed)
    , k(k)
    , pos(pos)
    , fwd_hash_nomonos(new HASH_TYPE[seeds.size()])
    , rev_hash_nomonos(new HASH_TYPE[seeds.size()])
    , fwd_hash(new HASH_TYPE[seeds.size()])
    , rev_hash(new HASH_TYPE[seeds.size()])
    , hash_arr(new HASH_TYPE[num_hashes_per_seed * seeds.size()])
  {
    if (k == 0)
      throw std::invalid_argument("BlindSeedNtHash: k must be greater than 0");
    if (seeds.empty())
      throw std::invalid_argument("BlindSeedNtHash: empty seeds");

    size_t pos_n = 0;
    if (utils::is_invalid_kmer(seq_ptr + pos, k, pos_n)) {
      throw std::invalid_argument(
        "BlindSeedNtHash: initial k-mer contains invalid characters");
    }

    seed::parsed_seeds_to_blocks(seeds, k, blocks, monomers);
    init_shift_tables();

    seed::ntmsm64_init(seq_ptr + pos,
                       blocks,
                       monomers,
                       k,
                       static_cast<unsigned>(blocks.size()),
                       num_hashes_per_seed,
                       fwd_shift_table,
                       rev_shift_table,
                       fwd_hash_nomonos.get(),
                       rev_hash_nomonos.get(),
                       fwd_hash.get(),
                       rev_hash.get(),
                       hash_arr.get());
  }

  BlindSeedNtHash(const BlindSeedNtHash& obj)
    : seq(obj.seq)
    , num_hashes_per_seed(obj.num_hashes_per_seed)
    , k(obj.k)
    , pos(obj.pos)
    , blocks(obj.blocks)
    , monomers(obj.monomers)
    , fwd_shift_table(obj.fwd_shift_table)
    , rev_shift_table(obj.rev_shift_table)
    , fwd_hash_nomonos(new HASH_TYPE[obj.blocks.size()])
    , rev_hash_nomonos(new HASH_TYPE[obj.blocks.size()])
    , fwd_hash(new HASH_TYPE[obj.blocks.size()])
    , rev_hash(new HASH_TYPE[obj.blocks.size()])
    , hash_arr(new HASH_TYPE[obj.num_hashes_per_seed * obj.blocks.size()])
  {
    std::memcpy(fwd_hash_nomonos.get(),
                obj.fwd_hash_nomonos.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos.get(),
                obj.rev_hash_nomonos.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(fwd_hash.get(),
                obj.fwd_hash.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(rev_hash.get(),
                obj.rev_hash.get(),
                obj.blocks.size() * sizeof(HASH_TYPE));
    std::memcpy(hash_arr.get(),
                obj.hash_arr.get(),
                obj.num_hashes_per_seed * obj.blocks.size() *
                  sizeof(HASH_TYPE));
  }

  BlindSeedNtHash(BlindSeedNtHash&&) = default;

  void roll(char char_in)
  {
    seed::ntmsm64_forward_core(
      [this, char_in](unsigned idx) {
        if (idx == k)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(blocks.size()),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos.get(),
      rev_hash_nomonos.get(),
      fwd_hash.get(),
      rev_hash.get(),
      hash_arr.get());

    seq.pop_front();
    seq.push_back(char_in);
    ++pos;
  }

  void roll_back(char char_in)
  {
    seed::ntmsm64_backward_core(
      [this, char_in](unsigned idx) {
        if (idx == 0)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[idx - 1]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(blocks.size()),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos.get(),
      rev_hash_nomonos.get(),
      fwd_hash.get(),
      rev_hash.get(),
      hash_arr.get());

    seq.pop_back();
    seq.push_front(char_in);
    --pos;
  }

  void peek(char char_in)
  {
    const auto n = blocks.size();
    auto fwd_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto fwd_hash_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_cpy = std::make_unique<HASH_TYPE[]>(n);

    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));

    seed::ntmsm64_forward_core(
      [this, char_in](unsigned idx) {
        if (idx == k)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[idx]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(n),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos_cpy.get(),
      rev_hash_nomonos_cpy.get(),
      fwd_hash_cpy.get(),
      rev_hash_cpy.get(),
      hash_arr.get());
  }

  void peek_back(char char_in)
  {
    const auto n = blocks.size();
    auto fwd_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_nomonos_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto fwd_hash_cpy = std::make_unique<HASH_TYPE[]>(n);
    auto rev_hash_cpy = std::make_unique<HASH_TYPE[]>(n);

    std::memcpy(fwd_hash_nomonos_cpy.get(),
                fwd_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));
    std::memcpy(rev_hash_nomonos_cpy.get(),
                rev_hash_nomonos.get(),
                n * sizeof(HASH_TYPE));

    seed::ntmsm64_backward_core(
      [this, char_in](unsigned idx) {
        if (idx == 0)
          return static_cast<unsigned char>(char_in);
        return static_cast<unsigned char>(seq[idx - 1]);
      },
      blocks,
      monomers,
      k,
      static_cast<unsigned>(n),
      num_hashes_per_seed,
      fwd_shift_table,
      rev_shift_table,
      fwd_hash_nomonos_cpy.get(),
      rev_hash_nomonos_cpy.get(),
      fwd_hash_cpy.get(),
      rev_hash_cpy.get(),
      hash_arr.get());
  }

  const HASH_TYPE* hashes() const { return hash_arr.get(); }
  ssize_t get_pos() const { return pos; }
  unsigned get_hash_num() const { return num_hashes_per_seed * blocks.size(); }
  NUM_HASHES_TYPE get_hash_num_per_seed() const { return num_hashes_per_seed; }
  K_TYPE get_k() const { return static_cast<K_TYPE>(seq.size()); }
  const HASH_TYPE* get_forward_hash() const { return fwd_hash.get(); }
  const HASH_TYPE* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::deque<char> seq;
  NUM_HASHES_TYPE num_hashes_per_seed;
  K_TYPE k;
  ssize_t pos;
  std::vector<SpacedSeedBlocks> blocks;
  std::vector<SpacedSeedMonomers> monomers;

  std::vector<std::array<HASH_TYPE, 4>> fwd_shift_table;
  std::vector<std::array<HASH_TYPE, 4>> rev_shift_table;

  std::unique_ptr<HASH_TYPE[]> fwd_hash_nomonos;
  std::unique_ptr<HASH_TYPE[]> rev_hash_nomonos;
  std::unique_ptr<HASH_TYPE[]> fwd_hash;
  std::unique_ptr<HASH_TYPE[]> rev_hash;
  std::unique_ptr<HASH_TYPE[]> hash_arr;

  void init_shift_tables()
  {
    fwd_shift_table.resize(k + 1);
    rev_shift_table.resize(k + 1);
    for (size_t d = 0; d <= k; ++d) {
      for (unsigned i = 0; i < 4; ++i) {
        HASH_TYPE seed_f = seed::BASE_SEEDS[i];
        HASH_TYPE seed_r = seed::BASE_SEEDS[3 - i];
        fwd_shift_table[d][i] = utils::roll_next(seed_f, d);
        rev_shift_table[d][i] = utils::roll_next(seed_r, d);
      }
    }
  }
};

} // namespace nthash