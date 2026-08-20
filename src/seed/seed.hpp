#pragma once

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

namespace nthash::seed {

class SeedNtHash
{
public:
  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::string>& seeds,
             unsigned num_hashes_per_seed,
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
             unsigned num_hashes_per_seed,
             K_TYPE k,
             size_t pos = 0)
    : SeedNtHash(seq.data(), seq.size(), seeds, num_hashes_per_seed, k, pos)
  {
  }

  SeedNtHash(const char* seq,
             size_t seq_len,
             const std::vector<std::vector<unsigned>>& seeds,
             unsigned num_hashes_per_seed,
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
             unsigned num_hashes_per_seed,
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
    if (seed::is_invalid_kmer(seq.data() + pos + 1, k, pos_n)) {
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
    if (seed::is_invalid_kmer(seq.data() + pos - 1, k, pos_n)) {
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
  unsigned get_hash_num_per_seed() const { return num_hashes_per_seed; }
  K_TYPE get_k() const { return k; }
  const HASH_TYPE* get_forward_hash() const { return fwd_hash.get(); }
  const HASH_TYPE* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::string_view seq;
  unsigned num_hashes_per_seed;
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
        fwd_shift_table[d][i] = internal::roll_next(seed_f, d);
        rev_shift_table[d][i] = internal::roll_next(seed_r, d);
      }
    }
  }

  bool init()
  {
    size_t pos_n = 0;
    while (pos <= seq.size() - k) {
      if (seed::is_invalid_kmer(seq.data() + pos, k, pos_n)) {
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

} // namespace nthash::seed