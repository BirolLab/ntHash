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

class BlindSeedNtHash
{
public:
  /**
   * Construct a BlindSeedNtHash object for spaced seeds on-the-fly.
   */
  BlindSeedNtHash(const char* seq_ptr,
                  const std::vector<std::string>& seeds,
                  unsigned num_hashes_per_seed,
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
    if (k == 0) {
      throw std::invalid_argument("BlindSeedNtHash: k must be greater than 0");
    }
    if (seeds.empty()) {
      throw std::invalid_argument("BlindSeedNtHash: empty seeds");
    }
    size_t pos_n = 0;
    if (seed::is_invalid_kmer(seq_ptr + pos, k, pos_n)) {
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
                  unsigned num_hashes_per_seed,
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
    if (k == 0) {
      throw std::invalid_argument("BlindSeedNtHash: k must be greater than 0");
    }
    if (seeds.empty()) {
      throw std::invalid_argument("BlindSeedNtHash: empty seeds");
    }
    size_t pos_n = 0;
    if (seed::is_invalid_kmer(seq_ptr + pos, k, pos_n)) {
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

  BlindSeedNtHash(BlindSeedNtHash&&) noexcept = default;

  void roll(char char_in)
  {
    seed::ntmsm64_forward_core(
      [this, char_in](unsigned idx) {
        if (idx == k) {
          return static_cast<unsigned char>(char_in);
        }
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
        if (idx == 0) {
          return static_cast<unsigned char>(char_in);
        }
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
        if (idx == k) {
          return static_cast<unsigned char>(char_in);
        }
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
        if (idx == 0) {
          return static_cast<unsigned char>(char_in);
        }
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
  unsigned get_hash_num_per_seed() const { return num_hashes_per_seed; }
  K_TYPE get_k() const { return static_cast<K_TYPE>(seq.size()); }
  const HASH_TYPE* get_forward_hash() const { return fwd_hash.get(); }
  const HASH_TYPE* get_reverse_hash() const { return rev_hash.get(); }

private:
  std::deque<char> seq;
  unsigned num_hashes_per_seed;
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
        fwd_shift_table[d][i] = internal::roll_next(seed_f, d);
        rev_shift_table[d][i] = internal::roll_next(seed_r, d);
      }
    }
  }
};

} // namespace nthash::seed