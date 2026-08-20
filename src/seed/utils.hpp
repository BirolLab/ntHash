#pragma once

#include "../internal/internal.hpp"

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

using internal::HASH_TYPE;
using internal::K_TYPE;

using SpacedSeedBlocks = std::vector<std::array<unsigned, 2>>;
using SpacedSeedMonomers = std::vector<unsigned>;

constexpr HASH_TYPE BASE_SEEDS[4] = { internal::SEED_A,
                                      internal::SEED_C,
                                      internal::SEED_G,
                                      internal::SEED_T };

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
             unsigned m2,
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
        const auto idx = internal::CONVERT_TAB[c];
        fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
        rh_seed ^= rev_shift_table[pos][idx];
      }
    }

    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = static_cast<unsigned char>(kmer_seq[pos]);
      const auto idx = internal::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = internal::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t = h_val[base] *
               (i_hash ^ static_cast<HASH_TYPE>(k) * internal::MULTISEED);
      t ^= t >> internal::MULTISHIFT;
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
  unsigned m2,
  const std::vector<std::array<HASH_TYPE, 4>>& fwd_shift_table,
  const std::vector<std::array<HASH_TYPE, 4>>& rev_shift_table,
  HASH_TYPE* fh_nomonos,
  HASH_TYPE* rh_nomonos,
  HASH_TYPE* fh_val,
  HASH_TYPE* rh_val,
  HASH_TYPE* h_val)
{
  for (unsigned i_seed = 0; i_seed < m; i_seed++) {
    HASH_TYPE fh_seed = internal::roll_next(fh_nomonos[i_seed]);
    HASH_TYPE rh_seed = rh_nomonos[i_seed];

    for (const auto& block : seeds_blocks[i_seed]) {
      const auto i_in = block[1];
      const auto i_out = block[0];
      const auto char_in = char_at(i_in);
      const auto char_out = char_at(i_out);

      const auto idx_in = internal::CONVERT_TAB[char_in];
      const auto idx_out = internal::CONVERT_TAB[char_out];

      fh_seed ^= fwd_shift_table[k - i_out][idx_out];
      fh_seed ^= fwd_shift_table[k - i_in][idx_in];
      rh_seed ^= rev_shift_table[i_out][idx_out];
      rh_seed ^= rev_shift_table[i_in][idx_in];
    }

    rh_seed = internal::roll_back(rh_seed);
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = char_at(pos + 1);
      const auto idx = internal::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = internal::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t = h_val[base] *
               (i_hash ^ static_cast<HASH_TYPE>(k) * internal::MULTISEED);
      t ^= t >> internal::MULTISHIFT;
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
  unsigned m2,
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
    HASH_TYPE rh_seed = internal::roll_next(rh_nomonos[i_seed]);

    for (const auto& block : seeds_blocks[i_seed]) {
      const auto i_in = block[0];
      const auto i_out = block[1];
      const auto char_in = char_at(i_in);
      const auto char_out = char_at(i_out);

      const auto idx_in = internal::CONVERT_TAB[char_in];
      const auto idx_out = internal::CONVERT_TAB[char_out];

      fh_seed ^= fwd_shift_table[k - i_out][idx_out];
      fh_seed ^= fwd_shift_table[k - i_in][idx_in];
      rh_seed ^= rev_shift_table[i_out][idx_out];
      rh_seed ^= rev_shift_table[i_in][idx_in];
    }

    fh_seed = internal::roll_back(fh_seed);
    fh_nomonos[i_seed] = fh_seed;
    rh_nomonos[i_seed] = rh_seed;

    for (const auto pos : seeds_monomers[i_seed]) {
      const auto c = char_at(pos); // Fixed backward monomer index offset
      const auto idx = internal::CONVERT_TAB[c];
      fh_seed ^= fwd_shift_table[k - 1 - pos][idx];
      rh_seed ^= rev_shift_table[pos][idx];
    }

    fh_val[i_seed] = fh_seed;
    rh_val[i_seed] = rh_seed;

    const auto base = i_seed * m2;
    h_val[base] = internal::canonical(fh_seed, rh_seed);
    for (unsigned i_hash = 1; i_hash < m2; i_hash++) {
      auto t = h_val[base] *
               (i_hash ^ static_cast<HASH_TYPE>(k) * internal::MULTISEED);
      t ^= t >> internal::MULTISHIFT;
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

/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
[[nodiscard]] inline bool
is_invalid_kmer(const char* seq, unsigned k, size_t& pos_n)
{
  for (size_t i = k; i-- > 0;) {
    if (internal::SEED_TAB[(unsigned char)seq[i]] == internal::SEED_N) {
      pos_n = i;
      return true;
    }
  }
  return false;
}

} // namespace nthash::seed
