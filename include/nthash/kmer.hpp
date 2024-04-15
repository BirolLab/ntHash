#pragma once
#include "nthash/internal.hpp"
#include "nthash/nthash.hpp"

namespace {

using nthash::CONVERT_TAB;
using nthash::CP_OFF;
using nthash::DIMER_TAB;
using nthash::RC_CONVERT_TAB;
using nthash::SEED_N;
using nthash::SEED_TAB;
using nthash::srol;
using nthash::srol_table;
using nthash::sror;
using nthash::TETRAMER_TAB;
using nthash::TRIMER_TAB;

/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
inline bool is_invalid_kmer(const char *seq, unsigned k, size_t &pos_n) {
  for (int i = (int)k - 1; i >= 0; i--) {
    if (SEED_TAB[(unsigned char)seq[i]] == SEED_N) {
      pos_n = i;
      return true;
    }
  }
  return false;
}

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
inline uint64_t base_forward_hash(const char *seq, unsigned k) {
  uint64_t h_val = 0;
  for (unsigned i = 0; i < k - 3; i += 4) {
    h_val = srol(h_val, 4);
    uint8_t loc = 0;
    loc += 64 * CONVERT_TAB[(unsigned char)seq[i]];     // NOLINT
    loc += 16 * CONVERT_TAB[(unsigned char)seq[i + 1]]; // NOLINT
    loc += 4 * CONVERT_TAB[(unsigned char)seq[i + 2]];
    loc += CONVERT_TAB[(unsigned char)seq[i + 3]];
    h_val ^= TETRAMER_TAB[loc];
  }
  const unsigned remainder = k % 4;
  h_val = srol(h_val, remainder);
  if (remainder == 3) {
    uint8_t trimer_loc = 0;
    trimer_loc += 16 * CONVERT_TAB[(unsigned char)seq[k - 3]]; // NOLINT
    trimer_loc += 4 * CONVERT_TAB[(unsigned char)seq[k - 2]];
    trimer_loc += CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 0;
    dimer_loc += 4 * CONVERT_TAB[(unsigned char)seq[k - 2]];
    dimer_loc += CONVERT_TAB[(unsigned char)seq[k - 1]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1]];
  }
  return h_val;
}

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
inline uint64_t next_forward_hash(uint64_t fh_val, unsigned k,
                                  unsigned char char_out,
                                  unsigned char char_in) {
  uint64_t h_val = srol(fh_val);
  h_val ^= SEED_TAB[char_in];
  h_val ^= srol_table(char_out, k);
  return h_val;
}

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Forward hash value rolled back
 */
inline uint64_t prev_forward_hash(uint64_t fh_val, unsigned k,
                                  unsigned char char_out,
                                  unsigned char char_in) {
  uint64_t h_val = fh_val ^ srol_table(char_in, k);
  h_val ^= SEED_TAB[char_out];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
inline uint64_t base_reverse_hash(const char *seq, unsigned k) {
  uint64_t h_val = 0;
  const unsigned remainder = k % 4;
  if (remainder == 3) {
    uint8_t trimer_loc = 0;
    trimer_loc += 16 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]]; // NOLINT
    trimer_loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[k - 2]];
    trimer_loc += RC_CONVERT_TAB[(unsigned char)seq[k - 3]];
    h_val ^= TRIMER_TAB[trimer_loc];
  } else if (remainder == 2) {
    uint8_t dimer_loc = 0;
    dimer_loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[k - 1]];
    dimer_loc += RC_CONVERT_TAB[(unsigned char)seq[k - 2]];
    h_val ^= DIMER_TAB[dimer_loc];
  } else if (remainder == 1) {
    h_val ^= SEED_TAB[(unsigned char)seq[k - 1] & CP_OFF];
  }
  for (int i = (int)(k - remainder) - 1; i >= 3; i -= 4) {
    h_val = srol(h_val, 4);
    uint8_t loc = 0;
    loc += 64 * RC_CONVERT_TAB[(unsigned char)seq[i]];     // NOLINT
    loc += 16 * RC_CONVERT_TAB[(unsigned char)seq[i - 1]]; // NOLINT
    loc += 4 * RC_CONVERT_TAB[(unsigned char)seq[i - 2]];
    loc += RC_CONVERT_TAB[(unsigned char)seq[i - 3]];
    h_val ^= TETRAMER_TAB[loc];
  }
  return h_val;
}

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled hash value for the reverse-complement
 */
inline uint64_t next_reverse_hash(uint64_t rh_val, unsigned k,
                                  unsigned char char_out,
                                  unsigned char char_in) {
  uint64_t h_val = rh_val ^ srol_table(char_in & CP_OFF, k);
  h_val ^= SEED_TAB[char_out & CP_OFF];
  h_val = sror(h_val);
  return h_val;
}

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
inline uint64_t prev_reverse_hash(uint64_t rh_val, unsigned k,
                                  unsigned char char_out,
                                  unsigned char char_in) {
  uint64_t h_val = srol(rh_val);
  h_val ^= SEED_TAB[char_in & CP_OFF];
  h_val ^= srol_table(char_out & CP_OFF, k);
  return h_val;
}

} // namespace
