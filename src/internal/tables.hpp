#pragma once

#include <array>
#include <cstdint>

namespace nthash::internal {

// 64-bit random seeds corresponding to bases and their complements
// from the generate_seeds script with rng=42
constexpr uint64_t SEED_A = 0x3eb13b9046685257;
constexpr uint64_t SEED_C = 0x22310aefe5d92bca;
constexpr uint64_t SEED_G = 0x83677b6b400f4886;
constexpr uint64_t SEED_T = 0x9fe74a14e3be311b;
constexpr uint64_t SEED_N = 0x0000000000000000;

// offset for the complement base in the random seeds table
constexpr uint8_t CP_OFF = 0x07;

constexpr int ASCII_SIZE = 256;

[[nodiscard]] constexpr std::array<uint64_t, ASCII_SIZE>
generate_seed_table() noexcept
{
  std::array<uint64_t, ASCII_SIZE> tab{};
  tab['A'] = tab['a'] = SEED_A;
  tab['C'] = tab['c'] = SEED_C;
  tab['G'] = tab['g'] = SEED_G;
  tab['T'] = tab['t'] = tab['U'] = tab['u'] = SEED_T;
  tab['A' & CP_OFF] = SEED_T;
  tab['C' & CP_OFF] = SEED_G;
  tab['T' & CP_OFF] = SEED_A;
  tab['U' & CP_OFF] = SEED_A;
  tab['G' & CP_OFF] = SEED_C;
  return tab;
}

[[nodiscard]] constexpr std::array<uint8_t, ASCII_SIZE>
generate_convert_table() noexcept
{
  std::array<uint8_t, ASCII_SIZE> tab{};
  for (auto& val : tab) {
    val = ASCII_SIZE - 1;
  }
  tab['A'] = tab['a'] = 0;
  tab['C'] = tab['c'] = 1;
  tab['G'] = tab['g'] = 2;
  tab['T'] = tab['t'] = tab['U'] = tab['u'] = 3;
  return tab;
}

[[nodiscard]] constexpr std::array<uint8_t, ASCII_SIZE>
generate_rc_convert_table() noexcept
{
  std::array<uint8_t, ASCII_SIZE> tab{};
  for (auto& val : tab) {
    val = ASCII_SIZE - 1;
  }
  tab['A'] = tab['a'] = 3;
  tab['C'] = tab['c'] = 2;
  tab['G'] = tab['g'] = 1;
  tab['T'] = tab['t'] = tab['U'] = tab['u'] = 0;
  return tab;
}

alignas(64) inline constexpr auto SEED_TAB = generate_seed_table();
alignas(64) inline constexpr auto CONVERT_TAB = generate_convert_table();
alignas(64) inline constexpr auto RC_CONVERT_TAB = generate_rc_convert_table();

} // namespace nthash::internal