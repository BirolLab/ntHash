#include <nthash.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

std::vector<std::string>
get_data(unsigned n, unsigned l)
{
  const std::string chars = "ACGT";
  std::vector<std::string> data;
  while (n--) {
    std::string seq;
    seq.reserve(n);
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::uniform_int_distribution<size_t> dist(0, chars.size() - 1);
    for (size_t i = 0; i < l; i++) {
      seq.append(std::string(1, chars[dist(rng)]));
    }
    data.push_back(std::move(seq));
  }
  return data;
}

int
main()
{
  const auto data = get_data(1000000, 100);
  auto t1 = std::chrono::system_clock::now();
  uint64_t sum = 0;
  for (const auto& seq : data) {
    nthash::NtHash h(seq, 3, 64);
    while (h.roll()) {
      sum += h.hashes()[2];
    }
  }
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = (t2 - t1);
  std::cout << elapsed.count() << std::endl;
  std::cout << sum << std::endl;
  return 0;
}