#include <nthash.hpp>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

std::vector<std::string>
get_data(unsigned n, unsigned l)
{
  const char chars[] = "ACGT";
  std::vector<std::string> data;
  data.reserve(n);
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<size_t> dist(0, 3);
  while (n--) {
    std::string seq;
    seq.reserve(l);
    for (size_t i = 0; i < l; i++) {
      seq.push_back(chars[dist(rng)]);
    }
    data.push_back(std::move(seq));
  }
  return data;
}

void
print_usage(const char* prog_name)
{
  std::cerr << "Usage: " << prog_name << " [options]\n"
            << "Options:\n"
            << "  -n <int>    Number of sequences (default: 1000000)\n"
            << "  -l <int>    Sequence length (default: 100)\n"
            << "  -h <int>    Number of hashes (default: 3)\n"
            << "  -k <int>    K-mer size (default: 64)\n";
}

int
main(int argc, char** argv)
{
  // Default values
  unsigned n = 1000000;
  unsigned l = 100;
  unsigned num_hashes = 3;
  unsigned k = 64;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-n" && i + 1 < argc) {
      n = std::stoul(argv[++i]);
    } else if (arg == "-l" && i + 1 < argc) {
      l = std::stoul(argv[++i]);
    } else if (arg == "-h" && i + 1 < argc) {
      num_hashes = std::stoul(argv[++i]);
    } else if (arg == "-k" && i + 1 < argc) {
      k = std::stoul(argv[++i]);
    } else {
      print_usage(argv[0]);
      return 1;
    }
  }

  std::cout << "Generating " << n << " sequences of length " << l << "...\n";
  const auto data = get_data(n, l);

  std::cout << "Hashing (k=" << k << ", hashes=" << num_hashes << ")...\n";
  auto t1 = std::chrono::system_clock::now();
  uint64_t sum = 0;

  for (const auto& seq : data) {
    nthash::kmer::NtHash h(seq, num_hashes, k);
    while (h.roll()) {
      sum += h.hashes()[num_hashes - 1];
    }
  }

  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = (t2 - t1);

  std::cout << "Time: " << elapsed.count() << " s\n";
  std::cout << "Checksum: " << sum << std::endl;

  return 0;
}