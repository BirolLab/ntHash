#include <iostream>
#include <nthash.hpp>
#include <string>

int
main()
{
  const std::string seq = "ATCGTACGATGCATGCATGCTGACG";
  const unsigned num_hashes = 3;
  const unsigned kmer_size = 6;
  nthash::NtHash nth(seq, num_hashes, kmer_size);
  while (nth.roll()) {
    std::cout << seq.substr(nth.get_pos(), nth.get_k()) << '\t';
    for (unsigned i = 0; i < nth.get_hash_num(); i++) {
      std::cout << std::hex << "0x" << nth.hashes()[i] << '\t';
    }
    std::cout << std::endl;
  }
  return 0;
}