#ifndef WT_LARGE_ALPHABET
#define WT_LARGE_ALPHABET

#include <ctime>
#include <fstream>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

class WTLargeAlphabet {
  private:
    // Wt with alpahebt partitioning
    sdsl::wt_int<sdsl::bit_vector_il<512>> wt;

    // wt size
    uint64_t n;

  public:
    WTLargeAlphabet(sdsl::int_vector<> keys);
    WTLargeAlphabet(std::ifstream& in);
    ~WTLargeAlphabet();
    uint64_t access(size_t i);
    uint64_t select(size_t i, uint64_t c);
    uint64_t rank(size_t i, uint64_t c);
    void print();
    uint64_t size();
    int size_bytes();

    void serialize(std::ofstream& out);
};

#endif  // FULLY_INDEXABLE_DICTIONARY_HPP
