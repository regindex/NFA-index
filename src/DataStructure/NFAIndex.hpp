#ifndef NFA_INDEX_HPP_2
#define NFA_INDEX_HPP

#include <sdsl/bit_vectors.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "BitVector.hpp"
#include "FullyIndexableDictionary.hpp"
#include "WTLargeAlphabet.hpp"

struct coord {
    uint64_t l;
    uint64_t t;
};
template <typename B>
class NFA_Index {
  private:
    uint64_t n;      // number of nodes
    uint64_t m;      // number of edges
    uint64_t sigma;  // cardinality of the alphabet
    uint64_t p;      // width of the partial order
    bool debug;      // debug
    uint64_t k;

    std::unordered_map<char, uint64_t> alphabet_match;
    // std::vector<char> alfabeto;

    // alphabet chains mapping
    std::vector<std::unique_ptr<SuccinctSet>> ALPHABET_DICT;
    // chain for every node
    std::unique_ptr<Custom_BV<B>> CHAIN;
    // out degree for every node
    std::unique_ptr<Custom_BV<B>> OUT_DEG;
    // out edge using format (chain, lavbel)
    std::unique_ptr<WTLargeAlphabet> OUT;  // use perfect hashing function
    // in degree for every node
    std::unique_ptr<Custom_BV<B>> IN_DEG;
    // list of incoming label for edges ordered by target node
    std::unique_ptr<Custom_BV<B>> IN;
    // final states
    std::unique_ptr<Custom_BV<B>> FINAL;

    std::unordered_map<uint64_t, uint64_t> matching;

    uint64_t get_chain_cardinality(uint64_t i);
    uint64_t convert_p_sigma_pair(uint64_t chain, uint64_t label);
    std::pair<uint64_t, uint64_t> inverse_convert_p_sigma_pair(uint64_t encoded_value);

    uint64_t out_ijk(uint64_t i, uint64_t j, uint64_t k, uint64_t a);
    uint64_t in_leq_h(uint64_t i, uint64_t a, uint64_t h);
    uint64_t in_geq_z(uint64_t i, uint64_t a, uint64_t z);
    uint64_t max_lambda(uint64_t i, uint64_t a);

  public:
    NFA_Index(std::string input_file, size_t index, size_t build_type, bool debug);

    uint64_t count(std::string alpha);
    std::vector<std::pair<uint64_t, uint64_t>> locate(std::string alpha);
    bool membership(std::string alpha);

    void store(std::string input_file);
    void load(std::ifstream& input);

    void print_size(std::string);
    void print();
};

#endif  // NFA_INDEX_HPP
