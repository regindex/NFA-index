#ifndef SUCCINCT_SET_HPP
#define SUCCINCT_SET_HPP

#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <vector>

class SuccinctSet {
  private:
    // uint128_t
    // just does bittrick over this
    sdsl::bit_vector bitvec;
    sdsl::rank_support_v<1> rank_support;
    sdsl::select_support_mcl<1> select_support;

    size_t universe_size;
    size_t set_size;

  public:
    // Constructor
    SuccinctSet(size_t u, const std::vector<size_t>& elements);
    // Constructor to load from stream
    SuccinctSet(std::istream& in);
    // Rank: A.rank(x) = |{y ∈ A | y ≤ x}|
    size_t rank(size_t x) const;

    // Select: A.select(i) = x such that x ∈ A and A.rank(x) = i
    size_t select(size_t i) const;

    // Predecessor: Largest element ≤ x
    int pred(size_t x) const;

    // Strict-Successor: Smallest element > x
    int succ(size_t x) const;

    // Membership: Check if x ∈ A
    bool is_member(size_t x) const;

    void print() const;

    int size_bytes();

    void serialize(std::ostream& out) const;
};

#endif  // SUCCINCT_SET_HPP
