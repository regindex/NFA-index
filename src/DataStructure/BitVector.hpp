#ifndef INDEX_BITVECTOR_HPP
#define INDEX_BITVECTOR_HPP

#include <fstream>
#include <sdsl/bit_vectors.hpp>
template <typename B>
class Custom_BV {
  private:
    B bit_vec;
    typename B::rank_0_type rank0;
    typename B::rank_1_type rank1;
    typename B::select_0_type select0;
    typename B::select_1_type select1;

  public:
    Custom_BV(B b);
    Custom_BV(std::ifstream& in);
    size_t access(size_t i);
    size_t select(int i, size_t c);
    size_t rank(int i, size_t c);
    void print();
    size_t size();
    int size_bytes();
    void serialize(std::ofstream& out);
};

#endif  // INDEX_BITVECTOR_HPP
