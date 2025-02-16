#include "BitVector.hpp"

template <typename B>
Custom_BV<B>::Custom_BV(B b) {
    this->bit_vec = *new B(b);
    rank0 = typename B::rank_0_type(&bit_vec);
    rank1 = typename B::rank_1_type(&bit_vec);
    select0 = typename B::select_0_type(&bit_vec);
    select1 = typename B::select_1_type(&bit_vec);
}

template <typename B>
Custom_BV<B>::Custom_BV(std::ifstream& in) {
    bit_vec.load(in);
    rank0.load(in, &bit_vec);
    rank1.load(in, &bit_vec);
    select0.load(in, &bit_vec);
    select1.load(in, &bit_vec);
}

// access to bit at position i
template <typename B>
size_t Custom_BV<B>::access(size_t i) {
    return bit_vec[i];
}

// rank of c at position u
template <typename B>
size_t Custom_BV<B>::rank(int i, size_t c) {
    if (i == 0) return 0;
    if (c == 0)
        return this->rank0(i);
    else
        return this->rank1(i);
}

// get position of i-th c value
template <typename B>
size_t Custom_BV<B>::select(int i, size_t c) {
    if (i <= 0)
        return 0;
    else {
        if (i > this->rank(bit_vec.size(), c))
            return bit_vec.size();
        else if (c == 0)
            return select0(i);
        else
            return select1(i);
    }
}

template <typename B>
void Custom_BV<B>::print() {
    std::cout << bit_vec << std::endl;
}

template <typename B>
size_t Custom_BV<B>::size() {
    return bit_vec.size();
};

template <typename B>
void Custom_BV<B>::serialize(std::ofstream& out) {
    bit_vec.serialize(out);
    rank0.serialize(out);
    rank1.serialize(out);
    select0.serialize(out);
    select1.serialize(out);
}

template <typename B>
int Custom_BV<B>::size_bytes() {
    return sdsl::size_in_bytes(bit_vec) + sdsl::size_in_bytes(rank0) + sdsl::size_in_bytes(rank1) +
           sdsl::size_in_bytes(select0) + sdsl::size_in_bytes(select1);
}

template class Custom_BV<sdsl::bit_vector>;
template class Custom_BV<sdsl::bit_vector_il<64>>;
template class Custom_BV<sdsl::bit_vector_il<128>>;
template class Custom_BV<sdsl::bit_vector_il<256>>;
template class Custom_BV<sdsl::bit_vector_il<512>>;
template class Custom_BV<sdsl::bit_vector_il<1024>>;
template class Custom_BV<sdsl::bit_vector_il<2048>>;
template class Custom_BV<sdsl::rrr_vector<127>>;
template class Custom_BV<sdsl::rrr_vector<63>>;
template class Custom_BV<sdsl::rrr_vector<40>>;
template class Custom_BV<sdsl::rrr_vector<15>>;
template class Custom_BV<sdsl::rrr_vector<8>>;
template class Custom_BV<sdsl::sd_vector<sdsl::bit_vector>>;