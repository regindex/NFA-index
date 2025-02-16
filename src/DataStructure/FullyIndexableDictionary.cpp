#include "FullyIndexableDictionary.hpp"

// Constructor
SuccinctSet::SuccinctSet(size_t u, const std::vector<size_t>& elements) {
    this->universe_size = u;
    bitvec = sdsl::bit_vector(u + 1, 0);  // Initialize bit vector with all 0s

    for (auto x : elements) {
        bitvec[x] = 1;  // Mark positions of elements in the set
    }
    // Initialize rank and select structures
    rank_support = sdsl::rank_support_v<1>(&bitvec);
    select_support = sdsl::select_support_mcl<1>(&bitvec);
    this->set_size = rank_support.rank(u + 1);
}

SuccinctSet::SuccinctSet(std::istream& in) {
    bitvec.load(in);
    rank_support.load(in, &bitvec);
    select_support.load(in, &bitvec);

    in.read(reinterpret_cast<char*>(&universe_size), sizeof(universe_size));
    in.read(reinterpret_cast<char*>(&set_size), sizeof(set_size));
}

// Rank: A.rank(x) = |{y ∈ A | y ≤ x}|
size_t SuccinctSet::rank(size_t x) const {
    if (x > universe_size) return set_size;  // Out of bounds
    return rank_support(x + 1);
}

// Select: A.select(i) = x such that x ∈ A and A.rank(x) = i
size_t SuccinctSet::select(size_t i) const {
    if (i < 1 || i > set_size) throw std::out_of_range("Invalid rank index");
    return select_support(i);
}

// Predecessor: Largest element ≤ x
int SuccinctSet::pred(size_t x) const {
    size_t r = rank(x);
    return r > 0 ? select(r) : (-1);  // -1 indicates no predecessor
}

// Strict-Successor: Smallest element > x
int SuccinctSet::succ(size_t x) const {
    size_t r = rank(x);
    return r < set_size ? select(r + 1) : (-1);  // -1 indicates no successor
}

// Membership: Check if x ∈ A
bool SuccinctSet::is_member(size_t x) const {
    if (x > universe_size) return false;  // Out of bounds
    return bitvec[x] == 1;
}

void SuccinctSet::print() const {
    for (size_t i = 0; i <= universe_size; ++i) {
        if (bitvec[i]) std::cout << i << " ";
    }
    std::cout << std::endl;
}

int SuccinctSet::size_bytes() {
    return sdsl::size_in_bytes(bitvec) + sdsl::size_in_bytes(rank_support) + sdsl::size_in_bytes(select_support) +
           sizeof(universe_size) + sizeof(set_size);
}

// Serialize the SuccinctSet to a stream
void SuccinctSet::serialize(std::ostream& out) const {
    bitvec.serialize(out);
    rank_support.serialize(out);
    select_support.serialize(out);

    out.write(reinterpret_cast<const char*>(&universe_size), sizeof(universe_size));
    out.write(reinterpret_cast<const char*>(&set_size), sizeof(set_size));
}