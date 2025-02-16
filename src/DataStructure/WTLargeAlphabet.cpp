#include "WTLargeAlphabet.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

WTLargeAlphabet::WTLargeAlphabet(sdsl::int_vector<> keys) { sdsl::construct_im(wt, keys); }
WTLargeAlphabet::WTLargeAlphabet(std::ifstream& in) { wt.load(in); }
WTLargeAlphabet::~WTLargeAlphabet() {}

uint64_t WTLargeAlphabet::access(size_t i) { return this->wt[i]; }
uint64_t WTLargeAlphabet::select(size_t i, uint64_t c) { return this->wt.select(i, c); }
uint64_t WTLargeAlphabet::rank(size_t i, uint64_t c) { return this->wt.rank(i, c); }
void WTLargeAlphabet::print() { std::cout << this->wt << std::endl; }
uint64_t WTLargeAlphabet::size() { return this->wt.size(); }

int WTLargeAlphabet::size_bytes() { return sdsl::size_in_bytes(this->wt); }

void WTLargeAlphabet::serialize(std::ofstream& out) { this->wt.serialize(out); }
