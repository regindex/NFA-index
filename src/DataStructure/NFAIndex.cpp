#include "NFAIndex.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>  // For standard exceptions
#include <unordered_set>

bool write_map_to_ostream_u(const std::unordered_map<uint64_t, uint64_t>& my_map, std::ostream& os) {
    // Write the number of elements
    uint64_t map_size = my_map.size();
    os.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    if (!os.good()) {
        return false;
    }

    // Write each key-value pair
    for (const auto& pair : my_map) {
        os.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
        os.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second));
        if (!os.good()) {
            return false;
        }
    }

    return true;
}
bool read_map_from_istream_u(std::unordered_map<uint64_t, uint64_t>& my_map, std::istream& is) {
    // Read the number of elements
    uint64_t map_size = 0;
    is.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    if (!is.good()) {
        return false;
    }

    // Read each key-value pair and insert it into the map
    for (uint64_t i = 0; i < map_size; ++i) {
        uint64_t key = 0;
        uint64_t value = 0;
        is.read(reinterpret_cast<char*>(&key), sizeof(key));
        is.read(reinterpret_cast<char*>(&value), sizeof(value));
        if (!is.good()) {
            return false;
        }
        my_map[key] = value;
    }

    return true;
}

// Writes the map to an std::ostream (can be a file or a stringstream)
bool write_map_to_ostream(const std::unordered_map<char, uint64_t>& my_map, std::ostream& os) {
    if (!os) {
        std::cerr << "Error: Invalid ostream." << std::endl;
        return false;
    }

    uint64_t map_size = my_map.size();
    if (!os.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size))) {
        std::cerr << "Error writing map size to ostream." << std::endl;
        return false;
    }

    for (const auto& pair : my_map) {
        if (!os.write(&pair.first, sizeof(pair.first))) {
            std::cerr << "Error writing key to ostream." << std::endl;
            return false;
        }
        if (!os.write(reinterpret_cast<const char*>(&pair.second), sizeof(pair.second))) {
            std::cerr << "Error writing value to ostream." << std::endl;
            return false;
        }
    }

    return true;
}

// Reads the map from an std::istream
std::unordered_map<char, uint64_t> read_map_from_istream(std::istream& is) {
    std::unordered_map<char, uint64_t> my_map;

    if (!is) {
        std::cerr << "Error: Invalid istream." << std::endl;
        return my_map;
    }

    uint64_t map_size;
    if (!is.read(reinterpret_cast<char*>(&map_size), sizeof(map_size))) {
        std::cerr << "Error reading map size from istream. Stream might be empty or corrupted" << std::endl;
        return my_map;
    }

    my_map.reserve(map_size);

    for (uint64_t i = 0; i < map_size; ++i) {
        char key;
        uint64_t value;

        if (!is.read(&key, sizeof(key))) {
            std::cerr << "Error reading key from istream. Stream might be corrupted" << std::endl;
            return my_map;
        }
        if (!is.read(reinterpret_cast<char*>(&value), sizeof(value))) {
            std::cerr << "Error reading value from istream. Stream might be corrupted" << std::endl;
            return my_map;
        }
        my_map[key] = value;
    }

    return my_map;
}

// Funzioni per scrivere e leggere la struct coord
std::ostream& write_binary(std::ostream& os, const coord& c) {
    os.write(reinterpret_cast<const char*>(&c.l), sizeof(c.l));
    os.write(reinterpret_cast<const char*>(&c.t), sizeof(c.t));
    return os;
}

std::istream& read_binary(std::istream& is, coord& c) {
    is.read(reinterpret_cast<char*>(&c.l), sizeof(c.l));
    is.read(reinterpret_cast<char*>(&c.t), sizeof(c.t));
    return is;
}

// Funzioni per scrivere e leggere std::vector<coord>
std::ostream& write_binary(std::ostream& os, const std::vector<coord>& vec) {
    size_t size = vec.size();
    os.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto& c : vec) {
        write_binary(os, c);
    }
    return os;
}

std::istream& read_binary(std::istream& is, std::vector<coord>& vec) {
    size_t size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);  // Resize del vector prima della lettura

    for (size_t i = 0; i < size; ++i) {
        read_binary(is, vec[i]);
    }
    return is;
}

// Funzioni per scrivere e leggere la unordered_map
template <typename K>
std::ostream& write_binary(std::ostream& os, const std::unordered_map<K, std::vector<coord>>& map) {
    size_t size = map.size();
    os.write(reinterpret_cast<const char*>(&size), sizeof(size));

    for (const auto& pair : map) {
        size_t key_size = pair.first.size();
        os.write(reinterpret_cast<const char*>(&key_size), sizeof(key_size));
        os.write(pair.first.data(), key_size);

        write_binary(os, pair.second);  // Scrive il vector<coord>
    }
    return os;
}

template <typename K>
std::istream& read_binary(std::istream& is, std::unordered_map<K, std::vector<coord>>& map) {
    size_t size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));

    map.clear();
    map.reserve(size);

    for (size_t i = 0; i < size; ++i) {
        size_t key_size;
        is.read(reinterpret_cast<char*>(&key_size), sizeof(key_size));

        std::vector<char> key_buffer(key_size);
        is.read(key_buffer.data(), key_size);
        std::string key(key_buffer.begin(), key_buffer.end());

        std::vector<coord> vec;
        read_binary(is, vec);  // Legge il vector<coord>

        map[key] = vec;
    }
    return is;
}

struct VertexProperties {
    bool final;  // Store final_state as string to parse manually
    std::string positions;
    std::string chains;
};

struct EdgeProperties {
    std::string label;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProperties, EdgeProperties> Graph;

// Function to get the value from a comma-separated string at the specified index// Helper function to extract the value
// at a specified index from a comma-separated string and convert it to an integer
std::string get_value_at_index(const std::string& input, size_t index) {
    // Create a stringstream from the input string
    std::stringstream ss(input);
    std::string value;
    std::vector<std::string> values;

    // Split the input string by commas
    while (std::getline(ss, value, ',')) {
        // Trim any leading/trailing spaces from the value
        value.erase(0, value.find_first_not_of(' '));
        value.erase(value.find_last_not_of(' ') + 1);
        values.push_back(value);
    }

    // Check if the index is within bounds
    if (index > 0 && index <= values.size()) {
        return values[index - 1];  // return the value at the given index (1-based)
    } else {
        return "Invalid index";  // Handle out-of-bounds error
    }
}

// label in [0,\sigma-1], p in
template <typename B>
uint64_t NFA_Index<B>::convert_p_sigma_pair(uint64_t chain, uint64_t label) {
    uint64_t result = 0;

    // Shift the chain value left by sigma bits
    result = chain << sigma;

    // Combine the shifted chain with the label using bitwise OR
    result |= label;

    // return matching[result]
    return result;
}

// load index
template <typename B>
NFA_Index<B>::NFA_Index(std::string input_file, size_t index, size_t build_type, bool debug) {
    // Load NFA from file and construct index
    this->debug = debug;
    this->k = 0;
    if (build_type == 1) {
        std::ifstream input_stream(input_file, std::ios::binary);
        this->load(input_stream);
    } else if (build_type == 0) {
        try {

            std::ifstream input_stream(input_file);
            if (!input_stream) {
                std::cerr << "Error opening file." << std::endl;
            }
            // create temporary graph from file
            Graph g;
            boost::dynamic_properties dp;

            dp.property("label", get(&EdgeProperties::label, g));
            dp.property("final", get(&VertexProperties::final, g));  // Map final as a string

            dp.property("positions", get(&VertexProperties::positions, g));
            dp.property("chains", get(&VertexProperties::chains, g));
            try {
                boost::read_graphml(input_stream, g, dp);
            } catch (const std::exception& e) {
                std::cerr << "Error reading GraphML file: " << e.what() << std::endl;
            }

            // Fill the vertices_list with all the vertices
            std::vector<Graph::vertex_descriptor> vertices_list(boost::num_vertices(g));
            size_t i = 0;
            for (auto v : boost::make_iterator_range(boost::vertices(g))) {
                vertices_list[i++] = v;
            }

            // First, update the "chains" and "positions" properties by storing only the value at index `i` for each
            // vertex
            for (auto v : vertices_list) {
                // Update the "chains" and "positions" fields to store the value at the specified index
                g[v].chains = get_value_at_index(g[v].chains, index);        // Store the integer value as string
                g[v].positions = get_value_at_index(g[v].positions, index);  // Store the integer value as string
            }
            // After updating the properties, sort the vertices by the selected "chains" and "positions" values at the
            // specified index
            std::sort(vertices_list.begin(), vertices_list.end(),
                      [&g](Graph::vertex_descriptor v1, Graph::vertex_descriptor v2) {
                          // Extract the updated chain and position values as integers
                          int chain_value1 = std::stoi(g[v1].chains);
                          int chain_value2 = std::stoi(g[v2].chains);

                          // If the chain values are the same, compare by the positions value
                          if (chain_value1 == chain_value2) {
                              int pos_value1 = std::stoi(g[v1].positions);
                              int pos_value2 = std::stoi(g[v2].positions);
                              return pos_value1 < pos_value2;  // Compare by position if chains are the same
                          }

                          // Otherwise, compare by "chains"
                          return chain_value1 < chain_value2;
                      });
            // create index

            this->n = num_vertices(g);
            this->m = num_edges(g);
            sdsl::bit_vector tmp_out_deg(this->n + this->m, 0);
            sdsl::bit_vector tmp_in_deg(this->n + this->m, 0);
            sdsl::bit_vector tmp_final(this->n, 0);

            size_t c = 0, ot_d = 0, in_d = 0;
            std::string last_chain = "";
            sdsl::bit_vector tmp_chain(this->n);
            // iterate sorted vertices
            for (const auto& v : vertices_list) {
                // Update final
                if (g[v].final) tmp_final[c] = 1;
                // update chain
                if (g[v].chains != last_chain) {
                    tmp_chain[c] = 1;
                    last_chain = g[v].chains;
                }
                ++c;

                // out_degree
                size_t out_deg_v = out_degree(v, g);
                for (i = 0; i < out_deg_v; i++) ++ot_d;
                tmp_out_deg[ot_d++] = 1;

                // in_degree
                size_t in_deg_v = in_degree(v, g);
                for (i = 0; i < in_deg_v; i++) ++in_d;

                tmp_in_deg[in_d++] = 1;
            }
            this->p = std::stoi(last_chain);
            this->CHAIN = std::make_unique<Custom_BV<B>>(tmp_chain);
            this->IN_DEG = std::make_unique<Custom_BV<B>>(tmp_in_deg);
            this->OUT_DEG = std::make_unique<Custom_BV<B>>(tmp_out_deg);
            this->FINAL = std::make_unique<Custom_BV<B>>(tmp_final);

            // edges info for IN and OUT
            std::vector<std::tuple<int, int, char>> edge_list;

            std::set<char> alphabet;
            // Iterate over all edges in the graph and store (source, target, label_as_int)
            for (auto e : boost::make_iterator_range(boost::edges(g))) {
                // Get the source and target vertices of the edge
                auto source = boost::source(e, g);
                auto target = boost::target(e, g);

                // Get the edge label and convert it to integer
                char label = g[e].label[0];
                alphabet.insert(label);
                // Store the edge (source, target, label_as_int) in the edge list
                edge_list.push_back(std::make_tuple(source, target, label));
            }
            this->sigma = alphabet.size();
            std::vector<char> alphabet_vector(alphabet.begin(), alphabet.end());
            std::sort(alphabet_vector.begin(), alphabet_vector.end());

            for (size_t i = 0; i < this->sigma; ++i) {
                this->alphabet_match[alphabet_vector[i]] = i;
            }  // to store the matchings

            // compute OUT and dictionary
            // Sort by source ->label->target
            std::sort(edge_list.begin(), edge_list.end(), [](const auto& e1, const auto& e2) {
                if (std::get<0>(e1) == std::get<0>(e2)) {
                    if (std::get<2>(e1) == std::get<2>(e2)) {
                        return std::get<1>(e1) < std::get<1>(e2);  // Compare by target if source and label are equal
                    }
                    return std::get<2>(e1) < std::get<2>(e2);  // Compare by label if source is equal
                }
                return std::get<0>(e1) < std::get<0>(e2);  // Compare by source vertex
            });

            std::vector<std::vector<size_t>> tmp_dict(this->p);
            sdsl::int_vector<0> tmp_out(this->m, 0);
            i = 0;
            uint64_t value = 0;
            for (const auto& e : edge_list) {
                uint64_t entering_chain = this->CHAIN->rank(std::get<1>(e) + 1, 1);
                auto key = convert_p_sigma_pair(entering_chain, this->alphabet_match[std::get<2>(e)]);
                if (!(this->matching.find(key) != this->matching.end())) {
                    this->matching[key] = value;
                    ++value;
                }
            }

            for (const auto& e : edge_list) {
                uint64_t entering_chain = this->CHAIN->rank(std::get<1>(e) + 1, 1);
                tmp_out[i++] =
                    this->matching[convert_p_sigma_pair(entering_chain, this->alphabet_match[std::get<2>(e)])];
                tmp_dict.at(entering_chain - 1).push_back(this->alphabet_match[std::get<2>(e)]);
            }
            this->OUT = std::make_unique<WTLargeAlphabet>(tmp_out);  // OUT

            // alphabet dictionary
            for (const auto& v : tmp_dict) {
                this->ALPHABET_DICT.push_back(std::make_unique<SuccinctSet>(this->sigma, v));
            }
            // compute IN
            // Sort by target->label->source
            std::sort(edge_list.begin(), edge_list.end(), [](const auto& e1, const auto& e2) {
                if (std::get<1>(e1) == std::get<1>(e2)) {
                    if (std::get<2>(e1) == std::get<2>(e2)) {
                        return std::get<0>(e1) < std::get<0>(e2);  // Compare by source if source and target are equal
                    }
                    return std::get<2>(e1) < std::get<2>(e2);  // Compare by label if target is equal
                }
                return std::get<1>(e1) < std::get<1>(e2);  // Compare by target vertex
            });

            sdsl::bit_vector tmp_in(this->m, 0);
            tmp_in[0] = 1;

            std::string tmp_in_string = "" + get<2>(edge_list[0]);
            for (size_t i = 1; i < edge_list.size(); ++i) {
                // IN[k] != IN [k − 1] -> different label -> put 1
                if (std::get<2>(edge_list[i]) != std::get<2>(edge_list[i - 1])) tmp_in[i] = 1;
                // k and k-1 edge different chain
                if (this->CHAIN->rank(std::get<1>(edge_list[i]) + 1, 1) !=
                    this->CHAIN->rank(std::get<1>(edge_list[i - 1]) + 1, 1))
                    tmp_in[i] = 1;
            }

            this->IN = std::make_unique<Custom_BV<B>>(tmp_in);

        } catch (const std::runtime_error& e) {
            // Handle specific runtime errors (like file not found)
            std::cerr << "Runtime error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            // Handle any other exceptions derived from std::exception
            std::cerr << "Error: " << e.what() << std::endl;
        } catch (...) {
            // Catch any other exceptions not covered by the previous handlers
            std::cerr << "An unknown error occurred!" << std::endl;
        }
    }
    if (debug) this->print();
}
template <typename B>
void NFA_Index<B>::print() {
    if (debug) {
        std::cout << "p: " << this->p << std::endl;
        std::cout << "n: " << this->n << std::endl;
        std::cout << "m: " << this->m << std::endl;
        std::cout << "sigma: " << this->sigma << std::endl;
        std::cout << "Chain: ";
        this->CHAIN->print();
        std::cout << "In_degree: ";
        this->IN_DEG->print();
        std::cout << "Out_degree: ";
        this->OUT_DEG->print();
        std::cout << "OUT: ";
        this->OUT->print();
        std::cout << "OUT (as list of pair): ";
        std::cout << std::endl;
        std::cout << "IN: ";
        this->IN->print();

        std::cout << "FINAL: ";
        this->FINAL->print();
        std::cout << "Alphabet mapping:\n";
        for (const auto& pair : this->alphabet_match) {
            std::cout << "(" << pair.first << "->" << pair.second << ") ";
        }
        std::cout << std::endl;

        for (size_t i = 0; i < this->p; ++i) {
            std::cout << "Sigma_" << i + 1 << ": ";
            this->ALPHABET_DICT[i]->print();
        }
    }
}

// get cardinality of the i-th chain
template <typename B>
uint64_t NFA_Index<B>::get_chain_cardinality(uint64_t i) {
    return this->CHAIN->select(i + 1, 1) - this->CHAIN->select(i, 1);
}

// op_1: out(Q_j[1,k], i,a)
//  the number of edges labeled with 'a' that leaves states from Q_j[1,k]
//  and enter in the i-th chain
//  k is the number of nodes to be taken
template <typename B>
uint64_t NFA_Index<B>::out_ijk(uint64_t i, uint64_t j, uint64_t k, uint64_t a) {
    if (k == 0) return 0;
    uint64_t l = this->CHAIN->select(j, 1);
    // issue if j = 0 since j-1 = very big number????

    uint64_t r = l + k - 1;
    uint64_t x = this->OUT_DEG->rank(this->OUT_DEG->select(l, 1), 0);
    uint64_t y = this->OUT_DEG->rank(this->OUT_DEG->select(r + 1, 1), 0);

    // check x compared to y
    if (x == y) {
        if (debug)
            std::cout << "out_ijk("
                      << "i:" << i << " j:" << j << " k:" << k << " a:" << a << ")->" << 0 << std::endl;
        return 0;
    } else {
        uint64_t i_a = convert_p_sigma_pair(i, a);
        if (!(this->matching.find(i_a) != this->matching.end())) return 0;
        auto v = matching[i_a];
        uint64_t res_x = this->OUT->rank(x, v);
        uint64_t res_y = this->OUT->rank(y, v);
        if (debug)
            std::cout << "out_ijk("
                      << "i:" << i << " j:" << j << " k:" << k << " a:" << a << ")->" << res_y - res_x << std::endl;
        return res_y - res_x;
    }
}

/*  op2: in(Q_i[1,k], a) ≤ h
        return k, the biggest value such that the  number of edges labeled with 'a'
        that enter states in Q_j[1,k] is smaller or equal than 'h'

      count_op2: -> number of nodes that can be retrieved

    IMPORTANT: if return 0 means that we cannot take the first nodes since it has a edge labelled a (and h=0) for
            example. If return one just take the first node and like this so on

    ritorna l'indice da [1,card_i]
     */
template <typename B>
uint64_t NFA_Index<B>::in_leq_h(uint64_t i, uint64_t a, uint64_t h) {
    if (!this->ALPHABET_DICT.at(i - 1)->is_member(a)) {
        return get_chain_cardinality(i);
    }
    uint64_t l = this->CHAIN->select(i, 1);
    // uint64_t r = this->CHAIN->select(i + 1, 1) - 1;

    uint64_t x = this->IN_DEG->rank(this->IN_DEG->select(l, 1), 0);

    // f is
    uint64_t f = this->IN->rank(x, 1);
    uint64_t g = this->IN->select((f + this->ALPHABET_DICT.at(i - 1)->rank(a)), 1);

    /* std::cout << "in_leq_h->i:" << i << " h:" << h << " l:" << l << " r:" << r << " x:" << x << " y:" << y << " f:"
       << f
              << " g:" << g << " a:" << a << std::endl; */

    // FARE CHECK SU LABEL ENTRATE IN K SE PIÙ GRANDE DIMINUISCI
    // Case 1): h == 0
    if (h == 0) {
        uint64_t tmp = this->IN_DEG->rank(this->IN_DEG->select(g, 0), 1) + 1;
        return tmp - l;  // k = tmp-l
    } else if (h > 0) {  // case 2): h > 0
        uint64_t h_1 = this->IN->select((f + this->ALPHABET_DICT.at(i - 1)->rank(a) + 1), 1) - g;
        if (h_1 <= h) return this->get_chain_cardinality(i);
        // uint64_t tmp = this->IN_DEG->rank(this->IN_DEG->select(g + h, 0), 1) + 1;

        // controllo del massimo numero di archi
        uint64_t tmp_2 = this->IN_DEG->rank(this->IN_DEG->select(g + h + 1, 0), 1) + 1;
        return tmp_2 - l - 1;  // k = p-l
    } else {                   // case 3): h < 0 -> this isn't accepted
        throw std::invalid_argument("h inside function in_leq_h cannot be negative!");
    }
}

/* op3: in(Q_i[1,t], a) ≥ z
        return t such that the number of edges entering Q_i[1,t] with label 'a'
        is greater that z if it exist

        if it return a value it exist if it returns 0 doesn't exist (doesn't exist also means all set it think?)

        TODO ISSUE?-> if last node has edge entering in it return |Q_i| so it's a problem?

*/
template <typename B>
uint64_t NFA_Index<B>::in_geq_z(uint64_t i, uint64_t a, uint64_t z) {
    // check if z>=1
    if (z < 1) throw std::invalid_argument("z inside function in_geq_z should be greater or equal than 1!");
    if (!this->ALPHABET_DICT.at(i - 1)->is_member(a)) return 0;

    uint64_t l = this->CHAIN->select(i, 1);  // inizio catena
    // uint64_t r = this->CHAIN->select(i + 1, 1) - 1;  // fine catena

    uint64_t x = this->IN_DEG->rank(this->IN_DEG->select(l, 1), 0);  // archi prima della catena

    uint64_t f = this->IN->rank(x, 1);

    // smallest edge labelled 'a' reaching Q_i
    uint64_t g = this->IN->select((f + this->ALPHABET_DICT.at(i - 1)->rank(a)), 1);

    // number of edges labelled 'a' reaching Q_i
    uint64_t h_1 = this->IN->select((f + this->ALPHABET_DICT.at(i - 1)->rank(a) + 1), 1) - g;

    // case where we have less edges than x -> 0
    if (h_1 < z) {
        return 0;  // k = tmp-l
    }
    // case where we have at least z edges

    // return node reached by z-th edge
    uint64_t tmp = this->IN_DEG->rank(this->IN_DEG->select(g + z, 0), 1) + 1;
    return tmp - l;  // k = p-l
}

/*op4: max_lambda(Q_i[h]) < a
        the largest h s.t. the max lambda entering in Q_i[h] is < a

 */
template <typename B>
uint64_t NFA_Index<B>::max_lambda(uint64_t i, uint64_t a) {
    int succ = this->ALPHABET_DICT.at(i - 1)->succ(a);
    if (succ != -1) {
        // succ is defined
        return this->in_leq_h(i, succ, 0);
    } else {
        return this->get_chain_cardinality(i);
    }
}

// return the set of all states reached by alpha
template <typename B>
uint64_t NFA_Index<B>::count(std::string alpha) {
    // iterate over all labels
    if (debug) std::cout << "--------------searching: " << alpha << "------------" << std::endl;
    // salvo gli indici qua dentro da [1,card

    std::vector<uint64_t> l, l_1, t, t_1;

    for (uint64_t i = 1; i <= this->p; ++i) {
        l.push_back(0);
        t.push_back(0);
        l_1.push_back(0);
        t_1.push_back(this->get_chain_cardinality(i));
    }
    for (uint64_t idx = 0; idx < alpha.size(); ++idx) {
        uint64_t a = this->alphabet_match[alpha[idx]];
        if (debug) std::cout << "#---------step: " << alpha[idx] << "----------#" << std::endl;

        // iterate over all chain
        for (uint64_t i = 1; i <= this->p; ++i) {
            if (debug)
                std::cout << "-chain: " << i << " l'_" << i << ": " << l_1[i - 1] << " t'_" << i << ": " << t_1[i - 1]
                          << std::endl;

            // compute c and d
            uint64_t c = 0, d = 0;
            for (uint64_t j = 1; j <= this->p; ++j) {
                c += this->out_ijk(i, j, l_1[j - 1], a);
                d += this->out_ijk(i, j, t_1[j - 1], a);
            }
            if (debug) std::cout << "  c: " << c << " d: " << d << std::endl;

            // issue
            if (c == d) {
                // doesn't exist
                // calcola $l_1 e mettili uguali
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                t[i - 1] = l[i - 1];
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;

            } else if (d > c) {
                // exist
                t[i - 1] = this->in_geq_z(i, a, d);
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;
            } else if (c > d) {
                throw std::invalid_argument("c > d inside count function");
            }
        }
        bool end = true;
        for (uint64_t i = 1; i <= this->p; ++i) {
            // if count = 0 we can stop
            if (t[i - 1] != l[i - 1]) end = false;
            t_1[i - 1] = t[i - 1];
            l_1[i - 1] = l[i - 1];
        }
        if (end) return 0;
    }
    uint64_t total = 0;
    for (uint64_t i = 1; i <= this->p; ++i) {
        total += t[i - 1] - l[i - 1];
    }
    return total;
}

// locate state reached by pattern alpha
template <typename B>
std::vector<std::pair<uint64_t, uint64_t>> NFA_Index<B>::locate(std::string alpha) {
    // iterate over all labels
    if (debug) std::cout << "--------------searching: " << alpha << "------------" << std::endl;
    // salvo gli indici qua dentro da [1,card

    std::vector<uint64_t> l, l_1, t, t_1;

    for (uint64_t i = 1; i <= this->p; ++i) {
        l.push_back(0);
        t.push_back(0);
        l_1.push_back(0);
        t_1.push_back(this->get_chain_cardinality(i));
    }
    for (uint64_t idx = 0; idx < alpha.size(); ++idx) {
        uint64_t a = this->alphabet_match[alpha[idx]];
        if (debug) std::cout << "#---------step: " << alpha[idx] << "----------#" << std::endl;

        // iterate over all chain
        for (uint64_t i = 1; i <= this->p; ++i) {
            if (debug)
                std::cout << "-chain: " << i << " l'_" << i << ": " << l_1[i - 1] << " t'_" << i << ": " << t_1[i - 1]
                          << std::endl;

            // compute c and d
            uint64_t c = 0, d = 0;
            for (uint64_t j = 1; j <= this->p; ++j) {
                c += this->out_ijk(i, j, l_1[j - 1], a);
                d += this->out_ijk(i, j, t_1[j - 1], a);
            }
            if (debug) std::cout << "  c: " << c << " d: " << d << std::endl;

            // issue
            if (c == d) {
                // doesn't exist
                // calcola $l_1 e mettili uguali
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                t[i - 1] = l[i - 1];
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;

            } else if (d > c) {
                // exist
                t[i - 1] = this->in_geq_z(i, a, d);
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;
            } else if (c > d) {
                throw std::invalid_argument("c > d inside count function");
            }
        }
        bool end = true;
        for (uint64_t i = 1; i <= this->p; ++i) {
            // if count = 0 we can stop
            if (t[i - 1] != l[i - 1]) end = false;
            t_1[i - 1] = t[i - 1];
            l_1[i - 1] = l[i - 1];
        }
        if (end) return std::vector<std::pair<uint64_t, uint64_t>>(0);
    }
    std::vector<std::pair<uint64_t, uint64_t>> res;
    for (uint64_t i = 0; i < this->p; ++i) {
        res.push_back({t_1[i], l_1[i] + 1});
    }
    return res;
}

// check if alpha belong to the language of N
template <typename B>
bool NFA_Index<B>::membership(std::string alpha) {
    // iterate over all labels
    if (debug) std::cout << "--------------searching: " << alpha << "------------" << std::endl;
    // salvo gli indici qua dentro da [1,card

    std::vector<uint64_t> l, l_1, t, t_1;

    for (uint64_t i = 1; i <= this->p; ++i) {
        l.push_back(0);
        t.push_back(0);
        l_1.push_back(0);
        t_1.push_back(this->get_chain_cardinality(i));
    }
    for (uint64_t idx = 0; idx < alpha.size(); ++idx) {
        uint64_t a = this->alphabet_match[alpha[idx]];
        if (debug) std::cout << "#---------step: " << alpha[idx] << "----------#" << std::endl;

        // iterate over all chain
        for (uint64_t i = 1; i <= this->p; ++i) {
            if (debug)
                std::cout << "-chain: " << i << " l'_" << i << ": " << l_1[i - 1] << " t'_" << i << ": " << t_1[i - 1]
                          << std::endl;

            // compute c and d
            uint64_t c = 0, d = 0;
            for (uint64_t j = 1; j <= this->p; ++j) {
                c += this->out_ijk(i, j, l_1[j - 1], a);
                d += this->out_ijk(i, j, t_1[j - 1], a);
            }
            if (debug) std::cout << "  c: " << c << " d: " << d << std::endl;

            // issue
            if (c == d) {
                // doesn't exist
                // calcola $l_1 e mettili uguali
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                t[i - 1] = l[i - 1];
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;

            } else if (d > c) {
                // exist
                t[i - 1] = this->in_geq_z(i, a, d);
                l[i - 1] = this->in_leq_h(i, a, c);
                if (l[i - 1] >= 1) {
                    l[i - 1] = std::min(l[i - 1], this->max_lambda(i, a));
                }
                if (debug) std::cout << "t_" << i << ":" << t[i - 1] << "  l_" << i << ":" << l[i - 1] << std::endl;
            } else if (c > d) {
                throw std::invalid_argument("c > d inside count function");
            }
        }
        bool end = true;
        for (uint64_t i = 1; i <= this->p; ++i) {
            // if count = 0 we can stop
            if (t[i - 1] != l[i - 1]) end = false;
            t_1[i - 1] = t[i - 1];
            l_1[i - 1] = l[i - 1];
        }
        if (end) return false;
    }
    return true;
}

template <typename B>
void NFA_Index<B>::store(std::string input_file) {
    std::ofstream out_stream(input_file, std::ios::binary);
    if (!out_stream) {
        std::cerr << "Error opening file for writing!" << std::endl;
    }

    // TODO final states
    out_stream.write(reinterpret_cast<const char*>(&n), sizeof(n));
    out_stream.write(reinterpret_cast<const char*>(&m), sizeof(m));
    out_stream.write(reinterpret_cast<const char*>(&sigma), sizeof(sigma));
    out_stream.write(reinterpret_cast<const char*>(&p), sizeof(p));
    out_stream.write(reinterpret_cast<const char*>(&k), sizeof(k));

    // Serialize each structure
    this->CHAIN->serialize(out_stream);
    this->OUT_DEG->serialize(out_stream);
    this->IN_DEG->serialize(out_stream);
    this->IN->serialize(out_stream);
    this->OUT->serialize(out_stream);
    this->FINAL->serialize(out_stream);

    for (uint64_t i = 0; i < this->p; ++i) this->ALPHABET_DICT[i]->serialize(out_stream);

    // write_binary(out_stream, this->prefix_table);
    write_map_to_ostream(this->alphabet_match, out_stream);
    write_map_to_ostream_u(this->matching, out_stream);
    out_stream.close();
}

template <typename B>
void NFA_Index<B>::load(std::ifstream& input) {
    input.read(reinterpret_cast<char*>(&n), sizeof(n));
    input.read(reinterpret_cast<char*>(&m), sizeof(m));
    input.read(reinterpret_cast<char*>(&sigma), sizeof(sigma));
    input.read(reinterpret_cast<char*>(&p), sizeof(p));
    input.read(reinterpret_cast<char*>(&k), sizeof(k));

    this->CHAIN = std::make_unique<Custom_BV<B>>(input);
    this->OUT_DEG = std::make_unique<Custom_BV<B>>(input);
    this->IN_DEG = std::make_unique<Custom_BV<B>>(input);
    this->IN = std::make_unique<Custom_BV<B>>(input);
    this->OUT = std::make_unique<WTLargeAlphabet>(input);
    this->FINAL = std::make_unique<Custom_BV<B>>(input);

    for (uint64_t i = 0; i < this->p; ++i) {
        this->ALPHABET_DICT.push_back(std::make_unique<SuccinctSet>(input));
    }

    // read_binary(input, this->prefix_table);
    this->alphabet_match = read_map_from_istream(input);
    read_map_from_istream_u(this->matching, input);
}

template <typename B>
void NFA_Index<B>::print_size(std::string file_path) {
    std::ofstream ofs1(file_path);  // Opens in output mode (truncates if exists)
    uint64_t chain_size = this->CHAIN->size_bytes();
    ofs1 << "CHAIN: " << chain_size << " bytes" << std::endl;

    uint64_t final_size = this->FINAL->size_bytes();
    ofs1 << "FINAL: " << final_size << " bytes" << std::endl;

    uint64_t out_deg_size = this->OUT_DEG->size_bytes();
    ofs1 << "OUT_DEG: " << out_deg_size << " bytes" << std::endl;

    uint64_t in_deg_size = this->IN_DEG->size_bytes();
    ofs1 << "IN_DEG: " << in_deg_size << " bytes" << std::endl;

    uint64_t in_size = this->IN->size_bytes();
    ofs1 << "IN: " << in_size << " bytes" << std::endl;

    uint64_t out_size = this->OUT->size_bytes();
    ofs1 << "OUT: " << out_size << " bytes" << std::endl;

    uint64_t alphabet_dict_size = 0;
    for (const auto& dict : this->ALPHABET_DICT) {
        alphabet_dict_size += dict->size_bytes();
    }
    ofs1 << "ALPHABET_DICT: " << alphabet_dict_size << " bytes" << std::endl;
    uint64_t variable_size = sizeof(n) + sizeof(m) + sizeof(sigma) + sizeof(p) + sizeof(k);
    ofs1 << "Variable size: " << variable_size << " bytes" << std::endl;

    ofs1 << "Index size: "
         << chain_size + final_size + out_deg_size + in_deg_size + in_size + out_size + alphabet_dict_size +
                variable_size
         << " bytes" << std::endl;

    ofs1.close();
}

template class NFA_Index<sdsl::bit_vector>;
template class NFA_Index<sdsl::bit_vector_il<64>>;
template class NFA_Index<sdsl::bit_vector_il<128>>;
template class NFA_Index<sdsl::bit_vector_il<256>>;
template class NFA_Index<sdsl::bit_vector_il<512>>;
template class NFA_Index<sdsl::bit_vector_il<1024>>;
template class NFA_Index<sdsl::bit_vector_il<2048>>;
template class NFA_Index<sdsl::rrr_vector<15>>;
template class NFA_Index<sdsl::rrr_vector<63>>;
template class NFA_Index<sdsl::rrr_vector<127>>;
template class NFA_Index<sdsl::rrr_vector<40>>;
template class NFA_Index<sdsl::rrr_vector<8>>;
template class NFA_Index<sdsl::sd_vector<sdsl::bit_vector>>;