#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/program_options.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <fstream>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <string>
#include <vector>

#include "DataStructure/FullyIndexableDictionary.hpp"
#include "DataStructure/NFAIndex.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    std::string input_file = "", operation = "", out_file = "index.sdsl";
    size_t build_type = 0, chain_index = 1;
    bool debug = false;
    try {
        po::options_description desc("Options");
        desc.add_options()("help,h", "Help screen")(
            "build,b", po::value<size_t>()->value_name("BUILD_TYPE"),
            "Build construction type: 0 read a NFA a build index, 1 load the index from file")(
            "input,i", po::value<std::string>()->value_name("INPUT_PATH"), "Input path of the NFA or Index")(
            "output,o", po::value<std::string>()->value_name("OUTPUT_PATH"),
            "Output path for the 0 construction type (default: index.sdsl)")(
            "operation,p", po::value<std::string>()->value_name("OPERATION"),
            "Operation to be performed: count, locate or membership")(
            "chain,c", po::value<std::size_t>()->value_name("CHAIN"),
            "Use Chain partition (CHAIN) to build the Index (if loaded from NFA, default 1)")(
            "debug,d", "Debug options to add some debug output")(
            "size,s", po::value<std::string>()->value_name("SIZE_PATH"),
            "Print size of the index structures in Bytes inside the SIZE_PATH file");
        po::variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("input")) {
            input_file = vm["input"].as<std::string>();
        }
        if (vm.count("build")) {
            build_type = vm["build"].as<size_t>();
        }
        if (vm.count("operation")) {
            operation = vm["operation"].as<std::string>();
        }
        if (vm.count("chain")) {
            chain_index = vm["chain"].as<size_t>();
        }
        if (vm.count("output")) {
            out_file = vm["output"].as<std::string>();
        }
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        if (vm.count("debug")) {
            debug = true;
        }
        // sds
        NFA_Index<sdsl::bit_vector_il<512>> NFA(input_file, chain_index, build_type, debug);
        if (build_type == 0) NFA.store(out_file);

        if (vm.count("size")) {
            NFA.print_size(vm["size"].as<std::string>());
        }

        // TODO read values since input operations is present
        if (operation == "count") {
            std::string alpha = "";
            while (std::cin >> alpha) {
                std::cout << NFA.count(alpha) << std::endl;
            }
        } else if (operation == "locate") {
            std::string alpha = "";
            while (true) {
                std::cin >> alpha;
                auto res = NFA.locate(alpha);
                for (auto e : res) {
                    std::cout << e.first << " " << e.second << " ";
                }
                std::cout << std::endl;
            }
        } else if (operation == "membership") {
            std::string alpha = "";
            while (true) {
                std::cin >> alpha;
                auto res = NFA.membership(alpha);
                std::cout << res << std::endl;
            }
        }
    } catch (const po::error& ex) {
        std::cerr << ex.what() << '\n';
    }
    return 0;
}