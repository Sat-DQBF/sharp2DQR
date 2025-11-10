#include "algo/algo.hpp"
#include "common/dqcir.hpp"
#include "common/utils.hpp"

#include "abc_interface/abc_interface.hpp"

#include <cxxopts.hpp>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <readline/readline.h>
#include <readline/history.h>


int main(int argc, char** argv) {
    cxxopts::Options options("sharp2DQR", "sharp2DQR");

    options.add_options()   ("i,input", "Input File", cxxopts::value<std::string>())
                            ("r,rename", "Rename variables", cxxopts::value<bool>()->default_value("false"))
                            ("no_reorder", "Disable BDD reordering", cxxopts::value<bool>()->default_value("false"))
                            ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
                            ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Cleaning
    {
        std::vector<std::string> files = {"./twoDQR.v", "./twoDQR_inv.v", "./twoDQR_inv.pla", "./check_skolem.v", "./skolem.v", "./skolem.aig", "./result.txt"};
        for (auto& file : files) {
            if (file_exists(file)) {
                std::filesystem::remove(file);
            }
        }
    }

    print_debug("DEBUG mode enabled");
    print_info("Algorithm = 2DQR");

    if (!result.count("input")) {
        print_error("No input file specified");
    }
    std::string input_file = result["input"].as<std::string>();
    print_info(("file = " + input_file).c_str());

    dqcir p;
    if (split_string(input_file, ".").back() == "dqcir") {
        p.from_file(input_file, result["rename"].as<bool>());
    } else {
        print_error("File extension must be .dicir");
    }
    Algo algo(p);
    algo.abc.verbose = result["verbose"].as<bool>();
    algo.args.reorder = !result["no_reorder"].as<bool>();
    algo.run();
}
