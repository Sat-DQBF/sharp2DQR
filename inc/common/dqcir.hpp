#ifndef DQCIR_HPP
#define DQCIR_HPP

#include <string>
#include <unordered_map>
#include <vector>

#include "common/verilog.hpp"

// Data structure for DQBF in circuit form
class dqcir {
   public:
    // Number of variables
    uint64_t var_cnt;

    dqcir();

    // Initialize from dqcir file
    void from_file(std::string path, bool rename = false);
    std::unordered_map<std::string, std::string> rename_map;

    // Print problem info
    void print_stat(bool detailed = false);

    // Universal variables, (name)
    std::vector<std::string> u_vars;

    // Existential variables, (name, dependency set)
    std::vector<std::pair<std::string, std::vector<std::string>>> e_vars;

    VerilogModule phi;
};
#endif