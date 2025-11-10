#ifndef VERILOG_HPP
#define VERILOG_HPP

#include "common/circuit.hpp"

#include <string>
#include <unordered_map>
#include <vector>

class Wire {
  public:
    std::string name;
    Circuit assign;

    Wire();
    Wire(std::string name);
    Wire(std::string name, Circuit assign);
};

class VerilogModule {
  public:
    std::string name;
    std::vector<std::string> PI;
    std::unordered_map<std::string, int> PI_size;
    std::vector<std::string> PO;
    std::vector<Wire> wires;
    struct submodule {
        VerilogModule& module;
        std::string name;
        std::vector<Circuit> inputs;
    };

    std::vector<submodule> submodules;

    VerilogModule();
    VerilogModule(std::string name);
    static VerilogModule from_file(std::string filename);
    std::string to_string();

  private:
    bool _from_file = false;
    std::string filename;
};

#endif