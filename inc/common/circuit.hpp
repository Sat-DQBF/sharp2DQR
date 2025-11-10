#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include <string>
#include <vector>

struct Input {
  std::string name;
  bool sign;
};

struct Gate {
  std::string name;
  std::string operation;
  std::vector<Input> inputs;
};

class Circuit {
  public:
    std::string name;
    std::vector<Circuit> children;
    std::string operation;

    Circuit();
    Circuit(std::string name);

    bool is_leaf();

    Circuit operator~();
    Circuit operator&(Circuit other);
    Circuit operator|(Circuit other);
    Circuit operator^(Circuit other);
    Circuit operator==(Circuit other);
    Circuit implies(Circuit other);

    std::string to_string();
};

#endif