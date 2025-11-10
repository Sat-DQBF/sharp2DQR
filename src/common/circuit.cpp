#include "common/circuit.hpp"
#include "common/utils.hpp"

#include <assert.h>

Circuit::Circuit() {
    name = "";
    children = std::vector<Circuit>();
    operation = "";
}

Circuit::Circuit(std::string name) {
    this->name = name;
    children = std::vector<Circuit>();
    operation = "";
}

bool Circuit::is_leaf() {
    return children.size() == 0;
}

Circuit Circuit::operator~() {
    Circuit cir = Circuit();
    cir.operation = "~";
    cir.children.push_back(*this);
    return cir;
}

Circuit Circuit::operator&(Circuit other) {
    Circuit cir = Circuit();
    cir.operation = "&";
    cir.children.push_back(*this);
    cir.children.push_back(other);
    return cir;
}

Circuit Circuit::operator|(Circuit other) {
    Circuit cir = Circuit();
    cir.operation = "|";
    cir.children.push_back(*this);
    cir.children.push_back(other);
    return cir;
}

Circuit Circuit::operator^(Circuit other) {
    Circuit cir = Circuit();
    cir.operation = "^";
    cir.children.push_back(*this);
    cir.children.push_back(other);
    return cir;
}

Circuit Circuit::operator==(Circuit other) {
    Circuit cir = Circuit();
    cir.operation = "==";
    cir.children.push_back(*this);
    cir.children.push_back(other);
    return cir;
}

Circuit Circuit::implies(Circuit other) {
    return (~(*this)) | other;
}

std::string Circuit::to_string() {
    if (is_leaf()){
        return name;
    } else {
        if (operation == "~") {
            assert(children.size() == 1);
            if (children[0].is_leaf()) {
                return "~" + children[0].name;
            } else {
                if (children[0].operation == "~") {
                    return children[0].children[0].to_string();
                }
                return "~(" + children[0].to_string() + ")";
            }
        } else if (operation == "==") {
            assert(children.size() == 2);
            // std::string tmp;
            // if (children[0].is_leaf()) {
            //     tmp = children[0].name;
            // } else {
            //     tmp = "(" + children[0].to_string() + ")";
            // }
            // if (children[1].is_leaf()) {
            //     return tmp + " == " + children[1].name;
            // } else {
            //     return tmp + " == (" + children[1].to_string() + ")";
            // }
            // return tmp;
            // return ((~children[0] | children[1]) & (~children[1] | children[0])).to_string();
            return (~(children[0] ^ children[1])).to_string();
        } else if (operation == "&" || operation == "|" || operation == "^") {
            assert(children.size() > 0);
            std::vector<std::string> child_strs;
            std::string tmp;
            for (auto child : children) {
                tmp = child.to_string();
                if (child.is_leaf()) {
                    child_strs.push_back(tmp);
                } else {
                    child_strs.push_back("(" + tmp + ")");
                }
            }
            return join(child_strs, operation);
        } else if (operation == "~&" || operation == "~|") {
            assert(children.size() > 0);
            std::vector<std::string> child_strs;
            std::string tmp;
            for (auto child : children) {
                tmp = child.to_string();
                if (child.is_leaf()) {
                    child_strs.push_back(tmp);
                } else {
                    child_strs.push_back("(" + tmp + ")");
                }
            }
            return "~(" + join(child_strs, operation.substr(1, 1)) + ")";
        } else {
            print_error(("Unsupported operation: " + operation).c_str());
            return "";
        }
    }
}