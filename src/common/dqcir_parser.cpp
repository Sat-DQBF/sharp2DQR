#include "common/dqcir.hpp"
#include "common/utils.hpp"

#include <assert.h>
#include <fstream>
#include <iostream>
#include <set>

// Read from dqcir file
void dqcir::from_file(std::string path, bool rename) {
    if (!path.size()) {
        print_error("No file given");
    }
    if (split_string(path, ".").back() != "dqcir") {
        print_warning("File does not ends in .dqcir");
    }
    std::ifstream file(path);
    if (!file) {
        print_error("Error on opening file");
    }

    uint line_cnt = 0;
    std::string line;
    std::string name;
    std::vector<std::string> parts;
    phi.name = "phi";

    var_cnt = 0;

    // Read U/E variables
    while (getline(file, line)) {
        line_cnt++;
        if (line[0] == '#') {
            continue;
        }
        parts = split_string(line, "= (),\n\r");
        if (parts[0] == "forall") {
            for (auto it = parts.begin() + 1; it < parts.end(); it++) {
                if (rename) {
                    name = "u" + std::to_string(u_vars.size() + 1);
                } else {
                    name = *it;
                }
                rename_map[*it] = name;
                u_vars.push_back(name);
                phi.PI.push_back(name);
            }
        } else if (parts[0] == "depend") {
            if (rename) {
                name = "e" + std::to_string(e_vars.size() + 1);
            } else {
                name = parts[1];
            }
            rename_map[parts[1]] = name;
            e_vars.emplace_back(name, std::vector<std::string>());
            phi.PI.push_back(name);
            for (auto it = parts.begin() + 2; it < parts.end(); it++) {
                e_vars.back().second.push_back(rename_map[*it]);
            }
        } else if (parts[0] == "exists") {
            for (auto it = parts.begin() + 1; it < parts.end(); it++) {
                if (rename) {
                    name = "e" + std::to_string(e_vars.size() + 1);
                } else {
                    name = *it;
                }
                rename_map[*it] = name;
                e_vars.emplace_back(name, std::vector<std::string>(u_vars));
                phi.PI.push_back(name);
            }
        } else if (parts[0] == "output") {
            if (rename) {
                name = "out";
            } else {
                name = parts[1];
            }
            rename_map[parts[1]] = name;
            phi.PO.push_back(name);
            break;
        } else {
            parse_err_msg(line_cnt, "Expected output variable before circuit");
        }
    }

    var_cnt = u_vars.size() + e_vars.size();

    std::set<std::string> allowed_oprs = {"and", "or", "not", "nand", "nor", "xor"};
    std::unordered_map<std::string, std::string> opr_map = {
        {"and", "&"}, {"or", "|"}, {"not", "~"}, {"nand", "~&"}, {"nor", "~|"}, {"xor", "^"}};

    int wire_cnt = 0;
    while (getline(file, line)) {
        line_cnt++;
        if (line.size() > 0 && line[0] != '#' && line[0] != '\n' && line[0] != '\r' && line[0] != ' ') {
            parts = split_string(line, "= (),\n\r");
            if (allowed_oprs.find(parts[1]) != allowed_oprs.end()) {
                if (rename) {
                    if (rename_map.find(parts[0]) == rename_map.end()) {
                        name = "w" + std::to_string(++wire_cnt);
                    } else {
                        name = rename_map[parts[0]];
                    }
                } else {
                    name = parts[0];
                }
                rename_map[parts[0]] = name;
                phi.wires.push_back(Wire(name));
                phi.wires.back().assign.operation = opr_map[parts[1]];
                for (int i = 2; i < parts.size(); i++) {
                    name = parts[i];
                    if (name[0] != '-') {
                        assert(rename_map.find(name) != rename_map.end());
                        name = rename_map[name];
                        phi.wires.back().assign.children.push_back(Circuit(name));
                    } else {
                        name = name.substr(1, name.size() - 1);
                        assert(rename_map.find(name) != rename_map.end());
                        name = rename_map[name];
                        phi.wires.back().assign.children.push_back(~Circuit(name));
                    }
                }
            } else {
                parse_err_msg(line_cnt, "Unsupported operator");
            }
        }
    }
};
