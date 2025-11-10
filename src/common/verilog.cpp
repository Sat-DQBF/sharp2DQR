#include "common/verilog.hpp"
#include "common/utils.hpp"

#include <assert.h>
#include <algorithm>
#include <fstream>
#include <set>
#include <vector>

Wire::Wire() {
    name = "";
    assign = Circuit();
}

Wire::Wire(std::string name) {
    this->name = name;
    this->assign = Circuit();
}

Wire::Wire(std::string name, Circuit assign) {
    this->name = name;
    this->assign = assign;
}

VerilogModule::VerilogModule() {
    name = "";
    PI = std::vector<std::string>();
    PO = std::vector<std::string>();
    wires = std::vector<Wire>();
}

VerilogModule::VerilogModule(std::string name) {
    this->name = name;
    PI = std::vector<std::string>();
    PO = std::vector<std::string>();
    wires = std::vector<Wire>();
}

std::string VerilogModule::to_string() {
    if (_from_file) {
        return read_file(filename);
    }

    std::string s;
    s += "module " + name + " (\n";
    s += "    " + join(PI, ",\n    ") + ",\n";
    s += "    " + join(PO, ",\n    ") + "\n";
    s += ");\n\n";
    for (auto t: PI) {
        if (PI_size.find(t) != PI_size.end() && PI_size[t] > 1) {
            s += "input [" + std::to_string(PI_size[t] - 1) + ":0] " + t + ";\n";
        } else {
            s += "input " + t + ";\n";
        }
    }
    for (auto t: PO) {
        s += "output " + t + ";\n";
    }

    for (auto& submodule : submodules) {
        for (auto& output : submodule.module.PO) {
            s += "wire " + submodule.name + "_" + output + ";\n";
        }
    }

    for (auto wire : wires) {
        if (std::find(PO.begin(), PO.end(), wire.name) != PO.end())
            continue;
        s += "wire " + wire.name + ";\n";
    }
    for (auto wire : wires) {
        s += "assign " + wire.name + " = " + wire.assign.to_string() + ";\n";
    }

    std::vector<VerilogModule> unique_submodules;
    std::set<std::string> unique_submodule_names;

    for (auto& submodule : submodules) {
        if (unique_submodule_names.find(submodule.module.name) == unique_submodule_names.end()) {
            unique_submodule_names.insert(submodule.module.name);
            unique_submodules.push_back(submodule.module);
        }

        s += submodule.module.name + " " + submodule.name + " (\n";
        assert(submodule.inputs.size() == submodule.module.PI.size());
        for (int i = 0; i < submodule.inputs.size(); i++) {
            s += "    ." + submodule.module.PI[i] + " (" + submodule.inputs[i].to_string() + "),\n";
        }
        for (auto& output : submodule.module.PO) {
            s += "    ." + output + " (" + submodule.name + "_" + output + ")";
            if (output != submodule.module.PO.back()) {
                s += ",\n";
            } else {
                s += "\n";
            }
        }
        s += ");\n";
    }
    s += "endmodule\n";

    for (auto& submodule : unique_submodules) {
        s += submodule.to_string();
    }
    return s;
}

// Warning: This function is not robust and may not work for all verilog files
VerilogModule VerilogModule::from_file(std::string filename) {
    VerilogModule module;
    module._from_file = true;
    if (!file_exists(filename)) {
        print_error(("File " + filename + " does not exist").c_str());
    }
    module.filename = filename;
    std::ifstream file(filename);
    read_until(file, "module ");
    module.name = read_until(file, " ", 0);
    read_until(file, "input ");
    module.PI = split_string(read_until(file, ";", 0, " \n\t"), ",");

    read_until(file, "output ");
    module.PO = split_string(read_until(file, ";", 0, " \n\t"), ",");
    file.close();
    
    return module;
}