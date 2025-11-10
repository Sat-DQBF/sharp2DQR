#include "algo/algo.hpp"
#include "common/utils.hpp"

#include <assert.h>

Circuit Algo::bv_at(std::string s, int i) {
    return Circuit(s + std::to_string(i));
}

Circuit Algo::r_eq(int i, int j) {
    return bv_at("R", i) == bv_at("R", j);
}

Circuit Algo::r_eq(int low1, int low2, int diff) {
    Circuit tmp;
    tmp.operation = "&";
    for (int i = 0; i <= diff; i++) {
        tmp.children.push_back(bv_at("R", low1 + i) == bv_at("R", low2 + i));
    }
    return tmp;
}

Circuit Algo::rn_eq(int i, int j) {
    return bv_at("Rn", i) == bv_at("Rn", j);
}

Circuit Algo::rn_eq(int low1, int low2, int diff) {
    Circuit tmp;
    tmp.operation = "&";
    for (int i = 0; i <= diff; i++) {
        tmp.children.push_back(bv_at("Rn", low1 + i) == bv_at("Rn", low2 + i));
    }
    return tmp;
}

Circuit Algo::r_eq_rn(int high, int low) {
    assert (low <= high);
    Circuit tmp("");
    tmp.operation = "&";
    for (int i = low; i <= high; i++) {
        tmp.children.push_back(bv_at("R", i) == bv_at("Rn", i));
    }
    return tmp;
}