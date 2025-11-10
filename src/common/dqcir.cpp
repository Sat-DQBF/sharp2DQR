#include "common/dqcir.hpp"

dqcir::dqcir() {
    
}

void dqcir::print_stat(bool detailed) {
    printf("-------- Stat --------\n%lu universal var(s):\n", u_vars.size());
    if (detailed) {
        for (auto& x : u_vars) {
            printf("%s ", x.c_str());
        }
        putchar_unlocked('\n');
    }
    printf("%zu existential var(s):\n", e_vars.size());
    if (detailed) {
        for (auto& y : e_vars) {
            printf("%s: ", y.first.c_str());
            for (auto& x : y.second) {
                printf("%s ", x.c_str());
            }
            putchar_unlocked('\n');
        }
    }
}
