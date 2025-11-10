#include "algo/algo.hpp"
#include "common/utils.hpp"

#include <assert.h>
#include <fstream>
#include <fmt/format.h>
#include <set>

extern "C"
{
    // Functions in base/abci
    Aig_Man_t * Abc_NtkToDar( Abc_Ntk_t * pNtk, int fExors, int fRegisters );
    DdNode * Bbr_NodeGlobalBdds_rec( DdManager * dd, Aig_Obj_t * pNode, int nBddSizeMax, int fDropInternal, ProgressBar * pProgress, int * pCounter, int fVerbose );
    DdNode * Extra_TransferPermute( DdManager * ddSource, DdManager * ddDestination, DdNode * f, int * Permute );
}

DdNode* Algo::get_candidate() {
    DdNode* bCand_f1_dd_graph_closure = dd_graph_closure.reached; Cudd_Ref(bCand_f1_dd_graph_closure);
    DdManager* dd = dd_graph_closure.dd;

    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);
    auto assign = [&](int i, bool val) {
        // https://stackoverflow.com/questions/55246590/cudd-manipulation-of-bdds
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
        bCand_f1_dd_graph_closure = Cudd_bddAnd(dd, bCand_f1_dd_graph_closure, val ? src : Cudd_Not(src)); Cudd_Ref(bCand_f1_dd_graph_closure);
        exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
        Cudd_RecursiveDeref(dd, src);
    };

    auto eq = [&](int i, int j) {
        DdNode* l = Cudd_bddIthVar(dd, i); Cudd_Ref(l);
        DdNode* r = Cudd_bddIthVar(dd, j); Cudd_Ref(r);
        DdNode* eq = Cudd_bddXnor(dd, l, r); Cudd_Ref(eq);
        bCand_f1_dd_graph_closure = Cudd_bddAnd(dd, bCand_f1_dd_graph_closure, eq); Cudd_Ref(bCand_f1_dd_graph_closure);
        exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, l); Cudd_Ref(exist_abs_cube);
        Cudd_RecursiveDeref(dd, l);
        Cudd_RecursiveDeref(dd, r);
        Cudd_RecursiveDeref(dd, eq);
    };

    std::set<std::string> z1(p.e_vars[0].second.begin(), p.e_vars[0].second.end());

    assign(dd_graph_closure.idx_lookup["k"], 0); // k = 0
    assign(dd_graph_closure.idx_lookup["y_k"], 1); // y_k = 1

    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z1.find(p.u_vars[i]) == z1.end()) {
            assign(dd_graph_closure.idx_lookup["x"] + i, 0); // u_i = 0
        }
    }
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        eq(dd_graph_closure.idx_lookup[p.e_vars[0].second[i]], dd_graph_closure.idx_lookup["target z_k"] + i);
    }
    assign(dd_graph_closure.idx_lookup["target y_k"], 0); // target y_k = 0
    bCand_f1_dd_graph_closure = Cudd_bddExistAbstract(dd, bCand_f1_dd_graph_closure, exist_abs_cube); Cudd_Ref(bCand_f1_dd_graph_closure);

    DdNode* bCand_f1 = Extra_TransferPermute(dd, dd_phi.dd, bCand_f1_dd_graph_closure, dd_graph_closure_to_dd_phi); Cudd_Ref(bCand_f1);
    Cudd_RecursiveDeref(dd, bCand_f1_dd_graph_closure);
    return bCand_f1;
}

int Algo::check_candidate(DdNode * candidate, char * cex) {
    DdManager* dd = dd_phi.dd;
    DdNode* y2 = Cudd_bddIthVar(dd, p.u_vars.size() + 1); Cudd_Ref(y2);

    DdNode* bNPhi_x_f1_y2 = Cudd_bddIte(dd, candidate, dd_phi.bNPhi_x_1_y2, dd_phi.bNPhi_x_0_y2); Cudd_Ref(bNPhi_x_f1_y2);

    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);

    DdNode* bNPhi_x_f1_1 = Cudd_bddAnd(dd, bNPhi_x_f1_y2, y2); Cudd_Ref(bNPhi_x_f1_1);
    DdNode* bNPhi_x_f1_0 = Cudd_bddAnd(dd, bNPhi_x_f1_y2, Cudd_Not(y2)); Cudd_Ref(bNPhi_x_f1_0);

    exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, y2); Cudd_Ref(exist_abs_cube);
    bNPhi_x_f1_0 = Cudd_bddExistAbstract(dd, bNPhi_x_f1_0, exist_abs_cube); Cudd_Ref(bNPhi_x_f1_0);

    std::set<std::string> z2(p.e_vars[1].second.begin(), p.e_vars[1].second.end());
    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z2.find(p.u_vars[i]) == z2.end()) {
            DdNode* z = Cudd_bddIthVar(dd, u_var_idx_to_dd_phi_idx[i]); Cudd_Ref(z);
            exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, z); Cudd_Ref(exist_abs_cube);
            Cudd_RecursiveDeref(dd, z);
        }
    }

    bNPhi_x_f1_1 = Cudd_bddExistAbstract(dd, bNPhi_x_f1_1, exist_abs_cube); Cudd_Ref(bNPhi_x_f1_1);

    DdNode* bNPhi_z_f1_star = Cudd_bddAnd(dd, bNPhi_x_f1_0, bNPhi_x_f1_1); Cudd_Ref(bNPhi_z_f1_star);

    Cudd_RecursiveDeref(dd, y2);
    Cudd_RecursiveDeref(dd, bNPhi_x_f1_y2);
    Cudd_RecursiveDeref(dd, bNPhi_x_f1_0);
    Cudd_RecursiveDeref(dd, bNPhi_x_f1_1);
    Cudd_RecursiveDeref(dd, exist_abs_cube);

    if (Cudd_bddPickOneCube(dd, bNPhi_z_f1_star, cex)) {
        // printf("CEX: ");
        // for (int i = 0; i < dd->size; i++) {
        //     printf("%i", cex[i]);
        // }
        // printf("\n");
        Cudd_RecursiveDeref(dd, bNPhi_z_f1_star);
        return 0;
    }
    Cudd_RecursiveDeref(dd, bNPhi_z_f1_star);
    return 1;
}

int Algo::evaluate_candidate(DdNode * candidate, char * assignment) {
    DdManager* dd = dd_phi.dd;
    DdNode* candidate_ = candidate; Cudd_Ref(candidate_);
    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        DdNode* src = Cudd_bddIthVar(dd, u_var_idx_to_dd_phi_idx[dd_graph_closure.idx_lookup[p.e_vars[0].second[i]] - dd_graph_closure.idx_lookup["x"]]); Cudd_Ref(src);
        candidate_ = Cudd_bddAnd(dd, candidate_, assignment[i] == 1 ? src : Cudd_Not(src)); Cudd_Ref(candidate_);
        exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
        Cudd_RecursiveDeref(dd, src);
    }

    DdNode* result = Cudd_bddExistAbstract(dd, candidate_, exist_abs_cube); Cudd_Ref(result);
    int result_val = !Cudd_bddLeq(dd, result, Cudd_Not(dd->one));
    Cudd_RecursiveDeref(dd, result);
    Cudd_RecursiveDeref(dd, candidate_);
    Cudd_RecursiveDeref(dd, exist_abs_cube);
    return result_val;
}

DdNode* Algo::force_assignment_cube(char * z, bool val) {
    DdManager* dd = dd_graph_closure.dd;

    std::string z_str = "";
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        z_str += z[i] == 1 ? "1" : "0";
    }
    // print_info(fmt::format("Forcing z = {} to be {}", z_str, val).c_str());
    if (dd_graph_closure.forced_val.find(z_str) != dd_graph_closure.forced_val.end()) {
        print_error("Patch already exists");
    }
    dd_graph_closure.forced_val[z_str] = val;
    DdNode* cube = dd->one; Cudd_Ref(cube);

    auto assign = [&](int i, bool val) {
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
        cube = Cudd_bddAnd(dd, cube, val ? src : Cudd_Not(src)); Cudd_Ref(cube);
        Cudd_RecursiveDeref(dd, src);
    };

    auto fix = [&](int i) {
        DdNode* r_i = Cudd_bddIthVar(dd, i); Cudd_Ref(r_i);
        DdNode* r_i_next = Cudd_bddIthVar(dd, dd_graph_closure.register_size + i); Cudd_Ref(r_i_next);
        DdNode* r_i_eq_r_i_next = Cudd_bddXnor(dd, r_i, r_i_next); Cudd_Ref(r_i_eq_r_i_next);
        cube = Cudd_bddAnd(dd, cube, r_i_eq_r_i_next); Cudd_Ref(cube);
        Cudd_RecursiveDeref(dd, r_i);
        Cudd_RecursiveDeref(dd, r_i_next);
        Cudd_RecursiveDeref(dd, r_i_eq_r_i_next);
    };

    assign(dd_graph_closure.idx_lookup["k"], 0);

    if (val == 0) {
        assign(dd_graph_closure.idx_lookup["y_k"], 1);
        assign(dd_graph_closure.idx_lookup["y_k"] + dd_graph_closure.register_size, 0);
    } else {
        assign(dd_graph_closure.idx_lookup["y_k"], 0);
        assign(dd_graph_closure.idx_lookup["y_k"] + dd_graph_closure.register_size, 1);
    }

    std::set<std::string> z1(p.e_vars[0].second.begin(), p.e_vars[0].second.end());

    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        assign(dd_graph_closure.idx_lookup[p.e_vars[0].second[i]], z[i] == 1);
    }

    // printf(" to be %d\n", !(val == 0));
    // Everything except y_k is fixed
    fix(dd_graph_closure.idx_lookup["k"]);
    for (int i = dd_graph_closure.idx_lookup["y_k"] + 1; i < dd_graph_closure.register_size; i++) {
        fix(i);
    }

    return cube;
}

DdNode* Algo::patch_transition_with_cex(DdNode * curr_transition, DdNode* candidate, char * cex) {
    char * patch = (char *) malloc(sizeof(char) * (p.e_vars[0].second.size()));
    char * patch_ptr = patch;
    std::set<std::string> z0(p.e_vars[0].second.begin(), p.e_vars[0].second.end());
    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z0.find(p.u_vars[i]) != z0.end()) {
            *(patch_ptr) = (cex[u_var_idx_to_dd_phi_idx[i]] == 1);
            patch_ptr++;
        }
    }

    DdNode * patch_cube = force_assignment_cube(patch, evaluate_candidate(candidate, patch)); Cudd_Ref(patch_cube);
    DdNode * next_transition = Cudd_bddOr(dd_graph_closure.dd, curr_transition, patch_cube); Cudd_Ref(next_transition);

    Cudd_RecursiveDeref(dd_graph_closure.dd, patch_cube);
    free(patch);
    return next_transition;
}


DdNode* Algo::get_valid_candidate(DdNode * curr_transition) {
    dd_graph_closure.forced_val.clear();
    DdManager* dd = dd_graph_closure.dd;
    DdNode* candidate = get_candidate(); Cudd_Ref(candidate);
    if (curr_transition == nullptr) {
        curr_transition = dd_graph_closure.transition; Cudd_Ref(curr_transition);
    }
    DdNode* next_transition;
    int result;
    char * cex = (char *) malloc(sizeof(char) * dd_phi.dd->size);
    DdNode* initial = dd_graph_closure.reached; Cudd_Ref(initial);
    while (!check_candidate(candidate, cex)) {
        next_transition = patch_transition_with_cex(curr_transition, candidate, cex); Cudd_Ref(next_transition);
        // Cudd_RecursiveDeref(dd, curr_transition);
        curr_transition = next_transition; Cudd_Ref(curr_transition);
        Cudd_RecursiveDeref(dd, next_transition);
        
        Cudd_RecursiveDeref(dd, dd_graph_closure.reached);
        result = reachTransRelation(dd, initial, curr_transition, Cudd_Not(dd->one), pPars, dd_graph_closure.reached, dd_graph_closure.register_size);
        Cudd_AutodynEnable(dd, CUDD_REORDER_SYMM_SIFT);
        initial = dd_graph_closure.reached; Cudd_Ref(initial);
        if (result == -1) {
            print_info("2DQR UNKNOWN (Valid)");
            return nullptr;
        }
        assert (result == 1);

        Cudd_RecursiveDeref(dd_phi.dd, candidate);
        candidate = get_candidate(); Cudd_Ref(candidate);
    }
    return candidate;
}