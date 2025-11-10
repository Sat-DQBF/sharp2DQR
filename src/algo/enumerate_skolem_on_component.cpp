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

DdNode* Cudd_bddSubstitude(DdManager * dd, DdNode * f, DdNode * x, DdNode * g) {
    // Returns f with x substituted with g

    DdNode* f_0 = Cudd_bddAnd(dd, f, Cudd_Not(x)); Cudd_Ref(f_0);
    DdNode* f_0_abs = Cudd_bddExistAbstract(dd, f_0, x); Cudd_Ref(f_0_abs);
    DdNode* f_1 = Cudd_bddAnd(dd, f, x); Cudd_Ref(f_1);
    DdNode* f_1_abs = Cudd_bddExistAbstract(dd, f_1, x); Cudd_Ref(f_1_abs);

    DdNode* f_g = Cudd_bddIte(dd, g, f_1_abs, f_0_abs); Cudd_Ref(f_g);

    Cudd_RecursiveDeref(dd, f_0);
    Cudd_RecursiveDeref(dd, f_1);
    Cudd_RecursiveDeref(dd, f_0_abs);
    Cudd_RecursiveDeref(dd, f_1_abs);
    return f_g;
};

SparseInt Algo::count_1DQBF(DdNode* f1, DdNode* component_y2) {
    DdManager * dd = dd_phi.dd;
    DdNode* bNPhi_x_f1_y2 = Cudd_bddIte(dd, f1, dd_phi.bNPhi_x_1_y2, dd_phi.bNPhi_x_0_y2); Cudd_Ref(bNPhi_x_f1_y2);
    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);

    std::set<std::string> z2(p.e_vars[1].second.begin(), p.e_vars[1].second.end());
    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z2.find(p.u_vars[i]) == z2.end()) {
            DdNode* z = Cudd_bddIthVar(dd, u_var_idx_to_dd_phi_idx[i]); Cudd_Ref(z);
            exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, z); Cudd_Ref(exist_abs_cube);
            Cudd_RecursiveDeref(dd, z);
        }
    }

    DdNode* y2 = Cudd_bddIthVar(dd, p.u_vars.size() + 1); Cudd_Ref(y2);
    exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, y2); Cudd_Ref(exist_abs_cube);
    bNPhi_x_f1_y2 = Cudd_bddExistAbstract(dd, bNPhi_x_f1_y2, exist_abs_cube); Cudd_Ref(bNPhi_x_f1_y2);
    DdNode* part = Cudd_Not(Cudd_bddOr(dd, bNPhi_x_f1_y2, Cudd_Not(component_y2))); Cudd_Ref(part);
    long long cnt = (long long) Cudd_CountMinterm(dd, part, p.e_vars[1].second.size());
    Cudd_RecursiveDeref(dd, y2);
    Cudd_RecursiveDeref(dd, exist_abs_cube);
    Cudd_RecursiveDeref(dd, part);
    Cudd_RecursiveDeref(dd, bNPhi_x_f1_y2);
    return SparseInt(1) << cnt;
}

int Algo::evaluate_candidate_dd_projected_graph(DdNode * candidate, char * assignment) {
    DdManager* dd = dd_projected_graph.dd;
    DdNode* candidate_ = candidate; Cudd_Ref(candidate_);
    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
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

SparseInt Algo::enumerate_skolem_on_component(DdNode* component_y1, DdNode* component_y2) {
    DdManager * dd_A = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    DdNode * component_y1_dd_A = nullptr;
    DdNode * constraints = dd_A->one; Cudd_Ref(constraints);
    std::vector<DdNode *> candidate_list_dd_projected_graph;
    SparseInt solution_count = SparseInt(0);
    int it_cnt = 0;
    while (true) {
        it_cnt++;
        bool invalid = false;
        char* assignment = (char *) malloc(sizeof(char) * dd_A->size);
        if (Cudd_bddPickOneCube(dd_A, constraints, assignment) == 0) {
            break;
        }
        DdNode* curr_transition = dd_graph_closure.transition; Cudd_Ref(curr_transition);
        for (int i = 0; i < dd_A->size; i++) {
            assignment[i] = (assignment[i] == 1);
        }

        // Add the "different from previous candidates" constraint
        DdNode * accumulated_patches = Cudd_Not(dd_graph_closure.dd->one); Cudd_Ref(accumulated_patches);
        for (int i = 0; i < candidate_list_dd_projected_graph.size(); i++) {
            char * assignment_i = assignment + i * p.e_vars[0].second.size();
            std::string z_str = "";
            for (int j = 0; j < p.e_vars[0].second.size(); j++) {
                z_str += assignment_i[j] == 1 ? "1" : "0";
            }
            bool val = !evaluate_candidate_dd_projected_graph(candidate_list_dd_projected_graph[i], assignment_i);
            if (dd_graph_closure.forced_val.find(z_str) != dd_graph_closure.forced_val.end()) {
                if (val == dd_graph_closure.forced_val[z_str]) {
                    continue;
                } else {
                    invalid = true;
                    break;
                }
            }
            DdNode * patch_cube = force_assignment_cube(assignment_i, val); Cudd_Ref(patch_cube);
            accumulated_patches = Cudd_bddOr(dd_graph_closure.dd, accumulated_patches, patch_cube); Cudd_Ref(accumulated_patches);
            Cudd_RecursiveDeref(dd_graph_closure.dd, patch_cube);
        }

        // Check satisfiability under the constraint
        if (!invalid) {
            if (candidate_list_dd_projected_graph.size() != 0) {
                curr_transition = Cudd_bddOr(dd_graph_closure.dd, curr_transition, accumulated_patches); Cudd_Ref(curr_transition);
            }
            int result = reachTransRelation(dd_graph_closure.dd, dd_graph_closure.base_initial, curr_transition, Cudd_Not(dd_graph_closure.dd->one), pPars, dd_graph_closure.reached, dd_graph_closure.register_size);
            Cudd_AutodynEnable(dd_graph_closure.dd, CUDD_REORDER_SYMM_SIFT);
            if (result == -1) {
                print_info("2DQR UNKNOWN (Valid)");
                return SparseInt();
            } else if (check_prop() == 0) {
                invalid = true;
            }
        }

        // Block if the constraint is unsatisfiable
        if (invalid) {
            DdNode* blocking_cube = dd_A->one; Cudd_Ref(blocking_cube);
            for (int i = 0; i < dd_A->size; i++) {
                DdNode* src = Cudd_bddIthVar(dd_A, i); Cudd_Ref(src);
                blocking_cube = Cudd_bddAnd(dd_A, blocking_cube, assignment[i] == 1 ? src : Cudd_Not(src)); Cudd_Ref(blocking_cube);
                Cudd_RecursiveDeref(dd_A, src);
            }
            constraints = Cudd_bddAnd(dd_A, constraints, Cudd_Not(blocking_cube)); Cudd_Ref(constraints);
            Cudd_RecursiveDeref(dd_A, blocking_cube);
            Cudd_RecursiveDeref(dd_graph_closure.dd, accumulated_patches);

            dd_graph_closure.forced_val.clear();
            free(assignment);
            continue;
        }
        free(assignment);
        
        // Generate a valid candidate
        DdNode* candidate = get_valid_candidate(curr_transition); Cudd_Ref(candidate);
        // Count the number of Skolem functions with f1 = candidate
        solution_count += count_1DQBF(candidate, component_y2);
        Cudd_RecursiveDeref(dd_graph_closure.dd, curr_transition);

        // Add the candidate to the list
        DdNode* candidate_dd_projected_graph = Extra_TransferPermute(dd_phi.dd, dd_projected_graph.dd, candidate, dd_phi_to_dd_projected_graph); Cudd_Ref(candidate_dd_projected_graph);
        candidate_list_dd_projected_graph.push_back(candidate_dd_projected_graph);
        print_info(fmt::format("it: {}, |Cs|: {}, #: {}", it_cnt, candidate_list_dd_projected_graph.size(), solution_count.to_string()).c_str());
        Cudd_RecursiveDeref(dd_graph_closure.dd, candidate);
        
        // Add consistency constraints
        DdNode* src_y_k = Cudd_bddIthVar(dd_projected_graph.dd, p.e_vars[0].second.size()); Cudd_Ref(src_y_k);
        DdNode* dst_y_k = Cudd_bddIthVar(dd_projected_graph.dd, p.e_vars[0].second.size() * 2 + 1); Cudd_Ref(dst_y_k);
        int idx = candidate_list_dd_projected_graph.size() - 1;
        
        int * dd_projected_graph_to_dd_A = (int*) calloc(dd_projected_graph.dd->size, sizeof(int));
        for (int j = 0; j < p.e_vars[0].second.size(); j++) {
            Cudd_bddNewVar(dd_A);
            dd_projected_graph_to_dd_A[j] = j + p.e_vars[0].second.size() * idx;
        }
        DdNode* component_y1_dd_A_shifted;
        if (component_y1_dd_A == nullptr) {
            component_y1_dd_A = Extra_TransferPermute(dd_phi.dd, dd_A, component_y1, dd_phi_to_dd_A); Cudd_Ref(component_y1_dd_A);
            component_y1_dd_A_shifted = component_y1_dd_A; Cudd_Ref(component_y1_dd_A_shifted);
        } else {
            component_y1_dd_A_shifted = Cudd_bddSwapVariables(dd_A, component_y1_dd_A, dd_A->vars, dd_A->vars + (p.e_vars[0].second.size()) * idx, p.e_vars[0].second.size()); Cudd_Ref(component_y1_dd_A_shifted);
        }
        constraints = Cudd_bddAnd(dd_A, constraints, component_y1_dd_A_shifted); Cudd_Ref(constraints);
        Cudd_RecursiveDeref(dd_A, component_y1_dd_A_shifted);

        // If X_a^b -> X_{a'}^{b'} where b = ~f(a) and b' = f'(a'), then if we choose to differ from f at a, we must not differ from f' at a'
        // In otherwords, if we choose to differ from f at a and to differ from f' at a', then we must not have X_a^b -> X_{a'}^{b'} where b = ~f(a) and b' = f'(a')
        // f : candidate_list_dd_projected_graph[idx]

        DdNode* constraint = Cudd_bddSubstitude(dd_projected_graph.dd, dd_projected_graph.projected_graph, src_y_k, Cudd_Not(candidate_list_dd_projected_graph[idx])); Cudd_Ref(constraint);
        for (int i = 0; i < idx; i++) {
            // f' : candidate_list_dd_projected_graph[i]
            DdNode* f_p = Cudd_bddSwapVariables(dd_projected_graph.dd, candidate_list_dd_projected_graph[i], dd_projected_graph.dd->vars, dd_projected_graph.dd->vars + p.e_vars[0].second.size() + 1, p.e_vars[0].second.size()); Cudd_Ref(f_p);
            DdNode* constraint_p = Cudd_bddSubstitude(dd_projected_graph.dd, constraint, dst_y_k, f_p); Cudd_Ref(constraint_p);

            for (int j = 0; j < p.e_vars[0].second.size(); j++) {
                dd_projected_graph_to_dd_A[p.e_vars[0].second.size() + 1 + j] = j + p.e_vars[0].second.size() * i;
            }
            DdNode* constraint_p_ = Extra_TransferPermute(dd_projected_graph.dd, dd_A, constraint_p, dd_projected_graph_to_dd_A); Cudd_Ref(constraint_p_);
            constraints = Cudd_bddAnd(dd_A, constraints, Cudd_Not(constraint_p_)); Cudd_Ref(constraints);
            Cudd_RecursiveDeref(dd_projected_graph.dd, f_p);
            Cudd_RecursiveDeref(dd_projected_graph.dd, constraint_p);
            Cudd_RecursiveDeref(dd_A, constraint_p_);
        }
        Cudd_RecursiveDeref(dd_projected_graph.dd, src_y_k);
        Cudd_RecursiveDeref(dd_projected_graph.dd, dst_y_k);
        Cudd_RecursiveDeref(dd_projected_graph.dd, constraint);

        free(dd_projected_graph_to_dd_A);
        dd_graph_closure.forced_val.clear();
    }
    Cudd_RecursiveDeref(dd_A, constraints);
    for (int i = 0; i < candidate_list_dd_projected_graph.size(); i++) {
        Cudd_RecursiveDeref(dd_projected_graph.dd, candidate_list_dd_projected_graph[i]);
    }
    Cudd_RecursiveDeref(dd_A, component_y1_dd_A);
    Cudd_Quit(dd_A);

    return solution_count;
}