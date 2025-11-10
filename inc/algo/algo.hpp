#ifndef TWODQR_HPP
#define TWODQR_HPP

#include "abc_interface/abc_interface.hpp"
#include "common/dqcir.hpp"
#include "common/sparse_int.hpp"

#include <unordered_set>

class Algo {
public:
    dqcir& p;

    Saig_ParBbr_t Pars, * pPars = &Pars;
    ABC abc;
    
    struct args {
        bool reorder;
    };
    
    args args;
    Algo(dqcir& instance) : p(instance) {};
    
    Circuit bv_at(std::string s, int i);
    Circuit r_eq(int i, int j);
    Circuit r_eq(int low1, int low2, int diff);
    Circuit rn_eq(int i, int j);
    Circuit rn_eq(int low1, int low2, int diff);
    Circuit r_eq_rn(int high, int low);
    
    struct dd_graph_closure {
        DdManager * dd;
        DdNode * initial;
        DdNode * transition;
        DdNode * reached;
        DdNode * base_initial;
        
        VerilogModule verilog;
        size_t register_size;
        std::unordered_map<std::string, int> idx_lookup;

        std::unordered_map<std::string, bool> forced_val;
    };
    dd_graph_closure dd_graph_closure;
    void generate_graph_closure_transition_system();
    int setup_dd_graph_closure();
    int check_prop();
    int get_free_vars();

    DdNode* get_candidate();
    int check_candidate(DdNode * candidate, char * cex);
    int evaluate_candidate(DdNode * candidate, char * assignment);
    DdNode* force_assignment_cube(char * z, bool val);
    DdNode* patch_transition_with_cex(DdNode * curr_transition, DdNode* candidate, char * cex);
    DdNode* get_valid_candidate(DdNode * curr_transition = nullptr);
    
    struct dd_component {
        DdManager * dd;
        DdNode * transition;
        
        VerilogModule verilog;
        size_t register_size;
        std::unordered_map<std::string, int> idx_lookup;
    };
    dd_component dd_component;
    void generate_component_transition_system();
    int setup_dd_component();
    
    int* u_var_idx_to_dd_phi_idx;
    
    int component_cnt;
    SparseInt lower_bound;
    DdNode* base_candidate;
    DdNode* base_initial;
    
    struct dd_phi {
        DdManager * dd; // Vars: [x, y_1, y_2]
        DdNode * bPhi;
        DdNode * bNPhi;
        DdNode * bNPhi_x_0_y2;
        DdNode * bNPhi_x_1_y2;
        DdNode * bFreeVars_y1;
        DdNode * bDontCareVars_y1;
        DdNode * bDontCareVars_y2;
    };
    dd_phi dd_phi;
    int setup_dd_phi();
    SparseInt count_1DQBF(DdNode* f1, DdNode* component_y2);
    
    struct dd_projected_graph {
        DdManager * dd; // Vars: [z_1, y_1, z_1', y_1']
        // Path from (1, z_1, y_1) -> (1, z_1', y_1')
        DdNode * projected_graph;
        std::vector<DdNode*> candidate_list;
    };
    dd_projected_graph dd_projected_graph;
    int setup_dd_projected_graph(DdNode* graph_closure);
    int evaluate_candidate_dd_projected_graph(DdNode * candidate, char * assignment);
    

    // Transfering between dds
    int * dd_component_to_dd_phi;
    int * dd_graph_closure_to_dd_phi;
    int * dd_phi_to_dd_A;
    int * dd_phi_to_dd_projected_graph;

    int get_dont_care_vars();
    int force_dont_care_vars();

    SparseInt enumerate_skolem_on_component(DdNode* comopnent_y1, DdNode* component_y2);
    void run();
};

#endif