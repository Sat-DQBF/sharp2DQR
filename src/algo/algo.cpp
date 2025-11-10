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

void write_to_log(std::string content, std::string filename = "./result.txt", bool append = false) {
    std::ofstream out;
    if (append) {
        out.open(filename, std::ios_base::app);
    } else {
        out.open(filename);
    }
    out << content;
    out.close();
}


void Algo::run() {
    Bbr_ManSetDefaultParams(pPars);
    pPars->fTransRel = 1;
    pPars->fVerbose = abc.verbose;
    // pPars->nBddMax = 1000000; // Default: 50000
    pPars->nBddMax = INT_MAX / 4; // Default: 50000
    pPars->fSilent = abc.verbose ? 0 : 1;
    pPars->fReorder = args.reorder ? 0 : 1;

    if (setup_dd_phi() == -1) {
        print_info("2DQR UNKNOWN");
        return;
    };

    // Compute the set of dont care vars
    get_dont_care_vars();
    unsigned long long dont_care_count_y1 = (unsigned long long) Cudd_CountMinterm(dd_phi.dd, dd_phi.bDontCareVars_y1, p.e_vars[0].second.size());
    unsigned long long dont_care_count_y2 = (unsigned long long) Cudd_CountMinterm(dd_phi.dd, dd_phi.bDontCareVars_y2, p.e_vars[1].second.size());
    unsigned long long dont_care_count = dont_care_count_y1 + dont_care_count_y2;
    print_info((fmt::format("Number of dont care vars: {} + {}", dont_care_count_y1, dont_care_count_y2)).c_str());
    
    // Compute the set of free variables
    generate_graph_closure_transition_system();
    print_info("Generating BDD for Transition Relation");
    if (setup_dd_graph_closure() == -1) {
        print_info("2DQR UNKNOWN");
        return;
    };
    print_info("Finished BDD for Transition Relation");

    Cudd_AutodynEnable(dd_graph_closure.dd, CUDD_REORDER_SYMM_SIFT);
    Cudd_AutodynEnable(dd_phi.dd, CUDD_REORDER_SYMM_SIFT);
    // force_dont_care_vars();
    
    // dd_graph_closure_to_dd_phi
    dd_graph_closure_to_dd_phi = (int *) calloc(dd_graph_closure.dd->size, sizeof(int));
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        dd_graph_closure_to_dd_phi[dd_graph_closure.idx_lookup["target z_k"] + i] = u_var_idx_to_dd_phi_idx[dd_graph_closure.idx_lookup[p.e_vars[0].second[i]] - dd_graph_closure.idx_lookup["x"]];
    }

    // dd_phi_to_dd_A
    dd_phi_to_dd_A = (int *) calloc(dd_phi.dd->size, sizeof(int));
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        dd_phi_to_dd_A[u_var_idx_to_dd_phi_idx[dd_graph_closure.idx_lookup[p.e_vars[0].second[i]] - dd_graph_closure.idx_lookup["x"]]] = i;
    }

    // dd_phi_to_dd_projected_graph
    dd_phi_to_dd_projected_graph = (int *) calloc(dd_phi.dd->size, sizeof(int));
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        dd_phi_to_dd_projected_graph[u_var_idx_to_dd_phi_idx[dd_graph_closure.idx_lookup[p.e_vars[0].second[i]] - dd_graph_closure.idx_lookup["x"]]] = i;
    }

    // Check satisfiability
    if (reachTransRelation(dd_graph_closure.dd, dd_graph_closure.initial, dd_graph_closure.transition, Cudd_Not(dd_graph_closure.dd->one), pPars, dd_graph_closure.reached, dd_graph_closure.register_size) == -1) {
        print_info("2DQR UNSAT");
        return;
    };
    Cudd_AutodynEnable(dd_graph_closure.dd, CUDD_REORDER_SYMM_SIFT);

    if (!check_prop()) {
        print_info("2DQR UNKNOWN");
        return;
    }
    get_free_vars();
    dd_phi.bFreeVars_y1 = Cudd_bddAnd(dd_phi.dd, dd_phi.bFreeVars_y1, Cudd_Not(dd_phi.bDontCareVars_y1)); Cudd_Ref(dd_phi.bFreeVars_y1);
    print_info(("Number of free vars: " + std::to_string((unsigned long long) Cudd_CountMinterm(dd_phi.dd, dd_phi.bFreeVars_y1, p.e_vars[0].second.size()))).c_str());

    // Setup for constraints solver
    setup_dd_projected_graph(dd_graph_closure.reached);

    // Get a better initial
    dd_graph_closure.base_initial = dd_graph_closure.reached; Cudd_Ref(dd_graph_closure.base_initial);

    component_cnt = 0;
    lower_bound = SparseInt(1);
    
    // Compute the number of components
    generate_component_transition_system();
    print_info("Generating BDD for Component Transition Relation");
    if (setup_dd_component() == -1) {
        print_info("2DQR UNKNOWN");
        return;
    };
    print_info("Finished BDD for Component Transition Relation");

    // dd_component to dd_phi
    dd_component_to_dd_phi = (int *) calloc(dd_component.dd->size, sizeof(int));
    for (int i = 0; i < p.u_vars.size(); i++) {
        dd_component_to_dd_phi[dd_component.idx_lookup["x"] + i] = i;
    }
    
    
    char * assignment = (char *) malloc(dd_phi.dd->size * sizeof(char));
    DdNode * bRemaining;

    bRemaining = Cudd_Not(dd_phi.bDontCareVars_y1); Cudd_Ref(bRemaining);
    DdNode * bForcedVars = Cudd_bddAnd(dd_phi.dd, Cudd_Not(dd_phi.bDontCareVars_y1), Cudd_Not(dd_phi.bFreeVars_y1)); Cudd_Ref(bForcedVars);

    if (!Cudd_bddLeq(dd_phi.dd, bForcedVars, Cudd_Not(dd_phi.dd->one))) {
        int * dd_phi_to_dd_component = (int *) calloc(dd_phi.dd->size, sizeof(int));
        for (int i = 0; i < p.u_vars.size(); i++) {
            dd_phi_to_dd_component[i] = dd_component.idx_lookup["x"] + i;
        }
        DdNode* bInitial = Extra_TransferPermute(dd_phi.dd, dd_component.dd, bForcedVars, dd_phi_to_dd_component); Cudd_Ref(bInitial);
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["k"]); Cudd_Ref(src);
            bInitial = Cudd_bddAnd(dd_component.dd, bInitial, Cudd_Not(src)); Cudd_Ref(bInitial);
        }

        // Get the component reachable from the initial state
        DdNode * bReached, * component_y1, * component_y2;
        int result = reachTransRelation(dd_component.dd, bInitial, dd_component.transition, Cudd_Not(dd_component.dd->one), pPars, bReached, dd_component.register_size);
        Cudd_AutodynEnable(dd_component.dd, CUDD_REORDER_SYMM_SIFT);

        // Remove component from remaining
        DdNode * exist_abs_cube = dd_component.dd->one; Cudd_Ref(exist_abs_cube);
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["y_k"]); Cudd_Ref(src);
            exist_abs_cube = Cudd_bddAnd(dd_component.dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
        }
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["k"]); Cudd_Ref(src);
            exist_abs_cube = Cudd_bddAnd(dd_component.dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
            component_y1 = Cudd_bddAnd(dd_component.dd, bReached, Cudd_Not(src)); Cudd_Ref(component_y1);
            component_y2 = Cudd_bddAnd(dd_component.dd, bReached, src); Cudd_Ref(component_y2);
        }

        component_y1 = Cudd_bddExistAbstract(dd_component.dd, component_y1, exist_abs_cube); Cudd_Ref(component_y1);
        component_y2 = Cudd_bddExistAbstract(dd_component.dd, component_y2, exist_abs_cube); Cudd_Ref(component_y2);
        Cudd_RecursiveDeref(dd_component.dd, exist_abs_cube);

        DdNode * component_y1_dd_phi = Extra_TransferPermute(dd_component.dd, dd_phi.dd, component_y1, dd_component_to_dd_phi); Cudd_Ref(component_y1_dd_phi);
        DdNode * component_y2_dd_phi = Extra_TransferPermute(dd_component.dd, dd_phi.dd, component_y2, dd_component_to_dd_phi); Cudd_Ref(component_y2_dd_phi);
        bRemaining = Cudd_bddAnd(dd_phi.dd, bRemaining, Cudd_Not(component_y1_dd_phi)); Cudd_Ref(bRemaining);
        
        // Enumerate skolem on component
        component_y1_dd_phi = Cudd_bddAnd(dd_phi.dd, component_y1_dd_phi, dd_phi.bFreeVars_y1); Cudd_Ref(component_y1_dd_phi);
        print_info(fmt::format("Component size: {} {}", (long long) Cudd_CountMinterm(dd_phi.dd, component_y1_dd_phi, p.e_vars[0].second.size()), (long long) Cudd_CountMinterm(dd_phi.dd, component_y2_dd_phi, p.e_vars[1].second.size())).c_str());
        SparseInt component_subcount = enumerate_skolem_on_component(component_y1_dd_phi, component_y2_dd_phi);
        lower_bound = lower_bound * component_subcount;
        write_to_log("SAT LB_" + lower_bound.to_string() + " << " + std::to_string(dont_care_count) + "_" + std::to_string(component_cnt));
        print_info(fmt::format("Component {}: {}", component_cnt, component_subcount.to_string()).c_str());

        Cudd_RecursiveDeref(dd_component.dd, bReached);
        Cudd_RecursiveDeref(dd_component.dd, bInitial);
        Cudd_RecursiveDeref(dd_component.dd, component_y1);
        Cudd_RecursiveDeref(dd_phi.dd, component_y1_dd_phi);
        Cudd_RecursiveDeref(dd_component.dd, component_y2);
        Cudd_RecursiveDeref(dd_phi.dd, component_y2_dd_phi);
    }

    while (Cudd_bddPickOneCube(dd_phi.dd, bRemaining, assignment)) {
        component_cnt++;
        // Set initial state to a unexplored vertex
        DdNode* bInitial = dd_component.dd->one; Cudd_Ref(bInitial);
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["k"]); Cudd_Ref(src);
            bInitial = Cudd_bddAnd(dd_component.dd, bInitial, Cudd_Not(src)); Cudd_Ref(bInitial);
        }
        for (int i = 0; i < p.e_vars[0].second.size(); i++) {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup[p.e_vars[0].second[i]]); Cudd_Ref(src);
            bInitial = Cudd_bddAnd(dd_component.dd, bInitial, assignment[dd_component.idx_lookup[p.e_vars[0].second[i]] - dd_component.idx_lookup["x"]] ? src : Cudd_Not(src)); Cudd_Ref(bInitial);
        }
        
        // Get the component reachable from the initial state
        DdNode * bReached, * component_y1, * component_y2;
        int result = reachTransRelation(dd_component.dd, bInitial, dd_component.transition, Cudd_Not(dd_component.dd->one), pPars, bReached, dd_component.register_size);
        Cudd_AutodynEnable(dd_component.dd, CUDD_REORDER_SYMM_SIFT);
        
        // Remove component from remaining
        DdNode * exist_abs_cube = dd_component.dd->one; Cudd_Ref(exist_abs_cube);
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["y_k"]); Cudd_Ref(src);
            exist_abs_cube = Cudd_bddAnd(dd_component.dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
        }
        {
            DdNode * src = Cudd_bddIthVar(dd_component.dd, dd_component.idx_lookup["k"]); Cudd_Ref(src);
            exist_abs_cube = Cudd_bddAnd(dd_component.dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
            component_y1 = Cudd_bddAnd(dd_component.dd, bReached, Cudd_Not(src)); Cudd_Ref(component_y1);
            component_y2 = Cudd_bddAnd(dd_component.dd, bReached, src); Cudd_Ref(component_y2);
        }

        component_y1 = Cudd_bddExistAbstract(dd_component.dd, component_y1, exist_abs_cube); Cudd_Ref(component_y1);
        component_y2 = Cudd_bddExistAbstract(dd_component.dd, component_y2, exist_abs_cube); Cudd_Ref(component_y2);
        Cudd_RecursiveDeref(dd_component.dd, exist_abs_cube);

        DdNode * component_y1_dd_phi = Extra_TransferPermute(dd_component.dd, dd_phi.dd, component_y1, dd_component_to_dd_phi); Cudd_Ref(component_y1_dd_phi);
        DdNode * component_y2_dd_phi = Extra_TransferPermute(dd_component.dd, dd_phi.dd, component_y2, dd_component_to_dd_phi); Cudd_Ref(component_y2_dd_phi);

        bRemaining = Cudd_bddAnd(dd_phi.dd, bRemaining, Cudd_Not(component_y1_dd_phi)); Cudd_Ref(bRemaining);
        
        // Enumerate skolem on component
        component_y1_dd_phi = Cudd_bddAnd(dd_phi.dd, component_y1_dd_phi, dd_phi.bFreeVars_y1); Cudd_Ref(component_y1_dd_phi);
        print_info(fmt::format("Component size: {} {}", (long long) Cudd_CountMinterm(dd_phi.dd, component_y1_dd_phi, p.e_vars[0].second.size()), (long long) Cudd_CountMinterm(dd_phi.dd, component_y2_dd_phi, p.e_vars[1].second.size())).c_str());
        SparseInt component_subcount = enumerate_skolem_on_component(component_y1_dd_phi, component_y2_dd_phi);
        lower_bound = lower_bound * component_subcount;
        write_to_log("SAT LB_" + lower_bound.to_string() + " << " + std::to_string(dont_care_count) + "_" + std::to_string(component_cnt));
        print_info(fmt::format("Component {}: {}", component_cnt, component_subcount.to_string()).c_str());

        Cudd_RecursiveDeref(dd_component.dd, bReached);
        Cudd_RecursiveDeref(dd_component.dd, bInitial);
        Cudd_RecursiveDeref(dd_component.dd, component_y1);
        Cudd_RecursiveDeref(dd_phi.dd, component_y1_dd_phi);
        Cudd_RecursiveDeref(dd_component.dd, component_y2);
        Cudd_RecursiveDeref(dd_phi.dd, component_y2_dd_phi);
    }
    print_info(fmt::format("2DQR: {} component(s)", component_cnt).c_str());
    write_to_log("SAT EXACT_" + lower_bound.to_string() + " << " + std::to_string(dont_care_count) + "_" + std::to_string(component_cnt));
    print_info(fmt::format("Exact bound: {} << {}", lower_bound.to_string(), std::to_string(dont_care_count)).c_str());

    Cudd_RecursiveDeref(dd_phi.dd, bRemaining);
}