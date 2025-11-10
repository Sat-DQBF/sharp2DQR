#include "algo/algo.hpp"
#include "common/utils.hpp"

#include <algorithm>
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

void Algo::generate_graph_closure_transition_system() {
    dd_graph_closure.verilog.name = "graph_closure";
    
    dd_graph_closure.verilog.wires.resize(0);
    int max_dep_size = std::max(p.e_vars[0].second.size(), p.e_vars[1].second.size());
    dd_graph_closure.register_size = 2 + p.u_vars.size() + 1 + max_dep_size;
    for (int i = 0; i < dd_graph_closure.register_size; i++) {
        dd_graph_closure.verilog.PI.push_back("R" + std::to_string(i));
    }
    for (int i = 0; i < dd_graph_closure.register_size; i++) {
        dd_graph_closure.verilog.PI.push_back("Rn" + std::to_string(i));
    }
    dd_graph_closure.verilog.PO.push_back("I");
    dd_graph_closure.verilog.PO.push_back("T");

    // Register
    // -------------------------------------------------
    // | k | y_k |     x     | target y_k | target z_k |
    // -------------------------------------------------
    
    // Map register names to their index
    dd_graph_closure.idx_lookup["k"] = 0;
    dd_graph_closure.idx_lookup["y_k"] = 1;
    dd_graph_closure.idx_lookup["x"] = 2;
    for (int i = 0; i < p.u_vars.size(); i++) {
        dd_graph_closure.idx_lookup[p.u_vars[i]] = 2 + i;
    }
    dd_graph_closure.idx_lookup["target y_k"] = 2 + p.u_vars.size();
    dd_graph_closure.idx_lookup["target z_k"] = 2 + p.u_vars.size() + 1;

    // Precalculations
    std::set<std::string> x(p.u_vars.begin(), p.u_vars.end());
    std::set<std::string> z0(p.e_vars[0].second.begin(), p.e_vars[0].second.end());
    std::set<std::string> z1(p.e_vars[1].second.begin(), p.e_vars[1].second.end());
    std::set<std::string> z0_intersect_z1;
    std::set_intersection(z0.begin(), z0.end(), z1.begin(), z1.end(), std::inserter(z0_intersect_z1, z0_intersect_z1.end()));
    std::set<std::string> z0_minus_z1;
    std::set_difference(z0.begin(), z0.end(), z0_intersect_z1.begin(), z0_intersect_z1.end(), std::inserter(z0_minus_z1, z0_minus_z1.end()));
    std::set<std::string> z1_minus_z0;
    std::set_difference(z1.begin(), z1.end(), z0_intersect_z1.begin(), z0_intersect_z1.end(), std::inserter(z1_minus_z0, z1_minus_z0.end()));

    VerilogModule::submodule dest_0_1 = {.module = p.phi, .name = "phi_0_1", .inputs = std::vector<Circuit>()}; // y_k = 0, y_k' = 1
    VerilogModule::submodule dest_1_0 = {.module = p.phi, .name = "phi_1_0", .inputs = std::vector<Circuit>()}; // y_k = 1, y_k' = 0

    for (auto u : p.u_vars) {
        if (z0_minus_z1.find(u) != z0_minus_z1.end()) {
            dest_0_1.inputs.push_back(bv_at("R", dd_graph_closure.idx_lookup[u]));
        } else {
            dest_0_1.inputs.push_back(bv_at("Rn", dd_graph_closure.idx_lookup[u]));
        }
        if (z1_minus_z0.find(u) != z1_minus_z0.end()) {
            dest_1_0.inputs.push_back(bv_at("R", dd_graph_closure.idx_lookup[u]));
        } else {
            dest_1_0.inputs.push_back(bv_at("Rn", dd_graph_closure.idx_lookup[u]));
        }
    }
    dd_graph_closure.verilog.wires.push_back(Wire("neg_Rn_y_k"));
    dd_graph_closure.verilog.wires.back().assign = ~bv_at("Rn", dd_graph_closure.idx_lookup["y_k"]);
    dest_0_1.inputs.push_back(bv_at("R", dd_graph_closure.idx_lookup["y_k"]));
    dest_0_1.inputs.push_back(Circuit("neg_Rn_y_k"));
    dest_1_0.inputs.push_back(Circuit("neg_Rn_y_k"));
    dest_1_0.inputs.push_back(bv_at("R", dd_graph_closure.idx_lookup["y_k"]));

    dd_graph_closure.verilog.submodules.push_back(dest_0_1);
    dd_graph_closure.verilog.submodules.push_back(dest_1_0);

    // Initial:
    // (k, target k, z_k) and (y_k, target y_k, target z_k) are the same
    {
        Wire tmp("I");
        tmp.assign = Circuit("");
        tmp.assign.operation = "&";
        tmp.assign.children.push_back(~bv_at("R", dd_graph_closure.idx_lookup["k"])); // k' = 0
        tmp.assign.children.push_back(r_eq(dd_graph_closure.idx_lookup["y_k"], dd_graph_closure.idx_lookup["target y_k"])); // y_k' = target y_k'
        // target z_k' = z_k' \subseteq x'
        {
            if (p.e_vars[0].second.size() > 0) {
                int i = 0;
                for (auto& v : p.e_vars[0].second) {
                    tmp.assign.children.push_back(r_eq(dd_graph_closure.idx_lookup["target z_k"] + i, dd_graph_closure.idx_lookup[v]));
                    i++;
                }
            }
        }
        dd_graph_closure.verilog.wires.push_back(tmp);
    }

    // Transition

    std::vector<Wire> transition_vector;
    // Change unrelated bits
    {
        Wire tmp("change_repr");
        tmp.assign = Circuit("");
        tmp.assign.operation = "&";
        tmp.assign.children.push_back(r_eq_rn(dd_graph_closure.idx_lookup["y_k"], dd_graph_closure.idx_lookup["k"])); 
        tmp.assign.children.push_back(r_eq_rn(dd_graph_closure.register_size - 1, dd_graph_closure.idx_lookup["target y_k"]));
        {
            if (p.e_vars[0].second.size() > 0) {
                Circuit tmp_2("");
                tmp_2.operation = "&";
                for (auto& v : p.e_vars[0].second) {
                    tmp_2.children.push_back(r_eq_rn(dd_graph_closure.idx_lookup[v], dd_graph_closure.idx_lookup[v]));
                }
                tmp.assign.children.push_back((~bv_at("R", dd_graph_closure.idx_lookup["k"])).implies(tmp_2));
            }
        }
        {
            if (p.e_vars[1].second.size() > 0) {
                Circuit tmp_2("");
                tmp_2.operation = "&";
                for (auto& v : p.e_vars[1].second) {
                    tmp_2.children.push_back(r_eq_rn(dd_graph_closure.idx_lookup[v], dd_graph_closure.idx_lookup[v]));
                }
                tmp.assign.children.push_back((bv_at("R", dd_graph_closure.idx_lookup["k"])).implies(tmp_2));
            }
        }
        transition_vector.push_back(tmp);
    }

    // Stepping
    {
        Wire tmp("step");
        tmp.assign = Circuit("");
        tmp.assign.operation = "&";
        tmp.assign.children.push_back(bv_at("R", dd_graph_closure.idx_lookup["k"]) == (~bv_at("Rn", dd_graph_closure.idx_lookup["k"])));  // k = ~k'
        tmp.assign.children.push_back(r_eq_rn(dd_graph_closure.register_size - 1, dd_graph_closure.idx_lookup["target y_k"]));
        
        // Intersection is unchanged
        for (auto& v : z0_intersect_z1) {
            tmp.assign.children.push_back(r_eq_rn(dd_graph_closure.idx_lookup[v], dd_graph_closure.idx_lookup[v]));
        }
        tmp.assign.children.push_back((~bv_at("R", dd_graph_closure.idx_lookup["k"])).implies(~Circuit("phi_0_1_" + p.phi.PO[0])));
        tmp.assign.children.push_back((bv_at("R", dd_graph_closure.idx_lookup["k"])).implies(~Circuit("phi_1_0_" + p.phi.PO[0])));
        transition_vector.push_back(tmp);
    }

    Wire transition = Wire("T");
    transition.assign = Circuit("");
    transition.assign.operation = "|";
    for (auto& v : transition_vector) {
        transition.assign.children.push_back(Circuit(v.name));
        dd_graph_closure.verilog.wires.push_back(v);
    }
    dd_graph_closure.verilog.wires.push_back(transition);
}

int Algo::setup_dd_graph_closure() {
    dd_graph_closure.dd = Cudd_Init(dd_graph_closure.register_size * 2, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    DdManager * dd = dd_graph_closure.dd;

    int * perm = (int *) malloc(dd_graph_closure.register_size * 2 * sizeof(int));
    for (int i = 0; i < dd_graph_closure.register_size; i++) {
        perm[2 * i] = i;
        perm[2 * i + 1] = i + dd_graph_closure.register_size;
    }

    Cudd_ShuffleHeap(dd, perm);

    std::ofstream out("./graph_closure.v");
    out << dd_graph_closure.verilog.to_string();
    out.close();
    if (abc.read("./graph_closure.v")) {
        print_error("Failed to read file");
    }
    Abc_Ntk_t * pNtk = Abc_FrameReadNtk(abc.pAbc);
    Aig_Man_t * pMan = Abc_NtkToDar(pNtk, 0, 1);
    ProgressBar * pProgress = NULL;
    Aig_Obj_t * pObj;
    DdNode * bFunc;
    int i, Counter;

    if (args.reorder) {
        Cudd_AutodynEnable(dd, CUDD_REORDER_SYMM_SIFT);
    }
    Aig_ManCleanData(pMan);
    
    (Aig_ManConst1(pMan))->pData = dd->one; Cudd_Ref(dd->one);
    Aig_ManForEachCi(pMan, pObj, i) {
        pObj->pData = dd->vars[i]; Cudd_Ref(dd->vars[i]);
    }
    Counter = 0;

    Aig_ManForEachCo( pMan, pObj, i )
    {
        bFunc = Bbr_NodeGlobalBdds_rec( dd, Aig_ObjFanin0(pObj), pPars->nBddMax, 1, pProgress, &Counter, abc.verbose );
        if ( bFunc == NULL )
        {
            if ( abc.verbose )
            printf( "Constructing global BDDs is aborted.\n" );
            Aig_ManFreeGlobalBdds( pMan, dd );
            Cudd_Quit( dd ); 
            // reset references
            Aig_ManResetRefs( pMan );
            return -1;
        }
        bFunc = Cudd_NotCond( bFunc, Aig_ObjFaninC0(pObj) );  Cudd_Ref( bFunc ); 
        pObj->pData = bFunc;
    }
    Aig_ManResetRefs(pMan);
    if (args.reorder) {
        Cudd_ReduceHeap(dd, CUDD_REORDER_SYMM_SIFT, 1);
    }
    // Cudd_AutodynDisable(dd);

    // Modified from Aig_ManInitStateVarMap @ bbrReach.c:124
    DdNode ** pbVarsX, ** pbVarsY;
    
    // set the variable mapping for Cudd_bddVarMap()
    pbVarsX = ABC_ALLOC(DdNode *, dd_graph_closure.register_size);
    pbVarsY = ABC_ALLOC(DdNode *, dd_graph_closure.register_size);
    for ( i = 0; i < dd_graph_closure.register_size ; i++ )
    {
        pbVarsX[i] = dd->vars[i];
        pbVarsY[i] = dd->vars[dd_graph_closure.register_size + i];
    }
    Cudd_SetVarMap(dd, pbVarsX, pbVarsY, dd_graph_closure.register_size);
    ABC_FREE(pbVarsX);
    ABC_FREE(pbVarsY);

    dd_graph_closure.initial = Aig_ObjGlobalBdd((Aig_Obj_t *) Vec_PtrEntry(pMan->vCos, 0)); Cudd_Ref(dd_graph_closure.initial);
    dd_graph_closure.transition = Aig_ObjGlobalBdd((Aig_Obj_t *) Vec_PtrEntry(pMan->vCos, 1)); Cudd_Ref(dd_graph_closure.transition);
    return 0;
}

int Algo::check_prop() {
    DdNode* bS = dd_graph_closure.reached;
    DdManager* dd = dd_graph_closure.dd;

    DdNode* exist_abs_cube = dd->one; Cudd_Ref(exist_abs_cube);
    auto assign = [&](int i, bool val) {
        // https://stackoverflow.com/questions/55246590/cudd-manipulation-of-bdds
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
        bS = Cudd_bddAnd(dd, bS, val ? src : Cudd_Not(src)); Cudd_Ref(bS);
        exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, src); Cudd_Ref(exist_abs_cube);
        Cudd_RecursiveDeref(dd, src);
    };

    auto eq = [&](int i, int j, bool neg) {
        DdNode* l = Cudd_bddIthVar(dd, i); Cudd_Ref(l);
        DdNode* r = Cudd_bddIthVar(dd, j); Cudd_Ref(r);
        DdNode* eq = Cudd_bddXnor(dd, l, neg ? Cudd_Not(r) : r); Cudd_Ref(eq);
        bS = Cudd_bddAnd(dd, bS, eq); Cudd_Ref(bS);
        exist_abs_cube = Cudd_bddAnd(dd, exist_abs_cube, l); Cudd_Ref(exist_abs_cube);
        Cudd_RecursiveDeref(dd, l);
        Cudd_RecursiveDeref(dd, r);
        Cudd_RecursiveDeref(dd, eq);
    };

    std::set<std::string> z1(p.e_vars[0].second.begin(), p.e_vars[0].second.end());

    assign(dd_graph_closure.idx_lookup["k"], 0); // k = 0
    eq(dd_graph_closure.idx_lookup["y_k"], dd_graph_closure.idx_lookup["target y_k"], 1); // y_k = ~target y_k

    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z1.find(p.u_vars[i]) == z1.end()) {
            assign(dd_graph_closure.idx_lookup["x"] + i, 0); // u_i = 0
        }
    }
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        eq(dd_graph_closure.idx_lookup[p.e_vars[0].second[i]], dd_graph_closure.idx_lookup["target z_k"] + i, 0);
    }
    bS = Cudd_bddExistAbstract(dd, bS, exist_abs_cube); Cudd_Ref(bS);

    DdNode* target_y_k = Cudd_bddIthVar(dd, dd_graph_closure.idx_lookup["target y_k"]); Cudd_Ref(target_y_k);
    DdNode* bS_0 = Cudd_bddAnd(dd, bS, Cudd_Not(target_y_k)); Cudd_Ref(bS_0);
    bS_0 = Cudd_bddExistAbstract(dd, bS_0, target_y_k); Cudd_Ref(bS_0);
    DdNode* bS_1 = Cudd_bddAnd(dd, bS, target_y_k); Cudd_Ref(bS_1);
    bS_1 = Cudd_bddExistAbstract(dd, bS_1, target_y_k); Cudd_Ref(bS_1);

    DdNode* bCheck = Cudd_bddAnd(dd, bS_0, bS_1); Cudd_Ref(bCheck);
    return Cudd_bddLeq(dd, bCheck, Cudd_Not(dd->one));
}
