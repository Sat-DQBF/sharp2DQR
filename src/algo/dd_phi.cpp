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

int Algo::setup_dd_phi() {
    dd_phi.dd = Cudd_Init(p.u_vars.size() + p.e_vars.size(), 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    DdManager * dd = dd_phi.dd;

    // Cudd_MakeTreeNode(dd_phi.dd, 0, p.e_vars[1].second.size(), MTR_DEFAULT);
    // Cudd_MakeTreeNode(dd_phi.dd, p.e_vars[1].second.size(), p.u_vars.size() + p.e_vars.size() - p.e_vars[1].second.size(), MTR_DEFAULT);

    std::ofstream out("./phi.v");
    out << p.phi.to_string();
    out.close();
    if (abc.read("./phi.v")) {
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

    std::set<std::string> z_2(p.e_vars[1].second.begin(), p.e_vars[1].second.end());
    u_var_idx_to_dd_phi_idx = (int *) malloc(p.u_vars.size() * sizeof(int));
    // int grp_idx = 0;
    // int ungrp_idx = p.e_vars[1].second.size();
    // for (i = 0; i < p.u_vars.size(); i++) {
    //     pObj = (Aig_Obj_t *) Vec_PtrEntry(pMan->vCis, i);
    //     if (z_2.find(p.u_vars[i]) != z_2.end()) {
    //         pObj->pData = dd->vars[grp_idx]; Cudd_Ref(dd->vars[grp_idx]);
    //         u_var_idx_to_dd_phi_idx[i] = grp_idx;
    //         grp_idx++;
    //     } else {
    //         pObj->pData = dd->vars[ungrp_idx]; Cudd_Ref(dd->vars[ungrp_idx]);
    //         u_var_idx_to_dd_phi_idx[i] = ungrp_idx;
    //         ungrp_idx++;
    //     }
    // }
    // for (i = 0; i < p.e_vars.size(); i++) {
    //     pObj = (Aig_Obj_t *) Vec_PtrEntry(pMan->vCis, p.u_vars.size() + i);
    //     pObj->pData = dd->vars[ungrp_idx]; Cudd_Ref(dd->vars[ungrp_idx]);
    //     ungrp_idx++;
    // }
    for (i = 0; i < p.u_vars.size(); i++) {
        pObj = (Aig_Obj_t *) Vec_PtrEntry(pMan->vCis, i);
        pObj->pData = dd->vars[i]; Cudd_Ref(dd->vars[i]);
        u_var_idx_to_dd_phi_idx[i] = i;
    }
    for (i = 0; i < p.e_vars.size(); i++) {
        pObj = (Aig_Obj_t *) Vec_PtrEntry(pMan->vCis, p.u_vars.size() + i);
        pObj->pData = dd->vars[p.u_vars.size() + i]; Cudd_Ref(dd->vars[p.u_vars.size() + i]);
    }

    assert (Saig_ManCoNum(pMan) == 1);
    Counter = 0;
    pObj = (Aig_Obj_t *) Vec_PtrEntry(pMan->vCos, 0);
    dd_phi.bPhi = Bbr_NodeGlobalBdds_rec(dd, Aig_ObjFanin0(pObj), pPars->nBddMax, 1, pProgress, &Counter, abc.verbose);

    if ( dd_phi.bPhi == NULL )
    {
        if ( abc.verbose )
        printf( "Constructing global BDDs is aborted.\n" );
        Aig_ManFreeGlobalBdds( pMan, dd );
        Cudd_Quit( dd ); 
        // reset references
        Aig_ManResetRefs( pMan );
        return -1;
    }
    dd_phi.bPhi = Cudd_NotCond( dd_phi.bPhi, Aig_ObjFaninC0(pObj) );  Cudd_Ref( dd_phi.bPhi ); 
    pObj->pData = dd_phi.bPhi;
    Aig_ManResetRefs( pMan );
    
    if (args.reorder) {
        Cudd_ReduceHeap(dd, CUDD_REORDER_SYMM_SIFT, 1);
    }
    // Cudd_AutodynDisable(dd);

    dd_phi.bNPhi = Cudd_Not(dd_phi.bPhi); Cudd_Ref(dd_phi.bNPhi);

    DdNode* y1 = Cudd_bddIthVar(dd, p.u_vars.size()); Cudd_Ref(y1);
    dd_phi.bNPhi_x_0_y2 = Cudd_bddAnd(dd, dd_phi.bNPhi, Cudd_Not(y1)); Cudd_Ref(dd_phi.bNPhi_x_0_y2);
    dd_phi.bNPhi_x_0_y2 = Cudd_bddExistAbstract(dd, dd_phi.bNPhi_x_0_y2, y1); Cudd_Ref(dd_phi.bNPhi_x_0_y2);

    dd_phi.bNPhi_x_1_y2 = Cudd_bddAnd(dd, dd_phi.bNPhi, y1); Cudd_Ref(dd_phi.bNPhi_x_1_y2);
    dd_phi.bNPhi_x_1_y2 = Cudd_bddExistAbstract(dd, dd_phi.bNPhi_x_1_y2, y1); Cudd_Ref(dd_phi.bNPhi_x_1_y2);
    return 0;
}

int Algo::get_dont_care_vars() {
    DdNode* out_deg_nonzero = dd_phi.bNPhi; Cudd_Ref(out_deg_nonzero);
    DdNode* exists_abs_cube_y1 = dd_phi.dd->one; Cudd_Ref(exists_abs_cube_y1);
    DdNode* exists_abs_cube_y2 = dd_phi.dd->one; Cudd_Ref(exists_abs_cube_y2);
    std::set<std::string> z_1(p.e_vars[0].second.begin(), p.e_vars[0].second.end());
    std::set<std::string> z_2(p.e_vars[1].second.begin(), p.e_vars[1].second.end());
    for (int i = 0; i < p.u_vars.size(); i++) {
        if (z_1.find(p.u_vars[i]) == z_1.end()) {
            DdNode* src = Cudd_bddIthVar(dd_phi.dd, u_var_idx_to_dd_phi_idx[i]); Cudd_Ref(src);
            exists_abs_cube_y1 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y1, src); Cudd_Ref(exists_abs_cube_y1);
            Cudd_RecursiveDeref(dd_phi.dd, src);
        }
        if (z_2.find(p.u_vars[i]) == z_2.end()) {
            DdNode* src = Cudd_bddIthVar(dd_phi.dd, u_var_idx_to_dd_phi_idx[i]); Cudd_Ref(src);
            exists_abs_cube_y2 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y2, src); Cudd_Ref(exists_abs_cube_y2);
            Cudd_RecursiveDeref(dd_phi.dd, src);
        }
    }
    {
        DdNode* y1 = Cudd_bddIthVar(dd_phi.dd, p.u_vars.size()); Cudd_Ref(y1);
        exists_abs_cube_y1 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y1, y1); Cudd_Ref(exists_abs_cube_y1);
        exists_abs_cube_y2 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y2, y1); Cudd_Ref(exists_abs_cube_y2);
        Cudd_RecursiveDeref(dd_phi.dd, y1);
    }
    {
        DdNode* y2 = Cudd_bddIthVar(dd_phi.dd, p.u_vars.size() + 1); Cudd_Ref(y2);
        exists_abs_cube_y1 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y1, y2); Cudd_Ref(exists_abs_cube_y1);
        exists_abs_cube_y2 = Cudd_bddAnd(dd_phi.dd, exists_abs_cube_y2, y2); Cudd_Ref(exists_abs_cube_y2);
        Cudd_RecursiveDeref(dd_phi.dd, y2);
    }
    dd_phi.bDontCareVars_y1 = Cudd_Not(Cudd_bddExistAbstract(dd_phi.dd, out_deg_nonzero, exists_abs_cube_y1)); Cudd_Ref(dd_phi.bDontCareVars_y1);
    dd_phi.bDontCareVars_y2 = Cudd_Not(Cudd_bddExistAbstract(dd_phi.dd, out_deg_nonzero, exists_abs_cube_y2)); Cudd_Ref(dd_phi.bDontCareVars_y2);
    Cudd_RecursiveDeref(dd_phi.dd, out_deg_nonzero);
    Cudd_RecursiveDeref(dd_phi.dd, exists_abs_cube_y1);
    Cudd_RecursiveDeref(dd_phi.dd, exists_abs_cube_y2);
    return 0;
}

int Algo::force_dont_care_vars() {
    DdManager* dd = dd_graph_closure.dd;
    int * dd_phi_to_dd_graph_closure = (int *) calloc(dd_phi.dd->size, sizeof(int));
    for (int i = 0; i < p.u_vars.size(); i++) {
        dd_phi_to_dd_graph_closure[u_var_idx_to_dd_phi_idx[i]] = dd_graph_closure.idx_lookup[p.u_vars[i]];
    }

    DdNode* bDontCareVars_y1_dd_graph_closure = Extra_TransferPermute(dd_phi.dd, dd_graph_closure.dd, dd_phi.bDontCareVars_y1, dd_phi_to_dd_graph_closure); Cudd_Ref(bDontCareVars_y1_dd_graph_closure);
    DdNode* bDontCareVars_y2_dd_graph_closure = Extra_TransferPermute(dd_phi.dd, dd_graph_closure.dd, dd_phi.bDontCareVars_y2, dd_phi_to_dd_graph_closure); Cudd_Ref(bDontCareVars_y2_dd_graph_closure);

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
    // assign(dd_graph_closure.idx_lookup["k"], 0);
    assign(dd_graph_closure.idx_lookup["y_k"], 1);
    assign(dd_graph_closure.idx_lookup["y_k"] + dd_graph_closure.register_size, 0);


    DdNode* k = Cudd_bddIthVar(dd, dd_graph_closure.idx_lookup["k"]); Cudd_Ref(k);
    DdNode* cube_y1 = Cudd_bddAnd(dd, Cudd_Not(k) , bDontCareVars_y1_dd_graph_closure); Cudd_Ref(cube_y1);
    DdNode* cube_y2 = Cudd_bddAnd(dd, k , bDontCareVars_y2_dd_graph_closure); Cudd_Ref(cube_y2);
    DdNode* cube_y1_y2 = Cudd_bddOr(dd, cube_y1, cube_y2); Cudd_Ref(cube_y1_y2);
    cube = Cudd_bddAnd(dd, cube, cube_y1_y2); Cudd_Ref(cube);
    
    
    // Everything except y_k is fixed
    fix(dd_graph_closure.idx_lookup["k"]);
    for (int i = dd_graph_closure.idx_lookup["y_k"] + 1; i < dd_graph_closure.register_size; i++) {
        fix(i);
    }
    
    dd_graph_closure.transition = Cudd_bddOr(dd, dd_graph_closure.transition, cube); Cudd_Ref(dd_graph_closure.transition);

    Cudd_RecursiveDeref(dd, cube_y1_y2);
    Cudd_RecursiveDeref(dd, cube_y1);
    Cudd_RecursiveDeref(dd, cube_y2);
    Cudd_RecursiveDeref(dd, bDontCareVars_y1_dd_graph_closure);
    Cudd_RecursiveDeref(dd, bDontCareVars_y2_dd_graph_closure);
    Cudd_RecursiveDeref(dd, cube);
    return 0;
}

int Algo::get_free_vars() {
    DdNode* bS = dd_graph_closure.reached; Cudd_Ref(bS);
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

    DdNode* bFreeVars_dd_graph_closure = Cudd_bddAnd(dd, Cudd_Not(bS_0), Cudd_Not(bS_1)); Cudd_Ref(bFreeVars_dd_graph_closure);

    dd_phi.bFreeVars_y1 = Extra_TransferPermute(dd, dd_phi.dd, bFreeVars_dd_graph_closure, dd_graph_closure_to_dd_phi); Cudd_Ref(dd_phi.bFreeVars_y1);
    
    Cudd_RecursiveDeref(dd, bFreeVars_dd_graph_closure);
    Cudd_RecursiveDeref(dd, bS);
    Cudd_RecursiveDeref(dd, bS_0);
    Cudd_RecursiveDeref(dd, bS_1);
    Cudd_RecursiveDeref(dd, exist_abs_cube);
    return 0;
}