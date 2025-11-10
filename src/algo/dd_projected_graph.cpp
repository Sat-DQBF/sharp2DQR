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

int Algo::setup_dd_projected_graph(DdNode* graph_closure) {
    dd_projected_graph.dd =  Cudd_Init((p.e_vars[0].second.size() + 1) * 2, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

    int * dd_graph_closure_to_dd_projected_graph = (int *) calloc(dd_graph_closure.dd->size, sizeof(int));

    std::set<std::string> z1(p.e_vars[0].second.begin(), p.e_vars[0].second.end());

    DdManager* dd = dd_graph_closure.dd;
    DdNode* abs_cube = dd->one;

    auto assign = [&](int i, bool val) {
        // https://stackoverflow.com/questions/55246590/cudd-manipulation-of-bdds
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
        graph_closure = Cudd_bddAnd(dd, graph_closure, val ? src : Cudd_Not(src)); Cudd_Ref(graph_closure);
        abs_cube = Cudd_bddAnd(dd, abs_cube, src); Cudd_Ref(abs_cube);
        Cudd_RecursiveDeref(dd, src);
    };
    auto exists = [&](int i) {
        DdNode* src = Cudd_bddIthVar(dd, i); Cudd_Ref(src);
        abs_cube = Cudd_bddAnd(dd, abs_cube, src); Cudd_Ref(abs_cube);
        Cudd_RecursiveDeref(dd, src);
    };

    assign(dd_graph_closure.idx_lookup["k"], 0);

    for (auto& u : p.u_vars) {
        auto it = std::find(p.e_vars[0].second.begin(), p.e_vars[0].second.end(), u);
        if (it == p.e_vars[0].second.end()) {
            exists(dd_graph_closure.idx_lookup[u]);
        } else {
            dd_graph_closure_to_dd_projected_graph[dd_graph_closure.idx_lookup[u]] = p.e_vars[0].second.size() + 1 + std::distance(p.e_vars[0].second.begin(), it);
        }
    }
    dd_graph_closure_to_dd_projected_graph[dd_graph_closure.idx_lookup["y_k"]] = 2 * p.e_vars[0].second.size() + 1;
    for (int i = 0; i < p.e_vars[0].second.size(); i++) {
        dd_graph_closure_to_dd_projected_graph[dd_graph_closure.idx_lookup["target z_k"] + i] = i;
    }
    dd_graph_closure_to_dd_projected_graph[dd_graph_closure.idx_lookup["target y_k"]] = p.e_vars[0].second.size();

    graph_closure = Cudd_bddExistAbstract(dd_graph_closure.dd, graph_closure, abs_cube); Cudd_Ref(graph_closure);
    dd_projected_graph.projected_graph = Extra_TransferPermute(dd_graph_closure.dd, dd_projected_graph.dd, graph_closure, dd_graph_closure_to_dd_projected_graph); Cudd_Ref(dd_projected_graph.projected_graph);
    free(dd_graph_closure_to_dd_projected_graph);
    return 0;
}