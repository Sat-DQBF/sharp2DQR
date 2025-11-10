#ifndef ABC_INTERFACE_HPP
#define ABC_INTERFACE_HPP

#include <string>

extern "C"
{
    #define ABC_USE_STDINT_H
    #include "misc/util/abc_global.h"
    #include "misc/util/abc_namespaces.h"
    #include "misc/vec/vec.h"
    #include "base/abc/abc.h"
    #include "base/main/abcapis.h"
    #include "base/main/main.h"
    #include "base/main/mainInt.h"
    #include "bdd/bbr/bbr.h"
    #include "bdd/cudd/cudd.h"
    #include "bdd/dddmp/dddmp.h"

}

class ABC {
    public:
        bool verbose = false;
        Abc_Frame_t * pAbc;

        ABC() {
            pAbc = Abc_FrameGetGlobalFrame();
        }

        ~ABC() {
            Abc_Stop();
        }

        int exec(std::string command) {
            return Cmd_CommandExecute(pAbc, command.c_str());
        }

        int read(std::string filename, bool fraig = true) {
            int status = exec(("read " + filename).c_str());
            if (status != 0) {
                return status;
            }
            if (fraig) {
                return exec("fraig");
            } else {
                return exec("strash");
            }
        }

        int reach();

};

int Aig_ManVerifyUsingBdds_int_dd( Aig_Man_t * p, DdManager * dd, Saig_ParBbr_t * pPars, DdNode *& bReached );

int Aig_ManComputeReachable_dd( DdManager * dd, Aig_Man_t * p, DdNode ** pbParts, DdNode * bInitial, DdNode ** pbOutputs, Saig_ParBbr_t * pPars, int fCheckOutputs, DdNode *& bReached );

int reachTransRelation(DdManager * dd, DdNode * bInitial, DdNode * bTransition, DdNode * bNeg_property, Saig_ParBbr_t * pPars, DdNode *& bReached, int nStateSize);
int reachTransRelation_computeReachable(DdManager * dd, DdNode * bInitial, DdNode * bTransition, DdNode * bNeg_property, Saig_ParBbr_t * pPars, DdNode *& bReached, int nStateSize);


#endif 