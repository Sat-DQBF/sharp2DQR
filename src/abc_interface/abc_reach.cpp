#include "abc_interface/abc_interface.hpp"

extern "C"
{
    // Functions in base/abci
    Aig_Man_t * Abc_NtkToDar( Abc_Ntk_t * pNtk, int fExors, int fRegisters );
    int Abc_NtkDarReach( Abc_Ntk_t * pNtk, Saig_ParBbr_t * pPars );
    int Aig_ManComputeReachable( DdManager * dd, Aig_Man_t * p, DdNode ** pbParts, DdNode * bInitial, DdNode ** pbOutputs, Saig_ParBbr_t * pPars, int fCheckOutputs );
    DdNode * Bbr_bddComputeRangeCube( DdManager * dd, int iStart, int iStop );
    Abc_Cex_t * Aig_ManVerifyUsingBddsCountExample( Aig_Man_t * p, DdManager * dd, DdNode ** pbParts, Vec_Ptr_t * vOnionRings, DdNode * bCubeFirst, int iOutput, int fVerbose, int fSilent);
}

// Modified from bdd/bbr/bbrReach.c:501 Aig_ManVerifyUsingBdds_int
//          from bdd/bbr/bbrReach.c:238 Aig_ManComputeReachable
// Modifications:
// Assume precomputed global BDDs
// Assume p is transition relation
// Keep dd after computation
// Give the reached states after computation

int Aig_ManVerifyUsingBdds_int_dd( Aig_Man_t * p, DdManager * dd, Saig_ParBbr_t * pPars, DdNode *& bReached ) {
    int fCheckOutputs = !pPars->fSkipOutCheck;
    DdNode ** pbParts, ** pbOutputs;
    DdNode * bInitial, * bTemp;
    int RetValue, i;
    abctime clk = Abc_Clock();

    assert ( Saig_ManRegNum(p) == 0 );
    assert ( Saig_ManPiNum(p) % 2 == 0 );

    if ( pPars->fVerbose )
        printf( "Shared BDD size is %6d nodes.\n", Cudd_ReadKeys(dd) - Cudd_ReadDead(dd) );

    // check the runtime limit
    if ( pPars->TimeLimit && pPars->TimeLimit <= (Abc_Clock()-clk)/CLOCKS_PER_SEC )
    {
        printf( "Reached timeout after constructing global BDDs (%d seconds).\n",  pPars->TimeLimit );
        Cudd_Quit( dd );
        return -1;
    }

    if (Saig_ManPoNum(p) != 3) {
        printf( "The circuit should have exactly 2 POs for initial condition, transition relation and property.\n" );
        Cudd_Quit( dd );
        return -1;
    }
    // Property is the third PO
    pbOutputs = ABC_ALLOC( DdNode *, 1);
    pbOutputs[0] = Aig_ObjGlobalBdd((Aig_Obj_t *) Vec_PtrEntry(p->vCos, 2)); Cudd_Ref(pbOutputs[0]);

    // Transition relation is the second PO
    pbParts = ABC_ALLOC( DdNode *, 1);
    pbParts[0] = Aig_ObjGlobalBdd((Aig_Obj_t *) Vec_PtrEntry(p->vCos, 1)); Cudd_Ref(pbParts[0]);

    // Initial state is the first PO
    { // Modified from Aig_ManInitStateVarMap @ bbrReach.c:124
        DdNode ** pbVarsX, ** pbVarsY;
        int i;
        
        // set the variable mapping for Cudd_bddVarMap()
        pbVarsX = ABC_ALLOC(DdNode *, dd->size / 2);
        pbVarsY = ABC_ALLOC(DdNode *, dd->size / 2);
        for ( i = 0; i < Saig_ManPiNum(p) / 2 ; i++ )
        {
            pbVarsX[i] = dd->vars[i];
            pbVarsY[i] = dd->vars[Saig_ManPiNum(p) / 2 + i];
        }
        Cudd_SetVarMap(dd, pbVarsX, pbVarsY, Saig_ManPiNum(p) / 2);
        ABC_FREE(pbVarsX);
        ABC_FREE(pbVarsY);
    }
    bInitial = Aig_ObjGlobalBdd((Aig_Obj_t *) Vec_PtrEntry(p->vCos, 0)); Cudd_Ref(bInitial);

    // set reordering
    if ( pPars->fReorderImage )
        Cudd_AutodynEnable( dd, CUDD_REORDER_SYMM_SIFT );

    // check the result
    RetValue = -1;
    i = 2; // Property is the third PO
    if ( fCheckOutputs && !Cudd_bddLeq( dd, bInitial, Cudd_Not(pbOutputs[0]) ) )
    {
        DdNode * bIntersect;
        bIntersect = Cudd_bddIntersect( dd, bInitial, pbOutputs[0] );  Cudd_Ref( bIntersect );
        assert( p->pSeqModel == NULL );
        Abc_Print(1, "Counter-example derivation is not implemented for transition relation.\n");
        Cudd_RecursiveDeref( dd, bIntersect );
        if ( !pPars->fSilent )
            Abc_Print( 1, "Output %d of miter \"%s\" was asserted in frame %d. ", i, p->pName, -1 );
        RetValue = 0;
    }
    // explore reachable states
    if ( RetValue == -1 )
        RetValue = Aig_ManComputeReachable_dd( dd, p, pbParts, bInitial, pbOutputs, pPars, fCheckOutputs, bReached ); 

    // cleanup
    Cudd_RecursiveDeref( dd, bInitial );
    Cudd_RecursiveDeref( dd, pbParts[0] );
    Cudd_RecursiveDeref( dd, pbOutputs[0] );
    ABC_FREE( pbParts );
    ABC_FREE( pbOutputs );

    // report the runtime
    if ( !pPars->fSilent )
    {
    ABC_PRT( "Time", Abc_Clock() - clk );
    fflush( stdout );
    }
    return RetValue;
}

int Aig_ManComputeReachable_dd( DdManager * dd, Aig_Man_t * p, DdNode ** pbParts, DdNode * bInitial, DdNode ** pbOutputs, Saig_ParBbr_t * pPars, int fCheckOutputs , DdNode *& bReached)
{
    int fInternalReorder = 0;
    Bbr_ImageTree_t * pTree = NULL; // Suppress "might be used uninitialized"
    Bbr_ImageTree2_t * pTree2 = NULL; // Supprses "might be used uninitialized"
    DdNode * bCubeCs;
    DdNode * bCurrent;
    DdNode * bNext = NULL; // Suppress "might be used uninitialized"
    DdNode * bTemp;
    Cudd_ReorderingType method;
    int i, nIters, nBddSize = 0, status;
    int nThreshold = 10000;
    abctime clk = Abc_Clock();
    Vec_Ptr_t * vOnionRings;
    int fixedPoint = 0;

    status = Cudd_ReorderingStatus( dd, &method );
    if ( status )
        Cudd_AutodynDisable( dd );

    // start the image computation
    bCubeCs = Bbr_bddComputeRangeCube( dd, 0, Saig_ManPiNum(p) / 2);    Cudd_Ref( bCubeCs );
    if ( pPars->fPartition )
        pTree = Bbr_bddImageStart( dd, bCubeCs, 1, pbParts, Saig_ManPiNum(p) / 2, dd->vars + Saig_ManPiNum(p) / 2, pPars->nBddMax, pPars->fVerbose );
    else
        pTree2 = Bbr_bddImageStart2( dd, bCubeCs, 1, pbParts, Saig_ManPiNum(p) / 2, dd->vars + Saig_ManPiNum(p) / 2, pPars->fVerbose );
    Cudd_RecursiveDeref( dd, bCubeCs );
    if ( pTree == NULL )
    {
        if ( !pPars->fSilent )
            printf( "BDDs blew up during qualitification scheduling.  " );
        return -1;
    }

    if ( status )
        Cudd_AutodynEnable( dd, method );

    // start the onion rings
    vOnionRings = Vec_PtrAlloc( 1000 );

    // perform reachability analysis
    bCurrent = bInitial;   Cudd_Ref( bCurrent );
    bReached = bInitial;   Cudd_Ref( bReached );
    Vec_PtrPush( vOnionRings, bCurrent );  Cudd_Ref( bCurrent );
    for ( nIters = 0; nIters < pPars->nIterMax; nIters++ )
    { 
        // check the runtime limit
        if ( pPars->TimeLimit && pPars->TimeLimit <= (Abc_Clock()-clk)/CLOCKS_PER_SEC )
        {
            printf( "Reached timeout after image computation (%d seconds).\n",  pPars->TimeLimit );
            Vec_PtrFree( vOnionRings );
            // undo the image tree
            if ( pPars->fPartition )
                Bbr_bddImageTreeDelete( pTree );
            else
                Bbr_bddImageTreeDelete2( pTree2 );
            pPars->iFrame = nIters - 1;
            return -1;
        }

        // compute the next states
        if ( pPars->fPartition )
            bNext = Bbr_bddImageCompute( pTree, bCurrent );           
        else
            bNext = Bbr_bddImageCompute2( pTree2, bCurrent );  
        if ( bNext == NULL )
        {
            if ( !pPars->fSilent )
                printf( "BDDs blew up during image computation.  " );
            if ( pPars->fPartition )
                Bbr_bddImageTreeDelete( pTree );
            else
                Bbr_bddImageTreeDelete2( pTree2 );
            Vec_PtrFree( vOnionRings );
            pPars->iFrame = nIters - 1;
            return -1;
        }
        Cudd_Ref( bNext );
        Cudd_RecursiveDeref( dd, bCurrent );

        // remap these states into the current state vars
        bNext = Cudd_bddVarMap( dd, bTemp = bNext );                    Cudd_Ref( bNext );
        Cudd_RecursiveDeref( dd, bTemp );
        // check if there are any new states
        if ( Cudd_bddLeq( dd, bNext, bReached ) ) {
            fixedPoint = 1;
            break;
        }
        // check the BDD size
        nBddSize = Cudd_DagSize(bNext);
        if ( nBddSize > pPars->nBddMax )
            break;
        // check the result
        i = 2;
        if ( fCheckOutputs && !Cudd_bddLeq( dd, bNext, Cudd_Not(pbOutputs[0]) ) )
        {
            DdNode * bIntersect;
            bIntersect = Cudd_bddIntersect( dd, bNext, pbOutputs[0] );  Cudd_Ref( bIntersect );
            assert( p->pSeqModel == NULL );
            Abc_Print(1, "Counter-example derivation is not implemented for transition relation.\n");
            Cudd_RecursiveDeref( dd, bIntersect );
            if ( !pPars->fSilent )
                Abc_Print( 1, "Output %d of miter \"%s\" was asserted in frame %d. ", i, p->pName, Vec_PtrSize(vOnionRings) );
            Cudd_RecursiveDeref( dd, bReached );
            bReached = NULL;
            pPars->iFrame = nIters;
            break;
        }
        // get the new states
        bCurrent = Cudd_bddAnd( dd, bNext, Cudd_Not(bReached) );        Cudd_Ref( bCurrent );
        Vec_PtrPush( vOnionRings, bCurrent );  Cudd_Ref( bCurrent );
        // minimize the new states with the reached states
//        bCurrent = Cudd_bddConstrain( dd, bTemp = bCurrent, Cudd_Not(bReached) ); Cudd_Ref( bCurrent );
//        Cudd_RecursiveDeref( dd, bTemp );
        // add to the reached states
        bReached = Cudd_bddOr( dd, bTemp = bReached, bNext );           Cudd_Ref( bReached );
        Cudd_RecursiveDeref( dd, bTemp );
        Cudd_RecursiveDeref( dd, bNext );
        if ( pPars->fVerbose )
            fprintf( stdout, "Frame = %3d. BDD = %5d. ", nIters, nBddSize );
        if ( fInternalReorder && pPars->fReorder && nBddSize > nThreshold )
        {
            if ( pPars->fVerbose )
                fprintf( stdout, "Reordering... Before = %5d. ", Cudd_DagSize(bReached) );
            Cudd_ReduceHeap( dd, CUDD_REORDER_SYMM_SIFT, 100 );
            Cudd_AutodynDisable( dd );
            if ( pPars->fVerbose )
                fprintf( stdout, "After = %5d.\r", Cudd_DagSize(bReached) );
            nThreshold *= 2;
        }
        if ( pPars->fVerbose )
//            fprintf( stdout, "\r" );
            fprintf( stdout, "\n" );

        if ( pPars->fVerbose )
        {
            double nMints;
            nMints = Cudd_CountMinterm(dd, bReached, Saig_ManPiNum(p) / 2 );
            fprintf( stdout, "Reachable states = %.0f. (Ratio = %.4f %%)\n", nMints, 100.0*nMints/pow(2.0, Saig_ManPiNum(p) / 2) );
            fflush( stdout ); 
        }

    }
    Cudd_RecursiveDeref( dd, bNext );

    // free the onion rings
    Vec_PtrForEachEntry( DdNode *, vOnionRings, bTemp, i )
        Cudd_RecursiveDeref( dd, bTemp );
    Vec_PtrFree( vOnionRings );
    // undo the image tree
    if ( pPars->fPartition )
        Bbr_bddImageTreeDelete( pTree );
    else
        Bbr_bddImageTreeDelete2( pTree2 );
    if ( bReached == NULL )
        return 0; // proved reachable
    // report the stats
    if ( pPars->fVerbose )
    {
        double nMints;
        nMints = Cudd_CountMinterm(dd, bReached, Saig_ManPiNum(p) / 2 );
        if ( nIters > pPars->nIterMax || nBddSize > pPars->nBddMax )
            fprintf( stdout, "Reachability analysis is stopped after %d frames.\n", nIters );
        else
            fprintf( stdout, "Reachability analysis completed after %d frames.\n", nIters );
        fprintf( stdout, "Reachable states = %.0f. (Ratio = %.4f %%)\n", nMints, 100.0*nMints/pow(2.0, Saig_ManPiNum(p) / 2) );
        fflush( stdout );
    }
//ABC_PRB( dd, bReached );
    // Cudd_RecursiveDeref( dd, bReached );
    // SPG
    //
    if ( fixedPoint ) {
      if ( !pPars->fSilent ) {
        printf( "The miter is proved unreachable after %d iterations.  ", nIters );
      }
      pPars->iFrame = nIters - 1;
      return 1;
    }
    assert(nIters >= pPars->nIterMax || nBddSize >= pPars->nBddMax);
    if ( !pPars->fSilent )
      printf( "Verified only for states reachable in %d frames.  ", nIters );
    pPars->iFrame = nIters - 1;
    return -1; // undecided
}


int ABC::reach() {
    Saig_ParBbr_t Pars, * pPars = &Pars;
    Bbr_ManSetDefaultParams(pPars);
    pPars->fTransRel = 1;
    pPars->fVerbose = verbose;

    Abc_Ntk_t * pNtk = Abc_FrameReadNtk(pAbc);
    
    // Internal of Abc_NtkDarReach
    Aig_Man_t * pMan = Abc_NtkToDar(pNtk, 0, 1);

    // Internal of Aig_ManVerifyUsingBdds
    // Internal of Aig_ManVerifyUsingBdds_int
    printf("Number of Pi: %d\n", Saig_ManPiNum(pMan));
    DdManager * dd = Aig_ManComputeGlobalBdds(pMan, pPars->nBddMax, 1, 1, pPars->fVerbose);
    printf("Number of vars: %d\n", dd->size);
    printf("Order: ");
    for (int i = 0; i < dd->size; i++) {
        printf("%d ", dd->invperm[i]);
    }
    printf("\n");
    if (dd == NULL) {
        return -1;
    }

    DdNode * bReached;

    return Aig_ManVerifyUsingBdds_int_dd(pMan, dd, pPars, bReached);
}

int reachTransRelation(DdManager * dd, DdNode * bInitial, DdNode * bTransition, DdNode * bNeg_property, Saig_ParBbr_t * pPars, DdNode *& bReached, int nStateSize) {
    int fCheckOutputs = !pPars->fSkipOutCheck;
    DdNode * bTemp;
    int RetValue, i;
    abctime clk = Abc_Clock();

    if (pPars->fReorderImage)
        Cudd_AutodynEnable(dd, CUDD_REORDER_SYMM_SIFT);

    RetValue = -1;
    if ( fCheckOutputs && !Cudd_bddLeq( dd, bInitial, Cudd_Not(bNeg_property) ) )
    {
        DdNode * bIntersect;
        bIntersect = Cudd_bddIntersect( dd, bInitial, bNeg_property );  Cudd_Ref( bIntersect );
        if ( !pPars->fSilent ) 
            Abc_Print(1, "Counter-example derivation is not implemented for transition relation.\n");
        Cudd_RecursiveDeref( dd, bIntersect );
        if ( !pPars->fSilent )
            Abc_Print( 1, "Property violation was asserted in frame %d. ", -1 );
        RetValue = 0;
    }

    if ( RetValue == -1 ) {
        RetValue = reachTransRelation_computeReachable(dd, bInitial, bTransition, bNeg_property, pPars, bReached, nStateSize);
    }

    if ( !pPars->fSilent )
    {
    ABC_PRT( "Time", Abc_Clock() - clk );
    fflush( stdout );
    }
    return RetValue;
};



int reachTransRelation_computeReachable(DdManager * dd, DdNode * bInitial, DdNode * bTransition, DdNode * bNeg_property, Saig_ParBbr_t * pPars, DdNode *& bReached, int nStateSize) {
    int fCheckOutputs = !pPars->fSkipOutCheck;
    int fInternalReorder = 0;
    Bbr_ImageTree_t * pTree = NULL; // Suppress "might be used uninitialized"
    Bbr_ImageTree2_t * pTree2 = NULL; // Supprses "might be used uninitialized"
    DdNode * bCubeCs;
    DdNode * bCurrent;
    DdNode * bNext = NULL; // Suppress "might be used uninitialized"
    DdNode * bTemp;
    Cudd_ReorderingType method;
    int i, nIters, nBddSize = 0, status;
    int nThreshold = 10000;
    abctime clk = Abc_Clock();
    Vec_Ptr_t * vOnionRings;
    int fixedPoint = 0;

    status = Cudd_ReorderingStatus( dd, &method );
    if ( status )
        Cudd_AutodynDisable( dd );

    // start the image computation
    bCubeCs = Bbr_bddComputeRangeCube( dd, 0, nStateSize);    Cudd_Ref( bCubeCs );
    if ( pPars->fPartition )
        pTree = Bbr_bddImageStart( dd, bCubeCs, 1, &bTransition, nStateSize, dd->vars + nStateSize, pPars->nBddMax, pPars->fVerbose );
    else
        pTree2 = Bbr_bddImageStart2( dd, bCubeCs, 1, &bTransition, nStateSize, dd->vars + nStateSize, pPars->fVerbose );
    Cudd_RecursiveDeref( dd, bCubeCs );

    if ( pTree == NULL)
    {
        if ( !pPars->fSilent )
            printf( "BDDs blew up during qualitification scheduling.  " );
        return -1;
    }

    if ( status )
        Cudd_AutodynEnable( dd, method );

    // start the onion rings
    vOnionRings = Vec_PtrAlloc( 1000 );

    // perform reachability analysis
    bCurrent = bInitial;   Cudd_Ref( bCurrent );
    bReached = bInitial;   Cudd_Ref( bReached );
    Vec_PtrPush( vOnionRings, bCurrent );  Cudd_Ref( bCurrent );
    for ( nIters = 0; nIters < pPars->nIterMax; nIters++ )
    {
        // check the runtime limit
        if ( pPars->TimeLimit && pPars->TimeLimit <= (Abc_Clock()-clk)/CLOCKS_PER_SEC )
        {
            printf( "Reached timeout after image computation (%d seconds).\n",  pPars->TimeLimit );
            Vec_PtrFree( vOnionRings );
            // undo the image tree
            if ( pPars->fPartition )
                Bbr_bddImageTreeDelete( pTree );
            else
                Bbr_bddImageTreeDelete2( pTree2 );
            pPars->iFrame = nIters - 1;
            return -1;
        }

        // compute the next states
        if ( pPars->fPartition )
            bNext = Bbr_bddImageCompute( pTree, bCurrent );
        else
            bNext = Bbr_bddImageCompute2( pTree2, bCurrent );
        if ( bNext == NULL )
        {
            if ( !pPars->fSilent )
                printf( "BDDs blew up during image computation.  " );
            if ( pPars->fPartition )
                Bbr_bddImageTreeDelete( pTree );
            else
                Bbr_bddImageTreeDelete2( pTree2 );
            Vec_PtrFree( vOnionRings );
            pPars->iFrame = nIters - 1;
            return -1;
        }
        Cudd_Ref( bNext );
        Cudd_RecursiveDeref( dd, bCurrent );

        // remap these states into the current state vars
        bNext = Cudd_bddVarMap( dd, bTemp = bNext );                    Cudd_Ref( bNext ); 
        Cudd_RecursiveDeref( dd, bTemp );
        // check if there are any new states
        if ( Cudd_bddLeq( dd, bNext, bReached ) ) {
            fixedPoint = 1;
            break;
        }
        // check the BDD size
        nBddSize = Cudd_DagSize(bNext);
        if ( nBddSize > pPars->nBddMax )
            break;
        // check the result
        if ( fCheckOutputs && !Cudd_bddLeq( dd, bNext, Cudd_Not(bNeg_property) ) )
        {
            DdNode * bIntersect;
            bIntersect = Cudd_bddIntersect( dd, bNext, bNeg_property );  Cudd_Ref( bIntersect );
            if ( !pPars->fSilent )
                Abc_Print(1, "Counter-example derivation is not implemented for transition relation.\n");
            Cudd_RecursiveDeref( dd, bIntersect );
            if ( !pPars->fSilent )
                Abc_Print( 1, "Property violation was asserted in frame %d. ", nIters );
            Cudd_RecursiveDeref( dd, bReached );
            bReached = NULL;
            pPars->iFrame = nIters;
            break;
        }

        // get the new states
        bCurrent = Cudd_bddAnd( dd, bNext, Cudd_Not(bReached) );        Cudd_Ref( bCurrent );
        Vec_PtrPush( vOnionRings, bCurrent );  Cudd_Ref( bCurrent );
        // add to the reached states
        bReached = Cudd_bddOr( dd, bTemp = bReached, bNext );           Cudd_Ref( bReached );
        Cudd_RecursiveDeref( dd, bTemp );
        Cudd_RecursiveDeref( dd, bNext );
        if ( pPars->fVerbose )
            fprintf( stdout, "Frame = %3d. BDD = %5d. ", nIters, nBddSize );
        if ( fInternalReorder && pPars->fReorder && nBddSize > nThreshold )
        {
            if ( pPars->fVerbose )
                fprintf( stdout, "Reordering... Before = %5d. ", Cudd_DagSize(bReached) );
            Cudd_ReduceHeap( dd, CUDD_REORDER_SYMM_SIFT, 100 );
            Cudd_AutodynDisable( dd );
            if ( pPars->fVerbose )
                fprintf( stdout, "After = %5d.\r", Cudd_DagSize(bReached) );
            nThreshold *= 2;
        }
        if ( pPars->fVerbose )
        {
            double nMints;
            nMints = Cudd_CountMinterm(dd, bReached, nStateSize );
            fprintf( stdout, "Reachable states = %.0f. (Ratio = %.4f %%)\n", nMints, 100.0*nMints/pow(2.0, nStateSize) );
            fflush( stdout ); 
        }
    }
    Cudd_RecursiveDeref( dd, bNext );

    // free the onion rings
    Vec_PtrForEachEntry( DdNode *, vOnionRings, bTemp, i )
        Cudd_RecursiveDeref( dd, bTemp );
    Vec_PtrFree( vOnionRings );
    // undo the image tree
    if ( pPars->fPartition )
        Bbr_bddImageTreeDelete( pTree );
    else
        Bbr_bddImageTreeDelete2( pTree2 );
    if ( bReached == NULL )
        return 0; // proved reachable
    // report the stats
    if ( pPars->fVerbose )
    {
        double nMints;
        nMints = Cudd_CountMinterm(dd, bReached, nStateSize );
        if ( nIters > pPars->nIterMax || nBddSize > pPars->nBddMax )
            fprintf( stdout, "Reachability analysis is stopped after %d frames.\n", nIters );
        else
            fprintf( stdout, "Reachability analysis completed after %d frames.\n", nIters );
        fprintf( stdout, "Reachable states = %.0f. (Ratio = %.4f %%)\n", nMints, 100.0*nMints/pow(2.0, nStateSize) );
        fflush( stdout );
    }
//ABC_PRB( dd, bReached );
    // Cudd_RecursiveDeref( dd, bReached );
    // SPG
    //
    if ( fixedPoint ) {
      if ( !pPars->fSilent ) {
        printf( "The miter is proved unreachable after %d iterations.  ", nIters );
      }
      pPars->iFrame = nIters - 1;
      return 1;
    }
    assert(nIters >= pPars->nIterMax || nBddSize >= pPars->nBddMax);
    if ( !pPars->fSilent )
      printf( "Verified only for states reachable in %d frames.  ", nIters );
    pPars->iFrame = nIters - 1;
    return -1; // undecided
    
}