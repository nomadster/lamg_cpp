/*
 * File:   LAMGLSSolver.h
 * Author: nomadsoul
 *
 * Created on 16 agosto 2013, 10.08
 */
#include <vector>
#include <iostream>

#include "Settings.h"
#include "MDefs.h"
#include "MCFLSSolver.h"
#include "Levels.h"



#ifndef LAMGLSSOLVER_H
#define	LAMGLSSOLVER_H

/*#if( OPT_USE_NAMESPACES )
namespace MCFClass_di_unipi_it {
    using namespace OPTtypes_di_unipi_it;
#endif
 */
using namespace std;
using namespace blaze;

class LAMGLSSolver : public MCFLSSolver {
public:

    LAMGLSSolver(istream *iStrm = NULL);


    // LAMGLSSolver(const LAMGLSSolver& orig);

    virtual ~LAMGLSSolver();

    /** Da implementare */
    MCFLSSolver *clone(void);
    void SetGraph(cIndex tn, cIndex tm, cIndex_Set tSn, cIndex_Set tEn);
    BOOL SetD(cHpRow tD);

    void SetRHS(cHpRow tRHS);
    /* Errore massimo che mi è consentito fare nel sistema.
     * Se x è soluzione del sistema ottenuta con LAMG
     * allora deve essere che |Ax - b| <= tPrcsn
     * Non importa fare l'override perché tanto mi basta salvarmi il vettore
     * che è la cosa che fa già da solo */
    void SetPrcsn(cHpRow tPrcsn);

    LSSStatus SolveADAT(HpRow Res); //Res è il risultato.

    /* Viene usato durante la fase di build per indicare
     * lo stato del coarsening in cui ci si trova.  */
    enum CoarseningState {
        STATE_FINEST,
        STATE_ELIMINATION,
        STATE_AGG,
        STATE_DONE_COARSENING
    };

private:

    /* Elimina le strutture dati*/
    void MemDeAlloc();
    /* [true] se la memoria è da deallocare.
     * Serve perché il costruttore non è l'unico posto in cui si dealloca.
     * vedi setGraph()
     */
    bool _doDealloc;

    /* Mi dice se devo usare le SpMat o le DMat */
    bool useSparse;

    /* Per valutare la precisione del risultato */
    DynamicVector<double>* mPrcsn;
    /* Servono per essere compliant con MCFLSSolver */
    size_t MGITCnt;
    size_t MGReITCnt;
    size_t MGItr; // max number of iterations per step
    size_t MGReItr; // max number of ri-iterations with the same RHS or D.
    size_t DirSolvDim;

    /******************************************/

    template<class M>
    bool isSparseOK(M& m);

    /* Costruisce la laplaciana con MDM^t se D!=NULL, MM^t altrimenti
     * La Laplaciana corrisponde al livello-1 del multigrid
     */
    void buildLaplacian();


    /** Costruisce il multigrid */
    void buildMultigrid();

    /* Crea un nuovo livello di tipo LevelElimination del multigrid utilizzando
     * la low-degree elimination. Può tornare NULL in caso il coarsening non sia
     * riuscito */
    template <class Matrix>
    LevelElimination* coarsenElimination(const Matrix& finerMatrix);

    /* Crea un nuovo livello del multigrid utilizzando la aggregation */
    template <class Matrix>
    LevelAggregation* coarsenAggregation(const Matrix& A, SpMat& C, const DMat& tv);

    /* Ci dice se è possibile andare avanti con la creazione del multigrid
     * [TRUE] si, [FALSE] no.*/
    bool canCoarsen(Level* level);

    /* una passata di multigrid su multiGrid*/
    void solveCycle(void);
    DynamicVector<double>* x_0; //Soluzione
    DynamicVector<double>* r; //Residuo

    /** La matrice Laplaciana */
    SpMat* Lapl;

    /* Il termine noto (è una DMat perché potrebbero essere più di un vettore) */
    DynamicVector<double>* RHS;

    /* La gerarchia dei livelli del Multigrid */
    std::vector<Level*> multiGrid;

    /* Il numero di livelli di tipo AGGREGATION */
    size_t numAGGLvls; //?? Questo c'è in Matlab, ma serve?

    double aggrCycleIndex;

    //Qua mi devo calcolare quanti archi ha il FinestLevel
    size_t numFinestEdges;

    /* Mi torna il gamma dell'l-esimo livello*/
    double gamma(size_t l) {

        double gam = -inf;
        Level* lv = multiGrid[l];
        if (lv->isAboveElimination == true) {
            gam = 1.0;
        } else {
            double ee = (double) lv->A->numEdges();
            if (ee > (0.1 * (double) numFinestEdges)) {
                //numEdges[l] > 0.1numFinestEdges gamma = IlGammaGiusto

                if (lv->type == Level::AGGREGATION) {
                    LevelAggregation* agg = dynamic_cast<LevelAggregation*> (lv);
                    gam = agg->cycleIndex;
                } else
                    gam = 1.5;
            } else {
                //numEdges[l] non è abbastanza grosso!
                ee = (double) multiGrid[l + 1]->A->numEdges();
                gam = min(2.0, (ee / (double) numFinestEdges));
            }
        }
        return gam;
    }

};
/*
#if( OPT_USE_NAMESPACES )
}; // end( namespace MCFClass_di_unipi_it )
#endif
 */
#endif	/* LAMGLSSOLVER_H */