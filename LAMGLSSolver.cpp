
#include <blaze/Blaze.h>
#include <limits.h>
#include <iostream>
#include "LAMGLSSolver.h"
#include "MtxOps.h"

using namespace blaze;

LAMGLSSolver::LAMGLSSolver(istream* iStrm) {
cout << "Costruttore di LAMG" << endl;

    /** ATTENZIONE
     *  Commento queste variabili perché per ora mi è più comodo
     *  gestire tutto a compile-time. La comodità sta nel fatto che posso
     *  introdurre macro al volo dentro il file giusto nel momento che
     *  mi servono, sena dovermi ricordare di inizializzarle al default
     *  nel costruttore.
     * Senza contare poi che alcuni di questi parametri saranno di run-time
     * ed altri di compile-time.
     * Rimando questa decisione a progetto concluso.
     * Le variabili (ora macro) si trovano dunque nel file settings.h

    DfltdSfInpt(iStrm, buildStaticMultigrid, FALSE);
    DfltdSfInpt(iStrm, dirSolverSize, size_t(256));
     */
    /*
    uint32_t seed;
    DfltdSfInpt(iStr, seed, 0);
     */
    /* In order to reproduce certain series of random numbers, the seed of the
     * random number generator has to be set explicitly via the setSeed()
     * function. Otherwise a random seed is used for the random
     * number generation. */
    /*
    if (seed != 0)
        blaze::setSeed(seed); //Imposto il seed del RNG

     */
    useSparse = true;
    numAGGLvls = 0;
    this->RHS = NULL;
    this->mPrcsn = NULL;
    this->x_0 = NULL;
    this->r = NULL;
    /* Per "efficienza" si alloca già il posto per il numero massimo possibile
     * di puntatori
     */
    multiGrid.reserve(SETUP_MAX_LEVELS);
    /* Controlla il funzionamento del distruttore. Questo perché MemDeAlloc()
     * viene chiamata anche in SetGraph() e potrebbe accadere che il distruttore
     * venga chiamato su un oggetto senza elimenti da eliminare generando errori
     */
    _doDealloc = true;
    this->aggrCycleIndex = SETUP_CYCLE_INDEX;

    /****Compliancy con MCFLSSolver*/
    MGItr = 10;
    MGReItr = 5; // max number of ri-iterations with the same RHS or D.
    DirSolvDim = 256;
    /*************************/

cout << "Fine costruttore di LAMG" << endl;
}

void LAMGLSSolver::MemDeAlloc() {
cout << "MemDealloc()" << endl;
    if (_doDealloc) {
        //        if (M != NULL) {
        //            delete M;
        //            M = NULL;
        //        }
        //        if (Diag != NULL) {
        //            delete Diag;
        //            Diag = NULL;
        //        }
        if (RHS != NULL) {
            delete RHS;
            RHS = NULL;
        }
        if (mPrcsn != NULL) {
            delete mPrcsn;
            mPrcsn = NULL;
        }
        if (x_0 != NULL) {
            delete x_0;
            x_0 = NULL;
        }
        if (r != NULL) {
            delete r;
            r = NULL;
        }



        for (std::vector<Level*>::iterator it = multiGrid.begin(); it != multiGrid.end(); it++)
            delete (*it);
        numAGGLvls = 0;
        useSparse = true;
        _doDealloc = false;
    }
cout << "MemDealloc() fine" << endl;
    return;
}

void LAMGLSSolver::SetGraph(cIndex tn, cIndex tm, cIndex_Set tSn, cIndex_Set tEn) {

    cout << "SetGraph() inizio" << endl;

    if (n || m) { // this is not the first call
        n = m = 0;
        MemDeAlloc();
    } else if (((tn == tm) && (tm == 0))) {
        /* From Spec:
         * "Passing tn == tm == 0 is intended as a signal to the solver
         * to deallocate everything and wait for new orders;
         * in this case, all the other parameters are ignored.
         */
        n = m = 0;
        Sn = NULL;
        En = NULL;
        MemDeAlloc();
        return;
    }
    // call the method of the base class
    MCFLSSolver::SetGraph(tn, tm, tSn, tEn);
    _doDealloc = true;
    //tm è il numero di archi nel "problema originale"
    this->numFinestEdges = tm;

    mPrcsn = new DynamicVector<double>(n);
 cout << "SetGraph() fine" << endl;
}

BOOL LAMGLSSolver::SetD(cHpRow tD) {
cout << "SetD()" << endl;   
 MCFLSSolver::SetD(tD);
    //Ora ho D che è il vettore dei pesi degli archi, mi posso costruire la laplaciana.
    buildLaplacian();
    buildMultigrid();
cout << "SetD() finito" << endl;    
return true;
}

void LAMGLSSolver::buildLaplacian() {
cout << "buildLaplacian()" << endl;
    Lapl = new SpMat(n, n);

    /**************************************************************
     ********* TODO TODO TODO TODO TODO TODO TODO TODO ************
     **************************************************************
     * Importare da MCFMGLSSolver la struttura MyLess. In questo  *
     * modo sarà possibile ordinare gli archi in modo da poter    *
     * utilizzare append() e finalize() nella creazione della     *
     * Laplaciana. Più efficiente perché così usiamo il TMP!!!!   *
     * Senza questo accorgimento rischi di fare append() su righe *
     * non consecutive, pagando con UNDEFINED BEHAVIOUR!!!!       *
     * ----->(Forse mi conviene farlo in SetGraph()?)<-----       *
     * ----> Occhio che posso crearmi la laplaciana in SetGraph() *
     *  e successivamente in SetD() imposatre i valori dei nnZ.   *
     * Se faccio l'ordinamento devo esser sicuro di ordinare      *
     * anche i D[*].                                               *
     **************************************************************/


    for (size_t i = 0; i < m; i++) {
        (*Lapl)(Sn[i], En[i]) = -(D[i]);
        (*Lapl)(En[i], Sn[i]) = -(D[i]);
    }
    for (size_t i = 0; i < n; i++) {
        double sum_i = 0.0;
        for (SpMat::Iterator it = Lapl->begin(i); it != Lapl->end(i); it++) {

            sum_i += it->value();
        }
        (*Lapl)(i, i) = (sum_i);
    }
cout << "buildLaplacian() finito" << endl;
}

void LAMGLSSolver::buildMultigrid() {
cout << "buildMultigrid()" << endl;
    /* Quando siamo qui i livelli non ci sono ancora.
     * Quello che abbiamo è la matrice Laplaciana, e basta.*/
    CoarseningState cs = STATE_FINEST;
    Level* newLvl = NULL, *finerLvl = NULL;
    while (cs != STATE_DONE_COARSENING /* && non ci conviene risolvero direttamente */) {
        if (cs == STATE_FINEST) {
            // Costruisco il finestLevel
            newLvl = new LevelFinest(this->Lapl);
            numAGGLvls++;
            cs = STATE_ELIMINATION;
        } else if (cs == STATE_ELIMINATION) {
            // Fai CoarseningStrategyElimination(), se ci riesci crea un LevelElimination
            if (finerLvl->A->isSparse)
                newLvl = coarsenElimination(finerLvl->A->spM());
            else
                newLvl = coarsenElimination(finerLvl->A->dsM());
            //In ogni caso passa a STATE_AGG
            cs = STATE_AGG;
		cout << "Son qua" <<endl;
        } else if (cs == STATE_AGG) {

            if (!finerLvl->isRelaxationFast(multiGrid.size() - 1))
                cs = STATE_DONE_COARSENING;
            else {
                SpMat C;
                finerLvl->A->adjacency(C);
                if (finerLvl->A->isSparse)
                    newLvl = coarsenAggregation(finerLvl->A->spM(), C, finerLvl->TVsRefConst());
                else
                    newLvl = coarsenAggregation(finerLvl->A->dsM(), C, finerLvl->TVsRefConst());
            }
            if (newLvl == NULL)
                cs = STATE_DONE_COARSENING;
        }
cout << "Esco catena if su stato coarsening" << endl;
        /* COMPORTAMENTO: se non si è costruito nessuno livello, e lo stato era
         * ELIMINATION allora prova comunque a fare un AGGREGATION. Se invece
         * lo stato era già AGGREGATION allora bisogna terminare il coarsening.
         */



        //Se hai costruito un nuovo livello allora:
        if (newLvl != NULL) {
cout << "newlvl non nullo" << endl;
            //(1) controlla che si siano diminuiti i nodi, altrimenti dai errore!
            if (finerLvl != NULL) {
		cout << "Controllo coarsening" << endl;
                if (finerLvl->A->rows() == newLvl->A->rows())
                    cerr << "Errore!! Non è diminuito il numero di nodi durante il coarsening!\n";
                //TODO^^^^^^Emettere un vero errore!
		cout << "Controllo coarsening finito" << endl;
            }

            //(2) Aggiungi il nuovo livello al multigrid
            multiGrid.push_back(newLvl);

            if (newLvl->type == Level::AGGREGATION) {
                //(3.a) Se il nuovo livello è un LevelAggregation allora
                //    incrementa numAGGLvl
                numAGGLvls++;
            } else if (newLvl->type == Level::ELIMINATION) {
                cout << "setto ad aboveElimination" << endl;
		//(3.b) Se invece è un elimination allora imposta il
                // livello superiore ad aboveElimination()
                finerLvl->setAboveEliminationLevel(dynamic_cast<LevelElimination*> (newLvl)->lastStage());
            	cout << "setto ad aboveElimination ok" << endl;
            }
            //(4) Per finire se non si può andare avanti col coarsening
            // terminiamo.

            if (!canCoarsen(newLvl))
                cs = STATE_DONE_COARSENING;
            //Notare come se newLvl == NULL io non vengo qua e non controllo

            //(5) In ogni caso aggiorniamo i puntatori
            finerLvl = newLvl;
            newLvl = NULL;
            /* L'altro condizione di canCoarsen() in LAMG era che se il livello
             * appena costruito fosse stato nullo si andava avanti nel coarsening.
             * Ed è quello che succede qui poiché canCoarsen() la chiamo solo su
             * un livello che sono sicuro esiste :) */
        } else {
            cout << "newlvl era nullo" << endl;
            cs = STATE_DONE_COARSENING;
	}
    }
cout << "buildMultigrid() finito" << endl;
}



/*
 *****************************************************************************
 *********************** COARSENING ELIMINATION ******************************
 *****************************************************************************
 */
/* Questa funzione si occupa di selezionare tra i nodi candidati quali siano
 * lowdegree e quali no, secondo le regole riportate nella relazione.
 * Questa dichiarazione rende la funzione una funzione template di cui esistono
 * due specializzazioni. Inoltre il fatto che sia inline dovrebbe garantirmi che
 * il corpo della funzione venga sostituito alla sua chiamata come avviene in C
 * con le macro. Questo dovrebbe rendere il codice più veloce. */
template <class c>
inline void lowdegreesweep(c& m, size_t i, DynamicVector<NodeEliminationStatus>& status);

/* La specializzo nel caso di DMat*/
template<>
inline void lowdegreesweep(DMat& m, size_t i, DynamicVector<NodeEliminationStatus>& status) {
    bool hasLowDegreeNeighbour = false;
    /* Cerco vicini low-degree*/
    for (size_t j = 0; j < m.rows(); ++j) {
        if (j != i && status[j] == LOW_DEGREE) {
            hasLowDegreeNeighbour = true;
            status[i] = NOT_ELIMINATED;
            break;
        }
    }
    /* Se non aveva un vicino low_degree allora è lui il nodo low_degree!*/
    if (!hasLowDegreeNeighbour) {
        status[i] = LOW_DEGREE;
        for (size_t j = 0; j < m.rows(); ++j) {

            if (j != i)
                status[j] = NOT_ELIMINATED;
        }
    }
}

/* La specializzo nel caso di SpMat*/
template<>
inline void lowdegreesweep(SpMat& m, size_t i, DynamicVector<NodeEliminationStatus>& status) {
    bool hasLowDegreeNeighbour = false;
    for (SpMat::Iterator it = m.begin(i); it != m.end(i); ++it) {
        /* Se la seguente condizione è vera significa che il
         * nodo i ha vicini low_degree e non deve essere eliminato
         */
        if (it->index() != i && status[it->index()] == LOW_DEGREE) {
            hasLowDegreeNeighbour = true;
            status[i] = NOT_ELIMINATED;
            break;
        }
    }
    /* Se non aveva un vicino low_degree allora è lui il nodo low_degree!*/
    if (!hasLowDegreeNeighbour) {
        status[i] = LOW_DEGREE;
        for (SpMat::Iterator it = m.begin(i); it != m.end(i); it++) {
            if (it->index() != i) {

                status[it->index()] = NOT_ELIMINATED;
            }
        }
    }
}

/* Questo metodo template calcola P e q correttamente*/
template <class M>
inline void eliminationOperators(M& A, DynamicVector<size_t>& Cset, size_t fnode,
        DynamicVector<double>& q, SpMat& P, size_t& P_col, size_t P_row);

/* Specializzazione per le SpMat */
template <>
inline void eliminationOperators(SpMat& A, DynamicVector<size_t>& Cset, size_t fnode,
        DynamicVector<double>& q, SpMat& P, size_t& P_col, size_t P_row) {
    /* Inizializziamo la riga P_row con A.nonZeros(fnode) - 1 non nulli*/
    double scalingFactor = 1.0;
    /* Riserviamo spazio in ciascuna colonna di P per un adeguato numero
     * di elementi. Il -1 c'è perché (il reciproco del)l'elemento in
     * diagonale lo mettiamo in q
     */
    P.reserve(P_row, A.nonZeros(fnode) - 1);
    /* Per ciascuna f-riga prendo gli elementi in ogni c-colonna */
    DynamicVector<size_t>::Iterator ccol = Cset.begin();
    for (SpMat::Iterator frow = A.begin(fnode); frow != A.end(fnode); ++frow) {
        if (frow->index() == fnode) { //Elemento della diagonale
            q[P_row] = (1.0 / frow->value()); //Q ha elementi pari al numero di righe di R
            scalingFactor = -(q[P_row]);
        } else if (ccol != Cset.end()) {
            break; //Non ha senso andare avanti se abbiamo finito i ccol
        } else if (frow->index() == (*ccol)) {
            /* Elemento fuori della diagonale ed è anche un c-colonna */
            P.append(P_row, P_col, frow->value());
            P_col++;
            ccol++;
        }
    }

    P.finalize(P_row); //Finalizziamo la riga P_row

    /* Non dimentichiamo di scalare gli elementi della riga corrente
     * per  -scalingFactor */
    for (SpMat::Iterator it = P.begin(P_row); it != P.end(P_row); it++)
        it->value() *= -scalingFactor;
}

/* Specializzazione per le DMat */
template <>
inline void eliminationOperators(DMat& A, DynamicVector<size_t>& Cset, size_t fnode,
        DynamicVector<double>& q, SpMat& P, size_t& P_col, size_t P_row) {
    double scalingFactor = 1.0;
    P.reserve(P_row, A.nonZeros(fnode) - 1);
    DynamicVector<size_t>::Iterator ccol = Cset.begin();
    for (size_t frow = 0; frow < A.rows(); ++frow) {
        if (frow == fnode) { //Elemento sulla diagonale
            q[P_row] = (1.0 / A(frow, fnode));
            scalingFactor = -(q[P_row]);
        } else if (ccol != Cset.end()) {
            break; //Non ha senso andare avanti se abbiamo finito i ccol
        } else if (frow == (*ccol)) {
            P.append(P_row, P_col, A(frow, fnode));
            P_col++;
            ccol++;
        }
    }
    P.finalize(P_row);

    for (SpMat::Iterator it = P.begin(P_row); it != P.end(P_row); it++)
        it->value() *= -scalingFactor;
}

template <class Matrix >
LevelElimination * LAMGLSSolver::coarsenElimination(const Matrix & finerMatrix) {
cout << "coarsenElimination()" << endl;
    //Per salvare i P ed i q
    std::vector<qPStage*> cStages;
    /* Prealloco spazio per al più SETUP_ELIMINATION_MAX_STAGES così evito
     * possibili riallocazioni */
    cStages.reserve(SETUP_ELIMINATION_MAX_STAGES);

    //Copio la matrice così posso modificarla senza modificare quella del finerLevel
    Matrix A(finerMatrix);
    Matrix Acc; //Per lo Schur Complement System
    Matrix Acf; //Per lo Schur Complement System
    // Il numero di stage di eliminazione eseguiti sulla matrice A
    size_t stageNum = 0;

    /* Vector che vengono riusati*/
    DynamicVector<size_t> Cset; //Insieme dei nodi da non eliminare
    DynamicVector<size_t> Fset; //Insieme dei nodi da eliminare
    DynamicVector<size_t> degree; //Grado di ciascun nodo
    DynamicVector<size_t> candidate; //Nodi candidati all'eliminazione
    DynamicVector<NodeEliminationStatus> status; //Stato dei nodi
    /* P è sempre sparsa anche quando A è densa */
    SpMat P;
    DynamicVector<double> q;

    while (stageNum < SETUP_ELIMINATION_MAX_STAGES) {

        size_t A_rows = A.rows();
        if (A_rows <= MAX_DIRECT_SOLVE_SIZE)
            break;
        /* (1) Calcola il vettore degree della matrice. */
        degree.resize(A_rows);

        for (size_t i = 0; i < A_rows; ++i) {
            degree[i] = A.nonZeros(i);
            if (A(i, i) != 0)//Tolgo l'elemento sulla diagonale
                degree[i]--;
            else
                cerr << "La diagionale i-esima c'ha lelemento a zero!!! ERRORE MORTALE!!!" << endl;
        }

        /* (2) [f c ] = lowDegreeNodes(A,degree,MaxDegree)  */
        candidate.resize(A_rows);
        candidate = 0;
        /* Individuo i nodi candidati (degree[i] <= MAX_DEGREE) */
        size_t cnnZ = 0; //Devo usare questa variabili perché i valori "settati" di candidate comprendono il valore "0" cioè il nodo 0
        for (size_t i = 0; i < A_rows; ++i) {
            if (degree[i] <= SETUP_ELIMINATION_MAX_DEGREE) {
                candidate[cnnZ++] = i;
            }
        }
        /* I primi cnnZ elementi sono stati inizializzati */
        candidate.resize(cnnZ, true);

        status.resize(A_rows);
        status = HIGH_DEGREE; //Tutti gli elementi prendono HIGH_DEGREE
        //    status(candidate)  = 0; % Reset all relevant nodes to "not visited"
        for (size_t i = 0; i < candidate.size(); ++i) {
            status[candidate[i]] = NOT_DECIDED;
        }

        for (size_t k = 0; k < candidate.size(); ++k) {
            lowdegreesweep(A, candidate[k], status); //Template call
        }

        /* Adesso devo creare i vettori F e C, inserendovi i nodi che verranno
         * o meno eliminati */
        size_t nf = 0; //|Fset|
        size_t nc = 0; //|Cset|
        Cset.resize(A_rows, false);
        Fset.resize(A_rows, false);
        for (size_t i = 0; i < A_rows; ++i) {
            if (status[i] == LOW_DEGREE)
                Fset[nf++] = i; //Inserisco il nodo i nell'insieme F
            else
                Cset[nc++] = i; //Lo inserisco invece in C
        }

        /* L'insieme C non può mai essere vuoto, dobbiamo lasciargli almeno
         * un elemento
         */
        if (nc == 0) {
            Cset[nc++] = Fset[--nf];
            Fset[nf] = 0;

        }
        Cset.resize(nc, true); //nc non è mai 0
        Fset.resize(nf, true);

        /* FINE di (2) [f c ] = lowDegreeNodes(A,degree,MaxDegree) */


        if ((nf <= SETUP_ELIMINATION_MIN_ELIM_FRACTION * A_rows)) {
            /* Il coarsening non è abbastanza efficace perché andiamo ad eliminare
             * un numero, nf, di nodi che è inferiore alla minima soglia accettabile
             * di eliminazione. Ci fermiamo senza eliminare.*/
            break;
        }


        /* (3) Una volta individuati i nodi da eliminare devo calcolare gli
         * operatori P e q che mi consentono di eliminare tutti questi nodi
         *      [R, q] = eliminationOperators(A, f, index);
         */


        /* In MEX si crea una matrice columnMajor che poi verrà trasposta perché
         * le matrici di Matlab sono columnMajor e se uno vuole sfruttare la
         * rappresentazione interna deve usarle così.
         * Noi invece possiamo già generare la matrice trasposta costruendo
         * direttamente la matrice rowMajor e "sostituendo" ai termini riga
         * quelli colonna.
         */
        /* P è sempre una matrice sparsa perché il numero dei suoi nonzero
         * e basso anche quando A è densa */
        P.resize(nf, nc, false);
        q.resize(nf, false);

        /* Quanti elemenenti abbiamo salvato in Q (e quante righe di R
         * abbiamo costruito)
         */
        size_t P_row = 0;
        size_t P_col = 0;
        /* Per ogni f-riga di A prendi ciascun c-elemento di questa riga e formaci
         * la i-esima riga di R. L'inverso dell'f-esimo elemento di ciascuna f-riga
         * va messo in q. Inoltre ciascuna riga i di R va scalata di un fattore -q[i]
         */
        for (DynamicVector<size_t>::Iterator dvit = Fset.begin(); dvit != Fset.end(); dvit++) {
            eliminationOperators(A, Cset, (*dvit), q, P, P_col, P_row);
            P_row++; //Passiamo alla prosssima riga
        }
        //Salvo i dati così creati nel cStages
        cStages.push_back(new qPStage(new SpMat(P), new DynamicVector<double>(q),
                new DynamicVector<size_t>(Fset), new DynamicVector<size_t>(Cset), (nf + nc)));
        /* FINE (03) [R, q] = eliminationOperators(A, f, index);*/




        /* (4) Adesso devo calcolare il sistema complementare di Schur dato da
         *  A = Ac,c + Ac,f*R^t         */

        /* Acc è la sottomatrice di A che ha come elementi tutti gli a_ij tali
         * che (ij) appartiene a Cset x Cset. Quindi per ogni riga crow devo
         * prenderci tutti gli elementi che hanno indice pari a ccol */
        MtxOps::submatrix(A, Cset, Cset, Acc);

        /* Acf invece ha nc righe e nf colonne*/
        MtxOps::submatrix(A, Cset, Fset, Acf);
        /* Finalmente posso aggiornare A*/
        A = Acc + Acf * P;
        /* Acc_[nc x nc] + (Acf * P)_[nc x nc].
         * Dunque La matrice A risultante diventa una nc x nc e perdiamo esattamente nf nodi.*/
    }//Fine while(stageNum...)

    /* Usciti dal while dobbiamo creare il livello, se è possibile farlo */
    LevelElimination* ret = NULL;
    if (stageNum != 0) {

        if (useSparse) {
            if (isSparseOK(A)) {
                ret = new LevelElimination(new SpMat(A), cStages, multiGrid.back());
            } else {
                useSparse = false;
                ret = new LevelElimination(new DMat(A), cStages, multiGrid.back());
            }
        } else {

            ret = new LevelElimination(new DMat(A), cStages, multiGrid.back());
        }
    }
cout << "coarsenElimination() finito" << endl;
    return ret;
}
/****************************************************************************/

/*
 *****************************************************************************
 *********************** COARSENING AGGREGATION ******************************
 *****************************************************************************
 */

/* HELPER FUNCTION che calcola la affinity tra i nodi i e j definita come la
 * formula (3.10) a pagina 10 del paper lungo */
inline double affinity_l2(const DMat& tv, size_t u, size_t v) {
    /* Qua mi faccio due view sulle righe di tv, poi calcolo e ritorno C_ij */

    /* Qua si calcola
     *          (X_u,X_v)^2
     *  --------------------------
     *  (X_u,X_u)^2 * (X_v,X_v)^2
     *
     * con      (X,Y) = sum_{i=1}^K x^(k)y^(k)
     * che è anche noto con il nome di INNER PRODUCT
     * ovvero PRODOTTO SCALARE!!!
     */
    //(X_u,X_v)^2
    double num = 0.0;
    num = (row(tv, u) * trans(row(tv, v)));
    num *= num;

    //(X_u,X_u)
    double den1 = (row(tv, u) * trans(row(tv, u)));
    //(X_v,X_v)
    double den2 = (row(tv, v) * trans(row(tv, v)));
    //(X_u,X_u)^2 * (X_v,X_v)^2
    double den = pow(den1, 2) * pow(den2, 2);

    BLAZE_USER_ASSERT((den != 0.0), "Non è divertente dividere per zero quando si calcola la AFFINITY_L2!!\n");

    return num / den;
}

/* La matrice 'c' diventa la matrice che ha lo stesso sparsity pattern di 'adj'
 * ma che contiene le affinity calcolate usando i dati dei test vectors 'tv'
 * HELPER FUNCTION */
//inline void affinityMatrix(SpMat& c, /*std::vector<double>& cmax,*/ const SpMat& adj, const DMat & tv) {
//
//    c = adj;
//    //cmax non è usata da nessuna parte, per ora.
//    //cmax.clear();
//    //cmax.reserve(c.rows());
//
//    //double cmax_i = 0.0;
//    for (size_t i = 0; i < adj.rows(); ++i) {
//        for (SpMat::Iterator j = c.begin(i); j != c.end(i); ++j) {
//            //Sfrutto la simmetria di C
//            if (j->index() < i)
//                continue;
//            else {
//                (*j) = affinity_l2(tv, i, j->index());
//                c(j->index(), i) = j->value();
//                //Mi salvo anche i massimi, che poi mi servono (forse)
//                //if (cmax_i < j->value()) cmax_i = j->value();
//            }
//        }
//        //cmax.push_back(cmax_i);
//        //cmax_i = 0.0;
//    }
//}

/* c è la matrice di adiacenza. Al posto dei suoi nonzero metto le
 * affinity dei nodi, calcolate usando i TVs*/
inline void affinityMatrix(SpMat& c, std::vector<double>& max_aff, const DMat & tv) {

    double cmax_i = -inf;
    for (size_t i = 0; i < c.rows(); ++i) {
        for (SpMat::Iterator j = c.begin(i); j != c.end(i); ++j) {
            //Sfrutto la simmetria di C
            if (j->index() < i)
                continue;
            else {
                // (X_u,X_v) == (X_v,X_u) Quindi sfrutto la simmetria
                j->value() = affinity_l2(tv, i, j->index());
                c(j->index(), i) = j->value();

                if (cmax_i < j->value()) cmax_i = j->value();
            }
        }
        max_aff[i] = cmax_i;
    }
}

//Mi genera l'insieme deglis trong neighbors

inline void computeStrongNeighbors(const double delta, const std::vector<double>& max_aff, const SpMat& C, std::set<size_t>& strongNeighbors) {
    strongNeighbors.clear();
    for (size_t u = 0; u < C.rows(); u++) {
        for (SpMat::ConstIterator v = C.begin(u); v != C.end(u); v++) {
            if (v->index() <= u) //guardo solo il triangolo superiore, senza diagonale
                continue;
            //v->value() è c_uv
            //max_aff[u] è max{s!=u} c_us
            //max_aff[v->index()] è max{s!=v} c_sv
            //Quindi se c_uv >= delta*max{ max{s!=u} c_us, max{s!=v} c_sv}
            if (v->value() >= delta * max(max_aff[u], max_aff[v->index()])) {
                // u e v sono strong neighbors, li metto nel set

                strongNeighbors.insert(u);
                strongNeighbors.insert(v->index());
            }
        }
    }
}

//Il miglior seed per il nodo u, optrà essere aggiornata e resa figa quanto vuoi.

inline size_t bestSeed(const size_t& u, const std::set<size_t>& strNbhrs, const std::vector<int>&stat, const SpMat& C, bool& NOT_FOUND) {

    size_t s = 0;
    double cuv_max = -inf;
    for (std::set<size_t>::const_iterator v = strNbhrs.begin(); v != strNbhrs.end(); v++) {
        if (stat[(*v)] == NODE_UNDECIDED || stat[(*v)] == NODE_SEED) {
            //Se v è undecided o seed, controlla se è vicino di u
            if (C(u, (*v)) != 0) {
                // Lo è (C + derivata dalla adiacenza e quindi se ha un nonzero vuol dire che c'è un arco)
                // quindi è uno strongNeighbour seedABLE :)
                if (C(u, (*v)) > cuv_max) {

                    cuv_max = C(u, (*v)); //Aggiorno max{c_uv}
                    s = (*v); //Mi salvo l'indice del suo vicino che per ora è il massimo affine
                    NOT_FOUND = false; //Segnalo che ho trovato almeno uno strongNeighbor
                }
            }

        }
    }
    return s;
}

template <class Matrix >
LevelAggregation * LAMGLSSolver::coarsenAggregation(const Matrix& A, SpMat& C, const DMat & tv) {
cout << "coarsenAggregation()" << endl;
    /* In questo metodo l'algoritmo che eseguo è quello del paper lungo, in cui
     * la procedura bestSeed() è però ancora quella non ottimizzata. Provvederò
     * successivamente a sfruttare ilc oncetto di Energy (magari dopo che avrò
     * capito dove sta l'errore che c'è nello stesso paper).
     */

    size_t n = A.rows(); //Quanti nodi nel finerLevel
    //Calcolo la matrice delle affinity.
    /* Qua ci si tiene lo status di aggregazione dei nodi. Tutti sono undecided*/
    std::vector<int> status(n, NODE_UNDECIDED);
    /* Un nodo Hub è un nodo u tale che
     *
     *                 8* sum{v in E_u} |w_uv|*|E_v|
     *       |E_u| >= -----------------------------
     *                     sum{v in E_u} |w_uv|
     * i nodi hub sono SEED.
     */

    for (size_t u = 0; u < C.rows(); ++u) {
        double num = 0.0;
        double denom = 0.0;
        size_t degree = C.nonZeros(u);
        for (SpMat::Iterator v = C.begin(u); v != C.end(u); ++v) {
            num += (std::abs(v->value()) * C.nonZeros(v->index())); //v->value() è w_uv, C.nonZeros(v->index()) è |E_v|)
            denom += std::abs(v->value());
        }
        if (degree >= 8 * (num / denom))
            status[u] = NODE_SEED;
    }

    //La affinity massima di ciascun nodo
    std::vector<double> max_aff(n, -inf);
    //Adesso C diventa la matrice delle affinity.
    affinityMatrix(C, max_aff, tv);

    /* Queste le uso per salvarmi i vari status, b e nc degli stage*/
    std::vector< std::vector<int> > S_i(SETUP_MAX_AGGREGATION_STAGES, status);
    std::vector<double> B_i(SETUP_MAX_AGGREGATION_STAGES, inf);
    std::vector<size_t> nc_i(SETUP_MAX_AGGREGATION_STAGES, n);
    double delta = .9;
    double alpha = 1.0;
    size_t stageNum = 0;
    size_t nc;

    //L'insieme dei nodi che hanno strong neighbor
    std::set<size_t> strongNeighbors;
    bool NO_RESULT = true; //Perché è anche possibile che non si riesca ad aggregare!
    double maxCoarseningRatio = (SETUP_COARSENING_WORK_GUARD / this->aggrCycleIndex);

    while (stageNum < SETUP_MIN_AGGREGATION_STAGES || (alpha >= maxCoarseningRatio && stageNum < SETUP_MAX_AGGREGATION_STAGES)) {
        computeStrongNeighbors(delta, max_aff, C, strongNeighbors);

        if (strongNeighbors.empty()) {
            //Nessuno ha strongNeighbors, abbasso delta!
            S_i[stageNum] = status;
            nc_i[stageNum] = n;
            stageNum++;
            delta *= .6;

        } else {
            NO_RESULT = false; //C'è almeno un risultato sensato.

            nc = n;
            //aggregationStage(...)
            /* L'algoritmo che verrà eseguito nel for è il seguente:
             *
             * Per ogni u UNDECIDED
             *  Se u non ha StrongNeighbors
             *          saltalo;
             *  Altrimenti se non è più UNDECIDED
             *          saltalo; //Qualcuno l'ha fatto diventare seed
             *  Altrimenti
             *          Scegli s come il suo strongNeighbor che sia SEED o UNDECIDED, e che abbia la con maggiore affinity c_us
             *          Se esiste s
             *                  s diventa un SEED, se non lo era
             *                  u diventa aggregato di s
             *          Altrimenti
             *                  u rimane undecided e si va avanti
             * FINE
             *  */

            size_t u = 0;
            for (std::vector<int>::iterator status_u = status.begin(); status_u != status.end(); status_u++, u++) {

                if ((*status_u) != NODE_UNDECIDED) {
                    //u non è più undecided perché è stato modificato in un precedente sweep
                    continue;
                } else if (strongNeighbors.find(u) == strongNeighbors.end()) {
                    //u non ha strong neighbors, rimando la sua aggregazione.
                    //Quando calerà delta potrà diventarlo, prima però aggrego quelli forti
                    continue;
                } else {
                    bool NOT_FOUND = true;
                    size_t s = bestSeed(u, strongNeighbors, status, C, NOT_FOUND);
                    if (!NOT_FOUND) {
                        //Ho trovato un seed adeguato lo associo
                        status[s] = NODE_SEED;
                        status[u] = s;
                        // Sono il numero di aggregati, ogni volta che
                        // aggrego cala di 1 perché 2 nodi diventeranno 1 nodo solo.
                        nc = nc - 1;
                        //Non mi serve aggiornare i TV perché ho bestSeed stupido
                        //                        Devo aggiornare i TV
                        //                        for (size_t k = 0; k < K; ++k)
                        //                            tv(u, k) = tv(s, k);
                        //aggregateSize[s]++;
                        //aggregateSize[u] = aggregateSize[s];
                    }
                }

            }
            //fine aggregationStage(...)

            alpha = nc / n;
            if (alpha <= maxCoarseningRatio)
                B_i[stageNum] = 1 - alpha;
            else
                B_i[stageNum] = 1 + alpha;

            //Salvo il corrente status
            S_i[stageNum] = status;
            nc_i[stageNum] = nc;
            stageNum++;
            delta *= .6;
        }
    }

    if (NO_RESULT) {
        //Non è stato possibile fare coarsening, è il caso di tornare NULL
        return NULL;
    }

    size_t bestAggregate = 0;
    double minB_i = B_i[0];
    for (size_t i = 0; i < stageNum; i++) {
        if (minB_i < B_i[i]) {
            bestAggregate = i;
            minB_i = B_i[i];
        }
    }

    //Il miglior risultato di aggregazione
    status = S_i[bestAggregate];
    nc = nc_i[bestAggregate];

    // Qui ora faccio CoarseSetAffinityEnergyMex.computeAggregateIndex()
    std::set<int> seeds;
    for (size_t v = 0; v < status.size(); v++) {
        /* 1. stat = obj.status();
         * 2. stat(stat < 0) = 0; % Convert all undecided seeds to their own aggregates
         * 3. seeds = find(stat == 0);
         * 4. stat(seeds) = seeds; %%Aggrega i seed con i seed
         */
        if (status[v] == NODE_UNDECIDED) {
            status[v] = v; //(2) e (4) insieme: converto gli undecided in seed e li associo con loro stessi
            seeds.insert(v);
        } else if (status[v] == NODE_SEED) {
            status[v] = v; //(3) e (4)insieme: i seed li associo con loro stessi
            seeds.insert(v);
        }
    }

    // Adesso seeds contiene gli indici (a livello di matrice del finerLvl )
    // dei nodi che sono stati marcati come seed.
    size_t new_name = 0;
    for (std::set<int>::iterator s = seeds.begin(); s != seeds.end(); s++) {
        /* Con questo for do un nuovo nome ad ogni nodo seed.
         * Poi associo ogni nodo con il nuovo nodo seed.*/
        for (size_t i = 0; i < status.size(); i++)
            if (status[i] == (*s))
                status[i] = new_name;
        new_name++;
    }

    BLAZE_USER_ASSERT(new_name == nc, "Se non ho capito male i nuovi nomi dei seed"
            "corrispondono esattamente al numero di aggregati che creo.\n"
            "Se questa asserzione è falsa devo rivedere il mio ragionamento.\n");
    /* A questo punto status dovrebbe essere equivalente ad aggregateIndex di
     * CoarseSetAffinityEnergyMex. Cioè ad ogni seed è stato dato un nuovo nome
     * crescente e ogni nodo che era associato con quel seed è ora associato con
     * il nuovo nome di quel seed */

    //Ora creo la matrice T, una matrice sparsa che ha "1" in tutti i nonzero
    // e dimensione nc * n.
    SpMat T(nc, n);
    T.reserve(n); //I non zero sono pari a n
    /* i = obj.computeAggregateIndex(); << è il mio status
     * n = numel(i);
     * nc = obj.numAggregates;
     * T = sparse(i, 1 : n, ones(1, n), nc, n); */
    for (size_t i = 0; i < status.size(); ++i) {
        //i funge sia da indice per status, che da conteggio che rappresenta il j di 1:n
        T(status[i], i) = 1;
        /* ^^Uso questo, operator(), perché non ho garanzia che status[i] < status[i+1]
         * e se usassi APPEND avrei una possibile UNDEFINED BEHAVIOUR !!! */
    }


    /* Attenzione prima di calcolare T = inv(diag(T))*T mi conviene utilizzare T così
     * com'è per calcolarmi il Galerkin Coarsening R*A*P altrimenti mi tocca
     * ricostruirmi di nuovo T poiché
     *          [i,j]   = find(T);
     *          R       = sparse(i,j,ones(numel(i),1))
     * cioè R è una matrice con lo sparsity pattern di T, ma con 1 come valore
     * di tutti i nonzero, cioè esattamente T prima di essere scalato per
     * l'inverso della sua diagonale
     */

    SpMat TA = T*A;
    CompressedMatrix<double, columnMajor> Tt = trans(T);
    Matrix crsnA = TA*Tt;
    LevelAggregation* ret = NULL;

    //Ora posso scalare T.
    /* % Scale T to unit row-sums, so that the coarse system
     * % represents a [zero-sum] graph Laplacian
     *   T = (diag(sum(T,2))) \ T;
     */
    SpMat diag(nc, nc);
    for (size_t i = 0; i < nc; i++) {
        size_t sum_row_i = 0;
        for (SpMat::Iterator row_i = T.begin(i); row_i != T.end(i); row_i++) {
            sum_row_i++; //Tanto T c'ha tutti "1"
        }
        diag.append(i, i, (1.0 / (double) sum_row_i));
    }
    if (useSparse) {
        if (isSparseOK(crsnA)) {
            ret = new LevelAggregation(new SpMat(crsnA), new SpMat(T), new SpMat(diag * T), status, multiGrid.back(), aggrCycleIndex);
        } else {
            useSparse = false;
            ret = new LevelAggregation(new DMat(crsnA), new SpMat(T), new SpMat(diag * T), status, multiGrid.back(), aggrCycleIndex);
        }
    } else {

        ret = new LevelAggregation(new DMat(crsnA), new SpMat(T), new SpMat(diag * T), status, multiGrid.back(), aggrCycleIndex);
    }

    /* Qui devo aggiornare aggrCycleIndex: */
    //if ( numEdges <= 0.1 * finestNumEdges )
    size_t numCoarseEdges = ret->A->numEdges();
    double alphaEdge = ((double) numFinestEdges / (double) numCoarseEdges);
    this->aggrCycleIndex = blaze::max(1.0, blaze::min(2.0, (SETUP_COARSENING_WORK_GUARD / alphaEdge)));
cout << "coarsenAggregation() finito" << endl;
    return ret;
}


/****************************************************************************/

/* Se anche solo una delle condizioni testate è vera
 * si deve tornare falso.
 */
bool LAMGLSSolver::canCoarsen(Level * level) {
cout << "canCoarsen()" << endl;
    /* Dai commenti di O. Livne:
     * Decide whether the current level (level) can be further
     * coarsened based on the input options.
     *  - #levels exceeded, can't further coarsen
     *  - A is completely zero, no need to further coarsen
     *  - level is coarse,  no need to further coarsen
     */
    return ((numAGGLvls < SETUP_MAX_AGG_LEVELS) ||
            (multiGrid.size() < SETUP_MAX_LEVELS) ||
            (level->canCoarsen()));
cout << "canCoarsen() finito" << endl;
}

LAMGLSSolver::~LAMGLSSolver() {

    MemDeAlloc();
}

void LAMGLSSolver::SetPrcsn(cHpRow tPrcsn) {
cout << "SetPrcsn()" << endl;
    MCFLSSolver::SetPrcsn(tPrcsn);
    if (mPrcsn != NULL) {
        if (mPrcsn->size() != n) {
            mPrcsn->resize(n, false);
        }
    } else {
        mPrcsn = new DynamicVector<double>(n);
    }

    for (size_t i = 0; i < n; ++i)
        (*mPrcsn)[i] = tPrcsn[i];
cout << "SetPrcsn() finito" << endl;
}

void LAMGLSSolver::SetRHS(cHpRow tRHS) {
    if (RHS != NULL) {
        if (RHS->size() != n)
            RHS->resize(n, false);
    } else {
        RHS = new DynamicVector<double>(n);
    }

    for (size_t i = 0; i < n; ++i)
        (*RHS)[i] = tRHS[i];
}

MCFLSSolver* LAMGLSSolver::clone(void) {

    /**** From MCFMGLSSolver */
    LAMGLSSolver *clone = new LAMGLSSolver();
    clone->MGItr = MGItr;
    clone->MGReItr = MGReItr;
    clone->DirSolvDim = DirSolvDim;
    MCFLSSolver::clone(clone);
    /**********************************/

    cerr << "Should I do SomeThing Here with LAMGLSSolver vars ?!?!\n";

    return clone;
}

MCFLSSolver::LSSStatus LAMGLSSolver::SolveADAT(HpRow Res) {
    cout << "SolveADAT()" << endl;
    cout << "Ci sono " << multiGrid.size() << " livelli nel multigrid" << endl;
    size_t _zz = 0;
    if (x_0 == NULL) { //Prima chiamata a SolveADAT!
        x_0 = new DynamicVector<double>(n);
        r = new DynamicVector<double>(n);
        //ASSUMO che Res sia l'initialGuess!
        //        for (DynamicVector<double>::Iterator i = (*x_0).begin(); i != (*x_0).end(); ++i, ++_zz)
        //            (*i) = Res[_zz];
    }
    //Non contiamo la copia da HpRow a DynamicVector
    if (MCFLSSolveTimer)
        MCFLSSolveTimer->Start();

    LSSStatus SStatus;
    MGITCnt = 0;

    if (MGReITCnt == MGReItr)
        return kNError;
    else {
        do {
            SStatus = kOK;
            MGITCnt++;

            this->solveCycle();

            /* Controllo sulla precisione! */
            DynamicVector<double>::ConstIterator i, j;
            i = mPrcsn->begin();
            j = r->begin();
            for (; i != mPrcsn->end(); ++i, ++j) {
                if (abs((*j)) > (*i)) {
                    SStatus = kLwPrcsn;
                    break;
                }
            }
        } while (SStatus != kOK && MGITCnt < MGItr);


    }

    /* Fai la roba di Cyle + ProcessorSolve */




    //Non contiamo la copia da DynamicVector a HpRow
    if (MCFLSSolveTimer)
        MCFLSSolveTimer->Stop();

    /* Copia i risultati in Res[]*/
    _zz = 0;
    for (DynamicVector<double>::Iterator i = (*x_0).begin(); i != (*x_0).end(); ++i, ++_zz)
        Res[_zz] = (*i);

cout << "SolveADAT() finito" << endl;
    return ( SStatus);

}

void LAMGLSSolver::solveCycle(void) {
cout << "solveCycle()" << endl;
    /* Quante visite si possono fare a ciascun livello. Ogni valore dipende
     * dal gamma di ogni singolo livello. Questo fatto determina l'andamento
     * del ciclo. Vedi Paper-SHORT pag. 13 */
    std::vector<double> numVisits(multiGrid.size(), 0.0);
    /* Il gamma di ciascun livello. Se lo salva così può riusarlo.
     * Un valore di -oo significa che è ancora da settare */
    std::vector<double> gammas(multiGrid.size(), -inf);
    size_t finest;
    size_t coarsest = multiGrid.size() - 1;
    DynamicVector<double> RHS_orig;

    if (multiGrid.size() >= 2){
        if(multiGrid[1]->type == Level::ELIMINATION) {
        	Level* coarseLev = multiGrid[1];
        	RHS_orig = (*RHS);
        	coarseLev->restrict((*RHS));
        	if (coarseLev->A->rows() == 1) {
            		//Finito la soluzione è 0!
            		(*x_0) = 0;
            		return;
        	} else {
            		coarseLev->coarseType((*x_0));
            		finest = 1; //Si parte da 1 e b è stato già ristretto!
        	}
	} else {
        //Si parte dal primo livello.
        finest = 0;
    	}
  	size_t currLvl = finest;
    	size_t k;
    	double maxVisits;
    	while (true) {
        	if (currLvl == coarsest) {
            		k = currLvl - 1;
            		Level* lvl = multiGrid[currLvl];
            		if (currLvl == 0)
                		lvl->relax((*x_0), (*r), (*RHS), 1, false);
            		else {
                		/* Faccio 400 sweepate di GaussSeidel a botte di 10.
                 		* Ogni 10 controllo la convergenza. Se il residuo
                 		* sta dentro la precisione rThreshold, allora smetto
                 		*/
                		double rNormInitial = MtxOps::norma2((*r));
                		double rThreshold = max(1e-13, 1e-5 * (rNormInitial + DBL_EPS));
                		double rNorm = +inf;
                		size_t ceil = std::ceil((double) SOLVE_RELAX_FREQUENCY / (double) SOLVE_MAX_COARSEST_SWEEP);
                		for (size_t block = 0; block < ceil; ++block) {
                    			lvl->relax((*x_0), (*r), (*RHS), SOLVE_RELAX_FREQUENCY, false);
                    			rNorm = MtxOps::norma2((*r));
                    			if (rNorm < rThreshold)
                 		       		break;
                		}
            		}

        	} else {
            		if (currLvl == finest)
                		maxVisits = 1.0;
            		else {
                		gammas[currLvl] = this->gamma(currLvl);
                		maxVisits = gammas[currLvl] * numVisits[currLvl - 1];
            		}
            		if (numVisits[currLvl] < maxVisits)
               	 		k = currLvl + 1;
            		else
                		k = currLvl - 1;
       	 	}

        	if (k < finest)
            		break;
        	Level* fineLevel = multiGrid[currLvl];
        	Level* coarseLevel = multiGrid[currLvl + 1];
        	if (k > currLvl) {
            		numVisits[currLvl] = numVisits[currLvl] + 1;
            	//preProcess(currLvl)
            	if (!(fineLevel->isAboveElimination)) {
                	//Se non siamo above elimination dobbiamo fare relax.
                	//nuPre=1
                	fineLevel->relax((*x_0), (*r), (*RHS), 1, false);
            	}
            	coarseLevel->restrict((*RHS));
            	coarseLevel->coarseType((*x_0));
        	} else {
            	//postProcess(k)
            		coarseLevel->interpolate((*x_0));
            		if (!(fineLevel->isAboveElimination)) {
                	//nuPost=2
                	fineLevel->relax((*x_0), (*r), (*RHS), 2, false);
            		}
        	}
        	currLvl = k;
    	}//While
    }//if
    else {
    	cout << "C'ho un solo livello" << endl;
     	/* Faccio 400 sweepate di GaussSeidel a botte di 10.
         * Ogni 10 controllo la convergenza. Se il residuo
         * sta dentro la precisione rThreshold, allora smetto
         */
         double rNormInitial = MtxOps::norma2((*r));
         double rThreshold = max(1e-13, 1e-5 * (rNormInitial + DBL_EPS));
         double rNorm = +inf;
         size_t ceil = std::ceil((double) SOLVE_RELAX_FREQUENCY / (double) SOLVE_MAX_COARSEST_SWEEP);
         for (size_t block = 0; block < ceil; ++block) {
         	multiGrid[0]->relax((*x_0), (*r), (*RHS), SOLVE_RELAX_FREQUENCY, false);
         	rNorm = MtxOps::norma2((*r));
         	if (rNorm < rThreshold)
         		break;
         }

   }
   

    //postCycle(finest);
    double mean = MtxOps::mean((*x_0));

    for (size_t i = 0; i != x_0->size(); ++i)
        (*x_0)[i] -= mean;

    //Adesso calcolo r e poi ho finito.
    if (finest == 1) {
        //Si è risolto a partire dal secondo livello. Dobbiamo interpolare!
        multiGrid[1]->interpolate((*x_0));
        (*r) = multiGrid[0]->A->spM() * (*x_0) - RHS_orig;
    } else {
        (*r) = multiGrid[0]->A->spM() * (*x_0) - (*RHS);
    }
}

template<>
bool LAMGLSSolver::isSparseOK(DMat& m) {
    return false;
}

template<>
bool LAMGLSSolver::isSparseOK(SpMat& m) {
    /* Sparsa è OK finchè non ha il 65%+ di nonZero */
    return (((double) m.nonZeros() / (double) (m.rows() * m.columns())) <= 0.65);
}
