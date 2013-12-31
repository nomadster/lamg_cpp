/*
 * File:   Levels.h
 * Author: nomadsoul
 *
 * Created on 26 settembre 2013, 11.39
 */
#include <math.h>
#include "MDefs.h"
#include "Settings.h"
#ifndef LEVELS_H
#define	LEVELS_H

/*
 **************************************************************************
 * ENCAPSULATOR. I due costruttori garantiscono che ogni istanza abbia un *
 * solo tipo di matrice al suo interno. isSparse è pubblico così si può   *
 * chiamare il Getter corretto, pena la dereferenziazione di un NULL      *
 * pointer. Se costruito utilizzando come argomento un puntatore a matrice*
 * l'oggetto si farà carico di chiamare il distruttore.                   *
 **************************************************************************
 */
class Encapsulator {
private:
    SpMat* sp;
    DMat* ds;
    /* Se true vuol dire che il distruttore si deve occupare
     * di eliminare la memoria */
    bool useHeap;

    /* Il degree della matrice contenuta. Se lo calcolo una volta mi conviene
     * cachearlo qui */
    DynamicVector<double, blaze::columnVector> degree;
    bool isSetDegree;

    size_t _numEdges;
    bool isSetnumEdges;

public:
    /* Consente al chiamante di capire quale matrice è contenuta nell'oggetto*/
    bool isSparse;

    /* Una eccezione. */
    class NonPuoiInizializzareUnOggettoEncapsulatorSenzaMatrice {
    };

    /* In questo modo non è possibile costruire oggetti senza matrice */
    Encapsulator() {
        throw NonPuoiInizializzareUnOggettoEncapsulatorSenzaMatrice();
    }

    /* Costruiscono un encapsulator a partire da un puntatore a SpMat.
     * L'encapsulator NON si occuperà di liberare memoria.*/
    Encapsulator(SpMat* mtx);
    /* Costruiscono un encapsulator a partire da un puntatore a SpMat.
     * L'encapsulator SI occuperà di liberare memoria.*/
    Encapsulator(SpMat& mtx);

    /* Costruiscono un encapsulator a partire da un puntatore a DMat.
     * L'encapsulator NON si occuperà di liberare memoria.*/
    Encapsulator(DMat* mtx);

    /* Costruiscono un encapsulator a partire da un puntatore a DMat.
     * L'encapsulator SI occuperà di liberare memoria.*/
    Encapsulator(DMat& mtx);

    /* Consente di conoscere il numero delle RIGHE della matrice contenuta
     * senza conoscerne il tipo. Senza questo metodo potrebbe essere necessario
     * inserire if "brutti" in situazioni complicate.
     */
    const size_t rows() const;

    /* Consente di conoscere il numero delle COLONNE della matrice contenuta
     * senza conoscerne il tipo. Senza questo metodo potrebbe essere necessario
     * inserire if "brutti" in situazioni complicate.
     */
    const size_t cols() const;


    /* Consente di sapere se la matrice contenuta è una matrice 0*/
    bool isZeroMatrix() const;

    /* Torna una reference alla matrice DENSA contenuta nel'oggetto.
     * La variabile isSparse deve valere FALSE, altrimenti si vince una bella
     * dereferenziazione di NULL pointer. */
    DMat& dsM() const;

    /* Come sopra solo che la reference tornata è const. */
    const DMat& dsMconst() const;

    /* Torna una reference alla matrice SPARSA contenuta nel'oggetto.
     * La variabile isSparse deve valere TRUE, altrimenti si vince una bella
     * dereferenziazione di NULL pointer. */
    SpMat& spM() const;

    /* Come sopra solo che la reference tornata è const. */
    const SpMat& spMconst() const;

    /* Il parametro 'd' diviene la matrice diagonale della matrice contenuta. */
    void diag(SpMat& d) const;

    /* Il parametro 'Adj' diviene la matrice di adiacenza, 'strongAdj' quella
     * filtrata.
     * La strongAdjacency è la matrice di adiacenza a cui sono stati tolti
     * tutti quegli archi con costo |w_ij| inferiore di una certa soglia.
     * Facciamo ciò perché:
     *  "Second, edges with very small |wuv| are discarded during aggregation.
     *   Since relaxation converges fast at the nodes that become
     *   disconnected, they need not be coarsened at all; however, to keep the
     *   coarse-level matrix a Laplacian, we aggregate all of them into a
     *   single (dummy) aggregate." cfr.PaperLungo, pag.9
     */
    void strongAdjacency(SpMat& Adj, SpMat& strongAdj) const;

    /* adj diviene la matrice di adiacenza*/
    void adjacency(SpMat& adj) const;

    /* Calcola il degree, se lo salva e mi torna una const reference.*/
    const DynamicVector<double, blaze::columnVector>& degreeVct();

    /* Mi dice quanti ARCHI stanno nella matrice contenuta.*/
    size_t numEdges();

    /* Se useHeap è TRUE elimina la memoria occupata dal puntatore */
    ~Encapsulator();

    /* Operatore<< per la stampa su ostream */
    friend std::ostream & operator<<(std::ostream & os, Encapsulator & enc);
};

/* Questa classe contiene una copia degli operatori e degli insiemi creati
 * in un dato stage di eliminazione. Poteva essere una struct invece di una classe
 * con tutti i membri public, solo che la classe fornisce costruttore e
 * distruttore ed è più comoda da usare.
 */
class qPStage {
public:

    SpMat* P;
    DynamicVector<double>* q;
    DynamicVector<size_t>* Fset;
    DynamicVector<size_t>* Cset;
    size_t n;

    /* Si copia i puntatori */
    qPStage(SpMat*P, DynamicVector<double>* q, DynamicVector<size_t>* Fset, DynamicVector<size_t>* Cset, size_t n);

    /* E si preoccupa di eliminare la memoria */
    ~qPStage();
};

/*
 **************************************************************************
 ****************** LEVEL. Superclasse di tutti i livelli *****************
 ** Definisce attributi comuni e metodi specializzabili dalle sottoclassi *
 **************************************************************************
 */
class Level {
public:

    /* Tutti i possibili tipi di livello */
    enum LevelType {
        FINEST, //Il problema originale
        ELIMINATION, //Creato con lowdegree elimination
        AGGREGATION, //Creato con aggregazione di nodi affini
        COARSEST //Ultimo livello del multigrid
    };

    /* Il tipo del livello */
    LevelType type;
    /* La matrice dei coefficienti */
    Encapsulator* A;

    /* Costruttori con puntatore a matrice.*/
    template <class Matrix >
    Level(LevelType tp, Matrix * mtx);

    /* Costruttori con reference di matrice.*/
    template <class Matrix >
    Level(LevelType tp, Matrix & mtx);

    /* Il distruttore deve essere virtual perché almeno al momento di una delete
     * su un Level* viene invocato anche il distruttore del tipo attuale. */
    virtual ~Level();

    /* Essere un livello "above elimination" influisce sul metodo di
     * rilassamento utilizzato. */
    void setAboveEliminationLevel(qPStage& lastStage) {
        /* Forse ha senso metterci un qPStage& dentro la classe, tanto il tempo di vita degli oggetti è lo stesso*/
        //        elimP = (*lastStage.P);
        //        elimQ = (*lastStage.q);
        //        elimFset = (*lastStage.Fset);
        //        elimCset = (*lastStage.Cset);
        isAboveElimination = true;
    }

    /* Stabilisce se l'Asymptotic Convergence Factor (Acf) è abbastanza veloce
     * (TRUE) o meno (FALSE).
     * Il parametro mgIdx, che corrisponde all'indice del livello al'interno
     * della gerarchia multigrid, viene utilizzato per stimare il numero (nu)
     * di passaggi di rilassamento utilizzati per rilassare i tv.
     */
    bool isRelaxationFast(size_t mgIdx);

    /* Ci dice se il livello attuale è in grado (TRUE) o meno (FALSE) di
     * effettuare ulteriore coarsening. */
    bool canCoarsen() const;

    /* Rilassa i vettori x ed r utilizzando nu passaggi della strategia di
     * rilassamento. Gli argomenti x ed r vengono modificati.
     */
    void tvRelax(DynamicVector<double>& x, DynamicVector<double>& r, size_t nu);

    /* Aggiunge K tv al livello utilizzando SETUP_TV_SWEEP passaggi di
     * rilassamento per generarli. Questo metodo non ha bisogno di parametri
     * perché le informazioni richieste sono dentro l'oggetto.
     */
    void generateTV();

    /* Funzione di rilassamento utilizzata nella fase di solve.
     * E' il metodo di Gauss-Seidel e viene usato anche per rilassare i testVectors.
     * Il parametro b ha senso solo se homogeneous == true e in tal caso
     * corrisponde al RHS */
    void relax(DynamicVector<double>& x, DynamicVector<double>& r, const DynamicVector<double>& b, size_t nu, bool homogeneous) const;

    /* Ritorna il numero di TV del livello corrente.*/
    size_t tvNum() const;

    /* Ritorna una const reference ai test-vector*/
    const DMat& TVsRefConst() const;

    /*Ritorna una reference ai test-vector. Possono essere modificati!*/
    /* Ritorna una const reference ai test-vector*/
    DMat& TVsRef();
protected:

    /* Il numero di test-vector presenti in/richiesti da questo livello */
    size_t K;
    /* Matrice di numNodes*K vettori colonna. Ciascuna colonna è un test-vector*/
    DMat* x;
    /* Matrice di numNodes*K vettori colonna. r=A*x */
    DMat* r;

    //Se commento questo metodo in LAMG funziona tutto lo stesso, quindi non viene mai usato
    //    /* Rilassa usando le informazioni sui nodi coarse calcolati nell'elimination
    //     * stage. Modifica i parametri X e B */
    //    virtual void eliminationrelax(DynamicVector<double>& x, DynamicVector<double>& b, size_t nu);

public:
    /* Mi dice se questo livello è al di sopra di uno di Elimination. Questo
     * influisce sul metodo di rilassamento utilizzato */
    bool isAboveElimination;
    /* Metodi della fase di solve */
    //Coarse Type Operator.
    virtual void coarseType(DynamicVector<double>& x) const = 0;

    // Operatore che finerLevel -> coarserLevel
    virtual void restrict(DynamicVector<double>& RHS) = 0;

    // Operatore coarserLevel -> finerLevel
    virtual void interpolate(DynamicVector<double>& xc) const = 0;

};

/*
 **************************************************************************
 ******* LEVEL FINEST. Questa classe contiene il sistema originario *******
 **************************************************************************
 */
class LevelFinest : public Level {
public:

    /* Costruttori con puntatore a matrice.*/
    template <class Matrix>
    LevelFinest(Matrix* mtx);

    /* Costruttori con reference di matrice.*/
    template <class Matrix>
    LevelFinest(Matrix& mtx);

    /* Metodi della fase di solve */
    //Coarse Type Operator.

    virtual void coarseType(DynamicVector<double>& x) const {
        cerr << "La classe LevelFinest non ha un operatore coarseType()!\n";
    }

    // Operatore che finerLevel -> coarserLevel

    virtual void restrict(DynamicVector<double>& RHS) {
        cerr << "La classe LevelFinest non ha un operatore restrict()!\n";
    }

    // Operatore coarserLevel -> finerLevel

    virtual void interpolate(DynamicVector<double>& xc) const {
        cerr << "La classe LevelFinest non ha un operatore interpolate()!\n";
    }

};

/*
 **************************************************************************
 ***** LEVEL ELIMINATION. Un livello creato con LowDegree Elimination *****
 **************************************************************************
 */
class LevelElimination : public Level {
private:

    /* Contiene tutte le info degli stage di elimination */
    std::vector<qPStage*> crsnStages;

    /* Il "nome" nel sistema originario dei nodi rimasti alla fine del
     * coarsening per low-degree elimination */
    std::vector<size_t> cNames;
    std::vector<size_t> cNames_sorted;

    /* I bStage sono cose che servono nella fase di solve */
    std::vector<DynamicVector<double>*> bStages;

    /* Inizializza. Devo ancora cpaire bene cosa faccia. */
    void initOps(std::vector<qPStage*>& cStages, Level* finer);

public:

    /* Costruttori con puntatore a matrice.*/
    template <class Matrix>
    LevelElimination(Matrix* mtx, std::vector<qPStage*>& cStages, Level* finer);

    /* Costruttori con reference di matrice.*/
    template <class Matrix>
    LevelElimination(Matrix& mtx, std::vector<qPStage*>& cStages, Level* finer);

    /* Ritorna un riferimento all'ultimo qPStage contenuto nel vector dei qPstage */
    qPStage& lastStage() const;

    /* Il distruttore deve eliminare le gli stage contenuti nei rispettivi vector */
    ~LevelElimination();

    /* Metodi della fase di solve */
    /* Restringe il parametro x del finerLevel a questo livello usando i cNames
     * calcolati nel costruttore. X viene quindi modificato perdendo
     * tutti gli elementi che non corrispondono a quegl'indici.
     */
    void coarseType(DynamicVector<double>& x) const;

    // Operatore che finerLevel -> coarserLevel

    void restrict(DynamicVector<double>& RHS);

    // Operatore coarserLevel -> finerLevel

    void interpolate(DynamicVector<double>& xc) const;

};

/*
 **************************************************************************
 ** LEVEL AGGREGATION. Un livello creato per aggregazione di nodi affini **
 **************************************************************************
 */
class LevelAggregation : public Level {
private:
    /* Inizializza. Devo ancora cpaire bene cosa faccia. */
//    void initOps(Level* finer);
    void initOps(Level * finer, SpMat* R, SpMat* T, std::vector<int>& aggregateIndex, double cycleIndex);
    SpMat* T;
    SpMat* R; //T prima di essere moltiplicata per inv(diag)
    //Ricorda P = trans(R)
    std::vector<int>* aggregateIndex;


public:
    double cycleIndex;

    /* Costruttori con puntatore a matrice.*/
    template <class Matrix>
    LevelAggregation(Matrix* mtx, SpMat* R, SpMat* T, std::vector<int>& aggregateIndex, Level* finer, double cycleIndex);


    /* Distruttore deve eliminare T*/
    ~LevelAggregation();

    /* Metodi della fase di solve */
    //Coarse Type Operator.

    void coarseType(DynamicVector<double>& x) const;

    // Operatore che finerLevel -> coarserLevel

    void restrict(DynamicVector<double>& RHS);

    // Operatore coarserLevel -> finerLevel

    void interpolate(DynamicVector<double>& xc) const;

};


#endif	/* LEVELS_H */


