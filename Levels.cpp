#include <math.h>
#include "MDefs.h"
#include "Settings.h"
#include "Levels.h"
#include "MtxOps.h"

/* ****************************************************************************
 * IMPLEMENTAZIONE DI ENCAPSULATOR
 * ****************************************************************************/

Encapsulator::Encapsulator(SpMat* mtx) {
    sp = mtx;
    ds = NULL;
    isSparse = true;
    useHeap = true;
    isSetDegree = false;
    isSetnumEdges = false;
}

Encapsulator::Encapsulator(SpMat& mtx) {
    sp = &mtx;
    ds = NULL;
    isSparse = true;
    useHeap = false;
    isSetDegree = false;
    isSetnumEdges = false;
}

Encapsulator::Encapsulator(DMat* mtx) {
    sp = NULL;
    ds = mtx;
    isSparse = false;
    useHeap = true;
    isSetDegree = false;
    isSetnumEdges = false;
}

Encapsulator::Encapsulator(DMat& mtx) {
    sp = NULL;
    ds = &mtx;
    isSparse = false;
    useHeap = false;
    isSetDegree = false;
    isSetnumEdges = false;
}

const size_t Encapsulator::rows() const {
    if (isSparse)
        return spM().rows();
    else
        return dsM().rows();
}

const size_t Encapsulator::cols() const {
    if (isSparse)
        return spM().columns();
    else
        return dsM().columns();
}

bool Encapsulator::isZeroMatrix() const {
cout << "isZeroMatrix()" << endl;
    if (isSparse)
        return (sp->nonZeros() == 0);
    else
        return (ds->nonZeros() == 0);
cout << "isZeroMatrix() finito" << endl;
}

DMat& Encapsulator::dsM() const {
    return (*ds);
}

const DMat& Encapsulator::dsMconst() const {
    return (*ds);
}

SpMat& Encapsulator::spM() const {
    return (*sp);
}

const SpMat& Encapsulator::spMconst() const {
    return (*sp);
}

void Encapsulator::diag(SpMat& here) const {

    size_t n = this->rows();
    if (here.rows() != n)
        here.resize(n, n, false);
    else if (here.nonZeros() > 0)
        here.reset();

    here.reserve(n);
    if (isSparse) {
        for (size_t i = 0; i < n; i++) {
            here.append(i, i, (*sp)(i, i));
            here.finalize(i);
        }
    } else {
        for (size_t i = 0; i < n; i++)
            here(i, i) = (*ds)(i, i);
    }
}

void Encapsulator::strongAdjacency(SpMat& Adj, SpMat& strongAdj) const {
    size_t n = this->rows();
    //Copio ds in Adj così poi posso togliere tutti gli elementi diagonali
    if (isSparse)
        Adj = -(*sp);
    else
        Adj = -(*ds);

    for (size_t i = 0; i < n; ++i)
        Adj(i, i) = 0;
    strongAdj = Adj;

    /* Qua filtro le smallEntries che sono quei valori della matrice tali che
     *  |a_ij| < SETUP_AGGR_ELIMINATION_THRESHOLD * min( max(B(i),B(j) )
     * dove B è il vettore che contiene i massimi di ciascuna colonna di A */
    //Mi genero il vettore dei MAX, che contiene il massimo elemento di ogni COLONNA di Adj
    DynamicVector<double> B(n);
    for (size_t i = 0; i < n; i++) {
        double max_i = 0.0;
        for (SpMat::ConstIterator it = strongAdj.begin(i); it != strongAdj.end(i); it++) {
            if (it->value() > max_i) max_i = it->value();
        }
        B[i] = max_i;
    }
    //Adesso rimuovo le small entry da Adj. Posso sfruttare il fatto che Adj è simmetrica
    // e risparmiare anche molti calcoli.
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            //Sfrutto la simmetria per non valutare il triangolo inferiore
            //Inoltre salto l'elemento a_ii che è 0
            if (abs(strongAdj(i, i)) < SETUP_AGGREGATION_WEAK_EDJE_THRESHOLD * max(B[i], B[j])) {
                strongAdj(i, j) = 0;
                strongAdj(j, i) = 0;
                //^^Qua sfrutto la simmetria per rimuovere due smallEntry con un compare solo
            }
        }
    }
    //Se smallEntries.m non mente ho finito. Adj contiene la strongAdjacency
}

void Encapsulator::adjacency(SpMat& adj) const {
    if (isSparse)
        adj = -(*sp);
    else
        adj = -(*ds);
    for (size_t i = 0; i<this->rows(); ++i)
        adj(i, i) = 0;
}

//Degree è di double ora mi sa che è una bella cazzata questa.

const DynamicVector<double, blaze::columnVector>& Encapsulator::degreeVct() {
    if (isSetDegree)
        return degree;

    if (isSparse) {
        degree.resize(sp->rows());
        for (size_t i = 0; i < sp->rows(); i++) {
            degree[i] = sp->nonZeros(i);
            if ((*sp)(i, i) != 0)
                degree[i]--;
        }
    } else {
        degree.resize(ds->rows());
        for (size_t i = 0; i < ds->rows(); i++) {
            degree[i] = ds->nonZeros(i);
            if ((*ds)(i, i) != 0)
                degree[i]--;
        }
    }
    return degree;



}

size_t Encapsulator::numEdges() {
    if (isSetnumEdges == true)
        return _numEdges;
    else {
        _numEdges = 0;

        if (isSparse) {
            for (size_t i = 0; i < sp->rows(); i++) {
                _numEdges += sp->nonZeros(i);
                if ((*sp)(i, i) != 0) //L'elemento sulla diagonale non è un arco.
                    _numEdges--;
            }
        } else {
            for (size_t i = 0; i < ds->rows(); i++) {
                _numEdges += ds->nonZeros(i);
                if ((*ds)(i, i) != 0)
                    _numEdges--;
            }
        }
        isSetnumEdges = true;
    }
    return _numEdges;
}

Encapsulator::~Encapsulator() {
    if (useHeap) {
        if (sp) delete sp;
        else delete ds;
    }
}

std::ostream & operator<<(std::ostream & os, Encapsulator & enc) {
    if (enc.isSparse)
        return os << enc.spM();
    return os << enc.dsM();
}

/* ****************************************************************************
 * IMPLEMENTAZIONE DI QPSTAGE
 * ****************************************************************************/

qPStage::qPStage(SpMat*P, DynamicVector<double>* q, DynamicVector<size_t>* Fset, DynamicVector<size_t>* Cset, size_t n) {
    P = P;
    q = q;
    Fset = Fset;
    Cset = Cset;
    n = n;
}

qPStage::~qPStage() {
    delete P;
    delete q;
    delete Fset;
    delete Cset;
}

/* ****************************************************************************
 * IMPLEMENTAZIONE DI LEVEL
 * ****************************************************************************/

template <class Matrix >
Level::Level(LevelType tp, Matrix * mtx) {
    type = tp;
    K = 0;
    isAboveElimination = false;
    x = NULL;
    r = NULL;
    A = new Encapsulator(mtx);
}

template <class Matrix >
Level::Level(LevelType tp, Matrix & mtx) {
    type = tp;
    K = 0;
    isAboveElimination = false;
    x = NULL;
    r = NULL;
    A = new Encapsulator(mtx);
}

Level::~Level() {
    delete A;
    if (x) delete x;
    if (r) delete r;
}

bool Level::isRelaxationFast(size_t mgIdx) {
    /* NOTA BENE:
     * In Matlab queto metodo viene utilizzato prima della creazione di un livello
     * Aggregation per determinare se ha senso o meno andare avanti e tentare una
     * aggregazione di nodi. Ha quindi come parametro il livello finer, cioè
     * l'ultimo generato. Nella mia implementazione invece questo è un metodo
     * di classe che viene invocato SUL livello superiore. In questo modo
     * è possibile utilizzare direttamente i membri dell'ultimo livello, perché
     * sono membri di quell'oggetto. Quindi dove in Matlab c'è finerLevel.* io
     * userò this->*
     */

    /* Naive estimation of the last nu-initial sweeps */
    size_t nu = SETUP_RELAX_ACF_MIN_SWEEPS + 2 * mgIdx;
    size_t tvNu = SETUP_TV_SWEEPS;
    size_t initial = 3; // Discard these many initial iterations from ACF estimate
    if (tvNu < initial) {
        cerr << "#TV relaxations < #initial relaxations, case not yet supported\n";
        //Rendilo un errore vero
    }

    DynamicVector<double> _x(A->rows());
    DynamicVector<double> _r;

    for (size_t i = 0; i < A->rows(); i++) {
        /* Anche la rand fornita da Blaze ritorna valori in [0,1) :)
         * Quindi posso generare valori casuali in (-1,1) allo stesso modo */
        _x[i] = (double) 2.0 * blaze::rand<double>() - (double) 1.0;
    }

    if (A->isSparse)
        _r = -(A->spM()) * _x;
    else
        _r = -(A->dsM()) * _x;

    /* Ogni chiamata a tvRelax modifica i vettori passati come argomento!! */
    this->tvRelax(_x, _r, initial);

    /* Qua calcolo il denominatore che mi serve per trovare relaxACF.
     * Questo perché mi serve __QUESTA__ versione di _x per calcolaro e se lo
     * faccio qui posso evitare di farmi inutili copie
     */
    //norm(x-mean(x))
    DynamicVector<double> vec_mean(_x.size());
    vec_mean = MtxOps::mean(_x);
    vec_mean = _x - vec_mean;
    double denom = MtxOps::norma2(vec_mean);
    BLAZE_USER_ASSERT((denom != 0), "Il denominatore di relaxACF non può essere zero!");

    /* Ora che ho usato questa "versione di _x ed _r posso andare tranquillamente
     * avanti a modificarli perché non mi serviranno più. Risparmio delle copie */
    //[tv, rtv] = level.tvRelax(x , r, tvNu-initial);
    if ((tvNu - initial) > 0)
        this->tvRelax(_x, _r, tvNu - initial);

    //Adesso me li copio nei miei tv
    if (this->x == NULL)
        /* Se la matrice x non è ancora stata costruita va costruita e tv diventa
         * il primo tv, corrispondente quindi alla prima colonna di this->x */
        this->x = new DMat(_x.size(), 1UL);

    else
        /* Altrimenti vanno liminati quelli che c'erano e questo
         * diventa il primo tv comunque */
        this->x->resize(_x.size(), 1, false);

    // level.x = tv; tv è il mio _x attuale
    blaze::column((*this->x), 1) = _x;

    /* Qui devo fare la stessa cosa per _r */
    if (this->r == NULL)
        this->r = new DMat(_r.size(), 1);
    else
        this->r->resize(_r.size(), 1, false);
    //level.r = rtv;
    blaze::column((*this->r), 1) = _r;

    /* A questo punto posso continuare ad usare _x ed _r per fare l'ultimo
     * rilassamento [y, r]  = level.tvRelax(tv, rtv, nu-tvNu);
     * Userò poi y (cioè il mio _x modificato) per calcolare il numeratore
     * di relaxACF
     */

    //[y, r]  = level.tvRelax(tv, rtv, nu-tvNu);
    this->tvRelax(_x, _r, nu - tvNu);

    //Riuso vec_mean per calcolare il numeratore. _x adesso è y di Matlab
    vec_mean.resize(_x.size(), false);
    vec_mean = MtxOps::mean(_x);
    vec_mean = _x - vec_mean;
    // norm(y-mean(y))
    double num = MtxOps::norma2(vec_mean);

    BLAZE_USER_ASSERT(((nu - initial) != 0), "Il denominatore dell'esponente di relaxACF non può essere zero perché dividere per zero è sbagliato!");

    //relaxAcf = (norm(y-mean(y))/norm(x - mean(x)))^(1/(nu-initial)) =
    // = (num /denom)^(1/(nu-initial))
    double relaxAcf = pow((num / denom), ((double) 1.0 / ((double) (nu - initial))));

    /* Se converge con velocità maggiore di SETUP_MAX_COARSE_RELAX_ACF
     * allora si terminerà il coarsening
     */
    return (relaxAcf > SETUP_MAX_COARSE_RELAX_ACF);
}

bool Level::canCoarsen() const {
cout << "canCoarsen() di Level" << endl;
    return ( (A->rows() >= MAX_DIRECT_SOLVE_SIZE) && !A->isZeroMatrix());
cout << "canCoarsen() di Level finito" << endl;
}

void Level::tvRelax(DynamicVector<double>& x, DynamicVector<double>& r, size_t nu) {
    /* Nel nostro caso (opzioni di LAMG di default) usa gsrelax(). */
    if (nu > 0)
        relax(x, r, x, nu, true); //Gli passo x come b tanto non viene usato
}

void Level::generateTV() {
    if (K <= 1) {
        cout << "K <= 1. Ho presunto che per questi valori di K non si cia da"
                " fare niente in generateTV ma potrebbe essere una cosa falsa."
                " Io comunque esco\n";
        return;
    }
    size_t offset = 0; //Questa sarà la prima colonna libera
    if (this->x == NULL) /* Se x è vuota creo la matrice */
        x = new DMat(A->rows(), K);
    else {
        offset = this->x->columns();
        /* Altrimenti la estendo affinché possa contenere fino a K tv */
        this->x->extend(0, K - offset, true);

    }

    /* Idem per r */
    if (r == NULL)
        r = new DMat(A->rows(), K);
    else
        r->extend(0, K - offset, true);

    /* Ora genero tutti i tv che mancano per arrivare a K e li aggiungo
     * ad x ed r come colonne di queste matrici.
     */
    DynamicVector<double> _x(this->x->rows());
    DynamicVector<double> _r;
    for (size_t i = offset; i < K; i++) {
        //Inizializzo il vettore _x
        for (size_t j = 0; j < _x.size(); j++)
            _x[j] = (double) 2.0 * blaze::rand<double>() - (double) 1.0;
        //Calcolo _r
        if (A->isSparse)
            _r = -(this->A->spM()) * _x;
        else
            _r = -(this->A->dsM()) * _x;
        //Rilasso _x ed _r SETUP_TV_SWEEPS volte */
        this->tvRelax(_x, _r, SETUP_TV_SWEEPS);
        //Li aggiungo alle rispettive matrici
        blaze::column((*this->x), i + 1) = _x;
        blaze::column((*this->r), i + 1) = _r;
    }

}


/* Il calcolo dell'elemento x_i della (k+1)-esima iterazione utilizza solo
 * gli elementi di x^(k+1) che sono già stati calcolati (quelli con indice
 * j<i) e solo gli elementi di x^(k) che ancora non sono stati aggiornati
 * nella (k+1)-esima iterazione (quelli con j>i). Questo fa si che, a
 * differenza del metodo di Jacobi, sia possibile implementare il metodo
 * utilizzando un unico vettore x che viene sovrascritto man mano che i nuovi
 * valori vengono calcolati.*/

/* Una passata di GaussSeidel per il sistema Ax=0 */
template <class Mx>
inline void gsrelaxHom(const Mx& A, DynamicVector<double>& x, DynamicVector<double>& r);

/* Template specialization per matrice sparsa */
template <>
inline void gsrelaxHom(const SpMat& A, DynamicVector<double>& x, DynamicVector<double>& r) {
    for (size_t i = 0; i < x.size(); ++i) {
        /* x_i^(k+1) = 1/a_ii * ( b_i -SUM_A - SUM_B ) dove
         * SUM_A = sum_{j<i}a_ij*x_j^{k+1} <<elementi aggiornati in questa iterazione
         * SUM_B = sum_{j>i}a_ij*x_j^{k} << elementi ancora da aggiornare in questa iterazione
         */
        double sum_a = 0.0;
        double sum_b = 0.0;
        for (SpMat::ConstIterator j = A.begin(i); j != A.end(i); ++j) {
            if (j->index() < i)
                sum_a += (j->value() * x[j->index()]); //j->value() è a_ij;
            else
                sum_b += (j->value() * x[j->index()]);
        }
        x[i] = (sum_a - sum_b) / A(i, i);
    }
    //Ricalcolo il residuo r
    r = -A*x;
}

inline void gsrelaxHom(const DMat& A, DynamicVector<double>& x, DynamicVector<double>& r) {
    for (size_t i = 0; i < x.size(); ++i) {
        /* x_i^(k+1) = 1/a_ii * ( b_i -SUM_A - SUM_B ) dove
         * SUM_A = sum_{j<i}a_ij*x_j^{k+1} <<elementi aggiornati in questa iterazione
         * SUM_B = sum_{j>i}a_ij*x_j^{k} << elementi ancora da aggiornare in questa iterazione
         */
        double sum_a = 0.0;
        double sum_b = 0.0;
        for (size_t j = 0; j < x.size(); ++j) {
            if (j < i)
                sum_a += (A(i, j) * x[j]); //j->value() è a_ij;
            else
                sum_b += (A(i, j) * x[j]);
        }
        x[i] = (sum_a - sum_b) / A(i, i);
    }
    //Ricalcolo il residuo r
    r = -A*x;
}

/* Una passata di GaussSeidel per il sistema Ax=b */
template <class Mx>
inline void gsrelax(const Mx& A, DynamicVector<double>& x, DynamicVector<double>& r, const DynamicVector<double>& b);

/* Template specialization del metodo di GS completo per A matrice sparsa */
template <>
inline void gsrelax(const SpMat& A, DynamicVector<double>& x, DynamicVector<double>& r, const DynamicVector<double>& b) {
    for (size_t i = 0; i < x.size(); ++i) {
        /* x_i^(k+1) = 1/a_ii * ( b_i -SUM_A - SUM_B ) dove
         * SUM_A = sum_{j<i}a_ij*x_j^{k+1} <<elementi aggiornati in questa iterazione
         * SUM_B = sum_{j>i}a_ij*x_j^{k} << elementi ancora da aggiornare in questa iterazione
         */
        double sum_a = 0.0;
        double sum_b = 0.0;
        for (SpMat::ConstIterator j = A.begin(i); j != A.end(i); ++j) {
            if (j->index() < i)
                sum_a += (j->value() * x[j->index()]); //j->value() è a_ij;
            else
                sum_b += (j->value() * x[j->index()]);
        }
        x[i] = (b[i] - sum_a - sum_b) / A(i, i);
    }
    //Ricalcolo il residuo r
    r = -A*x;
}

/* Template specialization del metodo di GS completo per A matrice densa */
template <>
inline void gsrelax(const DMat& A, DynamicVector<double>& x, DynamicVector<double>& r, const DynamicVector<double>& b) {

    for (size_t i = 0; i < x.size(); ++i) {
        /* x_i^(k+1) = 1/a_ii * ( b_i -SUM_A - SUM_B ) dove
         * SUM_A = sum_{j<i}a_ij*x_j^{k+1} <<elementi aggiornati in questa iterazione
         * SUM_B = sum_{j>i}a_ij*x_j^{k} << elementi ancora da aggiornare in questa iterazione
         */
        double sum_a = 0.0;
        double sum_b = 0.0;
        for (size_t j = 0; j < x.size(); ++j) {
            if (j < i)
                sum_a += (A(i, j) * x[j]); //j->value() è a_ij;
            else
                sum_b += (A(i, j) * x[j]);
        }
        x[i] = (b[i] - sum_a - sum_b) / A(i, i);
    }
    //Ricalcolo il residuo r
    r = -A*x;
}

void Level::relax(DynamicVector<double>& x, DynamicVector<double>& r, const DynamicVector<double>& b, size_t nu, bool homogeneous) const {
    /* EliminationRelax non viene mai chiamato in LAMG quindi lo tolgo */
    if (isAboveElimination) {
        return;
    } else if (nu == 0) {
        return;
    } else if (homogeneous) { //B non ci interessa
        if (A->isSparse) {
            for (size_t k = 0; k < nu; ++k)
                gsrelaxHom(A->spM(), x, r);
        } else {
            for (size_t k = 0; k < nu; ++k)
                gsrelaxHom(A->dsM(), x, r);
        }
    } else {
        if (A->isSparse)
            for (size_t k = 0; k < nu; ++k)
                gsrelax(A->spM(), x, r, b);
        else
            for (size_t k = 0; k < nu; ++k)
                gsrelax(A->dsM(), x, r, b);
    }
}

size_t Level::tvNum() const {
    return K;
}

const DMat & Level::TVsRefConst() const {
    if (this->x == NULL)
        cerr << "Errore, dereferenzio un NULL pointer in TVsRefConst!\n";
    return (*(this->x));
}

DMat & Level::TVsRef() {
    if (this->x == NULL)
        cerr << "Errore, dereferenzio un NULL pointer in TVsRef!\n";
    return (*(this->x));
}

/**** Metodi PRIVATE di LEVEL *****/

//void Level::eliminationrelax(DynamicVector<double>& x, DynamicVector<double>& b, size_t nu) {
//    cerr << "Questo tipo di livello non implementa il metodo eliminationrelax()!\n";
//    return;
//}

/* ****************************************************************************
 * IMPLEMENTAZIONE DI LEVELFINEST
 * ****************************************************************************/

/* Helper function per i costruttori di LevelFinest*/
inline void setDfltTV(size_t & K) {
    /* LevelFinest ha come default il minimo tra TV_NUM e TV_MAX */
#if (TV_NUM <= TV_MAX )
    K = TV_NUM;
#else
    K = TV_MAX;
#endif
}

template <class Matrix >
LevelFinest::LevelFinest(Matrix * mtx) : Level(FINEST, mtx) {
    setDfltTV(this->K);
}

template <>
LevelFinest::LevelFinest(SpMat* mtx) : Level(FINEST, mtx) {
    setDfltTV(this->K);
}


template <class Matrix >
LevelFinest::LevelFinest(Matrix & mtx) : Level(FINEST, mtx) {
    setDfltTV(this->K);
}

/* ****************************************************************************
 * IMPLEMENTAZIONE DI LEVELELIMINATION
 * ****************************************************************************/

void LevelElimination::initOps(std::vector<qPStage*>& cStages, Level * finer) {
    this->K = finer->tvNum();
    this->crsnStages = cStages;


    /*c = 1:size(args.A,2); %#ok */
    size_t ss = this->A->cols();
    this->cNames.reserve(ss);
    for (size_t i = 0; i < ss; ++i) {
        this->cNames[i] = i;
    }

    /* for i = obj.stage(obj.numStages:-1:1)
     *      s = i{:};
     */
    for (std::vector<qPStage*>::reverse_iterator s = this->crsnStages.rbegin(); s != this->crsnStages.rend(); ++s) {
        /*c = s.c(c);*/
        for (size_t i = 0; i < (*s)->Cset->size(); ++i) {
            this->cNames[i] = (*(*s)->Cset)[this->cNames[i]];
        }
    }/* end (for)*/

    /* A questo punto this->cNames contiene i nomi riferiti all'A originale dei
     * nodi rimasti alla fine del coarsening. Operazione un po' costosa ma
     * funzionante */
    this->cNames_sorted = cNames;
    std::sort(this->cNames_sorted.begin(), this->cNames_sorted.end());
}

//Commentato perché a quanto pare non è più necessario
//void Level::eliminationrelax(DynamicVector<double>& x, DynamicVector<double>& b, size_t nu) {
//    /* Corrisponde a runWithDynamicVector<double> di RelaxElimination */
//    if (nu > 0) {
//        // Mi conviene almeno costruirmi x(c) (e forse dopo anche b(f))
//        DynamicVector<double> xc(elimCset.size());
//        size_t i = 0; //Indice in xc
//        for (DynamicVector<size_t>::Iterator c = elimCset.begin();
//                c != elimCset.end(); c++) {
//            xc[i++] = x[(*c)];
//        }
//
//        DynamicVector<double> bf(elimFset.size());
//        i = 0;
//        for (DynamicVector<size_t>::Iterator f = elimFset.begin();
//                f != elimFset.end(); f++) {
//            bf[i++] = b[(*f)];
//        }
//
//        /* Adesso calcolo obj.P*x(obj.c) + obj.q.*b(f);
//         * L'operatore .* di Matlab e l'operator* in blaze applicato a due
//         * vettori con la medesima transpose flag sortiscono lo stesso
//         * effetto e cioè quello di fare il component-wise product */
//        DynamicVector<double> res = (elimP * xc) + (elimQ * bf);
//
//        /* Adesso x(f) = res; cioè gli elementi di x che hanno come indice
//         * gli elementi di Fset prendono ciascuno l'i-esimo elemento di res
//         */
//        i = 0;
//        for (DynamicVector<size_t>::Iterator f = elimFset.begin();
//                f != elimFset.end(); f++) {
//            x[(*f)] = res[i++];
//        }
//    }
//    return;
//}

template <>
LevelElimination::LevelElimination(DMat* mtx, std::vector<qPStage*>& cStages, Level * finer) : Level(ELIMINATION, mtx), crsnStages(cStages) {
    this->initOps(cStages, finer);
}

template <>
LevelElimination::LevelElimination(SpMat* mtx, std::vector<qPStage*>& cStages, Level * finer) : Level(ELIMINATION, mtx), crsnStages(cStages) {
    this->initOps(cStages, finer);
}


template <>
LevelElimination::LevelElimination(DMat& mtx, std::vector<qPStage*>& cStages, Level * finer) : Level(ELIMINATION, mtx), crsnStages(cStages) {
    this->initOps(cStages, finer);
}

template <>
LevelElimination::LevelElimination(SpMat& mtx, std::vector<qPStage*>& cStages, Level * finer) : Level(ELIMINATION, mtx), crsnStages(cStages) {
    this->initOps(cStages, finer);
}


qPStage & LevelElimination::lastStage() const {
    return (*(crsnStages.back()));
}

LevelElimination::~LevelElimination() {
    for (std::vector<qPStage*>::iterator it = crsnStages.begin(); it != crsnStages.end(); it++)
        delete (*it);
    // crsnStages.clear(); //Serve?!

    for (std::vector<DynamicVector<double>*>::iterator it = bStages.begin(); it != bStages.end(); it++)
        delete (*it);
    // bStages.clear(); //Serve?!
}

void LevelElimination::coarseType(DynamicVector<double>& x) const {
    // x = x(obj.c);
    MtxOps::subVecAssign(x, cNames_sorted);
}

void LevelElimination::restrict(DynamicVector<double>& RHS) {
    this->bStages.reserve(this->crsnStages.size() + 1);
    size_t curr_bstage = 0;
    bStages[curr_bstage] = new DynamicVector<double>(RHS);
    for (std::vector<qPStage*>::const_iterator s = this->crsnStages.begin(); s != this->crsnStages.end(); ++s) {
        //b = b(s.c) + s.PT * b(s.f);
        const DynamicVector<double>& RHS_orig = (*(bStages[curr_bstage]));

        MtxOps::subVecExtract(RHS, RHS_orig, (*(*s)->Cset)); //RHS = b(s.c)
        bStages[curr_bstage + 1] = new DynamicVector<double>(RHS); //Lo salvo

        MtxOps::subVecExtract(RHS, RHS_orig, (*(*s)->Fset)); //RHS = b(s.f)
        RHS = blaze::trans((*(*s)->P)) * RHS; //RHS = P^t * b(s.f)
        (*(bStages[curr_bstage + 1])) += RHS; //b = b(s.c) + P^t*b(s.f)
        curr_bstage++;
    }
    RHS = (*(bStages.back()));
}

void LevelElimination::interpolate(DynamicVector<double>& xc) const {
    DynamicVector<double> x;
    size_t q = (crsnStages.size() - 1);

    for (std::vector<qPStage*>::const_reverse_iterator sit = crsnStages.rbegin(); sit != crsnStages.rend(); ++sit, ++q) {
        const qPStage& s = (*(*sit));
        x.resize(s.n, false); //Qesto serve per aumentare x
        //x(s.f) = s.P * cx + s.q .* bStage{q}(s.f)
        DynamicVector<double> b_f;
        MtxOps::subVecExtract(b_f, (*bStages[q]), (*s.Fset));
        DynamicVector<double> x_f = (((*s.P) * xc) + ((*s.q) * b_f));
        size_t f = 0;
        for (DynamicVector<size_t>::ConstIterator fit = s.Fset->begin(); fit != s.Fset->end(); ++fit, ++f)
            x[(*fit)] = x_f[f];
        //x(s.c) = xc;
        size_t c = 0;
        for (DynamicVector<size_t>::ConstIterator cit = s.Cset->begin(); cit != s.Cset->end(); ++cit, ++c)
            x[(*cit)] = xc[c];
        xc = x;
    }
}

/* ****************************************************************************
 * IMPLEMENTAZIONE DI LEVELAGGREGATION
 * ****************************************************************************/

void LevelAggregation::initOps(Level * finer, SpMat* R, SpMat* T, std::vector<int>& aggregateIndex, double cycleIndex) {
    this->K = finer->tvNum() + TV_INC;

    if (this->K > TV_MAX)
        this->K = TV_MAX;
    this->R = R;
    this->T = T;
    this->aggregateIndex = new std::vector<int>(aggregateIndex.size());

    for (size_t i = 0; i < aggregateIndex.size(); i++)
        (*this->aggregateIndex)[i] = aggregateIndex[i];
    this->cycleIndex = cycleIndex;


}

template <>
LevelAggregation::LevelAggregation(DMat* mtx, SpMat* R, SpMat* T, std::vector<int>& aggregateIndex, Level * finer, double cycleIndex) : Level(AGGREGATION, mtx) {
    this->initOps(finer, R, T, aggregateIndex, cycleIndex);
}

template <>
LevelAggregation::LevelAggregation(SpMat* mtx, SpMat* R, SpMat* T, std::vector<int>& aggregateIndex, Level * finer, double cycleIndex) : Level(AGGREGATION, mtx) {
    this->initOps(finer, R, T, aggregateIndex, cycleIndex);
}


LevelAggregation::~LevelAggregation() {

    delete T;
    delete R;
}

void LevelAggregation::coarseType(DynamicVector<double>& x) const {
    x = (*T) * x;
}

void LevelAggregation::restrict(DynamicVector<double>& RHS) {
    RHS = (*R) * RHS;
}

void LevelAggregation::interpolate(DynamicVector<double>& xc) const {
    const SpMat& P = blaze::trans((*R));
    xc = P * xc;
}
