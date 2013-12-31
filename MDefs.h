/*
 * File:   MDefs.h
 * Author: nomadsoul
 *
 * Created on 17 agosto 2013, 18.07
 */

#include "blaze/Blaze.h"
#include "OPTtypes.h"

#ifndef MDEFS_H
#define	MDEFS_H

enum NodeEliminationStatus {
    NOT_DECIDED,
    LOW_DEGREE,
    HIGH_DEGREE,
    ZERO_DEGREE,
    NOT_ELIMINATED
};

#define NODE_UNDECIDED -2
#define NODE_SEED -1 //In Matlab usano 0 ma li gli indici di nodo partono da 1

//Un nodo aggregato prende invece l'indice del nodo con cui è stato aggregato


/* For defining Dense Matrices we use MxN double precision Dynamic Matrix
 * stored in rowMajor order
 */
using blaze::DynamicMatrix;
using blaze::rowMajor;
typedef DynamicMatrix<double, rowMajor> DMat;
typedef blaze::DenseRow<DMat, rowMajor> DRow;


/* For defining Sparse Matrices we use MxN double precision Compressed Matrix
 * stored in rowMajor order */
using blaze::CompressedMatrix;
typedef CompressedMatrix<double, rowMajor> SpMat;
//using blaze::SparseRow;
typedef blaze::SparseRow<SpMat, rowMajor> SpRow; // La view riga

/* Vettore denso per il RHS*/
using blaze::DynamicVector;
typedef DynamicVector<double> DVect;


#endif	/* MDEFS_H */

