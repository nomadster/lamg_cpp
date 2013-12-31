/*
 * File:   MtxOps.h
 * Author: nomadsoul
 *
 * Created on 18 novembre 2013, 22.19
 */
#include <blaze/util/Assert.h>
#include <math.h>
#include "MDefs.h"
#ifndef MTXOPS_H
#define	MTXOPS_H

class MtxOps {
public:

    /* subm = orig(row_idx,col_idx)*/
    static inline void submatrix(const SpMat& orig, const DynamicVector<size_t>& row_idx, const DynamicVector<size_t>& col_idx, SpMat & subm) {
        //Resize senza preservare gli eventuali elementi
        subm.resize(row_idx.size(), col_idx.size(), false);
        for (size_t i = 0; i < orig.rows(); i++) {
            subm.reserve(i, orig.nonZeros(i));
            for (DynamicVector<size_t>::ConstIterator crows = row_idx.begin();
                    crows != row_idx.end(); crows++) {
                size_t j = 0;

                for (DynamicVector<size_t>::ConstIterator fcol = col_idx.begin();
                        fcol != col_idx.end(); fcol++) {

                    subm.append(i, j++, orig((*crows), (*fcol)));
                }
            }
            subm.finalize(i);
        }
    };

    /* subm = orig(row_idx,col_idx)*/
    static inline void submatrix(const DMat& orig, const DynamicVector<size_t>& row_idx, const DynamicVector<size_t>& col_idx, DMat & subm) {
        //Resize senza preservare gli eventuali elementi
        subm.resize(row_idx.size(), col_idx.size(), false);
        for (size_t i = 0; i < orig.rows(); i++) {
            for (DynamicVector<size_t>::ConstIterator crows = row_idx.begin();
                    crows != row_idx.end(); crows++) {
                size_t j = 0;
                for (DynamicVector<size_t>::ConstIterator fcol = col_idx.begin();
                        fcol != col_idx.end(); fcol++) {

                    subm(i, j++) = orig((*crows), (*fcol));
                }
            }

        }
    };

    //    enum _a {
    //        USE_ALL_ROWS, // bla(:,set);
    //        USE_ALL_COLS // bla(set,:)
    //    };
    //
    //    /* subm = orig(idx,:) oppure subm = orig(:,idx)*/
    //    static inline void semicolon(const DMat& orig, const std::vector<size_t>& idx, _a d, DMat& subm) {
    //        if (d == USE_ALL_COLS) { /* subm = orig(idx,:) */
    //            //Ridimensiono la subm
    //            subm.resize(idx.size(), orig.columns(), false);
    //            //popolo subm prendendo i valori alle idx-righe di orig
    //            for (std::vector<size_t>::const_iterator i = idx.begin(); i != idx.end(); ++i) {
    //
    //                blaze::DenseRow<DMat> r_sub = blaze::row(subm, (*i));
    //                blaze::DenseRow<const DMat> r_orig = blaze::row(orig, (*i));
    //                r_sub = r_orig; //Ce stanno li SmartEspressioTempleits e usamoli no?
    //            }
    //        } else {
    //            subm.resize(orig.rows(), idx.size(), false);
    //            for (std::vector<size_t>::const_iterator j = idx.begin(); j != idx.end(); ++j) {
    //                blaze::DenseColumn<DMat> c_sub = blaze::column(subm, (*j));
    //                blaze::DenseColumn<const DMat> c_orig = blaze::column(orig, (*j));
    //                c_sub = c_orig; //Ammazzafigata li smartespressiontempleits!
    //            }
    //        }
    //    }

    /* Riduce il vettore v a dimensione idx.size() in modo tale che gli
     * elementi rimanenti abbiano i valori che precedentemente stavano nelle
     * posizioni corrispondenti agli elementi di idx. In Matlab orig = orig(idx)
     */
    static inline void subVecAssign(DynamicVector<double>& v, const std::vector<size_t>& idx) {
        BLAZE_USER_ASSERT((v.size() >= idx.size()), "In subVecAssign v.size() è MINORE di idx.size()\n");
        if (v.size() == idx.size())
            return;
        for (size_t i = 0; i < idx.size(); ++i)
            v[i] = v[idx[i]];
        v.resize(idx.size(), true);
    }

    /* subv = orig(set);*/
    static inline void subVecExtract(DynamicVector<double>& subv, const DynamicVector<double>& orig, DynamicVector<size_t>& set) {
        BLAZE_USER_ASSERT((orig.size() >= set.size()), "In subVectExtract orig.size() è MINORE di set.size()\n");
        if (orig.size() == set.size())
            return;
        subv.resize(set.size(), false);
        size_t i = 0;
        for (DynamicVector<size_t>::Iterator j = set.begin(); j != set.end(); ++j, ++i) {
            subv[i] = orig[(*j)];
        }
    }

    /* Norma2 a.k.a. Euclidean Norm a.k.a L2 distance a.k.a L2 norm:
     * sqrt( x[0]^2 + x[1]^2 + ... + x[n]^2 )*/
    static inline double norma2(DynamicVector<double>& v) {
        double ret = 0.0;
        for (size_t i = 0; i < v.size(); i++) {
            ret += (v[i] * v[i]);
        }
        return sqrt(ret);
    }

    /* Calcola la media degli elementi del vettore */
    static inline double mean(DynamicVector<double>& v) {
        double ret = 0.0;
        for (size_t i = 0; i < v.size(); i++)
            ret += v[i];
        return ret / (double) (v.size());
    }

    //    static inline void subVecAssign(const DynamicVector<double>& orig, const std::vector<size_t>& idx, DynamicVector<double>& subm) {
    //
    //        subm.resize(idx.size(), false);
    //        size_t i = 0;
    //        //popolo subm prendendo i valori alle idx-righe di orig
    //        for (std::vector<size_t>::const_iterator it = idx.begin(); it != idx.end(); ++it, ++i)
    //            subm[i] = orig[(*it)];
    //
    //    }
    //
    //    static inline void semicolon(const DMat& orig, const DynamicVector<size_t>& idx, _a d, DMat& subm) {
    //        if (d == USE_ALL_COLS) { /* subm = orig(idx,:) */
    //            //Ridimensiono la subm
    //            subm.resize(idx.size(), orig.columns(), false);
    //            //popolo subm prendendo i valori alle idx-righe di orig
    //            for (DynamicVector<size_t>::ConstIterator i = idx.begin(); i != idx.end(); ++i) {
    //
    //                blaze::DenseRow<DMat> r_sub = blaze::row(subm, (*i));
    //                blaze::DenseRow<const DMat> r_orig = blaze::row(orig, (*i));
    //                r_sub = r_orig; //Ce stanno li SmartEspressioTempleits e usamoli no?
    //            }
    //        } else {
    //            subm.resize(orig.rows(), idx.size(), false);
    //            for (DynamicVector<size_t>::ConstIterator j = idx.begin(); j != idx.end(); ++j) {
    //                blaze::DenseColumn<DMat> c_sub = blaze::column(subm, (*j));
    //                blaze::DenseColumn<const DMat> c_orig = blaze::column(orig, (*j));
    //                c_sub = c_orig; //Ammazzafigata li smartespressiontempleits!
    //            }
    //        }
    //    }
};
#endif	/* MTXOPS_H */

