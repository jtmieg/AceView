#ifndef _POLYNOME_H
#define _POLYNOME_H

#include "ac.h"
#include "matrix.h"

#define minAbs 0.0000001
#define INDEXMAX 1000
#define GMAX 256
#define SMAX (256*sizeof(short))
typedef struct termStruct {
  int type ;
  double complex z ;/* complex scalar multiplier. If zero, the whole TT is NULL */
  int sqrt1,sqrt2 ;
  short sigma[GMAX] ; /* sigma     matrices : non-commutative list of index "ab" means sigma_a sigma-bar_b */
  short sigB[GMAX] ;  /* sigma-bar matrices : non-commutative list of index "ab" means sigma-bar_a sigma_b */
  int N ; /* Taylor degree in x  symbol */
  char x[GMAX]  ; /* a,b,  i,j  x(meaning chi)  symbol to exponentiate */
  short sqrti ; /* exp(i pi/4) */
  char theta[GMAX] ; /* grassman indices, they anticommute ans quare to zero */
  short g[GMAX] ;     /* Lorentz metric */
  short gg[GMAX] ;    /* group metric */
  short eps[GMAX] ; /* espislon anti symmetric set of n times 4 indices */ 
  short mm[4][GMAX] ;      /*(k p q r)_mu momenta :  "1 ab" means the product p_a p_b */
                          /*  "1 a"  "2 b" means the product p_a q_b */
  int  denom[4] ;     /* number of terms of the form 1/k^, 1/(k+p)^2, 1/(k+p+q)^2, 1/(k+p+q+r)^2 */
  int Id2 ; /* Pauli identity matrix, needed its value is 2 when we trace */
  short freeIndex[GMAX] ;
} TT ;

typedef struct polynomeStruct *POLYNOME ;
struct polynomeStruct {
  AC_HANDLE h ;
  int id ;
  BOOL isFlat, isSum, isProduct ;
  POLYNOME p1, p2 ;
  TT tt ;
} ;

extern short firstDummyIndex ;
short newDummyIndex (void) ;


POLYNOME reduceIndices (POLYNOME pp) ;

POLYNOME newMultiSum (AC_HANDLE h, POLYNOME ppp[]) ;
POLYNOME expand (POLYNOME pp) ;
POLYNOME newPolynome (AC_HANDLE h) ;
void showPol (POLYNOME pp) ;
POLYNOME squareMomentaCleanUp (POLYNOME pp) ;
POLYNOME dimIntegralDo (POLYNOME pp, int pass) ;
BOOL freeIndex (POLYNOME pp) ;

/**** New Polynomes of all types ***/
POLYNOME newPolynome (AC_HANDLE h) ;
void freePolynome (POLYNOME pp) ;
POLYNOME copyPolynome (POLYNOME p1, AC_HANDLE h) ;
POLYNOME newScalar (complex float z, AC_HANDLE h) ;
POLYNOME newG (short mu, short nu, AC_HANDLE h) ;
POLYNOME newEpsilon (short a, short b, short c, short d, AC_HANDLE h) ;
/* newAG (... ,  0) antisymmetric link 1/2(ac bd - ad bc). Optionally adding the i epsilon 
 * newAG (... , +1) is the self-dual projector P+ = 1/4 ( ac bd - ad bc) + i/2 epsilon(abcd)
 * newAG (... , -1) is the self-dual projector P- = 1/4 ( ac bd - ad bc) - i/2 epsilon(abcd)
 *     WE have (P+)^2 = (P+),   (P-)^2 = (P-),   (P+)(P-) = (P-)(P+) = 0 
 */
POLYNOME newAG (short a, short b, short c, short d, int parity, AC_HANDLE h) ;
POLYNOME newAG6 (short a, short b, short c, short d, short e, short f, AC_HANDLE h) ;
POLYNOME newK (short cc, AC_HANDLE h) ;
POLYNOME newP (short cc, AC_HANDLE h) ;
POLYNOME newQ (short cc, AC_HANDLE h) ;
POLYNOME newR (short cc, AC_HANDLE h) ;
POLYNOME newPQR (int pqr, short mu, AC_HANDLE h) ;
POLYNOME newSigma (short cc, AC_HANDLE h) ;
POLYNOME newSigB (short cc, AC_HANDLE h) ;
POLYNOME polynomeScale (POLYNOME pp, double complex z) ;
POLYNOME newSum (POLYNOME p1, POLYNOME p2, AC_HANDLE h) ;
POLYNOME newProduct (POLYNOME p1, POLYNOME p2, AC_HANDLE h) ;
POLYNOME newMultiSum (AC_HANDLE h, POLYNOME ppp[]) ;
POLYNOME newMultiProduct (AC_HANDLE h, POLYNOME ppp[]) ;

#endif
