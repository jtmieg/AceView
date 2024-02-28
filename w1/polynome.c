 /*  File: polynome.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg 2024
 *-------------------------------------------------------------------
 * This file is part of a toolbox for formal calculus in QFT
 *
 * Description:
 *    Handle Polynomes: 
 *                 show, add, mutiply, exponentiate
 *    Implements t'Hooft dimensional regularized integrals
 *                 by defining momenta and denominators 
 *    The elments of the polynomes are
 *          complex numbers, 
 *          rational fractions (to avoid rounding errors)
 *          symbols (x,y,...), 
 *          Grassman variables (i,j,k...)
 *          Metric with dummy indices and epsilon symbols
 *          Pauli matrices sigma and sigmaBar SIG. SIGB
 *          Momenta with a vector index, 
 *          Specialized propagator denominators 1/k^2(k+p)^2(k+p+q)^2
 *    Handle matrices of polynomes PMX
 *                 show, add, mutiply, exponentiate, determinant.
 * Exported functions:
 *              See Header file: polynome.h
 * HISTORY:
 * Last edited:
 * Created: Feb 2024 (mieg)
 *-------------------------------------------------------------------
 */

#include "polynome.h"

/***********************************************************************************************************************************************/
/***************************************************************** Utilities  ******************************************************************/
/***********************************************************************************************************************************************/

short firstDummyIndex = 'a' ;

short newDummyIndex (void)
{
  short cc = firstDummyIndex++ ;
  
  if (cc >= INDEXMAX)
    messcrash ("too many dummy indices") ;
  return (short)cc  ;
}


/* Warning : we have no provision to store a big array in the database */
 /***********************************************************************************************************************************************/
/* look for all indices, 
 * edit = FALSE: count occurences of a,b,c...
 * edit = TRUE: replace by a continuous list liberating the high values
 */
static void reduceIndicesTtDo (short *cp, KEYSET ks, BOOL edit)
{  
  int i ;
  
  for (i = 0 ; *cp && i<GMAX ; cp++, i++)
    {
      if (! edit)
	keySet (ks, *(short *)cp) ++ ;
      else
	*cp = keySet (ks, *(short *)cp) ;
    }
  return ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static void reduceIndicesTt (POLYNOME pp, KEYSET ks, BOOL edit)
{
  reduceIndicesTtDo (pp->tt.g, ks, edit) ;
  reduceIndicesTtDo (pp->tt.gg, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigma, ks, edit) ;
  reduceIndicesTtDo (pp->tt.sigB, ks, edit) ;
  reduceIndicesTtDo (pp->tt.eps, ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[0], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[1], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[2], ks, edit) ;
  reduceIndicesTtDo (pp->tt.mm[3], ks, edit) ;
} /* reduceIndicesTt */

/***********************************************************************************************************************************************/

static void reduceIndicesDo (POLYNOME pp, KEYSET ks, BOOL edit)
{
  if (! pp) 
    return ;
  reduceIndicesDo (pp->p1, ks, edit) ;
  reduceIndicesDo (pp->p2, ks, edit) ;
  if (pp->tt.type)
    reduceIndicesTt (pp, ks, edit) ;
} /* reduceIndices */

/***********************************************************************************************************************************************/

static int reduceIndicesIsProduct (POLYNOME pp)
{
  int i, j, n ;
  KEYSET ks = keySetCreate () ;
  KEYSET ks2 = keySetCreate () ;
  KEYSET ks3 = keySetCreate () ;

  /* count the occurence of all indices */
  reduceIndicesDo (pp, ks, FALSE) ;
  /* contract the list by dropping unused values */
  for (i = j = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 1)    /* block and do not rename the free indices present once */
	{
	  keySet (ks2, 'a' + i)  = 'a' + i ;
	  keySet (ks3, 'a' + i)  = 1 ;       /* block it */ 
	  if (firstDummyIndex <= 'a' + i)
	    firstDummyIndex = 'a' + i + 1 ;
	}
      else if (n > 2)
	{
	  showPol (pp) ;
	  messcrash ("no index (%c) should be repeated %d > 2 times", 'a'+i, n) ;
	}
    }
  for (i = 0 ; 'a' + i < keySetMax (ks) ; i++)
    {
      n = keySet (ks, 'a' + i) ;
      if (n == 2)
	{
	  for (j = 0 ; ; j++)
	    if (keySet (ks3, 'a' + j)  == 0)  /* first unblocked index  */
	      {
		keySet (ks2, 'a' + i)  = 'a' + j ;
		keySet (ks3, 'a' + j)  = 1 ;
		if (firstDummyIndex <= 'a' + j)
		  firstDummyIndex = 'a' + j + 1 ;
		break ;
	      }
	}
    }
  /* rename all indices */
  reduceIndicesDo (pp, ks2, TRUE) ;

  keySetDestroy (ks) ; 
  keySetDestroy (ks2) ; 
  keySetDestroy (ks3) ; 
  return 'a' + j ;
} /* reduceIndicesIsProduct */

/***********************************************************************************************************************************************/

POLYNOME reduceIndices (POLYNOME pp)
{
  static int level = 0 ;

  level++ ;
  if (0 && level == 1)
    {
      pp = squareMomentaCleanUp (pp) ;
    }
  if (pp && pp->isSum)
    {
      if (pp->p1) 
	pp->p1 = reduceIndices (pp->p1) ;
      if (pp->p2) 
	pp->p2 = reduceIndices (pp->p2) ;
    }
  if (pp && ! pp->isSum) reduceIndicesIsProduct (pp) ;
  level-- ;
  if (level == 0 && pp)
    {
      if (! pp->isSum)
	reduceIndicesIsProduct (pp) ;
    }
  return pp ;
} /* reduceIndices */

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static int polOrder (const void *va, const void *vb)
{
  const POLYNOME p1 = *(POLYNOME *)va ;
  const POLYNOME p2 = *(POLYNOME *)vb ;
  float complex z1 = p1->tt.z ;
  float complex z2 = p2->tt.z ;
  short sqrti1 = (p1->tt.sqrti + 40000) % 4 ;
  short sqrti2 = (p2->tt.sqrti + 40000) % 4 ;
  int id1 = p1->tt.Id2 ;
  int id2 = p2->tt.Id2 ;
  int s ;
  p1->tt.z = 0 ;       
  p2->tt.z = 0 ;
  p1->tt.sqrt1 = 0 ;       
  p2->tt.sqrti = 0 ;
  p1->tt.Id2 = 0 ;
  p2->tt.Id2 = 0 ;
  memset (p1->tt.freeIndex, 0, SMAX) ;
  memset (p2->tt.freeIndex, 0, SMAX) ;
  s = memcmp (&(p1->tt), &(p2->tt), sizeof(TT)) ;
  p1->tt.z = z1 ;       
  p2->tt.z = z2 ;
  p1->tt.sqrti = sqrti1 ;       
  p2->tt.sqrti = sqrti2 ;
  p1->tt.Id2 = id1 ;
  p2->tt.Id2 = id2 ;       

  return s ;
}

/***********************************************************************************************************************************************/

static void sortPolGetSum (POLYNOME pp, Array aa)
{
  if (! pp)
    return ;
  if (pp->isSum)
    {
      sortPolGetSum (pp->p1, aa) ;
      sortPolGetSum (pp->p2, aa) ;
    }
  else if (pp)
    array (aa, arrayMax (aa), POLYNOME) = pp ;
  return ;
}

/***********************************************************************************************************************************************/

static POLYNOME sortReduceSum (Array aa)
{
  POLYNOME pp = 0 ;
  int ii, jj ;

  if (!aa)
    return 0 ;

  if (arrayMax (aa) > 1)
    arraySort (aa, polOrder) ;

  for (ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
    {	/* add identical terms */				
      int jj ;
      POLYNOME q1 = arr (aa, ii, POLYNOME) ; 
      if (! q1 || ! q1->tt.type)
	continue ;
      for (jj = ii + 1 ; jj < arrayMax (aa) ; jj++)
	{
	  POLYNOME q2 = arr (aa, jj, POLYNOME) ; 
	  if (q2)
	    {
	      int n = polOrder (&q1, &q2) ;
	      if (n == 0)
		{
		  q1->tt.z += q2->tt.z ;
		  arr (aa, jj, POLYNOME)  = 0 ;
		}
	      else if (n < 0)
		break ;
	    }
	}
    }
  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
    {  /* elimimate the killed terms */
      POLYNOME qq = arr (aa, ii, POLYNOME) ; 
      if (qq && (qq->isSum || qq->isProduct || qq->tt.z != 0))
	arr (aa, jj++, POLYNOME) = qq ;
    }
  arrayMax (aa) = jj ;
  
  if (arrayMax (aa) == 0)
    pp = 0 ;
  else if (arrayMax (aa) == 1)
    pp = arr (aa, 0, POLYNOME) ;
  else 
    {
      POLYNOME qq ;
      
      pp =  newPolynome (pp->h) ;
      for (qq = pp, ii = 0 ; ii + 1 < arrayMax (aa) ; ii++)
	{
	  qq->isSum = TRUE ;
	  qq->p1 = arr (aa, ii, POLYNOME) ;
	  if (ii < arrayMax (aa) - 2)
	    { qq->p2 = newPolynome (pp->h) ; qq = qq->p2 ; }
	}
      qq->p2 = arr (aa, ii, POLYNOME) ;
    }
  return pp ;
} /* sortReduceSum */
    
/***********************************************************************************************************************************************/

static POLYNOME sortPol (POLYNOME pp)
{
  if (! pp)
    return 0 ;
  else if (pp->isSum)
    {
      int ii = 0, jj = 0 ;
      Array aa = arrayCreate (32, POLYNOME) ;
      
      sortPolGetSum (pp, aa) ;
      if (0)
	{
	  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
	    {
	      POLYNOME qq = arr (aa, ii, POLYNOME) ; 
	      qq = sortPol (qq) ;
	      if (qq)
		arr (aa, jj++, POLYNOME) = qq ;
	    }
	  arrayMax (aa) = jj ;
	}
      pp = sortReduceSum (aa) ;
      arrayDestroy (aa) ;
    }
  else if (pp->isProduct)
    {
      pp->p1 = sortPol (pp->p1) ;
      pp->p2 = sortPol (pp->p2) ;
      if (pp->p1 && pp->p2)
	return pp ;
      return 0 ;
    }
  else if (! pp->tt.type ||  cabs (pp->tt.z) < minAbs)
    return 0 ;

  return pp ;
}

/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/
/***********************************************************************************************************************************************/

static BOOL ppIsNumber (POLYNOME pp)
{
  TT tt = pp->tt ;
  BOOL ok = tt.type ;
  if (tt.x[0])  ok = FALSE ;
  if (tt.mm[0][0])  ok = FALSE ;
  if (tt.mm[1][0])  ok = FALSE ;
  if (tt.mm[2][0])  ok = FALSE ;
  if (tt.mm[3][0])  ok = FALSE ;
  if (tt.g[0])  ok = FALSE ;
  if (tt.gg[0])  ok = FALSE ;
  if (tt.sigma[0])  ok = FALSE ;
  if (tt.sigB[0])  ok = FALSE ;
  if (tt.eps[0])  ok = FALSE ;
  return ok ;
} /* ppIsNumber */

/*******************************************************************************************/

static void cleanTtSqrt (TT *ttp)
{
  if (ttp->sqrt1 == 0) ttp->sqrt1 = 1 ;
  if (ttp->sqrt2 == 0) ttp->sqrt2 = 1 ;
  
  ttp->sqrti = (ttp->sqrti + 40000)%4 ;
  switch (ttp->sqrti)
    {
    case 0:  break ;
    case 1:  break ;
    case 2: ttp->sqrti = 0 ; ttp->z *= 1I ; break ;
    case 3: ttp->sqrti = 1 ; ttp->z *= 1I ; break ;
    }

  if (*ttp->theta)
    {
      BOOL more = TRUE ;
      int sign = 1 ;
      while (more)
	{
	  more = FALSE ;
	  for (char *cp = ttp->theta ; ! more && *cp ; cp++)
	    for (char *cq = cp + 1 ; !more && *cq ; cq++)
	      {
		if (*cp == *cq)
		  { ttp->z = 0 ; return ; }
		if (*cp > *cq)
		  { char cc = *cp ; *cp = *cq ; *cq = cc ; more = TRUE ; sign = -sign ; break ;}
	      }
	}
      ttp->z *= sign ;
    }

  if (ttp->sqrt1 * ttp->sqrt2 > 1)
    {
      int i ;
      /* reduce to common denominator */
      for (i = 2 ; i <= ttp->sqrt1 && i <= ttp->sqrt2 ; i++)
	{
	  if ((ttp->sqrt1 % i == 0) && (ttp->sqrt2 % i == 0))
	    {
	      ttp->sqrt1 /= i ;
	      ttp->sqrt2 /= i ;
	      i = 1 ;
	      continue ;
	    }
	}
      /* remove squares */
      for (i = 2 ; i*i <= ttp->sqrt1 ; i++)
	{
	  if (ttp->sqrt1 % (i*i) == 0)
	    {
	      ttp->sqrt1 /= (i*i)  ;
	      ttp->z *= i ;
	      i = 1 ;
	      continue ;
	    }
	}
      for (i = 2 ; i*i <= ttp->sqrt2 ; i++)
	{
	  if (ttp->sqrt2 % (i*i) == 0)
	    {
	      ttp->sqrt2 /= (i*i)  ;
	      ttp->z /= i ;
	      i = 1 ;
	      continue ;
	    }
	}
    }
} /* cleanTtSqrt */

/*******************************************************************************************/

static void gprintf (char *prefix, short *sp)
{
  printf (" %s_", prefix) ;
  sp-- ;
  while (*++sp)
    printf ("%c", (char)(*sp % 256)) ;
  printf (" ") ;
}

  static void showTT (POLYNOME pp)
{
  TT tt ;
  int i ; 
  
  if (!pp) 
    return ;
  tt = pp->tt ;  
  if (!tt.type)
    return ;
  
  cleanTtSqrt (&tt) ;

  if (1 && ! ppIsNumber (pp))
    {
      if (cabs (tt.z + 1) < minAbs)
	printf (" -") ;
      else if (cabs (tt.z - 1) > minAbs)
	{ tt.z = nicePrint ("", tt.z) ; }
    }
  else
    { tt.z = nicePrint ("", tt.z) ; }
  if (tt.sqrt1 > 1 || tt.sqrt2 > 1)
    {
      printf ("sqrt(%d",tt.sqrt1) ;
      if (tt.sqrt2 > 1) printf ("/%d",tt.sqrt2) ;
      printf (")") ;
    }
  switch (tt.sqrti)
    {
    case 0 :  break ;
    case 1:  printf(" eta ") ; break ;
    }


  if (*tt.x) printf (" %s ", tt.x) ;
  if (*tt.g) gprintf ("g", tt.g) ;
  if (*tt.gg) gprintf ("gg", tt.gg) ;
  if (*tt.eps) gprintf ("epsilon", tt.eps) ;

  if (*tt.sigma) gprintf ("s", tt.sigma) ;
  if (*tt.sigB) gprintf ("sB", tt.sigB) ;

  if (tt.mm[0][0]) gprintf ("k", tt.mm[0]) ;
  if (tt.mm[1][0]) gprintf ("p", tt.mm[1]) ;
  if (tt.mm[2][0]) gprintf ("q", tt.mm[2]) ;
  if (tt.mm[3][0]) gprintf ("r", tt.mm[3]) ;

  i = 0 ;
  if      (tt.denom[0]) printf (" %c k^%d ", i++ ? ' ' : '/', 2 * tt.denom[0]) ;
  if      (tt.denom[1]) printf (" %c (k+p)^%d ", i++ ? ' ' : '/', 2 * tt.denom[1]) ;
  if      (tt.denom[2]) printf (" %c (k+p+q)^%d ", i++ ? ' ' : '/', 2 * tt.denom[2]) ;
  if      (tt.denom[3]) printf (" %c (k+p+q+r)^%d ", i++ ? ' ' : '/', 2 * tt.denom[3]) ;

  if (*tt.theta) printf(" theta(%s) ", tt.theta) ;
  return ;
} /* showTT */

/*******************************************************************************************/

void showPol (POLYNOME pp)
{
  POLYNOME p1, p2 ;
  static int nn = 0, level = 0 ;

   if (!pp)
     {
       if (level == 0) printf ("(NULL) #### zero term ###\n") ;
       return ;
     }
   if (! level) nn = 0 ;
   level++ ;
   p1 = pp->p1 ;
   p2 = pp->p2 ;
   if (p1 == pp) messcrash ("pp == pp->p1 in showPol") ;
   if (p2 == pp) messcrash ("pp == pp->p2 in showPol") ;
   if (p1 && p1 == p2) messcrash ("pp->p1 == pp->p2 in showPol") ;

   if (pp->isProduct)
     { printf ("(");  showPol (p1) ; printf (")*(");  showPol(p2) ; printf(")") ; }
   else if (pp->isSum)
     { printf ("(");  showPol (p1) ; printf (")+") ; if (p2 && !p2->isSum) printf("(");  showPol(p2) ; if (p2 && !p2->isSum) printf(")") ; }
   else if (pp->tt.type)
     { showTT (pp) ; nn++ ; }
   level-- ;

   if (! level) printf ("  ### %d term%s ###\n", nn, nn > 1 ? "s" : "") ;
} /* showPol */

/***********************************************************************************************************************************************/
/******************************************* New Polynomes of all types ************************************************************************/
/***********************************************************************************************************************************************/

POLYNOME newPolynome (AC_HANDLE h)
{
  static int id = 0 ;
  int n = sizeof (struct polynomeStruct) ;
  POLYNOME pp = (POLYNOME) halloc (n, h) ;
  memset (pp, 0, n) ;
  pp->id = ++id ;
  pp->h = h ;
  return pp ;
}

/***********************************************************************************************************************************************/

void freePolynome (POLYNOME pp)
{
  if (pp && pp->id)
    {
      pp->id = 0 ;
      ac_free (pp->p1) ;
      ac_free (pp->p2) ;
      /* do NOT       ac_free (pp->h) ; pp may be a component of a larger pp with same h */
      ac_free (pp) ;
    }
  return ;
}

/*************************************************************************************************/

POLYNOME copyPolynome (POLYNOME p1, AC_HANDLE h)
{
  POLYNOME pp = 0 ; 
 
  if (p1)
    {
      int id ;
      pp = newPolynome (h) ;
      id = pp->id ;
      memcpy (pp, p1, sizeof (struct polynomeStruct)) ;
      pp->id = id ;
      if (pp->p1)
	pp->p1 =  copyPolynome (pp->p1, h) ;
      if (pp->p2)
	pp->p2 =  copyPolynome (pp->p2, h) ;
    }
  return pp ;
}

/*************************************************************************************************/

POLYNOME newScalar (complex float z, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type = 1 ;
  p->tt.z = z ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newG (short mu, short nu, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.g[0] = mu ;
  p->tt.g[1] = nu ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newEpsilon (short a, short b, short c, short d, AC_HANDLE h)
{
  POLYNOME pp ;
  
  pp = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0 ;
  pp->tt.eps[0] = a ;
  pp->tt.eps[1] = b ;
  pp->tt.eps[2] = c ;
  pp->tt.eps[3] = d ;
  return pp ;
}

/*************************************************************************************************/
/* newAG (... ,  0) antisymmetric link 1/2(ac bd - ad bc). Optionally adding the i epsilon 
 * newAG (... , +1) is the self-dual projector P+ = 1/4 ( ac bd - ad bc) + i/2 epsilon(abcd)
 * newAG (... , -1) is the self-dual projector P- = 1/4 ( ac bd - ad bc) - i/2 epsilon(abcd)
 *     WE have (P+)^2 = (P+),   (P-)^2 = (P-),   (P+)(P-) = (P-)(P+) = 0 
 */
POLYNOME newAG (short a, short b, short c, short d, int parity, AC_HANDLE h)
{
  POLYNOME pp, ppp[4] ;

  pp = ppp[0] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = c ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;

  pp = ppp[1] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1.0/2.0 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = c ;

  ppp[2] = ppp[3] = 0 ;

  if (parity)
    {
      ppp[0]->tt.z /= 2 ;
      ppp[1]->tt.z /= 2 ;

      pp = ppp[2] = newEpsilon (a, b, c, d, h) ;
      pp->tt.z = (parity ) / 4.0 ;             /* I * parity,  in Minkovski */
      if (parity == 2 || parity == -2)
	return pp ;
    }
  
  pp = newMultiSum (h, ppp) ;
  return pp ; 
} /* newAG */

/*************************************************************************************************/

POLYNOME newAG6 (short a, short b, short c, short d, short e, short f, AC_HANDLE h)
{
  POLYNOME pp, ppp[7] ;

  pp = ppp[0] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = f ;

  pp = ppp[1] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = d ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  pp = ppp[2] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = b ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = a ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  pp = ppp[3] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = -1 ;
  pp->tt.g[0] = c ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = a ;
  pp->tt.g[3] = e ;
  pp->tt.g[4] = b ;
  pp->tt.g[5] = d ;

  pp = ppp[4] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = e ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = f ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = d ;

  pp = ppp[5] = newPolynome (h) ;
  pp->tt.type = 1 ;
  pp->tt.z = 1 ;
  pp->tt.g[0] = a ;
  pp->tt.g[1] = f ;
  pp->tt.g[2] = b ;
  pp->tt.g[3] = d ;
  pp->tt.g[4] = c ;
  pp->tt.g[5] = e ;

  ppp[6] = 0 ;
  
  pp = newMultiSum (h, ppp) ;
  return pp ; 
} /* newAG6 */

/*************************************************************************************************/

POLYNOME newK (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[0][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newP (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[1][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newQ (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[2][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newR (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.mm[3][0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newPQR (int pqr, short mu, AC_HANDLE h)
{
  POLYNOME p0 = pqr >=0 ? newK (mu, h) : 0 ;
  POLYNOME p1 = pqr >=1 ? newP (mu, h) : 0 ;
  POLYNOME p2 = pqr >=2 ? newQ (mu, h) : 0 ;
  POLYNOME p3 = pqr >=3 ? newR (mu, h) : 0 ;
  POLYNOME ppp[5] = {p0, p1, p2, p3, 0} ;

  return newMultiSum (h, ppp) ;
}

/*************************************************************************************************/

POLYNOME newSigma (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type =1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigma[0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newSigB (short cc, AC_HANDLE h)
{
  POLYNOME p = newPolynome (h) ;
  p->tt.type = 1 ;
  p->tt.z = 1 ;
  p->tt.Id2 = 1 ;
  p->tt.sigB[0] = cc ;
  return p ;
}

/*************************************************************************************************/

POLYNOME newSum (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  if (p1 && ! p1->isSum && ! p1->isProduct && ! p1->tt.z) p1 = 0 ;
  if (p2 && ! p2->isSum && ! p2->isProduct && ! p2->tt.z) p2 = 0 ;
  if (p1 && p2)
    {
      POLYNOME p = newPolynome (h) ;
      p->p1 = copyPolynome (p1, h) ;
      p->p2 = copyPolynome (p2, h) ; 
      if (p1->tt.type && p2->tt.type && p1->tt.Id2 + p2->tt.Id2 == 1)
	messcrash ("Cannot add aPauli matrix to a number\n") ;;
      p->isSum = TRUE ;
      if (p1 == p2)
	messcrash ("p1 == p2 in newSum") ;
      return p ;
    }
  else if (p1)
    return copyPolynome (p1, h) ;
  else if (p2)
    return copyPolynome (p2, h) ;
  return 0 ;
}

/*************************************************************************************************/

POLYNOME newProduct (POLYNOME p1, POLYNOME p2, AC_HANDLE h)
{
  if (p1 && ! p1->isSum && ! p1->isProduct && ! p1->tt.z) p1 = 0 ;
  if (p2 && ! p2->isSum && ! p2->isProduct && ! p2->tt.z) p2 = 0 ;
  if (p1 && p2)
    {
      POLYNOME p = newPolynome (h) ;
      p->p1 = copyPolynome (p1, h) ;
      p->p2 = copyPolynome (p2, h) ; 
      p->isProduct = TRUE ;
      return p ;
    }
  return 0 ;
}

/*************************************************************************************************/

POLYNOME polynomeScale (POLYNOME pp, double complex z)
{
  if (! pp)
    return 0 ;

  pp = copyPolynome (pp, pp->h) ;
  if (pp->isSum)
    {
      pp->p1 = polynomeScale (pp->p1, z) ;
      pp->p2 = polynomeScale (pp->p2, z) ;
    }
  else if (pp->isProduct)
    {
      pp->p1 = polynomeScale (pp->p1, z) ;
    }
  else
    pp->tt.z *= z ;
  
  return pp ;
} /* polynomeScale */

/*************************************************************************************************/

POLYNOME newMultiSum (AC_HANDLE h, POLYNOME ppp[])
{
  POLYNOME pp, p1, p2 ;
  int i = 0 ;

  while (ppp[i]) i++ ;
  if (i <= 1) return copyPolynome(ppp[0], h) ;

  pp = ppp[--i] ;
  while (i > 0) 
    {
      p2 = pp ;
      p1 = ppp[--i] ;
      pp = newSum (p1, p2, h) ;
    }
  return pp ;
}

/*************************************************************************************************/

POLYNOME newMultiProduct (AC_HANDLE h, POLYNOME ppp[])
{
  POLYNOME pp, p1, p2 ;
  int i = 0 ;

  while (ppp[i]) i++ ;
  if (i <= 1) return copyPolynome (ppp[0], h) ;

  pp = ppp[--i] ;
  while (i > 0)
    {
      p2 = pp ;
      p1 = ppp[--i] ;
      pp = newProduct (p1, p2, h) ;
    }
  return pp ;
}

/*************************************************************************************************/
/*************************************************************************************************/

