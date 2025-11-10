/*
 * sa.introns.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to intron support and discovery
 */

#include "sa.h"

/**************************************************************/

static int confirmedIntronsOrder (const void *va, const void *vb)
{
  const INTRON *up = va ;
  const INTRON *vp = vb ;
  int n ;

  /* chrom order */
  n = up->chrom - vp->chrom ; if (n) return n ;
  /* pos order */
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  /* run order */
    n = up->run - vp->run ; if (n) return n ;
  return 0 ;
} /* confirmedIntronOrder */

/**************************************************************/

static long int confirmedIntronsCompress (BigArray aaa) 
{
  long int ii, jj = 0, iMax = bigArrayMax (aaa) ;
  bigArraySort (aaa, confirmedIntronsOrder) ;
  INTRON *up, *vp, *wp ;

  if (iMax)
    for (ii = 0, jj = 0, up = bigArrp (aaa, 0, INTRON), vp = up, wp = up ; ii < iMax ; ii++, up++)
      {
	if (vp && up->a1 == vp->a1 && up->a2 == vp->a2 && up->chrom == vp->chrom && up->run == vp->run)
	  vp->n += up->n ;
	else
	  {
	    vp = wp ;
	    if (wp < up) *wp = *up ;
	    wp++ ; jj++ ;
	  }
      }
  bigArrayMax (aaa) = jj ;
  return jj ;
} /* confirmedIntronsCompress */

/**************************************************************/

void saIntronsCumulate (BigArray aaa, Array aa) 
{
  int iMax = arrayMax (aa) ;
  
  if (iMax)
    {
      long int jj = bigArrayMax (aaa) ;
      bigArrayp (aaa, jj + iMax - 1, INTRON)->n = 0 ; /*  make room */
      memcpy (bigArrp (aaa, jj, INTRON), arrp (aa, 0, INTRON), iMax * sizeof (INTRON)) ;
      confirmedIntronsCompress (aaa) ;
    }
} /* confirmedIntronsCumulate */

/**************************************************************/

void saIntronsExport (const PP *pp, BigArray aaa)
{
  long int iMax = aaa ? confirmedIntronsCompress (aaa) : 0 ;
  
  if (iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (pp->outFileName, ".introns.tsf", 0, h) ;
      INTRON *up ;
      long int ii ;
      
      aceOutf (ao, "### Intron support in tsf format: chrom__a1_a2,  run, i (format for one integer), nb of supporting reads\n") ;
      aceOutf (ao, "### Call bin/tsf -i %s -I tsf -O table -o my_table.txt to reformat this file into an excell compatible tab delimited table\n",
	       aceOutFileName (ao)
	       ) ;
      for (ii = 0, up = bigArrp (aaa, ii, INTRON) ; ii < iMax ; ii++, up++)
	{
	  aceOutf (ao, "%s__%d_%d\t%s\ti\t%d\n"
		   , dictName (pp->bbG.dict, up->chrom >> 1) + 2, up->a1, up->a2
		   , dictName (pp->runDict, up->run)
		   , up->n
		   ) ;
	}
      
      ac_free (h) ;
    }
} /* saIntronsExport */

/**************************************************************/
/**************************************************************/
/* sliding intron, clip errors in overlap r636 r452 */
void saIntronsOptimize (BB *bb, ALIGN *vp, ALIGN *wp, Array dnaG)  
{
  A_ERR *epX, *epY ;
  int x2 = vp->x2 ;
  int y1 = wp->x1 ;
  int nEx = vp->nErr ;
  int nEy = wp->nErr ;
  int i, j, nE ;
  int bestN, bestI, bestJ ;
  int dy = vp->x2 - wp->x1 + 1 ;
  int da = vp->a2 - wp->a1 ;

  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
      memset (wp, 0, sizeof (ALIGN)) ;
      wp->chain = -1 ;
      return ;
    }
  if (vp->x1 > wp->x1)
    {
      int dx = vp->x1 - wp->x1 ;
      wp->x1 = vp->x1 ;
      if (wp->a1 < wp->a2) wp->a1 += dx ;
      else wp->a1 -= dx ;
    }
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && nEx + nEy > 0)
    {  
      for (nE = 0 ; nE < nEx ; nE++)
	{  /* count the vp errors that cannot be clipped */
	  epX = arrp (vp->errors, nE, A_ERR) ;
	  int zX = epX->iShort ; /* first bad base */
	  if (zX + 1 >= y1) /* bio coords */
	    break ;
	}
      /* if we clip both ali at y1, we keep nE errors
       * let us try to find a better position
       */
      i = nE ; j = 0 ; /* next error in vp and wp */
      nE = nE + nEy ;
      bestN = nE ; bestI = i ; bestJ = j ; 
      epX = epY = 0 ;
      while (bestN && (i < nEx || j < nEy))
	{
	  int zX = x2 ;
	  if (i < nEx)
	    {
	      epX = arrp (vp->errors, i, A_ERR) ;
	      zX = epX->iShort ; /* first bad base starting from the left */
	      zX = zX  < x2 - 1 ? zX : x2 - 1 ;
	    }
	  int zY = x2 ;
	  if (j < nEy)
	    {
	      epY = arrp (wp->errors, j, A_ERR) ;
	      zY = epY->iShort  ; /* last bad base starting from the right */
	      switch (epY->type)
 		{   /* we need to adjust the coordinates */
		case INSERTION:
		  break ;
		case INSERTION_DOUBLE:
		  zY+= 1 ;
		  break ;
		case INSERTION_TRIPLE:
		  zY += 2 ;
		  break ;
		case TROU:
		case TROU_DOUBLE:
		case TROU_TRIPLE:
		  zY -- ;
		  break ;
		default:
		  break ;
		}
	      zY = zY  < x2  ? zY : x2 ;
	    }
	  if (epX && epY && epX->iLong == epY->iLong && epX->iShort == epY->iShort)
	    {
	      /* merge the 2 alignments */
	      vp->a2 = wp->a2 ;
	      vp->x2 = wp->x2 ;
	      vp->nErr += wp->nErr ;
	      vp->errors = vp->errors ? vp->errors : wp->errors ;
	      if (vp->errors && vp->nErr)
		epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
	      memset (wp, 0, sizeof (ALIGN)) ;
	      wp->chain = -1 ;
	      return ;
	    }

	  if (zX < zY)
	    { nE++ ; i++ ; }
	  else if (zX >= zY && zY < x2)
	    {
	      nE-- ;
	      if (nE < bestN)
		{ bestN = nE ; bestI = i ; bestJ = j + 1 ; }
	      j++ ;
	    }
	  else
	    break ;
	}
      if (bestJ > 0) /* if bestJ == 0, wp->x1/a1 is well positioned */
	{ /* we must clip these errors */
	  epY = arrp (wp->errors, bestJ - 1, A_ERR) ;
	  int zA = epY->iLong + ((wp->chrom & 0x1) ? -1 : 1) ;
	  int zY = epY->iShort + 1 ; /* we need to find the first exact */
	  switch (epY->type)
	    {   
	    case INSERTION:
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_DOUBLE:
	      zY += 1 ;
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_TRIPLE:
	      zY += 2 ;
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case TROU:
	      zY -- ;
	      break ;
	    case TROU_DOUBLE:
	      zY -- ;
	      zA += ((wp->chrom & 0x1) ? -1: 1) ;
	      break ;
	    case TROU_TRIPLE:
	      zY -- ;
	      zA += ((wp->chrom & 0x1) ? -2 : 2) ;
	      break ;
	    default:
	      break ;
	    }
	  wp->x1 = zY + 1 ;  /* bio coords */
	  wp->a1 = (wp->chrom & 0x1 ?  arrayMax(dnaG) - zA : zA + 1) ;
	  wp->nErr -= bestJ ;
	  for (j = 0, epY = arrp (wp->errors, 0, A_ERR) ; j < wp->nErr ; epY++, j++)
	    *epY = *(epY + bestJ) ;
	  arrayMax (wp->errors) = wp->nErr ;
	}
      if (bestI < nEx) /* we can clip the vp errors */
	{
	  epX = arrp (vp->errors, bestI, A_ERR) ;
	  int zA = epX->iLong  + ((wp->chrom & 0x1) ? +1 : -1) ;
	  int zX = epX->iShort - 1 ; /* we need to find the last exact */
	  vp->nErr -= nEx - bestI ;
	  vp->x2 = zX + 1 ;  /* bio coords */
	  vp->a2 = (vp->chrom & 0x1 ?  arrayMax(dnaG) - zA : zA + 1) ;
	  arrayMax (vp->errors) = vp->nErr ;
	}
    }

  
  BOOL isIntronDown = TRUE ;
  BOOL foundDonor = FALSE ;
  BOOL foundAcceptor = FALSE ;
  int donor = 0, acceptor = 0 ;
  BOOL isReadDown = TRUE ;
  dy = vp->x2 - wp->x1 + 1 ;
  
  /* trim if possible the vp->x2 end of the first exon on a known 'donor' */
  if (vp->a1 < vp->a2)
    {
      isReadDown = TRUE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;
      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; isIntronDown = FALSE ; }
      else if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ;  isIntronDown = FALSE ; }
      if (donor > acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 >= donor && vp->a2 <= donor + dy - 1) 
	{  /* move back to the canonical donor site */
	  dy = vp->a2 - donor + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;
	}
      if (donor && vp->a2 == donor - 1)
	foundDonor = TRUE ;
    }
  else
    {
      isReadDown = FALSE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;

      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; isIntronDown = FALSE ; }
      if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ;  isIntronDown = FALSE ; }
      if (donor && donor < acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 <= donor && vp->a2 >= donor - dy + 1)
	
	{  /* move back to the canonical donor site */
	  dy = donor - vp->a2 + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	}      
      if (donor && vp->a2 == donor + 1)
	foundDonor = TRUE ;
    }

  /* alternativelly  trim if possible the wp->x2 start of the second  exon on a known 'acceptor' */
  dy = vp->x2 - wp->x1 + 1 ;
  if (isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 < acceptor + 1 && wp->a1 + dy >= acceptor + 1)
	{
	  dy = acceptor + 1 - wp->a1 ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	}
      if (wp->a1 == acceptor + 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 -= dy ;
	    }
	}
    }
  else if (!isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 > acceptor - 1 && wp->a1 - dy <= acceptor - 1)
	{
	  dy = wp->a1 - acceptor + 1 ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;
	}
      if (wp->a1 == acceptor - 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 += dy ;
	    }
	}
    }

  BOOL gt_ag = FALSE ;
  BOOL gc_ag = FALSE ;
  BOOL ct_ac = FALSE ;
  BOOL ct_gc = FALSE ;

  /* no annotated junction, look for gt_ag */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && isReadDown && ! foundDonor && !foundAcceptor  && wp->a1 > 3 && vp->a2 > dy)  
    {
      unsigned char *cp = arrp (dnaG, vp->a2 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, wp->a1 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;

      if (bb->nIntronSupportPlus >= bb->nIntronSupportMinus)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { gt_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { ct_ac = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { gc_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { ct_ac = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { gt_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { ct_gc = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	}
    ok1:
      if (bestI >= 0)
	{
	  bestI = (bestI ? bestI : (dy+1) / 2) ; /* avod zero to keep vp oriented */
	  dy -= bestI ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;

	  dy = bestI ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	}	    
    }
  else if (dy > 0 && !isReadDown && ! foundDonor && !foundAcceptor && wp->a1 > dy && vp->a2 > 3) 
    {
      /* move backwards on the genome */
      unsigned char *cp = arrp (dnaG, wp->a1 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, vp->a2 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;
      if (bb->nIntronSupportPlus < bb->nIntronSupportMinus)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	}
    ok2:
      if (bestI >= 0)
	{
	  dy -= bestI ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;

	  dy = bestI ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	}	    
    }
  
  /* In any case trim the second exon a bouts francs */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0)
    { 
      wp->x1 += dy ;
      wp->a1 += (wp->a1 < wp->a2 ? dy : -dy) ;
    }

  dy = vp->x2 - wp->x1 + 1 ;
  da = vp->a2 - wp->a1 ;

  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
      memset (wp, 0, sizeof (ALIGN)) ;
      wp->chain = -1 ;
      return ;
    }

  /* register exact intron support */
  if (isReadDown && donor == vp->a2 + 1 && acceptor == wp->a1 - 1)
    {  /* this intron is confirmed */
      INTRON *zp = arrayp (bb->confirmedIntrons, arrayMax (bb->confirmedIntrons), INTRON) ;
      zp->chrom = vp->chrom ;
      zp->run = bb->run ;
      zp->n = 1 ;
      zp->a1 = isIntronDown ? donor : acceptor ;
      zp->a2 = isIntronDown ? acceptor : donor ;
      bb->nIntronSupportPlus += (isIntronDown ? 1 : 0) ;
      bb->nIntronSupportMinus += (isIntronDown ? 0 : 1) ;
    }
  else if (! isReadDown && donor == vp->a2 - 1 && acceptor == wp->a1 + 1)
    {  /* this intron is confirmed */
      INTRON *zp = arrayp (bb->confirmedIntrons, arrayMax (bb->confirmedIntrons), INTRON) ;
      zp->chrom = vp->chrom ;
      zp->run = bb->run ;
      zp->n = 1 ;
      zp->a1 = isIntronDown ? acceptor : donor ;
      zp->a1 = isIntronDown ? donor : acceptor ;
      bb->nIntronSupportPlus += (isIntronDown ? 0 : 1) ;
      bb->nIntronSupportMinus += (isIntronDown ? 1 : 0) ;
    }
  gt_ag = gc_ag ^ gt_ag ^ ct_ac ^ ct_gc ; /* for compiler happiness */

} /* alignOptimizeIntron */

/**************************************************************/
/**************************************************************/
/**************************************************************/
