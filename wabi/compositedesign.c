#include "ac.h"
#include "cdna.h"
#include "makemrna.h"

#define ARRAY_CHECK
#define MALLOC_CHECK

/**********************************************************************************/
/**********************************************************************************/
/* This strategy is design to handle deep sequencing
 * We construct an oriented graph where each exon is a vertex
 * Each intronis an arc weighted by the number and samples of the supporting tags
 * The introns are ordered by weight
 * For each intron in decreasing order extend it into a path as follows
 *   If the next intron is not yet in a path, extend to the following exon connected by the heaviest intron
 *   If it is attach to the lowest path arriving at this exon from the other side
 *   unless the way i arrive on the exon is much smaller than the way the best path exits the exon
 *   or enough full paths have been constructed
 * Stop when all introns are used or enough paths have been constructed 
 */
typedef struct compositeDesignStruct { Array exons, introns, covering ; BitSet bbIn, bbOut ; KEYSET reads ; int c1, c2 ; AC_HANDLE h ; } DS ;
typedef struct dsVertexStruct { int a1, a2, donor, acceptor, cover, end, path, type, score, nn, clipable, flag ; KEY endKey ; } DSX ;

/**********************************************************************************/

static void showDsFlags (SC *sc, Array aa, int flag)
{
  DSX *up ;
  int ii ;
  int mask = (1 << 16) - 1 ;
  if (aa)
    {
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	{
	  if ((up->flag & flag) != flag)
	    continue ;

	  if (0) if (up->a1 > 6000 || up->a1 < 6300) continue ;
	  
	  fprintf( stderr, "%d:: %d %d\t %d %d\t d=%d a=%d c=%d s=%d e=%d path=%d\t"
		   , ii, up->a1, up->a2, sc->a1 + up->a1 - 1, sc->a1 + up->a2 - 1
		   , up->donor, up->acceptor, up->cover, up->score, up->end, up->path
		   ) ;
	  if (flag) fprintf( stderr, "flags %d:%d\t", up->flag & mask, up->flag >> 16) ;
	  if (gGene & up->type) fprintf (stderr, "Gene ") ;
	  if (gGap & up->type) fprintf (stderr, "Gap ") ;
	  if (gLink & up->type) fprintf (stderr, "Link ") ;
	  if (gGhost & up->type) fprintf (stderr, "Ghost ") ;
	  if (gFuseToGhost & up->type) fprintf (stderr, "FuseToGhost ") ;
	  if (gDroppedGene & up->type) fprintf (stderr, "DroppedGene ") ;
	  if (gSuspect & up->type) fprintf (stderr, "Suspect ") ;
	  if (gMicro & up->type) fprintf (stderr, "Micro-") ;
	  if (gX & up->type) fprintf (stderr, "Exon ") ;
	  if (gI & up->type) fprintf (stderr, "Intron ") ;
	  if (gJ & up->type) fprintf (stderr, "Fuzzy ") ;
	  if (gS & up->type) fprintf (stderr, "SL ") ;
	  if (gS0 & up->type) fprintf (stderr, "SL0 ") ;
	  if (gReal5p & up->type) fprintf (stderr, "r5p ") ;
	  if (g5 & up->type) fprintf (stderr, "X5 ") ;
	  if (g3 & up->type) fprintf (stderr, "X3 ") ;
	  if (gD & up->type) fprintf (stderr, "Debut ") ;
	  if (gDF & up->type) fprintf (stderr, "DebutFuzzy ") ;
	  if (gCompleteCDS & up->type) fprintf (stderr, "CDS ") ;
	  if (gF & up->type) fprintf (stderr, "Fin ") ;
	  if (gFF & up->type) fprintf (stderr, "FinFuzzy ") ;
	  if (gReal3p & up->type) fprintf (stderr, "r3p ") ;
	  if (gA & up->type) fprintf (stderr, "polyA ") ;
	  if (gB & up->type) fprintf (stderr, "Alter ") ;
	  if (gCompleteCDS & up->type) fprintf (stderr, "CDS ") ;
	  if (gStolen & up->type) fprintf (stderr, "Stolen ") ;
	  if (gPredicted & up->type) fprintf (stderr, "Predicted ") ;

	  if (up->clipable) fprintf (stderr, " *") ;
	  
	  fprintf (stderr, "\n") ;
	}
    }
  return ;
} /* showDsFlags */

static void showDs (SC *sc, Array aa)
{
  return showDsFlags (sc, aa, 0) ;
} /* showDs */

/**********************************************************************************/

static int dsScoreOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->score - b->score ; if (n) return -n ; /* large scores first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsScoreOrder */

/**********************************************************************************/

static int sFlagsOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->flag - b->flag ; if (n) return -n ; /* large scores first */
  n = a->a1 - b->a1 ; if (n) return n ;
  
  return 0 ;
} /* sFlagsOrder */

/**********************************************************************************/

static int dsCoverOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->cover - b->cover ; if (n) return -n ; /* large scores first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsCoverOrder */

/**********************************************************************************/

static int dsEndOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->end - b->end ; if (n) return -n ; /* large cover first */
  n = a->cover - b->cover ; if (n) return -n ; /* large cover first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsEndOrder */

/**********************************************************************************/

static int dsA1Order (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;
  
  return 0 ;
} /* dsA1rder */

/**********************************************************************************/

static int ssHappyFew (Array ss)
{
  int iMax = arrayMax (ss) ;
  if (iMax)
     {
       int ii, jj ;
       DSX *ssp, *ssp1 ;

       arraySort (ss, dsA1Order) ;
       for (ii = jj = 0, ssp1 = ssp = arrp (ss, 0, DSX) ; ii < iMax ; ii++, ssp++)
	 {
	   if (ssp->cover > 0)
	     {
	       if (jj < ii)
		 *ssp1 = *ssp ;
	       jj++ ; ssp1++ ;
	     }
	 }
       iMax = arrayMax (ss) = jj ;
     }

  arraySort (ss, dsA1Order) ;
  return iMax ;
} /* ssHappyFew */

/**********************************************************************************/
/* Rationalize the exon flags */
static void mrnaDesignCleanExons (SC *sc, Array ss)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;
  int ii ;
  DSX *up ;
  int iiMax = ssHappyFew (ss) ;


  /* create and boost some scores */
  for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
    {
      up->nn = ii ;
      up->score = up->cover ;
      if (up->type & gI) up->score *= 3 ;
      if (up->type & gA) up->score *= 3 ;
    }

  if (iiMax == 1)   
    {       /* eliminate weak short cloud exons */
      up = arrayp (ss,0, DSX) ;
      int da = up->a2 > up->a1 ? up->a2 - up->a1 : up->a1 - up->a2 ;
      if (da < 50 && up->score < 10)
	iiMax = arrayMax (ss) = 0 ;
    }

  /* if we have   2 introns with same acceptor
   * the exon bit going to the second donor
   * should have the support of the short intron
   */
  if (iiMax)
    for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
      {
	if (up->type & gI)
	  {
	    DSX *vp ;
	    int jj, ok2 = 0 ;
	    
	    for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->a2 < up->a2 && vp->a1 == vp[-1].a2 + 1   ; jj++, vp++)
	      if (vp->type & gI)
		{ ok2 = TRUE ; break ; }
	    if (ok2 && vp->a1 < up->a1 + 30)
	      for (int i = 1 ; ii + i < jj ; i++)
		{
		  if (up[i].score < vp->score)
		    up[i].score = vp->score ;
		}
	  }
      }

  /* polyA are saved */
  if (0 && iiMax)
    for (ii = iiMax - 1, up = arrp (ss, ii, DSX) ; ii >= 0 ; ii--, up--)
      {
	if (up->type & gA)
	  { 
	    DSX *vp ;
	    int i, jj, bestScore = up->score ;
	    int b2 = up->a1 - 1 ;
	    for (i = ii - 1, vp = arrp (ss, i, DSX) ; i >= 0 ; vp--, i--)
	      {
		if (vp->type & gI)
		  {
		    if (vp->a1 == b2 + 1)
		      break ;
		    continue ;
		  }
		if (vp->a2 != b2) break ;
		b2 = vp->a1 - 1 ;
		if (vp->score > bestScore) bestScore = vp->score ;
	      }
	    jj = i ;
	    for (i = jj + 1, vp = arrp (ss, i, DSX) ; i <= ii ; i++, vp++)
	      if (vp->type & gX)
		vp->score = bestScore ;
	  }
      }  
  
  /* overlapping elements are killed if below 1% */
  iiMax = ssHappyFew (ss) ;
  if (1)
    {
      int rho = 100 ; /* rejection rate */
      /* exact overlaps of exons introns */
      for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax - 1 ; ii++, up++)
	{
	  DSX *vp = up + 1 ;
	  if (vp->a1 == up->a1 && vp->a2 == up->a2)
	    {
	      if (rho * vp->cover < up->cover)
		vp->cover = vp->score = 0 ;
	      if (rho * up->cover < vp->cover)
		up->cover = up->score = 0 ;
	    }
	}

      iiMax = ssHappyFew (ss) ;
      /* generalized cassettes */
      for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax - 1 ; ii++, up++)
	{
	  int jj, kk ;
	  DSX *vp, *wp ;
	  for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->a1 == up->a1 ; jj++, vp++)
	    if (up->cover && vp->cover)
	      { /* vp is necesarilly longer than up try to extend up */
		int b2 = up->a2 ;
		int cover2 = up->cover ;
		int ns = 0 ;
		KEYSET cs = 0 ;
		for (kk = jj + 1, wp = vp + 1 ; kk  < iiMax && wp->a1 < b2 + 2 ; kk++, wp++)
		  {
		    if (wp->cover && wp->a1 == b2 + 1 && wp->a2 <= vp->a2)
		      {
			if (wp->cover > cover2 && ! (wp->type & (gReal5p | gReal3p)))
			  cover2 = wp->cover ;
			b2 = wp->a2 ;
			if (! cs) cs = keySetCreate () ;
			keySet (cs, ns++) = kk ;
		      }
		  }
		if (b2 == vp->a2)
		  { /* cassette found */
		    if (cover2 > rho * vp->cover)
		      vp->cover = vp->score = 0 ;
		    if (rho * cover2 < vp->cover)
		      {
			up->cover = up->score = 0 ;
			for (int m = 0 ; m < keySetMax (cs) ; m++)
			  {
			    wp = arrp (ss, keySet(cs, m), DSX) ;
			    wp->cover = wp->score = 0 ;
			  }
		      }
		  }
		keySetDestroy (cs) ;
	      }
	}
    }
  /* flag  gCompleteCDS introns: those linking gCompleteCDS exons */
  if (0)
    {
#ifdef JUNK
      for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
	{
	  if ((vp->type & gX) && (vp->type & gCompleteCDS))
	    {
	      int a21 = vp->a2 + 1 ;
	      for (jj = ii + 1, vp2 = vp +1 ; jj < eeMax && vp2->a1 <= a21 ; jj++, vp2++)
		if ((vp2->type & gI) && vp2->a1 == a21)
		  vp2->type |= gCompleteCDS ; /* flag, the donor is in a CDS */
	    }
	  if ((vp->type & gI) && (vp->type & gCompleteCDS))
	    {
	      int ok = 0, a21 = vp->a2 + 1 ;
	      for (jj = ii + 1, vp2 = vp +1 ; !ok && jj < eeMax && vp2->a1 <= a21 ; jj++, vp2++)
		if ((vp2->type & gX) && (vp2->type & gCompleteCDS) && vp2->a1 == a21)
		  ok = 1 ;
	      if (! ok) /* kill the flag, the intron acceptor does not hook to a CDS */
		vp2->type &= (~ gCompleteCDS) ;
	    }
	}
#endif
    }

  iiMax = ssHappyFew (ss) ;
  /* number the segments */
  for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
    up->nn = ii ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignCleanExons\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignCleanExons done\n") ;
    }
  ac_free (h) ;
   return ;
} /* mrnaDesignCleanExons */

/**********************************************************************************/

static void mrnaDesignGetElements (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii = 0, jje = 0, jji = 0, j, j2, j2Max, a1, *ip ;
  SMRNA *smrna ;
  Array introns, exons ;
  DSX *dsx, *dsy ;
  HIT *up ;
  Array ks ;
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;
  
  exons = ds->exons = arrayHandleCreate (100, DSX, ds->h) ;
  introns = ds->introns = arrayHandleCreate (100, DSX, ds->h) ;
  ds->reads = keySetHandleCreate (ds->h) ;

  if (1)
    {
      /* extract all exon fragments and all introns */
      if (smrnas && arrayMax(smrnas))
	for (ii = 0 ; ii < arrayMax(smrnas) ; ii++)
	  {
	    smrna = arrp (smrnas, ii, SMRNA) ;
	    if (! smrna->clones)
	      continue ;
	    for (j = 0, up = arrp(smrna->hits, 0, HIT) ; j < arrayMax(smrna->hits) ; up++, j++)
	      {
		if (up->a2 > up->a1 - 1)
		  if ((gX | gI) & up->type) /* exon */
		    {
		      if (gX & up->type) /* exon */
			dsx = arrayp (exons, jje++, DSX) ;
		      if (gI & up->type) /* intron */
			dsx = arrayp (introns, jji++, DSX) ;
		      dsx->a1 = up->a1 + smrna->a1 - 1 ;
		      dsx->a2 = up->a2 + smrna->a1 - 1 ;
		      dsx->type = up->type & (gX | gI) ;
		    }
	      }
	  }  
    }

  if (1)
    {
      /* include the Sl1 */

    }

    if (debug)
     {
       fprintf (stderr, "... getElements exons Z\n") ;
       showDs (sc, exons) ;
       fprintf (stderr, "... getElements exons Z done\n") ;
     }
   arraySort (exons, dsA1Order) ;
   arrayCompress (exons) ;
   if (debug)
     {
       fprintf (stderr, "... getElements exons A\n") ;
       showDs (sc, exons) ;
       fprintf (stderr, "... getElements exons A done\n") ;
     }
   /* now we split the exons in subparts */ 
   ks = arrayHandleCreate (2 * arrayMax (exons) + 1, int, h) ;
   for (j = j2 = 0, dsx = arrayp (exons, j, DSX) ; j < arrayMax (exons) ; j++, dsx++)
     {
      array (ks, j2++, int) = dsx->a1 ;
      array (ks, j2++, int) = dsx->a2 + 1 ;
     }
   arraySort (ks, intOrder) ;
   arrayCompress (ks) ;
   j2Max = keySetMax (ks) ;
   ds->exons = arrayHandleCreate (100, DSX, ds->h) ;
   for (jje = j = j2 = 0, dsx = arrayp (exons, j, DSX), ip = arrp (ks, 0, int) ; j < arrayMax (exons) ; j++, dsx++)
     {
       a1 = dsx->a1 ;
       while (j2 > 0 && *ip > a1) { j2-- ; ip-- ; }
       while (*ip < dsx->a1) { j2++; ip++ ;}
       while (++j2 < j2Max && *(++ip) <= dsx->a2 + 1)
	 {
	   dsy = arrayp (ds->exons, jje++, DSX) ;
	   dsy->a1 = a1 ; dsy->a2 = *ip - 1 ; dsy->type = gX ; a1 = *ip ;
	 }
     }

   arraySort (ds->exons, dsA1Order) ;
   arrayCompress (ds->exons) ;
   if (debug) 
     {
       fprintf (stderr, "... getElements exons final\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "... getElements exons final done\n") ;
     }

   arraySort (introns, dsA1Order) ;
   arrayCompress (introns) ;
   if (debug)
     {
       fprintf (stderr, "... getElements introns final\n") ;
       showDs (sc, introns) ;
       fprintf (stderr, "... getElements introns final done\n") ;
     }
   ac_free (h) ;
   return ;
} /* mrnaDesignGetElements */

/**********************************************************************************/
/* get the X Composite and their counts */
static BOOL mrnaDesignGetGraphElements (DS *ds, S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  BOOL debug = FALSE ; /* tyty TRUE ; */
  HIT *up, *up1 ;
  DSX *ssp ;
  int i1, ii, jj, iMax = 0, ir, iss = 0 ;
  int nXI = 0, nXE = 0, nXA = 0, nXSL = 0, nXends = 0 ;
  BOOL isDown ;
  Array units ;
  BSunit *uu ;
  OBJ Est = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  Array ss = arrayHandleCreate (256, DSX, ds->h) ;
  KEYSET ks1 = keySetHandleCreate (h) ;
  KEYSET ks2 = keySetHandleCreate (h) ;
  KEYSET ks5 = keySetHandleCreate (h) ;
  KEYSET ks3 = keySetHandleCreate (h) ;
  int iiMax = arrayMax ( s2m->plainHits) ;
  int ii1Max = arrayMax ( gmrna->hits) ;
  BOOL hasGoody = FALSE ;
  
  /* grab the introns in the plainHits EST alignments */
  isDown = (sc->a1 < sc->a2) ? TRUE : FALSE ;
  units = arrayHandleCreate (12, BSunit, h) ;

  for (ii = 0, up = arrp (s2m->plainHits, ii, HIT) ; ii < iiMax ; ii++, up++)
    for (i1 = 0, up1 = arrp (gmrna->hits, 0, HIT) ; i1 < ii1Max ; i1++, up1++)
      {
	int b1, b2 ;
	
	if (sc->a1 < sc->a2)
	  {
	    b1 = up1->a1 + sc->a1 - 1 ;
	    b2 = up1->a2 + sc->a1 - 1 ;
	  }
	else
	  {
	    b2 = - up1->a1 + sc->a1 + 1 ;
	    b1 = - up1->a2 + sc->a1 + 1 ;
	  }
	
	if (
	    ((up1->type & gI) && (up->type & gI) && up->a1 == b1 && up->a2 == b2) ||
	    ((up1->type & gX) && (up->type & gX) && up->a1 < b2 && up->a2 > b1)
	    )
	    {
	      BOOL flipped = up->x1 > up->x2 ? TRUE : FALSE ;
	      int a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
	      int a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
	      int xdx = up->x2 - up->x1 ;
	      int cdx = up->clipEnd - up->clipTop ;
	      if (cdx < 0) cdx = - cdx ;
	      if (xdx < 0) xdx = - xdx ;
	      
	      if (0 && debug)
		{
		  fprintf (stderr, "......getSS ii=%d a1=%d a2=%d %s \n", ii, a1, a2, name (up->est)) ;
		  showDs (sc, ss) ;
		  fprintf (stderr, "......getSS done\n") ;
		}
	      if ( ! strncmp (name(up->est), "XI_", 3) ||
		   ! strncmp (name(up->est), "XW_", 3) ||
		   ! strncmp (name(up->est), "XY_", 3)
		   )
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (up->type & gX)
		    { /* reclip the exon borders */
		      int dx = up->x2 - up->x1 - 18 ;
		      if (dx > 0)
			{
			  if (up->x1 == 1)
			    {
			      up->x1 += dx ;
			      if (up->a1 < up->a2) up->a1 += dx ;
			      else up->a1 -= dx ;
			    }
			  else if (up[1].est != up->est)
			    {
			      up->x2 -= dx ;
			      if (up->a1 < up->a2) up->a2 -= dx ;
			      else up->a2 += dx ;
			    }
			  a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
			  a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
			}
		    }
		    
		  if ((Est = bsCreate (up->est)))
		    {
		      int s = 0 ;
		      KEY intron = 0 ;
		      if (bsGetArray (Est, _Composite, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    s =  uu[0].i ; 
			  }
		      if (s && bsGetArray (Est, _Intron, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    intron = uu[0].k ;
			    if (keyFindTag (intron, _Other))
			      s /= 5 ;
			    if (! strncmp (name(up->est), "XW_", 3))
			      s *= 3 ;
			  }
		      
		      if (s > 0)
			{
			  DSX *ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ;
			  hasGoody = TRUE ;
			  s *= 1 ;
			  if (up->type & gI)
			    { nXI++; ssp->type = gI ; ssp->cover = s ; }
			  else if (up->type & gX)
			    { nXE++; ssp->type = gX ; ssp->cover = s ; ssp->clipable = 2 ; }
			}
		      
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	      
	      /* register the exon support and  shrink the higher expressed segments */
	      if ( ! strncmp (name(up->est), "XG_", 3) ||
		   ! strncmp (name(up->est), "XH_", 3)
		   )
		{
		  int iss0 = iss, b1, b2 ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if ((Est = bsCreate (up->est)))
		    {
		      if (bsGetArray (Est, _Composite, units, 3) && arrayMax (units) >= 3)
			{
			  if (up->reverse == flipped)
			    for (ir = 0 ; ir < arrayMax (units) ; ir += 3)
			      {
				uu = arrp (units, ir, BSunit) ;
				b1 = a1 + uu[0].i - 1 ; 	
				b2 = a1 + uu[1].i - 1 ;
				
				ssp = arrayp (ss, iss++, DSX) ;
				ssp->type = gX ;
				ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
				if (ssp->cover > 100)
				  hasGoody = TRUE ;
				if (ir == 0 || ir == arrayMax (units) - 3)
				  ssp->clipable = 1 ;
				if (0 && debug)
				  fprintf(stderr, "%s %d %d %d %s %s \n", name(up->est), ssp->a1, ssp->a2, ssp->cover, up->reverse ? "reverse" : "forward", isDown ? "Down":  " Up") ;
			      }
			  else
			    for (ir = arrayMax (units) - 3 ; ir >= 0 ; ir -= 3)
			      {
				uu = arrp (units, ir, BSunit) ;
				b1 = a2 - uu[1].i + 1 ; 	
				b2 = a2 - uu[0].i + 1 ;
				
				ssp = arrayp (ss, iss++, DSX) ;
				ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
				if (ssp->cover > 100)
				  hasGoody = TRUE ;
				if (ir == 0 || ir == arrayMax (units) - 3)
				  ssp->clipable = 1 ;
				if (0 && debug)
				  fprintf(stderr, "----%s %d %d %d %s %s \n", name(up->est), ssp->a1, ssp->a2, ssp->cover, up->reverse ? "reverse" : "forward", isDown ? "Down":  " Up") ;
			      }
			}
		      
		      bsDestroy (Est) ;
		      
		      /* shrink the higher expressed segments */
		      if ( ! strncmp (name(up->est), "XG_", 3) || ! strncmp (name(up->est), "XH_", 3))
			{
			  int ii ;
			  for (ii = iss0 ; ii < iss - 1 ; ii++)
			    {
			      DSX *ssq ;
			      ssp = arrayp (ss, ii, DSX) ;
			      ssq = ssp + 1 ;
			      if (ssp->a1 < ssp->a2 && ssp->a2 == ssq->a1 - 1)
				{
				  if (ssp->cover < ssq->cover && ssq->a2 > ssq->a1 + 10)
				    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
				  if (ssp->cover > ssq->cover && ssp->a2 > ssp->a1 + 10)
				    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
				}
			      if (ssp->a1 > ssp->a2 && ssp->a2 == ssq->a1 + 1)
				{
				  messcrash ("a1>a2 should never happen") ;
				  if (ssp->cover < ssq->cover && ssq->a1 > ssp->a1 + 10)
				    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
				  if (ssp->cover > ssq->cover && ssq->a2 > ssq->a1 + 10)
				    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
				}
			    }
			  if (0)
			    {
			      ssp = arrayp (ss, iss0, DSX) ;
			      ssp->a1 += (ssp->a1 < ssp->a2 ? +10 : -10) ;
			      ssp = arrayp (ss, iss - 1, DSX) ;
			      ssp->a2 += (ssp->a1 < ssp->a2 ? -10 : +10) ;
			      if (ssp->a1 > ssp->a2) ssp->cover = 0 ;
			    }
			  for (ii = iss0 ; ii < iss ; ii++)
			    {
			      ssp = arrayp (ss, ii, DSX) ;
			      if (0 && debug)
				fprintf(stderr, "... %s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
			    }
			}
		      continue ;
		    }
		}
	      
	      /* register the predicted cds support */
	      if (! strncmp (name(up->est), "Xcds_", 5))
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if ((Est = bsCreate (up->est)))
		    {
		      int s = 12 ;
		      nXE++ ;
		      ssp = arrayp (ss, iss++, DSX) ;
		      ssp->a1 = a1 ;
		      ssp->a2 = a2 ;
		      ssp->cover = s ;
		      ssp->type = gX | gCompleteCDS | gFF ;
		      ssp->clipable = 1 ;
		      
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	      
	      /* register the polyA */
	      if (! strncmp (name(up->est), "XA_", 3) ||
		  ! strncmp (name(up->est), "XSL", 3) ||
		  ! strncmp (name(up->est), "Xends_", 6)
		  )
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if ((Est = bsCreate (up->est)))
		    {
		      int s = 0 ;
		      if (bsGetArray (Est, _Composite, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    s =  uu[0].i ; 
			  }
		      if (s)
			{
			  if (! strncmp (name(up->est), "Xends_ELF_", 9) && ! isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ERF_", 9) && ! isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ERR_", 9) &&  isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ELR_", 9) &&  isDown)
			    continue ;
			  
			  hasGoody = TRUE ;
			  ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ;
			  ssp->cover = s ;
			  ssp->end = s ;
			  ssp->endKey = up->est ;
			  if (! strncmp (name(up->est), "XA_", 3))
			    {
			      nXA++ ;
			      ssp->a1 = ssp->a2 ;
			      ssp->type = gX | gA | g3 ;
			      ssp->clipable = FALSE ;
			    }
			  else if (! strncmp (name(up->est), "XSL", 3))
			    {
			      nXSL++ ;
			      ssp->a2 = ssp->a1 ;
			      ssp->type = gX | gS | g5 | gReal5p ;
			      ssp->clipable = FALSE ;
			    }
			  else
			    {
			      ssp->cover = 1 ; 
			      ssp->clipable = TRUE ;
			      
			      if (! strncmp (name(up->est), "Xends_ELF.", 10) && isDown) 
				{ nXends++ ; ssp->type = (gReal5p | gDF | gX) ; if (ssp->a2 > ssp->a1 + 30) { ssp->a1 += 0 ; ssp->a2 -= 0*30 ; }}
			      else if (! strncmp (name(up->est), "Xends_ERF.", 10) && isDown) 
				{ nXends++ ; ssp->type = (gReal3p | gFF | gX) ; if (ssp->a2 < ssp->a1 - 30) { ssp->a2 -= 0 ; ssp->a2 += 0*30 ; }} 
			      else if (! strncmp (name(up->est), "Xends_ERR.", 10) && ! isDown) 
				{ nXends++ ; ssp->type = (gReal5p | gDF | gX) ; if (ssp->a2 < ssp->a1 - 30) { ssp->a1 += 0 ; ssp->a2 += 0*30 ; }} 
			      else if (! strncmp (name(up->est), "Xends_ELR.", 10) && ! isDown) 
				{ nXends++ ; ssp->type = (gReal3p | gFF | gX) ; if (ssp->a2 > ssp->a2 + 30) { ssp->a1 -= 0 ; ssp->a2 -= 0*30 ; }}
			    }
			}
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	    }
      }

  if (debug)
    {
       fprintf (stderr, "mrnaDesignGetGraphElements\n") ;
       showDs (sc, ss) ;
       fprintf (stderr, "mrnaDesignGetGraphElements  done\n") ;
    }
  if (! hasGoody)
    { ac_free (h) ; return FALSE ;}
  arraySort (ss, dsA1Order) ;
  arrayCompress (ss) ;
  if (debug)
    {
       fprintf (stderr, "mrnaDesignGetGraphElements\n") ;
       showDs (sc, ss) ;
       fprintf (stderr, "mrnaDesignGetGraphElements  done\n") ;
       fprintf (stderr, "s2m->plainHits\n") ;
       showHits (s2m->plainHits) ;
       fprintf (stderr, "s2m->plainHits  done\n") ;
       fprintf (stderr, "gmrna->hits\n") ;
       showHits (gmrna->hits) ;
       fprintf (stderr, "gmrna->hits  done\n") ;
     }

  iMax = ssHappyFew (ss) ;
  if (iMax) /* kill overlapping introns below 1/1000 */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
         {
	   int i ;
	   DSX *ssp1 ;

	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       int b2 = ssp->a2 ;
	       for (i = 1, ssp1 = ssp + i ; i + jj < iMax && ssp1->a1 <= b2 ; i++, ssp1++)
		 {
		   if (ssp1->type & gI)
		     {
		       if (ssp->cover / 1000 > ssp1->cover)
			 ssp1->cover = 0 ;
		       if (ssp1->cover / 1000 > ssp->cover)
			 ssp->cover = 0 ;
		     }
		 }
	     }
	 }
     }
  
  iMax = ssHappyFew (ss) ;
  if (! iMax)
    goto done ;

  /* clip end elements overapping an intron */
  for (int ii = 0 ; ii < iMax ; ii++)
    {
      ssp = arrp (ss, ii, DSX) ;
      if (ssp->type & gReal5p)
	{
	  for (int jj = 0 ; jj < iMax ; jj++)
	    {
	      DSX *ssp1 = arrp (ss, jj, DSX) ;
	      if (ssp1->type & gI)
		if (ssp1->a1 > ssp->a1 + 20 && ssp1->a1 < ssp->a2)
		  ssp->a2 = ssp1->a1 - 1 ;
	    }
	}
      if (ssp->type & gReal3p)
	{
	  for (int jj = 0 ; jj < iMax ; jj++)
	    {
	      DSX *ssp1 = arrp (ss, jj, DSX) ;
	      if (ssp1->type & gI)
		if (ssp1->a2 > ssp->a1 && ssp1->a2 < ssp->a2 - 20)
		  ssp->a1 = ssp1->a2 + 1 ;
	    }
	}
    }
    
  /*  register all boundaries */
  keySet (ks1, 0) = 0 ;
  keySet (ks2, 0) = 0 ;
  for (int ii = 0, m = 1, n = 1 ; ii < iMax ; ii++)
    {
      ssp = arrp (ss, ii, DSX) ;
      keySet (ks1, m++) = ssp->a1 ;         /* begin */
      keySet (ks2, n++) = ssp->a2 ;         /* end */
      keySet (ks1, m++) = ssp->a2 + 1 ;     /* virtual begin */
      keySet (ks2, n++) = ssp->a1 - 1 ;     /* virtual end */
      if (ssp->type & gI)
	{
	  keySet (ks5, m++) = ssp->a2 + 1 ;     /* exon acceptor */
	  keySet (ks3, n++) = ssp->a1 - 1 ;     /* exon donor */
	}
    }
  keySetSort (ks1) ;
  keySetSort (ks2) ;
  keySetSort (ks5) ;
  keySetSort (ks3) ;
  keySetCompress (ks1) ;
  keySetCompress (ks2) ;
  keySetCompress (ks5) ;
  keySetCompress (ks3) ;
  
  /* split when needed */
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements before split\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements before split done\n") ;
    }
  if (1)
    {
      int m = 0, n = 0, b1 ;
      int mMax = keySetMax (ks1), nMax = keySetMax (ks2) ;
      Array ss1 = arrayHandleCreate (2*ii + 1 , DSX, ds->h ) ;
      int j5Max = keySetMax (ks5) ;
      int j3Max = keySetMax (ks3) ;
      
      for (int ii = 0, jj = 0 ; ii < iMax ; ii++)
	  {
	    int i4 = 0 ;
	    ssp = arrp (ss, ii, DSX) ;
	    if (ssp->type & gX)
	      {
		while (m > 0 && keySet (ks1, m) >= ssp->a1)
		  m-- ;
		while (n > 0 && keySet (ks2, n) >= ssp->a1)
		  n-- ;
		while (m < mMax  && keySet (ks1, m) < ssp->a1)
		  m++ ;
		while (n < nMax && keySet (ks2, n) < ssp->a1)
		  n++ ;
		b1 = ssp->a1 ;
		while (n < nMax && keySet (ks2, n) <= ssp->a2)
		  {
		    if (1 ||  /* otherwise in 0_923 je cree un trou avec l'exon au dessus */
			i4++ || (! ssp->clipable || keySet (ks2, n)  > b1 + 10))
		      {
			DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
			*ssp1 = *ssp ;
			ssp1->a1 = b1 ; ssp1->a2 = keySet (ks2, n) ;
		      }
		    b1 = keySet (ks2, n)  + 1 ;
		    n++ ;
		  }
		if (n == nMax) n-- ;
		if (b1 < ssp->a2 && (! ssp->clipable || b1 < ssp->a2 - 10))
		  {
		    DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
		    *ssp1 = *ssp ;
		    ssp1->a1 = b1 ; ssp1->a2 = ssp->a2 ;
		  }
	      }
	    else if (ssp->type & gI)
	      {
		DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
		*ssp1 = *ssp ;
	      }	    
	  }
      arrayDestroy (ss) ;
      ss = ss1 ;
      iMax = arrayMax (ss) ;
      
      for (int ii = 0, j3 = 0, j5 = 0 ; ii < iMax ; ii++)
	{
	  ssp = arrayp (ss1, ii, DSX) ;
	  if (ssp->type & gX)
	    {
	      int b1 = ssp->a1 ;
	      int b2 = ssp->a2 ;
	      
	      while (j5 > 0 && keySet (ks5, j5) > b1) j5-- ;
	      while (j5 < j5Max - 1 && keySet (ks5, j5) < b1) j5++ ;
	      if (keySet (ks5, j5) == b1) ssp->type |= g5 ;
	      
	      while (j3 > 0 && keySet (ks3, j3) > b2) j3-- ;
	      while (j3 < j3Max - 1 && keySet (ks3, j3) < b2) j3++ ;
	      if (keySet (ks3, j3) == b2) ssp->type |= g3 ;
	    }
	}
    }

  if (0 && debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after split\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after split done\n") ;
    }
	      
  arraySort (ss, dsA1Order) ;
	    
  if (0 && debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after split sort\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after split sort done\n") ;
    }
	      
  /* at each position choose the maximal value */
  if (1)
    {
      int jj = 0, gIX = gI | gX ;

      for (int ii = 0 ; ii < iMax ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  
	  int type = ssp->type & gIX ;
	  int b1 = ssp->a1 ;
	  int b2 = ssp->a2 ;
	  DSX *ssp1 ;
	  
	  for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax && (ssp1->type & gIX) == type && ssp1->a1 == b1 && ssp1->a2 == b2 ; jj++, ssp1++)
	    {
	      if ((ssp->type & gX) && (ssp1->type & (gDF | gFF)))
		ssp1->cover = ssp->cover ;
	      ssp->type |= ssp1->type ;
	      if (ssp->cover < ssp1->cover)
		ssp->cover = ssp1->cover ;
	      if (ssp->end < ssp1->end)
		{ ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
	      if (ssp->flag < ssp1->flag)
		ssp->flag = ssp1->flag ;
	      ssp1->cover = 0 ;
	      ssp1->flag = 0 ;
	    }
	  ii = jj - 1 ;
	}
    }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after top score\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after top score done\n") ;
    }

 
  if (iMax) /* merge continuous exons */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
	 {
	   if (ssp->a1 + 12 < ssp->a2 && (ssp->type & g5))
	     {
	       if ((ssp[1].type & gX) && ! (ssp[1].type & g5) && ssp->cover < ssp[1].cover && ssp[1].a1 == ssp->a2 + 1)
		 ssp->cover = ssp[1].cover ; 
	     }
	 }
       for (jj = 1, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
	 {
	   if (ssp->a1 +12 < ssp->a2 && (ssp->type & g3) && ssp[-1].a2 + 1 == ssp->a1)
	     {
	       if ((ssp[-1].type & gX) && ! (ssp[-1].type & g3) && ssp->cover < ssp[-1].cover)
		 ssp->cover = ssp[-1].cover ; 
	     }
	 }
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if ((ssp->type & gX) && ! (ssp->type & (gS | gA | g5 | g3)))
	     {
	       int i ;
	       int b1 = ssp->a2 + 1 ;
	       DSX *ssp1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; ssp1->a1 == b1 && ssp->type == ssp1->type && i < iMax ; ssp1++, i++)
		 {
		   if (ssp->cover < ssp1->cover)
		     ssp->cover = ssp1->cover ;
		   if (ssp->end < ssp1->end)
		     { ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
		   ssp->a2 = ssp1->a2 ;
		   ssp1->cover = 0 ;
		   b1 = ssp->a2 + 1 ;
		 }
	     }
	 }
     }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge exons \n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge exons done\n") ;
    }

 
  if (iMax) /* absorb in exons 2% of intron support */
     {
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       DSX *ssp1 ;
	       int i, c = ssp->cover/50 ;
	       int b1 = ssp->a1 ;
	       int b2 = ssp->a2 ;
	       if (c > 0)
		 {
		   for (i = -1, ssp1 = ssp - 1 ; ssp1->a1 == b1 && jj + i >= 0 ; ssp1--, i--)
		     if ((ssp1->type & gX) && ! (ssp1->type & (gS | gA | g5 | g3)))
		       ssp1->cover -= c ;
		   for (i = 1, ssp1 = ssp + 1 ; ssp1->a2 <= b2 && jj + i < iMax ; ssp1++, i++)
		     if ((ssp1->type & gX) && ! (ssp1->type & (gS | gA | g5 | g3)))
		       ssp1->cover -= c ;
		 }
	     }
	 }
     }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 3%% of intron score\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 3%% of intron score done\n") ;
    }


  if (iMax) /* absorb leaking exons at donor or acceptor site */
     {
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       DSX *ssp1 ;
	       int i ;
	       int b1 = ssp->a1 ;
	       int b2 = ssp->a2 ;

	       for (i = -1, ssp1 = ssp - 1 ; ssp1->a1 == b1 && jj + i >= 0 ; ssp1--, i--)
		 if ((ssp1->type & gX) && ssp1->a2 < b1 + 12 && ! (ssp1->type & (gA | g3)))
		   {
		     int k = 1 ;
		     while (ssp1[k].type & gI)
		       k++ ;
		     if ((ssp1[k].type & gX) && ssp1[k].a1 == ssp1->a2 + 1)
		       ssp1->cover = ssp1[k].cover ;
		     else
		       ssp1->cover = 0 ;
		   }
	       for (i = 1, ssp1 = ssp + 1 ; ssp1->a1 <= b2 && jj + i < iMax ; ssp1++, i++)
		 if ((ssp1->type & gX) && ssp1->a1 > b2 - 12 && ! (ssp1->type & (gS | g5) ))
		   {
		     int k = 1 ;
		     while (ssp1[-k].type & gI)
		       k++ ;
		     if ((ssp1[-k].type & gX) && ssp1[-k].a2 == ssp1->a1 - 1)
		       ssp1->cover = ssp1[-k].cover ;
		     else
		       ssp1->cover = 0 ;
		   }
	     }
	 }
     }

  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after absorb leaking exons\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after absorb leaking exons done\n") ;
    }
  
  
  /* regularize the real5p */
  if (1)
    {
      /* rehook the 1 bp SL to next exons */
      for (int ii = 0 ; ii < iMax  ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gS)
	    {
	      DSX *ssp1 ;
	      int jj, b1 = ssp->a2 + 1 ;
	      for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if ((ssp1->type & gX ))
		    {
		      if (ssp1->a1 < b1 + 30 && ssp1->a1 > ssp->a2 + 1)
			ssp->a2 = ssp1->a1 - 1 ;
		      break ;
		    }
		}
	    }
	}

      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{ /* clean up stretches of r5p */
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal5p)
	    {
	      DSX *ssp1 ;
	      int jj, b2, k ;
	      
	      /* discard gReal5p leaking above gS */
	      for (int pass = 0 ; pass < 2 ; pass++)
		for (jj = ii + 1, b2 = ssp->a2 + 1, ssp1 = ssp + 1, k = 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++, k++)
		  {
		    if (ssp1->type & gI)
		      {
			if (pass == 0)
			  b2 = ssp1->a2 + 1 ;
			continue ;
		      }
		    if (! (ssp1->type & gReal5p)) break ;
		    if (ssp1->type & gS)
		      {
			for (int m = 0 ; m < k ; m++)
			  if (! (ssp1[-m-1].type & gS))
			    ssp1[-m-1].type &= ~ gReal5p ;
			break ;
		      }
		    b2 = ssp1->a2 + 1 ;
		    if (0 && b2 > ssp->a1 + 30) break ;
		}
	    }
	}

      for (int pass = 0 ; pass < 2 ; pass++)
	for (int ii = 0 ; ii < iMax - 1 ; ii++)
	  { /* clean up stretches of r5p */
	    ssp = arrp (ss, ii, DSX) ;
	    if (ssp->type & gReal5p)
	      {
		DSX *ssp1 ;
		int jj, b2, k ;
		
		for (jj = ii + 1, b2 = ssp->a2 + 1, ssp1 = ssp + 1, k = 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++, k++)
		  {
		    if (ssp1->type & gI)
		      {
			if (pass == 0)
			  b2 = ssp1->a2 + 1 ;
			continue ;
		      }
		    if (! (ssp1->type & gReal5p)) break ;
		    if (ssp1->cover > ssp->cover)
		      for (int m = 0 ; m < k ; m++)
			ssp[m].cover = ssp1->cover ;
		    if (ssp1->end > ssp->end)
		      { ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
		    ssp1->type &= ~gReal5p ;
		    b2 = ssp1->a2 + 1 ;
		  }
	      }
	  }


      for (jj = 0 ; jj < iMax - 1 ; jj++)
	{ /* drop exons above real5p */
	  ssp = arrp (ss, jj, DSX) ;
	  if (ssp->type & (gReal5p | gS | gI))
	    break ;
	}
      if (ssp->type & (gReal5p | gS)) /* we are on top real5p with no intron above */
	for (int ii = 0 ; ii < jj ; ii++)
	  { /* if they are weak */
	    DSX* ssp1 = arrp (ss, ii, DSX) ;
	    
	    if (ssp1->cover < ssp->cover / 2)
	      ssp1->cover = 0 ;
	    else
	      break ;
	  }

      iMax = ssHappyFew (ss) ;
      for (jj = 0 ; jj < iMax ; jj++)
	{ /* drop internal exons above real5p */
	  ssp = arrp (ss, jj, DSX) ;
	  if (ssp->type & gReal5p)
	    {
	      DSX *ssp1 ;
	      int k = 0, jj, bestCover = ssp->cover, b2 = ssp->a1 - 1 ;
	      
	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) { if (k > 0) { k = 0; break ; } continue ;}
		  if (ssp1->type & (gA | gReal5p | gS | gReal3p)) {k = 0 ; break ; }
		  if (ssp1->cover > bestCover / 2) {k = 0 ; break ; }
		  b2 = ssp1->a1 - 1 ; k++ ;
		}
	      for (jj = ii - 1, ssp1 = ssp - 1 ; k > 0 && jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  b2 = ssp1->a1 - 1 ;
		  ssp1->cover = 0 ; k-- ;
		}
	    }
	}
      iMax = ssHappyFew (ss) ;
      for (int ii = iMax -1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (0 &&  (ssp->type & gReal5p))
	    {
	      DSX *ssp1, *ssp0 = ssp ;
	      int jj, bestCover = ssp->cover, b2 = ssp->a1 - 1 ;

	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gA | gReal3p)) break ;
		  if (ssp1->cover < bestCover / 100)
		    {
		      for ( ; jj >= 0 && ssp1->a2 == b2 && (ssp1->type & gReal5p) ; jj--, ssp1--)
			{
			  if (ssp1->type & gI) continue ;
			  b2 = ssp1->a1 - 1 ;
			  ssp1->type &= ~ gReal5p ;
			}
		      break ;
		    }
		  else
		    {
		      if (ssp1->cover > bestCover) bestCover = ssp1->cover ;
		      ssp0->type &= ~gReal5p ;
		      ssp1->type |= gReal5p ;
		      ssp0 = ssp1 ;
		      b2 = ssp1->a1 - 1 ;
		    }
		}
	    }
	}
    }
  
  /* regularize the real3p */
  if (1)
    {
      /* rehook the 1 bp polyA to previous exons */
      for (int ii = 0 ; ii < iMax  ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & (gA | gReal3p))
	    {
	      DSX *ssp1 ;
	      int jj, b2 = ssp->a1 - 1 ;
	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if ((ssp1->type & gX ))
		    {
		      if (ssp1->a2 > b2 - 30 && ssp->a1 > ssp1->a2 + 1)
			ssp->a1 = ssp1->a2 + 1 ;
		      break ;
		    }
		}
	    }
	}

      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal3p)
	    {
	      DSX *ssp1, *ssp2 = ssp ;
	      int jj, b2 = ssp->a2 + 1, jjA = ii ;

	      /* search for a quasi contiguous polyA */
	      for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax && ssp1->a1 <= b2+10 ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gS | gReal5p))
		    break ;
		  if (0 && (ssp1->type & gReal3p) && (ssp1->end > 3 * ssp->end)) /* tyty */
		    jjA = jj ;
		  if (ssp1->type & gA)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if (ssp1->type & gReal3p)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if ((ssp1->type & gCompleteCDS) && ssp1->a2 < b2 + 30)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if (ssp1->a2 > ssp2->a2 + 50) break ;
		  b2 = ssp1->a2 + 1 ;
		}
	      /* raise their cover */
	      b2 = ssp->a2 + 1 ;
	      int bestEnd = ssp->end ;
	      KEY bestKey = ssp->endKey ;
	      for (jj = ii , ssp1 = ssp ; jj <= jjA ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover < ssp->cover)
		    ssp1->cover = ssp->cover ;
		  if (ssp1->a1 > b2)
		    ssp1->a1 = b2 ;
		  if (ssp1->end < bestEnd)
		    { ssp1->end = bestEnd ; ssp1->endKey = bestKey ; }
		  if (ssp1->type & gReal3p)
		    if (ssp1->end > bestEnd)
		      { bestEnd = ssp1->end ; bestKey = ssp1->endKey ; }

		  b2 = ssp1->a2 + 1 ;
		  if (jj < jjA)
		    ssp1->type &= ~gReal3p ;
		  else
		    ssp1->type |= gReal3p ;
		}
	    }
	}
    }
    
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5p/r3p\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5p/r3p done\n") ;
    }

  if (iMax) /* merge SL in next exon */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if ((ssp->type & (gReal5p | gS)) && ! (ssp->type & (gA | gReal3p)))
	     {
	       int i ;
	       DSX *ssp1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; i < iMax ; ssp1++, i++)
		 {
		   if (ssp1->type & gI) continue ;
		   if (ssp1->cover <= 0) continue ;
		   if (ssp1->a1 == ssp->a2 + 1)
		     {
		       ssp1->type |= g5 ;
		     }
		   else
		     break ;
		 }
	     }
	 }
     }

  if (iMax) /* equalize the score of SL clusters */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gS)
	     {
	       int i, k = 0, cMax = ssp->cover ;
	       DSX *ssp1 ;
	       int b1 = ssp->a2 + 1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; i < iMax ; ssp1++, i++)
		 {
		   if (ssp1->type & gI)
		     continue ;
		   if (ssp1->a1 > b1)
		     break ;
		   if (ssp1->a1 > ssp->a1 + 40)
		     break ;
		   if (ssp1->type & gS)
		     {
		       k = i ;
		       if (ssp1->cover > cMax)
			 cMax = ssp1->cover ;
		     }
		   b1 = ssp1->a2 + 1 ;
		 }
	       for (i = jj ; i <= k ; i++)
		 {
		   ssp1 = arrp (ss, i, DSX) ;
		   if (ssp1->type & gI)
		     continue ;
		   ssp1->cover = cMax + i - jj ; /* by computing the lower one first we insure we do both */
		 }
	     }
	 }
     }

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL done\n") ;
    }

  if (0 && iMax) /* merge polyA in previous exon */
    {
      for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
	{
	  int jj1 = jj ;
	  if (ssp->cover <= 0)   continue ;
	  if ((ssp->type & gX) && (ssp->type & (gA)) && ! (ssp->type & (g5 | gS)))
	    {
	      int i ;
	      int c = ssp->cover ;
	      DSX *ssp1 ;

	      for (i = jj - 1, ssp1 = arrp (ss, i, DSX) ; i >= 0 ; i--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover <= 0) continue ;

		  if (ssp1->a2 >= ssp->a1 - 120)
		    {
		      jj1 = i ;
		      if (c < ssp1->cover)
			c = ssp1->cover ;
		    }
		  if (ssp1->type & gReal3p)
		    {
		      ssp->type |= gReal3p ;
		      ssp->cover = c ;
		      break ;
		    }
		}
	      for (i = jj - 1, ssp1 = arrp (ss, i, DSX) ; i >= jj1 ; i--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover <= c) ssp1->cover = c ;
		  if (ssp->type & gReal3p)
		    ssp1->type &= ~gReal3p ;
		}
	      jj = jj1 ;
	    }
	}
    }
  
  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge same, SL, pA\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge same, SL, pA done\n") ;
    }


  if (iMax) /* clip exons above first strong signal */
     {
       int bestCover = 0 ;
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     continue ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->type & (gS |gReal5p))
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	 }
       if (ssp->type & (gS | gReal5p))
	 {
	   /* I found a real5p, kill all exons above */
	   for (int i = 0 ; i < jj ; i++)
	     if (ssp[-i-1].a2 < ssp->a1 && (10*ssp[-i - 1].cover <= bestCover || (ssp[-i-1].type & (gDF | gFF))))
	       ssp[-i - 1].cover = 0 ;
	     else
	       break ;
	   if (! (ssp->type & (gS | g5)))
	     {
	       ssp->a1 -= 0 ;  /* was 15 move back to original position */
	       if (ssp->a1 < 1)
		 ssp->a1 = 1 ;
	     }
	 }
     }
  if (iMax) /* clip exons after last strong signal */

    {
       int bestCover = 0 ;
       for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     break ;
	   if (ssp->type & gReal5p)
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	 }
       if (ssp->type & gReal3p)
	 {
	   /* I found a real3p, kill all exons below */
	   DSX *ssp1 ;
	   int i ;
	   for (i = iMax - 1, ssp1 = arrp (ss, i, DSX) ; i > jj ; i--, ssp1--)
	     {
	       if (ssp1->type & gA)
		 break ;
	       else
		 ssp1->cover = 0 ;
	     }
	 }
     }
  if (iMax) /* clip exons after last polyA */
    {
       int bestCover = 0 ;
       for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     break ;
	   if (ssp->type & gA)
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	 }
       if (2 * ssp->cover >= bestCover)
	 if (ssp->type & (gA | gReal3p))
	   {
	     /* I found a polyA, kill all exons below */
	     for (int i = jj + 1 ; i < iMax ; i++)
	       ssp[i - jj].cover = 0 ;
	   }
    }
  
  iMax = ssHappyFew (ss) ;
  if (iMax)  /* kill intron ending below last exon */
    {
      int b2 = 0 ;
      for (ii = iMax - 1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gX)
	    { b2 = ssp->a2 ; break ; }
	}
      for (ii = iMax - 1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->a2 > b2 && (ssp->type & gI))
	    ssp->cover = 0 ;
	}
    }      
  
  iMax = ssHappyFew (ss) ;
  if (0 &&    /* bad idea check worm_WS294 gene chr3__0_692 */
      iMax) /* create a cut over the SL */
    {
       for (ii = 0, ssp = arrp (ss, 0, DSX) ; ii < iMax ; ii++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gS)
	     {
	       int k, b2 = ssp->a1 - 1 ;
	       for (k = 1 ; k < ii ; k++)
		 {
		   DSX *ssp1 = ssp - k ;
		   if (ssp1->type & gI) continue ;
		   if (ssp1->a2 != b2) break ;
		   if (! (ssp1->type & gS) && ssp1->cover < ssp->cover/30)
		     { ssp1->cover = 0 ; break ; }
		   b2 = ssp1->a1 - 1 ;
		 }
	     }
	 }
    }

 done:
  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  mrnaDesignCleanExons (sc, ss) ;
  ds->exons = ss ;
  
  if (0 || debug)
     {
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns done\n") ;
     }

   ac_free (h) ;
   return ds->exons && arrayMax (ds->exons) ? TRUE : FALSE ;
} /* mrnaDesignGetGraphElements */

/**********************************************************************************/
/* extend the path Down via the best supported route */
static int mrnaDesignExtendDownRaw (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxScore, int flag)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;
  
  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;

  for (vp = up, ii = nn0 ; ii >= 0 ; ii--, vp--)
    {
      if (vp->donor == path)
	{
	  if (vp->type & gI)
	    break ;
	  int score2 = ((vp->type & gI) ? vp->score/3 :  vp->score) ;
	  if (bestScore < score2)
	    bestScore = score2 ;
	}
    }

  for (nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
    {
      if (vp->a1 != a2) continue ;
      if ((vp->flag & flag) != flag) continue ;
      int score2 = ((vp->type & gI) ? vp->score/1 :  vp->score) ;
      if (score < score2) { score = score2 ; nn = vp->nn ; wp = vp ; }
    }

  if (! score)
    return length ;

  if (0 && nStop && 
      (wp->type & (gS | gS0 | gReal5p))
      )
    return length ;
  
  if (0 && ! (wp->type & (gA | gReal3p)))
    {
      if (nStop)
	{
	  if (up->score / 10 > score)
	    return length ;
	}
      if (up->score / 100 > score)
	return length ;

      if (score < 0 * bestScore / 100)
	return length ;
    }

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->donor = path ;
  
  if (wp->type & gI)
    { nIntron++ ; nStop = 0 ; }
  if (wp->type & gS)
    nStart++ ;
  if (wp->type & (gA | gReal3p))
    nStop++ ;
  score = ((wp->type & gI) ? wp->score/3 :  wp->score) ;
  if (score > bestScore)
    bestScore = score ;
  if (0) fprintf (stderr, "Path=%d  %d:%d >> %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
  length += mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore, flag) ;

  return length ;
} /* mrnaDesignExtendDownRaw */

/***********/

static int mrnaDesignExtendDownLoopCDS (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxScore, int flag)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;
  if (1) /* first try to chain the CDS */
    for (nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
      {
	if (vp->a1 != a2) continue ;
	if ((vp->flag & flag) != flag) continue ;
	if (! (vp->type &  (gCompleteCDS)))
	  continue ;
	if (score < vp->cover) { score = vp->cover ; nn = vp->nn ; wp = vp ; }
      }
  if (! score)
    return   mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn0, nIntron, nStart, nStop, bestScore, maxScore, flag) ;

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  wp->donor = path ;
  if (score > bestScore)
    bestScore = score ;
  length += mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore, flag) ;

  return length ;
} /* mrnaDesignExtendDownLoopCDS */

/***********/

static int mrnaDesignExtendDown (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn, int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxScore, int flag)
{
  if (useCDS)
    return mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore, flag) ;

  return mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore, flag) ;
} /* mrnaDesignExtendDown */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static int mrnaDesignExtendUp (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0
			       , int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxScore, int flag)
{ 
  int length = 0 ;
  int a1, ii, nn = 0, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  if (up->nn != nn0) messcrash ("nn != nn0 in  mrnaDesignExtendUp") ;
  a1 = up->a1 - 1 ;

  if (useCDS)
    for (nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
      {
	if (vp->a2 != a1) continue ;
	if ((vp->flag & flag) != flag) continue ;
	if (! (vp->type &  (gCompleteCDS)))
	  continue ;
	if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }
  if (! score)
    {
      useCDS = 0 ;
      for (nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
	{
	  if (vp->a2 != a1) continue ;
	  if ((vp->flag & flag) != flag) continue ;
	  if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
	}
    }

  if (! score)
    return length ;

  if (0 && (up->type & (gS | gReal5p)) && up->score > wp->score/2)
    return length ;

  if (0 && nStart &&
      (wp->type & (gA | gReal3p))
      )
    return length ;
    
  if (! (up->type & gI) && ! (wp->type & (gS | gS0 | gReal5p)))
    {
      if (up->type & gReal5p)
	{
	  if (up->score/ 10 > score)
	    return length ;
	}
      if (up->score / 100 > score)
	return length ;

      if (score < bestScore / 100)
	return length ;
    }

  /* find local exon score */
  if (0 && nStart && ! (wp->type & (gS | gReal5p | gCompleteCDS)))
    if (score < 0 * bestScore / 100)
      return length ;
    
  if (! useCDS)
    {
      if (0 && nStart && (wp->type & ( gA | gReal3p) && ! (wp->type & gCompleteCDS)))
	return 0 ;
      if (0 && nIntron && (wp->type & (gA | gReal3p)))
	return 0 ;
    }
  
  length = wp->a2 - wp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->donor = path ;
  if (! useCDS)
    {
      if (wp->type & gI)
	{ nIntron++ ; nStart = 0 ; }
      if (wp->type & (gS | gReal5p))
	nStart++ ;
      if (wp->type & (gA | gReal3p))
	nStop++ ;
    }
  if (wp->cover > bestScore)
    bestScore = wp->cover ;
  if (0) fprintf (stderr, "Path=%d  %d:%d << %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
  length += mrnaDesignExtendUp (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, useCDS, bestScore, maxScore, flag) ;

  return length ;
} /* mrnaDesignExtendUp */

/**********************************************************************************/

static int mrnaDesignExport (S2M *s2m, SC *sc, DS *ds, Array segs, KEYSET ks, int path, Array smrnas)
{
  int a1 = 0, ii,jj,  nn, maxScore = 0 ;
  DSX *up, *up2 ;
  SMRNA *smrna ;
  HIT *vp = 0 ;
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  
  Array exions = arrayHandleCreate (keySetMax(ks), DSX, h) ;
  smrna = arrayp(smrnas, path - 1, SMRNA) ;
  arrayDestroy (smrna->dnas) ;
  arrayMax (smrnas) = path ;
  if (! smrna->hits)
    smrna->hits = arrayHandleCreate (keySetMax(ks), HIT, s2m->h) ;
  else
    smrna->hits = arrayReCreate (smrna->hits, keySetMax(ks), HIT) ;
  if (! smrna->clones)
    smrna->clones = keySetHandleCreate (s2m->h) ; 
  keySetSort (ks) ; /* now the elements are in a1 order and correctly chained */
  for (ii = 0, maxScore = 100 ; ii < keySetMax (ks) ; ii++)
    {
      nn = keySet (ks, ii) ; 
      up = arrp (segs, nn, DSX) ;
      if (up->score > maxScore)
	maxScore = up->score ;
      up2 = arrayp (exions, ii, DSX) ;
      *up2 = *up ;
    }
   if (0 && debug)
     {
       fprintf (stderr, "mrnaDesignExport path %d exons/introns\n", path) ;
       showDs (sc, exions) ;
       fprintf (stderr, "mrnaDesignExport exons/introns done\n") ;
     }
#ifdef JUNK
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      up = arrp (exions, ii, DSX) ;
      if (up->type & (gS | gS0) ||
	  up->score > maxScore / 100)
	break ;
      if (0) up->score = 0 ;
    }
  if (0 && ii > 0 && (up->type & gI) && (up[-1].type & gX) && up[-1].score == 0)
    up[-1].score = up->score ;
  for (ii = keySetMax (ks) - 1 ; ii >= 0 ; ii--)
    {
      up = arrp (exions, ii, DSX) ;
      if (up->type & (gA) ||
	  up->score > maxScore / 100)
	break ;
      if (0) up->score = 0 ;
    }
#endif
  
   if (0 || debug)
     {
       fprintf (stderr, "### mrnaDesignExport path %d exons/introns\n", path) ;
       showDs (sc, exions) ;
       fprintf (stderr, "### mrnaDesignExport path %d exons/introns done\n", path) ;
     }

   for (vp = 0, ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      up = arrp (exions, ii, DSX) ;
      if (!ii && (up->type & gI))
	continue ;
      if (! up->score)
	continue ;
      if (!vp)
	a1 = smrna->a1 = up->a1 ; 
      smrna->a2 = up->a2 ; 
      if (!vp || ((up->type & (gX | gI)) != (vp->type & (gX | gI))))
	{	
	  vp = arrayp (smrna->hits, jj++, HIT) ;
	  vp->a1 = up->a1 - a1 + 1 ;
	}
      vp->a2 = up->a2 - a1 + 1 ;
      vp->type |= ( up->type  & (~ gFF)) ;
      if (jj > 1)
	vp->type &= (~ (gS | gS0 | gD | gReal5p | g5 | g3 | gDF)) ;
    }
  for (ii = 0 ; ii < jj - 1 ; ii++)
    {
      vp = arrayp (smrna->hits, ii, HIT) ;
      vp->type &= (~ (gA | gF | gReal3p | g5 | g3 | gFF)) ;
    }
  if (debug)
    showHits (smrna->hits) ;

  ac_free (h) ;
  return path ;
} /* mrnaDesignExport G_t_NC_003281.10_1_4561 */

/**********************************************************************************/

static BOOL mrnaDesignIsNewPath (Array ss, Array ksPaths, KEYSET ks0, int path, int xStart, int xStop, AC_HANDLE h)
{
  BOOL ok = TRUE ;
  KEYSET ks = keySetHandleCopy (ks0, h) ;
  int iMax = keySetMax (ks) ;
  keySetSort (ks) ;

  if (iMax && xStart)
    {
      int ii = keySet (ks, 0) ;
      DSX *ssp = arrp (ss, ii, DSX) ;
      if (xStart != ssp->a1)
	ok = FALSE ;
    }
  if (iMax && xStop)
    {
      int ii = keySet (ks, iMax - 1) ;
      DSX *ssp = arrp (ss, ii, DSX) ;
      if (xStop != ssp->a2)
	ok = FALSE ;
    }

  if (1)
    {
      int minCover = 0, maxCover = 0 ;
      int nI = 0 ;
      int nEnds = 0 ;
      
      for (int i = 0 ; i < iMax ; i++)
	{
	  int jj = keySet (ks, i) ;
	  DSX *ssp = arrp (ss, jj, DSX) ;
	  if (ssp->type & (gS | gReal3p | gReal5p))
	    nEnds++ ;
	  if (ssp->type & gI)
	    {
	      nI++ ;
	      if (! maxCover) minCover = ssp->cover ;
	      if (ssp->cover > maxCover) maxCover = ssp->cover ;
	      if (ssp->cover < minCover) minCover = ssp->cover ;
	    }
	}
      if (1000 * minCover < maxCover)
	ok = FALSE ;
      if (nI + nEnds == 0)
	ok = FALSE ;
    }
  
  for (int ii = 0 ; ok && ii < arrayMax (ksPaths) ; ii++)
    {
      KEYSET ks1 = array (ksPaths, ii, KEYSET) ;
      ok = FALSE ;
      if (! ks1 || keySetMax (ks1) != iMax)
	ok = TRUE ;
      else
	{
	  for (int i = 0 ; !ok && i < iMax ; i++)
	    if (keySet (ks, i) != keySet (ks1, i))
	      ok = TRUE ;
	}
    }

  /* flag the paths */
  for (int i = 0 ; i < iMax ; i++)
    {
      int ii = keySet (ks, i) ;
      DSX *up = arrp (ss, ii, DSX) ;
      if (ok) up->path = path ;
      up->donor = 0 ;
    }
  if (ok)
    array (ksPaths, arrayMax(ksPaths), KEYSET) = ks ;
  else
    keySetDestroy (ks) ;
  return ok ;
} /* mrnaDesignIsNewPath */

/**********************************************************************************/

static int mrnaDesignFlagDescendants (Array ss, int ii, int flag)
{
  int jj, n = 0 ;
  int iMax = arrayMax (ss) ;
  DSX *vp, *up = arrp (ss, ii, DSX) ;
  int a1 = up->a2 + 1 ;
  
  up->flag |= flag ;
  for (jj = ii + 1, vp = up + 1 ; jj < iMax ; vp++, jj++)
    {
      if (vp->a1 == a1)
	mrnaDesignFlagDescendants (ss, jj , flag) ;
      n++ ;
    }

  return n ;
} /* mrnaDesignFlagDescendants */

/**********************************************************************************/

static int mrnaDesignFlagAscendants (Array ss, int ii, int flag)
{
  int n = 0 ;
  int jj ;
  DSX *vp, *up = arrp (ss, ii, DSX) ;
  int a2 = up->a1 - 1 ;
  
  up->flag |= flag ;
  for (jj = ii - 1, vp = up - 1 ; jj >= 0 ; vp--, jj--)
    {
      if (vp->a2 == a2)
	mrnaDesignFlagAscendants (ss, jj , flag) ;
      n++ ;
    }
  return n ;
} /* mrnaDesignFlagAscendants */

/**********************************************************************************/

static void mrnaDesignSetEndShades (KEYSET light, KEYSET dark)
{
  int nLight = keySetMax (light) ;
  int nDark = keySetMax (dark) ;

  if (nLight + nDark)
    {
      AC_HANDLE h = ac_new_handle () ;
      vTXT txt = vtxtHandleCreate (h) ;
      const char *errors = 0 ;
      AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */
      
      keySetSort (light) ; keySetCompress (light) ;
      keySetSort (dark) ; keySetCompress (dark) ;
      for (int i = 0 ; i < nLight ; i++)
	{
	  const char *nam = name(keySet(light, i)) ;
	  vtxtPrintf (txt, "Sequence %s\n", nam) ;
	  if (! strncmp (nam, "Xends_ELF", 9) || ! strncmp (nam, "Xends_ERR", 9))
	    vtxtPrintf (txt, "Colour PALECYAN\n") ;
	  if (! strncmp (nam, "Xends_ERF", 9) || ! strncmp (nam, "Xends_ELR", 9))
	    vtxtPrintf (txt, "Colour PALEYELLOW\n") ;
	  vtxtPrintf (txt, "\n") ;
	}
      for (int i = 0 ; i < nDark ; i++)
	{
	  const char *nam = name(keySet(dark, i)) ;
	  vtxtPrintf (txt, "Sequence %s\n", nam) ;
	  if (! strncmp (nam, "Xends_ELF", 9) || ! strncmp (nam, "Xends_ERR", 9))
	    vtxtPrintf (txt, "Colour CYAN\n") ;
	  if (! strncmp (nam, "Xends_ERF", 9) || ! strncmp (nam, "Xends_ELR", 9))
	    vtxtPrintf (txt, "Colour YELLOW\n") ;
	  vtxtPrintf (txt, "\n") ;
	}
      ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;
      
      ac_free (h) ;
    }
  return ;
} /*  mrnaDesignSetEndShades */

/**********************************************************************************/

static int mrnaDesignFindStartEndPairs (Array ss, Array ss2, Array sFlags, Array starts, Array stops)
{
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  DSX *up, *up2, *vp = 0, *wp = 0, *sFlag = 0 ;
  int i, iMax = arrayMax (ss), iStart = 0, iStop = 0, startCover = 0, stopCover = 0 ;
  int nFlags = 0 ;
  KEYSET light = keySetHandleCreate (h) ;
  KEYSET dark = keySetHandleCreate (h) ;
  int nLight = 0, nDark = 0 ;
  
  for (i = 0, up = arrp (ss, 0, DSX) ; i < iMax ; i++, up++)
    {
      if (! (up->type & gS) && (up->type & gReal5p) && up->cover > 5 * up->end)
	{ up->type &= ~gReal5p ; keySet (light, nLight++) = up->endKey ; }
      if ((up->type & gReal3p) && up->cover > 5 * up->end)
      	{ up->type &= ~gReal3p ; keySet (light, nLight++) = up->endKey ; }
    }      
  
  for (i = 0, up2 = arrp (ss2, 0, DSX) ; i < iMax ; i++, up2++)
    {
      up = arrp (ss, up2->nn, DSX) ;
      up->path = i ; /* temporary */
    }

  for (i = 0, up = arrp (ss, 0, DSX), vp = 0 ; i < iMax ; i++, up++)
    {
      if ((up->type & gReal5p) && iStart < 15)
	{
	  if (!vp || vp->a2 < up->a1 + 1)
	    {
	      vp = arrayp (starts, iStart++, DSX) ;
	      *vp = *up ;
	      vp->nn = up->path ;
	      vp->path = i ;
	      vp->flag = 0 ;
	      vp->type = up->type & gS ;
	      startCover += vp->end ;
	    }
	  vp = up ;
	}
    }

  for (i = iMax - 1, up = arrp (ss, i, DSX), wp = 0 ; i >= 0 ; i--, up--)
    {
	if ((up->type & gReal3p) && iStop < 15)
	  {
	    if (! wp || wp->a1 > up->a2 + 1)
	      {
		wp = arrayp (stops, iStop++, DSX) ;
		*wp = *up ;
		wp->nn = up->path ;
		wp->path = i ;
		wp->flag = 0 ;
		wp->type = 0 ;
		stopCover += wp->end ;
	    }
	    wp = up ;
	  }
    }

  for (i = 0, up2 = arrp (ss2, 0, DSX), vp = 0 ; i < iMax ; i++, up++)
    up->path = up->flag = 0 ;

  /* normalize */
  if (startCover > 0)
    for (int i = 0 ; i < iStart ; i++)
      {
	DSX *vp = arrayp (starts, i, DSX) ;
	vp->type = vp->flag = 0 ;
	vp->end *=  10000.01 / startCover ;
	if (vp->type & gS)
	  vp->end += 6000 ;
	vp->type = 0 ;
	if ((vp->flag & gS) || vp->end > 2000)
	  keySet (dark, nDark++) = vp->endKey ;
	else
	  keySet (light, nLight++) = vp->endKey ;
      }
  if (stopCover > 0)
    for (int i = 0 ; i < iStop ; i++)
      {
	DSX *wp = arrayp (stops, i, DSX) ;
	wp->type = wp->flag = 0 ;
	wp->end *=  10000.01 / stopCover ;
	if (wp->end > 2000)
	  keySet (dark, nDark++) = wp->endKey ;
	else
	  keySet (light, nLight++) = wp->endKey ;
      }
  mrnaDesignSetEndShades (light, dark) ;
  /* sort */
  iStart = ssHappyFew (starts) ;
  iStop = ssHappyFew (stops) ;
  if (iStart > 1) arraySort (starts, dsEndOrder) ;
  if (iStop > 1) arraySort (stops, dsEndOrder) ;

  if (0 || debug)
    {
      for (int i = 0 ; i < iStart ; i++)
	{
	  vp = arrayp (starts, i, DSX) ;
	  fprintf (stderr, "..start %d : a1=%d startCover %d %% %ld\n"
		   , 1<<i, vp->a1, vp->end, (long) startCover * vp->end / 100 
		   ) ;
	}
      for (int i = 0 ; i < iStop ; i++)
	{
	  wp = arrayp (stops, i, DSX) ;
	  fprintf (stderr, "..stop %d : a2=%d stopCover %d %% %ld\n"
		   , 1<<i, wp->a2, wp->end, (long) stopCover * wp->end / 100
		   ) ;
	}
    }
  
  /* flag the segs and all their descendants or ascendants */
  for (int i = 0 ; i < iStart ; i++)
    {
      int flag = 1 << i ;
      vp = arrayp (starts, i, DSX) ;
      mrnaDesignFlagDescendants (ss, vp->path, flag) ;
    }
  for (int i = 0 ; i < iStop ; i++)
    {
      int flag = 1 << (i + 16) ;
      vp = arrayp (stops, i, DSX) ;
      mrnaDesignFlagAscendants (ss, vp->path, flag) ;
    }

  /* label the good starts (strongest or saves a new intron or a new stop) */
  int allFlag = 1 ;
  if (iStart)
    {
      int flag2 = 1 ;
      vp = arrp (starts, 0, DSX) ;
      vp->type = 1 ; allFlag |= flag2 ; vp->flag = flag2 ; 
      for (int ii = 1 ; ii < iStart ; ii++)
	{ /* regarder si on decouvre de nouveaux introns ou stops */
	  BOOL ok = FALSE ;

	  flag2 = 1 << ii ;
	  vp = arrp (starts, ii, DSX) ;
	  if (vp->end > 6000)  
	    ok = TRUE ;
	  else if (vp->end > 0)
	    ok = TRUE ;
	  else if (! (vp->flag & allFlag))
	    {
	      int i ;
	      for (i = 0, up = arrp (ss, 0, DSX) ; ! ok && i < iMax ; i++, up++)
		if ((up->flag & flag2) && ! (up->flag & allFlag))
		  {
		    if (up->type & gI)
		      { /* this starts yields a new intron */
			ok = TRUE ;
		      }
		    if (vp->end > 100 && (up->type & gX))
		      { /* this starts yields a new exon */
			ok = TRUE ;
		      }
		  }
	    }
	  if (ok)
	    { vp->type = 1 ; allFlag |= flag2 ; vp->flag = flag2 ; }
	}
    }
  if (iStop)
    {
      int flag2 = 1 << (16) ;
      wp = arrp (stops, 0, DSX) ;
      wp->type = 1 ; allFlag |= flag2 ; wp->flag = flag2 ;
      for (int ii = 1 ; ii < iStop ; ii++)
	{ /* regarder si on decouvre de nouveaux introns ou starts */
	  BOOL ok = FALSE ;

	  flag2 = 1 << (16 + ii) ;
	  wp = arrp (stops, ii, DSX) ;
	  if (wp->end > 6000)
	    ok = TRUE ;
	  else if (! (wp->flag & allFlag))
	    {
	      int i ;
	      for (i = 0, up = arrp (ss, 0, DSX) ; ! ok && i < iMax ; i++, up++)
		if ((up->flag & flag2) && ! (up->flag & allFlag))
		  {
		    if (up->type & gI)  /* this starts yields a new intron */
		      ok = TRUE ;
		    if (wp->end > 100 && up->type & (gX))  /* this starts yields a new exon */
		      ok = TRUE ;
		}
	    }
	  if (ok)
	    { wp->type = 1 ; allFlag |= flag2 ; wp->flag = flag2 ;}
	}
    }
  
  /* create the flag pairs */
  sFlag = arrayp (sFlags, 0, DSX) ; /* default, but do not set nFlags++ */
  int flag3 = 0 ;
  if (iStart)
    {
      for (int ii = 0 ; ii < iStart ; ii++)
	{
	  vp = arrp (starts, ii, DSX) ;
	  if (vp->type)
	    { /* foreach good start, select the highest matching stop */
	      BOOL ok = FALSE ;
	      int flag = 1 << (ii) ;
	      for (int j = 0 ; ! ok && j < iStop ; j++)
		{
		  wp = arrp (stops, j, DSX) ;
		  if (wp->a2 > vp->a1)
		    {
		      up = arrp (ss, wp->path, DSX) ;
		      if (up->flag & flag)
			{
			  int flag2 = 1 << (16 + j) ;
			  flag3 |= flag2 ;
			  if (debug)
			    fprintf (stderr, "Flag %d, a1 = %d, a2 = %d, score = %d / %d, flag = %d / %d\n"
				     , nFlags, vp->a1, wp->a2, vp->end, wp->end, flag, flag2 >> 16) ; 
			  sFlag = arrayp (sFlags, nFlags++, DSX) ;
			  sFlag->flag = flag | flag2 ;
			  ok = TRUE ;
			}
		    }
		}
	      if (! ok)
		{
		  if (debug)
		    fprintf (stderr, "Flag %d, a1 = %d, a2 = open, score = %d, flag = %d\n"
			     , nFlags, vp->a1, vp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag ;

		}
	    }
	  else if (vp->end > 10)
	    { /* we are in a descendant start, select its closest stop */
	      int b2 = 0 ;
	      int bestj = -1 ;
	      int flag = 1 << (ii) ;
	      for (int j = 0 ; j < iStop ; j++)
		{
		  wp = arrp (stops, j, DSX) ;
		  if (wp->a2 > vp->a1)
		    {
		      if (wp->end > vp->end/2)
			{
			  up = arrp (ss, wp->path, DSX) ;
			  if ((up->flag & flag) && (!b2 || up->a2 < b2))
			    {
			      b2 = up->a2 ;
			      bestj = j ;
			    }
			}
		    }
		}
	      if (bestj >= 0)
		{
		  int flag2 = 1 << (16 + bestj) ;
		  wp = arrp (stops, bestj, DSX) ;
		  flag3 |= flag2 ;
		  if (debug)
		    fprintf (stderr, "Descendant flag %d, a1 = %d, a2 = %d, score = %d / %d, flag = %d\n"
			     , nFlags, vp->a1, wp->a2, vp->end, wp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag | flag2 ;
		}
	      else
		{
		  if (debug)
		    fprintf (stderr, "Descendant flag %d, a1 = %d, a2 = open, score = %d, flag = %d\n"
			     , nFlags, vp->a1, vp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag ;
		}
	    }
	}
    }
  if (iStop)
    {
      for (int j = 0 ; j < iStop ; j++)
	{
	  wp = arrp (stops, j, DSX) ;
	  int flag2 = 1 << (16 + j) ;
	  if (! (flag2 & flag3)) /* new stop */
	    if (wp->type || wp->end > 10)
	      {
		/* select the closest start */
		int b1 = 0 ;
		int besti = -1 ;
		for (int i = 0 ; i < iStart ; i++)
		{
		  vp = arrp (starts, i, DSX) ;
		  if (vp->a1 < wp->a2)
		    {
		      if (vp->end > wp->end/2)
			{
			  up = arrp (ss, vp->path, DSX) ;
			  if ((up->flag & flag2) && (!b1 || up->a1 > b1))
			    {
			      b1 = up->a1 ;
			      besti = i ;
			    }
			}
		    }
		}
	      if (besti >= 0)
		{
		  int flag = 1 << (besti) ;
		  vp = arrp (starts, besti, DSX) ;
		  flag3 |= flag2 ;
		  if (debug)
		    fprintf (stderr, "New stop flag %d, a1 = %d, a2 = %d, score = %d, flags %d / %d\n"
			     , nFlags, vp->a1, wp->a2, wp->end, flag, flag2 >> 16) ;
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag | flag2 ;
		  sFlag->nn = wp->nn ;
		}
	      else
		{
		  if (debug)
		    fprintf (stderr, "New stop flag %d, a1 = open, a2 = %d, score = %d, flag = %d >> 16\n"
			     , nFlags, wp->a2, wp->end, flag2) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag2 ;
		  sFlag->nn = wp->nn ;
		}
	      }
	}
    }
  if (0) arraySort (sFlags, sFlagsOrder) ;
  if (nFlags) /* check if we captured all introns */
    {
      DSX *up, *up2 ;
      int ii, mask = (1 << 16) - 1 ;
      for (ii = 0, up2 = arrp (ss2, 0, DSX) ; ii < iMax ; ii++, up2++)
	{
	  up = arrp (ss, up2->nn, DSX) ;
	  if (up->type & gI) 
	    { /* is this intron is part of a flag pair */
	      BOOL ok = FALSE ;
	      for (int i = 0 ; ! ok && i < nFlags ; i++)
		{
		  sFlag = arrp (sFlags, i, DSX) ;
		  flag3 = sFlag->flag ;
		  if ((up->flag & flag3) == flag3)
		    ok = TRUE ;
		}
	      if (! ok)
		{ /* create a new flag pair */
		  int flag1 = 0, flag2 = 0 ;
		  if (up->flag & mask) /* take the best start */
		    {
		      flag1 = 1 ;
		      for (int i = 0 ; i < 8 ; i++)
		      if ((up->flag >> i) & 1)
			break ;
		      else
			flag1 <<= 1 ;
		    }
		  if (up->flag >> 16) /* take the best stop */
		    {
		      flag2 = 1 << 16 ;
		      for (int i = 0 ; i < 8 ; i++)
			if ((up->flag >> (16+i)) & 1)
			break ;
		      else
			flag2 <<= 1 ;
		    }
		  if (debug)
		    fprintf (stderr, "New pair flag %d : %d/%d inherited from the intron %d/%d\n"
			     , nFlags, flag1, flag2 >> 16, up->a1, up->a2
			     ) ;
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->nn = ii ;
		  sFlag->flag = flag1 | flag2 ;
		}
	    }
	}
      arrayCompress (sFlags) ;
      nFlags = arrayMax (sFlags) ;
    }
  
  return nFlags ; 
} /* mrnaDesignFindStartEndPairs */

/**********************************************************************************/
static int mrnaDesignFindPaths (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  BOOL debug = FALSE ;
  int iFlag, path, nIntron, nStart, nStop, useCDS ;
  int eeMax = ssHappyFew (ds->exons) ;
  Array segs, segs2 ;
  DSX *vp, *vp2, *sFlag ; ;
  KEYSET ks = keySetHandleCreate (ds->h) ;
  Array ksPaths = arrayHandleCreate (20, KEYSET, ds->h) ;
  int flag = 0 ;
  Array sFlags = arrayHandleCreate (16, DSX, ds->h) ;
  Array starts = arrayHandleCreate (16, DSX, ds->h) ;
  Array stops = arrayHandleCreate (16, DSX, ds->h) ;
  
  arraySort (ksPaths, dsScoreOrder) ; /* for computer happiness */

  
  if (0 && debug)
     {
       fprintf (stderr, "mrnaDesignFindPaths start\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignFindPaths start done\n") ;
     }
  
  /* sort a copy by cover */
  segs = ds->exons ;
  segs2 = arrayHandleCopy (segs, ds->h) ;
  arraySort (segs2, dsCoverOrder) ;

  if (arrayMax (segs))
    mrnaDesignFindStartEndPairs (segs, segs2, sFlags, starts, stops) ;

  /* construct all paths
   * recursion on start-end pairs specified by flags
   */
  for (path = 0, iFlag = 0 ; iFlag < arrayMax (sFlags) ; iFlag++)
    for (useCDS = 1, sFlag = arrp (sFlags, iFlag, DSX) ; useCDS >=0 ; useCDS--)
    {
      int ii2 = 0, bestIi2 = -1 ;
      flag = sFlag->flag ;
      vp2 = arrp (segs2, 0, DSX) ;
      int bigScore = 0 ;
      int maxScore = 0 ;
      int mask = (1 << 16) - 1 ;
      int flag1 = flag & mask ;
      int flag2 = flag ^ flag1 ;
      int iStart = -1, iStop = -1 ;
      int xStart = 0, xStop = 0 ;
      
      if (flag1)
	{
	  int i = 0, flag3 = flag1 ;
	  for (i = 0 ; flag3 ; flag3 >>=1, i++)
	    if (flag3 & 1)
	      break ;
	  DSX *wp = arrp (starts, i, DSX) ;
	  iStart = wp->path ;
	  xStart = wp->a1 ;
	}
      if (flag2)
	{
	  int i = 0, flag3 = flag2 >> 16 ;
	  for (i = 0 ; flag3 ; flag3 >>=1, i++)
	    if (flag3 & 1)
	      break ;
	  DSX *wp = arrp (stops, i, DSX) ;
	  iStop = wp->path ;
	  xStop = wp->a2 ;
	}

      if (1)
	{
	  fprintf (stderr, "====== mrnaDesignFindPaths start with iFlag %d/%d, flag %d %d\n"
		   , iFlag, arrayMax(sFlags), flag & mask, flag >> 16) ;
	  if (iFlag == -1)
	    {
	      showDsFlags (sc, ds->exons, flag) ;
	      fprintf (stderr, "mrnaDesignFindPaths start done\n") ;
	    }
	}
      
      for (ii2 = 0, vp2 = arrp (segs2, 0, DSX) ; ii2 < eeMax ; ii2++, vp2++)
	{
	  /* find highest score for that path */
	  vp = arrp (segs, vp2->nn, DSX) ;
	  if ((vp->flag & flag) != flag)
	    continue ;
	  if (vp->cover > bigScore)
	    bigScore = vp->cover ;
	}
      if (sFlag->nn)
	{ 
	  bestIi2 = sFlag->nn ;
	  vp2 = arrp (segs2, sFlag->nn, DSX) ;
	  vp = arrp (segs, vp2->nn, DSX) ;
	  maxScore = vp->score ;
	}
      else
	for (ii2 = 0, vp2 = arrp (segs2, 0, DSX) ; ii2 < eeMax ; ii2++, vp2++)
	  {
	    /* find highest scoring seg not yet incorporated in a path */
	    vp = arrp (segs, vp2->nn, DSX) ;
	    if ((vp->flag & flag) != flag)
	      continue ;
	    if (! vp->score)
	      continue ;
	    if (vp->path && vp2->nn != iStart && vp2->nn != iStop)
	      continue ;
	    if (vp->path)
	      continue ;
	    if ((vp->type & gX) && vp->score < bigScore / 100)
	      continue ;
	    if ((vp->type & gI) && vp->score < bigScore / 500)
	      continue ;
	    if (vp2->nn != iStart && (vp->type & gDF))
	      continue ;
	    if (vp2->nn != iStop && (vp->type & gFF))
	      continue ;
	    if (useCDS && ! (vp->type &  gCompleteCDS))
	      continue ;
	    /* do not start on a mini exon */
	    if (vp->a2 < vp->a1 + 10 && (vp->type & gX) && (!(vp->type & (gS | gA))))
	      continue ;
	    if (vp->cover > maxScore)
	      {
		maxScore = vp->cover ;
		bestIi2 = ii2 ;
	      }
	  }
      
      if (bestIi2 < 0)
	continue ;
      for (ii2 = bestIi2, vp2 = arrp (segs2, ii2, DSX) ; maxScore && ii2 < eeMax ; ii2++, vp2++)
	{
	  vp2 = arrp (segs2, ii2, DSX) ;
	  vp = arrp (segs, vp2->nn, DSX) ;
	  if ((vp->flag & flag) != flag)
	    continue ;
	  if (! vp->score)
	    continue ;
	  if (vp->path && (!vp2 || (vp2->nn != iStart && vp2->nn != iStop)))
	    continue ;
	  if ((vp->type & gI) && vp->score < bigScore / 1000)
	    continue ;
	  if (vp2->nn != iStart && vp2->nn != iStop)
	    if ((vp->type & gX) && vp->score < bigScore / 100)
	      continue ;
	  if ((vp->type & gI) && vp->score < bigScore / 1000)
	    continue ;
	  if (vp2 && vp2->nn != iStart && (vp->type & gDF))
	    continue ;
	  if (vp2 && vp2->nn != iStop && (vp->type & gFF))
	    continue ;
	  if (useCDS && ! (vp->type &  gCompleteCDS))
	    continue ;
	  /* do not start on a mini exon */
	  if (vp->a2 < vp->a1 + 10 && (vp->type & gX) && (!(vp->type & (gS | gA))))
	    continue ;

	  path++ ;
	  if (0 || debug)
	    {
	      fprintf (stderr, "+++New path %d, flag = %d/%d, start on %d %d score=%d maxScore=%d\n", path, flag & mask, flag >> 16, vp->a1, vp->a2, vp->score, maxScore) ;
	    }
	  ks = keySetReCreate (ks) ;
	  keySet (ks, 0) = vp->nn ;
	  vp->donor = path ;
	  if (1) /* ! (vp->type & gB) */
	    {
	      int bestScore = ((vp->type & gI) ? vp->score/3 :  vp->score) ;
	      nIntron = (vp->type & gI) ? 3 : 0 ;
	      nStart = 1 ;
	      nStop = (vp->type & (gA | gReal3p)) ? 1 : 0 ;
	      mrnaDesignExtendDown (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, bestScore, maxScore, flag) ;
	      nStart = (vp->type & gS) ? 1 : 0 ;
	      nStop = 1 ;
	      mrnaDesignExtendUp (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, bestScore, maxScore, flag) ;
	      if (mrnaDesignIsNewPath (segs, ksPaths, ks, path, xStart, xStop, ds->h))
		mrnaDesignExport (s2m, sc, ds, segs, ks, path, smrnas) ;
	      else
		path-- ;
	    }
	  if(useCDS) /* we only want the best CDS guy */
	    break ;
	}
    }
  return path ;
} /* mrnaDesignFindPaths */

/**********************************************************************************/
/* flag the reads contributing introns that were actually rejected as void or poor */
static void mrnaDesignCleanUp (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii ;
  KEYSET ks ;
  vTXT txt ;
  const char *errors = 0 ;
  AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */

  if (! keySetMax (ds->reads)) 
    return ;
  ks = query (ds->reads, "(IS XY_* || IS XW_*) && (ct_ac || other || small_deletion) && (ct_ac > 1 || other > 1 || (ct_ac && other)) && ! composite > 50)") ;
  if (keySetMax (ks))
    {
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "Sequence %s\n-D Is_read\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), &errors, 0, ds->h) ;
    } 
  keySetDestroy (ks) ;
  ac_free (db) ;
  return ;
} /* mrnaDesignCleanUp */

/**********************************************************************************/
/**********************************************************************************/
/* s2m->>compositeDesignCovering was created by  mrnaDesignUsingCompositeStrategy below
 * the final mRNAs have been contructed by makemrna.c
 * the task is to reintrocuce the SL. polyA etc flags
 * and the quantitative covering of all elements
 */
static void mrnaDesignSetOneCompletenessFlag (S2M *s2m, SC* sc, HIT *up, SMRNA *smrna, Array estHits, BOOL isTop)
{
  Array covering = s2m->compositeDesignCovering ;
  int iss, issMax = covering ? arrayMax (covering) : 0  ;
  DSX *ssp ;
  unsigned int topFlags = gS | gS0 | gReal5p ;
  unsigned int bottomFlags = gReal3p | gA ;
  unsigned int cFlags = isTop ? topFlags : bottomFlags ;
  int a1 = up->a1, a2 = up->a2 ;
  if (!(up->type & gX) || !(up->type & cFlags))
    return ;
  up->type &= ~cFlags ;

  if (issMax)
    for (iss = 0, ssp = arrp (covering, 0, DSX) ; iss < issMax ; iss++, ssp++)
      {
	if (ssp->a1 > a2)
	  break ;
	if (ssp->a2 < a1)
	  continue ;
	if (isTop && ssp->a1 >= a1 && ssp->a1 <= a1 + 10)
	  {
	    up->type |= (ssp->type & cFlags) ;
	  }
	if (!isTop && ssp->a2 >= a2 - 10 && ssp->a2 <= a2)
	  {
	    up->type |= (ssp->type & cFlags) ;
	  }
      }

  return ;
} /*  mrnaDesignSetOneCompletenessFlag */

/**********************************************************************************/
void mrnaDesignSetCompletenessFlags (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  int iii, j ;
  SMRNA *smrna ;

  for (iii = 0 ; iii >-1 && iii < arrayMax(smrnas) ; iii++) 
    {
      smrna = arrp (smrnas, iii, SMRNA) ;
      if (smrna->hits && arrayMax (smrna->hits))      
	{
	  /* search the flags of the first exon */
	  for (j = 0, up = arrp (smrna->hits, 0, HIT) ; j < 1 && j < arrayMax(smrna->hits);  up++, j++)
	    mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, TRUE) ;
	  
	  /* search the polyA flags of the last exon */
	  for (j = arrayMax(smrna->hits) - 1 ; j >= 0 && (up = arrp (smrna->hits, j, HIT)) ; j = -1)
	    mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, FALSE) ;
	}
    }
  
  return ;
} /*  mrnaDesignSetCompletenessFlags */

/**********************************************************************************/
/**********************************************************************************/

BOOL mrnaDesignUsingCompositeStrategy (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  AC_HANDLE h = 0 ;
  DS ds ;
  int i ;
  HIT *up ;
  BOOL ok = FALSE ;
  BOOL debug =  FALSE ;

  if (0 || arrayMax (smrnas) < 1) 
    return FALSE; 

  /* check if we are using composite reads */
  if (arrayMax (s2m->geneHits))
    for (i = 0, up = arrp (s2m->geneHits, 0, HIT), ok = FALSE ; !ok && i < arrayMax (s2m->geneHits) ; i++, up++)
      if (keyFindTag (up->est, _Composite))
	ok = TRUE ;
  if (! ok)
    return FALSE ;

  memset (&ds, 0, sizeof (ds)) ;
  ds.h = h = ac_new_handle () ;

  mrnaDesignGetElements (s2m, sc, &ds, smrnas) ;
  if (debug)
    {
      showDs (sc, ds.exons) ;
      showDs (sc, ds.introns) ; 
    }
  ok = FALSE ;
  if (mrnaDesignGetGraphElements (&ds, s2m, sc, gmrna, smrnas))
    {
      ok = TRUE ;
      smrnas = arrayReCreate (smrnas, 12, SMRNA) ;
      mrnaDesignFindPaths (s2m, sc, &ds, smrnas) ;
      if (0) mrnaDesignSetCompletenessFlags (s2m, sc, gmrna, smrnas) ;
      mrnaDesignCleanUp (s2m, sc, &ds, smrnas) ;
    }
  else
    arrayMax (smrnas) = 0 ;
  ac_free (h) ;
  return ok ;
} /*  mrnaDesignUsingCompositeStrategy */

/**********************************************************************************/
/**********************************************************************************/
#ifdef JUNK
void mrnaDesignAbsorbInIntron (void) 
{
  const char *errors = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_DB db = ac_open_db (0, 0) ; /* local database, cannot fail */
  const char *qqG = "select chr,a1,a2,x,x1,x2,c from x in ?Sequence where x like \"XH_*\" or x like \"XG_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], x1 in x->composite, x2 in x1[1], c in x2[1]" ;
  const char *qqI = "select chr,a1,a2,c from x in ?Sequence where x like \"XI_*\" or x like \"XJ_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], c in x->composite" ;
  const char *qqA = "select chr,a1,a2,x,c from x in ?Sequence where x like \"XH_*\" or x like \"XG_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], c in x->composite" ;

  AC_TABLE aI = ac_bql_table (db, qqI, 0, 0, &errors, h) ;
  AC_TABLE aG = ac_bql_table (db, qqG, 0, 0, &errors, h) ;
  AC_TABLE aA = ac_bql_table (db, qqA, 0, 0, &errors, h) ;

  int iI = 0, iIMax = aI->rows ;
  int iA = 0, iAMax = aA->rows ;
  int iG = 0, iGMax = aG->rows ;

  for (iI = 0 ; iI < iIMax ; iI++) 
    {
      KEY iChr = ac_table_key (aI, iI, 0, 0) ;
      int i1 = ac_table_int (aI, iI, 1, 0) ;
      int i2 = ac_table_int (aI, iI, 2, 0) ;
      int iC = ac_table_int (aI, iI, 3, 0) ;

      continue ;
      
      for (iA ; iA > 0 ; iA--)
	{
	  KEY aChr = ac_table_key (aA, iA, 0, 0) ;
	  int a1 = ac_table_int (aA, iA, 1, 0) ;
	  int a2 = ac_table_int (aA, iA, 2, 0) ;
	  if (aChr < iChr || a2 < i1)
	    break ;
	}
      for (iA ; iA < iAMax ; iA++)
	{
	  KEY aX = ac_table_key (aA, iA, 0, 0) ;
	  KEY aChr = ac_table_key (aA, iA, 1, 0) ;
	  int a1 = ac_table_int (aA, iA, 2, 0) ;
	  int a2 = ac_table_int (aA, iA, 3, 0) ;
	  int aC = ac_table_int (aA, iA, 4, 0) ;
	  if (aChr > iChr || a1 > i2)
	    break ;
	}
    }
  ac_free (h) ;
  return ;
} /* mrnaDesignAbsorbInIntron */
#endif
/**********************************************************************************/
/**********************************************************************************/

BOOL mrnaDesignCutOnePreMrna (KEY gene, KEY green, KEYSET ks)
{
  BOOL new = FALSE ;
  int nn = ks ? keySetMax (ks) : 0 ;
  const char *errors = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_DB db = ac_open_db (0, 0) ; /* local database, cannot fail */
  AC_OBJ Gene = ac_get_obj (db, "Transcribed_gene", name(gene), h) ;
  const char *qqMap = "select chr, a1, a2 from chr in  @->IntMap, a1 in chr[1], a2 in a1[1]" ;
  AC_TABLE map = ac_obj_bql_table (Gene, qqMap, 0, &errors, h) ;
  AC_OBJ Green = ac_get_obj (db, "Sequence", name(green), h) ;
  AC_TABLE greenMap = ac_obj_bql_table (Green, qqMap, 0, &errors, h) ;
  const char *qq2 = "select a1, a2, r, x1, x2 from a1 in @->Assembled_from, a2 in a1[1], r in a2[1], x1 in r[1], x2 in x1[1]" ;
  AC_TABLE afrom = ac_obj_bql_table (Gene, qq2, 0, &errors, h) ;
  AC_TABLE covers = 0 ;
  int ir, g1 = 0, g2 = 0, gr1 = 0, gr2 = 0, z1 = 0, z2 = 0, dz = 0, xStop = 0, xStart = 0 ;
  const char *rStart = "Xends_ELF" ;
  const char *rStop  = "Xends_ERF" ;
  int nBreaks = 0 ;
  int nGoodBreaks = 0 ;
  KEYSET breaks = keySetHandleCreate (h) ;
  KEYSET bCovers = keySetHandleCreate (h) ;
  
  for (ir = 0 ; map && ir < map->rows ; ir++)
    {
      g1 = ac_table_int (map, ir, 1, 0) ;
      g2 = ac_table_int (map, ir, 2, 0) ;
    }
  if (g1 > g2)
    {
      rStart = "Xends_ERR" ;
      rStop  = "Xends_ELR" ;
    }

  for (ir = 0 ; map && ir < map->rows ; ir++)
    {
      gr1 = ac_table_int (greenMap, ir, 1, 0) ;
      gr2 = ac_table_int (greenMap, ir, 2, 0) ;
    }

  /* get green offset in tg cooordinates */
  if (! afrom || ! map)
    goto done ;
  for (ir = 0 ; ir < afrom->rows ; ir++)
    if (! strcmp (name (green), ac_table_printable (afrom, ir, 2, "toto")))
      { /* green coords in tg */
	z1 = ac_table_int (afrom, ir, 0, 0) ;
	z2 = ac_table_int (afrom, ir, 1, 0) ;
	int z11 = ac_table_int (afrom, ir, 3, 0) ;
	dz = z1 - z11 ;
	break ;
      }
  if (! z1 || ! z2)
    goto done ;

  /* locate potential breakpoints in tg coordinates */
  int cStart = 0, cStop = 0 ;
  for (ir = 0 ; ir < afrom->rows ; ir++)
    {
      if (! strncmp (rStop, ac_table_printable (afrom, ir, 2, "toto"), 9))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStop = ac_table_int (afrom, ir, 1, 0) ;
	  cStop = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	}
      if (! strncmp ("XSL", ac_table_printable (afrom, ir, 2, "toto"), 3))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStart = ac_table_int (afrom, ir, 0, 0) ;
	  cStart = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	  if (xStop && xStop < xStart + 30)
	    {
	      int bb = (xStop + xStart) / 2 ;
	      if (bb > z1 + 50 && bb < z2 - 50)
		{
		  keySet (breaks, nBreaks) = bb - dz ;
		  keySet (bCovers, nBreaks) = -100 ;
		  nBreaks++ ;
		}
	      xStop = 0 ;
	    }
	}
      if (! strncmp (rStart, ac_table_printable (afrom, ir, 2, "toto"), 9))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStart = ac_table_int (afrom, ir, 0, 0) ;
	  cStart = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	  if (xStop && xStop < xStart + 30)
	    {
	      int bb = (xStop + xStart) / 2 ;
	      if (bb > z1 + 50 && bb < z2 - 50)
		{
		  keySet (breaks, nBreaks) = bb - dz ;
		  keySet (bCovers, nBreaks) = (cStop + cStart) / 2 ;
		  nBreaks++ ;
		}
	      xStop = 0 ;
	    }
	}
    }

  if (nBreaks)
    { /* check if the candidate breakpoint falls in an exon cover dip */
      const char *qq1 = "select x1, x2, c from x1 in @->Composite, x2 in x1[1], c in x2[1]" ;
      int ir = 0 ;

      keySet (breaks, nBreaks++) = z2 ;
      covers = ac_obj_bql_table (Green, qq1, 0, &errors, h) ;
      for (int ib = 0 ; ib < nBreaks ; ib++)
	{
	  int u1, u2, c, bb = keySet (breaks, ib) ;
	  for (; ir >= 0 ; ir--)
	    {
	      u1 = ac_table_int (covers, ir, 0, 0) ;
	      if (u1 <= bb)
		break ;
	    }
	  if (ir < 0) ir = 0 ;
	  for (; ir < covers->rows ; ir++)
	    {
	      u2 = ac_table_int (covers, ir, 1, 0) ;
	      c = ac_table_int (covers, ir, 2, 0) ;
	      if (u2 >= bb)
		break ;
	    }
	  if (ir >= covers->rows)
	    { ir-- ; continue ; }
	  u1 = ac_table_int (covers, ir, 0, 0) ;
	  c = ac_table_int (covers, ir, 2, 0) ;	      
	  if (u1 <= bb && u2 >= bb)
	    {
	      int ok = 0 ;
	      int bc = keySet (bCovers, ib) ;
	      if (bc == -100) ok = 3 ;
	      if (! ok || 4 * c < bc)
		{
		  for (int j = 1 ; j < 4  && j + ir < covers->rows ; j++)
		    if (1 * c < ac_table_int (covers, ir + j, 2, 0))
		      { ok = 1 ; break ; }
		  for (int j = 1 ; j < 4  && ir - j >= 0 ; j++)
		    if (1 * c < ac_table_int (covers, ir - j, 2, 0))
		      { ok += 2 ; break ; }
		}
	      if (ok != 3)
		keySet (breaks, ib) = 0 ;  /* drop it */
	      else
		nGoodBreaks++ ;
	    }
	}
    }

 done:
  if (nGoodBreaks)
    {
      char *dna = ac_obj_dna (Green, h) ;
      int d2 = strlen (dna) ;
      int d11 = 1, d12 ;
      vTXT txt = vtxtHandleCreate (h) ;

      vtxtPrintf (txt, "Sequence %s\n-D Is_read\n-D From_gene\n-D In_mRNA\n\n", name (green)) ;
      vtxtPrintf (txt, "cDNA_clone %s\n-D From_gene\n-D In_mRNA\n\n", name (green)) ;
      vtxtPrintf (txt, "-D DNA %s\n\n", name (green)) ;
      vtxtPrintf (txt, "-D cDNA_clone %s\n\n", name (green)) ;
      vtxtPrintf (txt, "-D Sequence %s\n\n", name (green)) ;
      /* compress */
      keySetSort (breaks) ;
      int jb = 0, bb0 = 50 ;
      for (int ib = 0 ; ib < nBreaks ; ib++)
	{
	  int bb = keySet (breaks, ib) ;
	  if (bb && bb > bb0 && bb < z2 - 50)
	    { keySet (breaks, jb++) = bb ; bb0 = bb + 50 ; }
	}
      keySet (breaks, jb++) = z2 ;
      nBreaks = keySetMax (breaks) = jb ;
      for (int ib = 0 ; ib < nBreaks ; ib++)
	{
	  char *nam = hprintf (h, "%s.%d", name (green), ib) ;
	  int bb = keySet (breaks, ib) ;

	  fprintf (stderr, "Breaking %s at position %d\n", name (green), bb) ;
	  
	  d12 = bb - 10 ;
	  if (d12 < 1)
	    d12 = 1 ;
	  if (d12 > d2)
	    d12 = d2 ;
	  dna = ac_zone_dna (Green, d11, d12, h) ;

	  vtxtPrintf (txt, "DNA %s\n%s\n\n", nam, dna) ;
	  vtxtPrintf (txt, "cDNA_clone %s\nRead %s\nFRom_gene %s\n\n", nam, nam, name (gene)) ;
	  vtxtPrintf (txt, "Sequence %s\nIs_read\nForward\n", nam) ;
	  if (map->rows)
	    vtxtPrintf (txt, "IntMap %s %d %d\n"
			, ac_table_printable (map, 0, 0, 0)
			, gr1 < gr2 ? gr1 + d11 - 1 : gr1 - d11 + 1
			, gr1 < gr2 ? gr1 + d12 - 1 : gr1 - d12 + 1
			) ;

	  for (int ir = 0; ir < covers->rows ; ir++)
	    {
	      int u1 = ac_table_int (covers, ir, 0, 0) ;
	      int u2 = ac_table_int (covers, ir, 1, 0) ;
	      int c = ac_table_int (covers, ir, 2, 0) ;
	      if (u1 > d12)
		break ;
	      if (u1 < d11) u1 = d11 ;
	      if (u2 > d12) u2 = d12 ;
	      if (u1 <= u2)
		vtxtPrintf (txt, "Composite %d %d %d\n", u1 - d11 + 1, u2 - d11 + 1, c) ;
	    }
	  vtxtPrintf (txt, "\n") ;
	  d11 = bb + 10 ;
	  if (d11 < 1)
	    d11 = 1 ;
	  if (d11 > d2)
	    d11 = d2 ;

	}
      if (vtxtPtr (txt))
	ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;
    }
  
  if (ks)
    keySet (ks, nn++) = green ;
  ac_free (h) ;
  ac_free (db) ;
  
  return new ;
} /* mrnaDesignCutOnePreMrna */

/**********************************************************************************/
/* cut the premra XG/XH exons reads on strong termination signals */

static void mrnaDesignCutGenePreMrna (KEY gene)
{
  KEYSET greens = queryKey (gene, ">Read ; Is_Read && (IS XG_* || IS XH_*)") ;
  for (int i = 0 ; i < keySetMax (greens) ; i++)
    {
      KEY green = keySet (greens, i) ;
      mrnaDesignCutOnePreMrna (gene, green, 0) ;
    }
  return ;
} /* mrnaDesignCutGenePreMrna */

/**********************************************************************************/
/* cut the premra XG/XH exons reads on strong termination signals */
void mrnaDesignCutAllPreMrna (KEYSET tgs)
{
  KEYSET genes = query (tgs, "CLASS Transcribed_genes") ;

  if (genes)
    for (int i = 0 ; i < keySetMax (genes) ; i++)
      mrnaDesignCutGenePreMrna (keySet (genes, i)) ;
  keySetDestroy (genes) ;
  
  return ;
} /* mrnaDesignCutAllPreMrna */
 
/**********************************************************************************/
/**********************************************************************************/

int mrnaDesignCleanEcho (void)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *errors = 0 ;
  AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */
  KEY _ps = str2tag ("ProcessStatus") ;  
  const char *qqxg = "select r, chr, a1, a2, s, x1, x2, c from r in ?Sequence where r like \"XG_*\" && !r#ProcessStatus , chr in r->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, x1 in r->Composite, x2 in x1[1], c in x2[1]" ;
  AC_TABLE xxg = ac_bql_table (db, qqxg, 0, 0, &errors, h) ;
  const char *qqxxi = "select r, chr, a1, a2, s, c from r in ?Sequence where r#composite && r like \"XI_*\", I in r->intron, chr in I->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, c in r->Composite" ;
  AC_TABLE xxi = ac_bql_table (db, qqxxi, 0, 0, &errors, h) ;
  const char *qqxxs = "select r, chr, a1, a2, s, c from r in ?Sequence where r#composite && !r#ProcessStatus && ! r like \"XG_*\" && ! r like \"XH_*\", chr in r->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, c in r->Composite" ;
  AC_TABLE xxs = ac_bql_table (db, qqxxs, 0, 0, &errors, h) ;
  KEYSET maps = keySetHandleCreate (h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  int ccc = 0, nKilled = 0 ;
  KEY oldxg = 0 ;
      
  if (xxg) for (int ir = 0 ; ir < xxg->rows ; ir++)
    {
      KEY map = ac_table_key (xxg, ir, 1, 0) ;
      if (map)
	keySetInsert (maps, map) ;
    }
  for (int ii = 0 ; ii < keySetMax (maps) ; ii++)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      KEY map = keySet (maps, ii) ;
      Array pm, plus = arrayHandleCreate (100000, int, h1) ;
      Array minus = arrayHandleCreate (100000, int, h1) ;

      /* construct a wiggle for each strand */
      for (int ir = 0 ; ir < xxi->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxi, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int a1 = (ac_table_int (xxi, ir, 2, 0) + 4) / 10 ;
	      int a2 = (ac_table_int (xxi, ir, 3, 0) + 4) / 10 ;
	      KEY strand = ac_table_key (xxi, ir, 4, 0) ;
	      int c = ac_table_int (xxi, ir, 5, 0) ;

	      if (strand == _Reverse) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      pm = a1 < a2 ? plus : minus ;
	      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      for (int a = a1 + 1 ; a < a2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  array (plus, a, int)  = c + c1 ; /* introns are cumulative */
		  array (minus, a, int)  = c + c1 ; /* introns are cumulative */
		}
	    }
	}

      for (int ir = 0 ; ir < xxg->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxg, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int b1, b2, a1 = ac_table_int (xxg, ir, 2, 0) ;
	      int a2 = ac_table_int (xxg, ir, 3, 0) ;
	      int x1 = ac_table_int (xxg, ir, 5, 0) ;
	      int x2 = ac_table_int (xxg, ir, 6, 0) ;
	      int c = ac_table_int (xxg, ir, 7, 0) ;

	      if (x1 > x2) continue ;
	      if (a1 <= a2) { b1 = a1 + x1 - 1 ; b2 = a1 + x2 - 1 ; }
	      else { b1 = a1 - x1 + 1 ; b2 = a1 - x2 + 1 ; }
	      if (b1 > b2) { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }
	      b1 = (b1 + 4) / 10 ; b2 = (b2 + 4) / 10 ;
	      pm = a1 < a2 ? plus : minus ;
	      for (int a = b1 + 1 ; a < b2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  if (c1 < c) array (pm, a, int)  = c ;
		}
	    }
	}

      /* now clean up */
      for (int ir = 0 ; xxs && ir < xxs->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxs, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int a1 = (ac_table_int (xxs, ir, 2, 0) + 4) / 10 ;
	      int a2 = (ac_table_int (xxs, ir, 3, 0) + 4) / 10 ;
	      KEY strand = ac_table_key (xxs, ir, 4, 0) ;
	      int n = 0, c = ac_table_int (xxs, ir, 5, 0) ;
	      long int ccL = 0 ;
	      
	      if (strand == _Reverse) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      pm = a1 > a2 ? plus : minus ; /* opposite strand wiggle */
	      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      for (int a = a1 + 1 ; a < a2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  ccL += c1 ; n++ ;
		}
	      if (n) ccL /= n ; /* average premessenger of other strand */
	      ccL /= 100 ;
	      if (ccL > c)
		ccc = 0 ;
	      else
		ccc = c - ccL ;
	      if (ccc < 0) ccc = 0 ;
	      vtxtPrintf (txt, "Sequence %s\nProcessStatus  %d\n-D Composite\n", ac_table_printable (xxs, ir, 0, "toto"), c) ;
	      if (ccc > 0)
		vtxtPrintf (txt, "Composite %d\n\n", ccc) ;
	      else
		{
		  nKilled++ ;
		  vtxtPrintf (txt, "-D Is_read\nColour BROWN\n\n", c) ;
		}
	    }
	}

      ccc = 0 ;
      oldxg = 0 ;
      for (int ir = 0 ; xxg && ir < xxg->rows ; ir++)
	{
	  KEY xg = ac_table_key (xxg, ir, 0, 0) ;
	  KEY map1 = ac_table_key (xxg, ir, 1, 0) ;
	  int cc = 0 ;
	  
	  if (xg && map == map1 && ! keyFindTag (xg, _ps))
 	    {
	      int b1, b2, n = 0, a1 = ac_table_int (xxg, ir, 2, 0) ;
	      int a2 = ac_table_int (xxg, ir, 3, 0) ;
	      int x1 = ac_table_int (xxg, ir, 5, 0) ;
	      int x2 = ac_table_int (xxg, ir, 6, 0) ;
	      int c = ac_table_int (xxg, ir, 7, 0) ;

	      if (!x1 || x1 > x2) continue ;
	      if (a1 < a2) { b1 = a1 + x1 - 1 ; b2 = a1 + x2 - 1 ; }
	      else { b1 = a1 - x1 + 1 ; b2 = a1 - x2 + 1 ; }
	      if (b1 > b2) { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }
	      b1 = (b1 + 4) / 10 ; b2 = (b2 + 4) / 10 ;
	      pm = a1 > a2 ? plus : minus ;
	      for (int a = b1 + 1 ; a < b2 ; a++)
		{
		  if (a < 0) invokeDebugger () ;
		  int c1 = array (pm, a, int) ;
		  cc += c1 ; n++ ;
		}
	      if (n) cc /= n ; /* average premessenger of other strand */
	      cc = c - cc / 100 ; /* new value */
	      if (cc < 0) cc = 0 ;
	      if (xg != oldxg)
		{
		  if (oldxg && ! ccc)
		    {
		      nKilled++ ;
		      vtxtPrintf (txt, "-D Is_read\nColour BROWN\n") ;
		    }
		  oldxg = xg ;
		  ccc = 0 ;
		  vtxtPrintf (txt, "\nSequence %s\nProcessStatus\n-D Composite\n", name (xg)) ;
		}
	      if (cc > 0 && x1 <= x2)
		{
		  ccc += cc ;
		  vtxtPrintf (txt, "Composite %d %d %d\n", x1, x2, cc) ;
		}
	    }
	}
      if (oldxg && ! ccc)
	{
	  nKilled++ ;
	  vtxtPrintf (txt, "-D is_read\nColour BROWN\n") ;
	}
      vtxtPrintf (txt, "\n") ;
    }

  if (vtxtPtr (txt))
    ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;

  fprintf (stderr, "mrnaDesignCleanEcho kill %d reads in %d maps\n", nKilled, keySetMax (maps)) ;
  ac_free (h) ;

  return nKilled ;
} /* mrnaDesignCleanEcho */

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
