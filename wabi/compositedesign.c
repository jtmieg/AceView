#include "ac.h"
#include "cdna.h"
#include "makemrna.h"

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
typedef struct dsVertexStruct { int a1, a2, donor, acceptor, cover, path, type, score, nn, clipable ; } DSX ;

/**********************************************************************************/

static void showDs (SC *sc, Array aa)
{
  DSX *up ;
  int ii ;
  
  if (aa)
    {
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	{
	  fprintf( stderr, "%d:: %d %d\t %d %d\t d=%d a=%d c=%d s=%d path=%d\t"
		   , ii, up->a1, up->a2, sc->a1 + up->a1 - 1, sc->a1 + up->a2 - 1
		   , up->donor, up->acceptor, up->cover, up->score, up->path
		   ) ;
	  
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
	  if (gStolen & up->type) fprintf (stderr, "Stolen ") ;
	  if (gPredicted & up->type) fprintf (stderr, "Predicted ") ;

	  if (up->clipable) fprintf (stderr, " *") ;
	  
	  fprintf (stderr, "\n") ;
	}
      showDs (0, 0) ;
    }
  return ;
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

static int dsCoverOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->cover - b->cover ; if (n) return -n ; /* large cover first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsScoreOrder */

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
static void mrnaDesignGetGraphElements (DS *ds, S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  BOOL debug = TRUE ;
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
			  }
		      
		      if (s > 0)
			{
			  DSX *ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ;
			  s *= 5 ;
			  if (up->type & gI)
			    { nXI++; ssp->type = gI ; ssp->cover = s ; }
			  else if (up->type & gX)
			    { nXE++; ssp->type = gX ; ssp->cover = s ; }
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
		      ssp->type = gX | gCompleteCDS | gFF ; /* useless */
		      ssp->type = gX | gFF ;
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
			  
			  ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ;
			  ssp->cover = s ;
			  ssp->clipable = TRUE ;
			  if (! strncmp (name(up->est), "XA_", 3))
			    {
			      nXA++ ;
			      ssp->a1 = ssp->a2 ;
			      ssp->type = gX | gA ;
			      ssp->clipable = FALSE ;
			    }
			  else if (! strncmp (name(up->est), "XSL", 3))
			    {
			      nXSL++ ;
			      ssp->a2 = ssp->a1 ;
			      ssp->type = gX | gS ;
			      ssp->clipable = FALSE ;
			    }
			  else if (! strncmp (name(up->est), "Xends_ELF_", 9) && isDown) 
			    { nXends++ ; ssp->type = (gReal5p | gDF | gX) ; if (0) { ssp->a1 += 1 ; ssp->a2 = ssp->a1 + 30 ; }}
			  else if (! strncmp (name(up->est), "Xends_ERF_", 9) && isDown) 
			    { nXends++ ; ssp->type = (gReal3p | gFF | gX) ; if (0) { ssp->a2 -= 1 ; ssp->a1 = ssp->a2 - 30 ; }} 
			  else if (! strncmp (name(up->est), "Xends_ERR_", 9) && ! isDown) 
			    { nXends++ ; ssp->type = (gReal5p | gDF | gX) ; if (0) { ssp->a1 += 1 ; ssp->a2 = ssp->a1 + 30 ; }} 
			  else if (! strncmp (name(up->est), "Xends_ELR_", 9) && ! isDown) 
			    { nXends++ ; ssp->type = (gReal3p | gFF | gX) ; if (0) { ssp->a2 -= 1 ; ssp->a1 = ssp->a2 - 30 ; }}
			  
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
  arraySort (ss, dsA1Order) ;
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
		       if (ssp->cover > 1000 * ssp1->cover)
			 ssp1->cover = 0 ;
		       if (ssp1->cover > 1000 * ssp->cover)
			 ssp->cover = 0 ;
		     }
		 }
	     }
	 }
     }
  
  iMax = ssHappyFew (ss) ;
  if (! iMax)
    goto done ;

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
		{
		  ssp1->cover = ssp->cover ;
		  ssp1->type &= (~(gDF | gFF)) ;
		}
	      ssp->type |= ssp1->type ;
	      if (ssp->cover < ssp1->cover)
		ssp->cover = ssp1->cover ;
	      ssp1->cover = 0 ;
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
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if ((ssp->type & gX) && ! (ssp->type & (gA | g3 | gReal3p)))
	     {
	       int i ;
	       DSX *ssp1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; i < iMax ; ssp1++, i++)
		 {
		   if (! ssp1->cover) continue ;
		   if (ssp1->type & gI) continue ;
		   if (ssp1->type & (g5 | gS | gS0)) break ;
		   if (ssp1->cover != ssp->cover) break ;
		   if (ssp1->a1 != ssp->a2 + 1) break ;
		   ssp->a2 = ssp1->a2 ; ssp->type |= ssp1->type ;
		   ssp1->cover = 0 ;
		   if (ssp->type & (gA | g3 | gReal3p)) break ;
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

 
  if (iMax) /* absorb in exons 10% of intron support */
     {
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       DSX *ssp1 ;
	       int i, c = ssp->cover/10 ;
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
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 10%% of intron score\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 10%% of intron score done\n") ;
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
		 if ((ssp1->type & gX) && ssp1->a2 < b1 + 10 && ! (ssp1->type & (gA | g3)))
		   { if ((ssp1[1].type & gX) && ssp1[1].a1 == ssp1->a2 + 1) ssp1->cover = ssp1[1].cover ; else ssp1->cover = 0 ;}
	       for (i = 1, ssp1 = ssp + 1 ; ssp1->a1 <= b2 && jj + i < iMax ; ssp1++, i++)
		 if ((ssp1->type & gX) && ssp1->a1 > b2 - 10 && ! (ssp1->type & (gS | g5) ))
		   { if ((ssp1[-1].type & gX) && ssp1[-1].a2 == ssp1->a1 - 1) ssp1->cover = ssp1[-1].cover ; else ssp1->cover = 0 ;}
	     }
	 }
     }

  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after absorb leaking exons\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after absrb leaking exons done\n") ;
    }
  
  
  /* regularize the real5p */
  if (1)
    {
      for (int ii = iMax -1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal5p)
	    {
	      DSX *ssp1, *ssp0 = ssp ;
	      int jj, bestCover = ssp->cover, b1 = ssp->a1 - 1 ;

	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 && ssp1->a2 == b1 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gA | gReal3p)) break ;
		  if (100 * ssp1->cover < bestCover)
		    {
		      for ( ; jj >= 0 && ssp1->a2 == b1 && (ssp1->type & gReal5p) ; jj--, ssp1--)
			{
			  if (ssp1->type & gI) continue ;
			  b1 = ssp1->a1 - 1 ;
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
		      b1 = ssp1->a1 - 1 ;
		    }
		}
	    }
	}
    }
  
  /* regularize the real3p */
  if (1)
    {
      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal3p)
	    {
	      DSX *ssp1, *ssp0 = ssp ;
	      int jj, bestCover = ssp->cover, b2 = ssp->a2 + 1 ;

	      for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gS | gS0 | gReal5p)) break ;
		  if (100 * ssp1->cover < bestCover) break ;
		  if (ssp1->cover > bestCover) bestCover = ssp1->cover ;
		  ssp0->type &= ~gReal3p ;
		  ssp1->type |= gReal3p ;
		  ssp0 = ssp1 ;
		  b2 = ssp1->a2 + 1 ;
		}
	    }
	}
    }
  
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5p\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5 done\n") ;
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

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL done\n") ;
    }

  if (iMax) /* merge polyA in previous exon */
    {
      for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
	{
	  if (ssp->cover <= 0)   continue ;
	  if ((ssp->type & gX) && (ssp->type & (gA | g3)) && ! (ssp->type & (g5 | gS)))
	    {
	      int i ;
	      DSX *ssp1 ;
	      for (i = jj - 1, ssp1 = arrp (ss, i, DSX) ; i >= 0 ; i--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover <= 0) continue ;
		  if (ssp1->a2 >= ssp->a1 - 120)
		    {
		      ssp1->type |=g3 ;
		      if (1 && ssp->a1 > ssp1->a2 +1)
			ssp->a1 = ssp1->a2 + 1 ;
		    }
		  break ;
		}
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
	   if (ssp->type & gReal5p)
	     {
	       break ;
	     }
	   if (0 && 2 * ssp->cover < bestCover)
	     break ;
	 }
       if (ssp->type & gReal5p)
	 {
	   /* I found a real5p, kill all exons above */
	   for (int i = 0 ; i < jj ; i++)
	     if (ssp[-i-1].a2 < ssp->a1 && (1*ssp[-i - 1].cover <= bestCover || (ssp[-i-1].type & (gDF | gFF))))
	       ssp[-i - 1].cover = 0 ;
	     else
	       break ;
	   if (! (ssp->type & g5))
	     {
	       ssp->a1 -= 15 ; /* move back to original position */
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
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	 }
       if (ssp->type & gA)
	 {
	   /* I found a polyA, kill all exons below */
	   for (int i = jj + 1 ; i < iMax ; i++)
	     ssp[i - jj].cover = 0 ;
	 }
     }

 done:
  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  ds->exons = ss ;
  
  if (debug)
     {
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns done\n") ;
     }

   ac_free (h) ;
   return ;
} /* mrnaDesignGetGraphElements */

/**********************************************************************************/
/* Rationalize the exon flags */
static void mrnaDesignCleanExons (S2M *s2m, SC *sc, DS *ds)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = TRUE ;
  int ii ;
  DSX *up ;
  int iiMax = ssHappyFew (ds->exons) ;


  /* create and boost some scores */
  for (ii = 0, up = arrp (ds->exons, 0, DSX) ; ii < iiMax ; ii++, up++)
    {
      up->nn = ii ;
      up->score = up->cover ;
      if (up->type & gI) up->score *= 3 ;
      if (up->type & gA) up->score *= 3 ;
    }

  if (iiMax == 1)   
    {       /* eliminate weak short cloud exons */
      up = arrayp (ds->exons,0, DSX) ;
      int da = up->a2 > up->a1 ? up->a2 - up->a1 : up->a1 - up->a2 ;
      if (da < 50 && up->score < 10)
	iiMax = arrayMax (ds->exons) = 0 ;
    }

  if (iiMax)
    {  /* erase consecutive gReal5p */
      for (ii = 0, up = arrp (ds->exons, 0, DSX) ; ii < iiMax ; ii++, up++)
	if (up->type & gReal5p)
	  { 
	    DSX *vp ;
	    int i ;
	    
	    for (i = ii + 1, vp = up+1 ; i < iiMax && (vp->type & gReal5p)  ; vp++, i++)
	      vp->type &= ~gReal5p ;
	    
	    /* if there is a  previous exon, shift the flag vp */
	    if (ii > 0 && (up[-1].type & gX))
	      {
		if (((up[-1].type & gS) || 10*up[-1].score > up->score)&& up[-1].a1 + 25 > up->a1)
		  {
		    up->type &= ~gReal5p ;
		    up[-1].type |= gReal5p ;
		    if (up[-1].score < up->score)
		      up[-1].score = up->score ;
		  }
	      }
	  }
    }

  if (iiMax)
    {  /* erase consecutive gReal3p */
      for (ii = iiMax - 1, up = arrp (ds->exons, ii, DSX) ; ii > 0 ; ii--, up--)
	if (up->type & gReal3p)
	  { 
	    DSX *vp ;
	    int i ;
	    
	    for (i = ii - 1, vp = up - 1 ; i >= 0 && (vp->type & gReal3p)  ; vp--, i--)
	      vp->type &= ~gReal3p ;
	    
	    /* if there is a following exon, shift the flag down */
	    if (ii < iiMax -1 && (up[1].type & gX))
	      {
		if (((up[1].type & gA) || 10*up[1].score > up->score)&& up[1].a2 - 25 < up->a2)
		  {
		    up->type &= ~gReal3p ;
		    up[-1].type |= gReal3p ;
		    if (up[1].score < up->score)
		      up[1].score = up->score ;
		  }
	      }
	  }
    }
  
  /* if we have   2 introns with same acceptor
   * the exon bit going to the second donor
   * should have the support of the short intron
   */
  if (iiMax)
    for (ii = 0, up = arrp (ds->exons, 0, DSX) ; ii < iiMax ; ii++, up++)
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
  if (iiMax)
    for (ii = iiMax - 1, up = arrp (ds->exons, ii, DSX) ; ii >= 0 ; ii--, up--)
      {
	if (up->type & gA)
	  { 
	    DSX *vp ;
	    int i, jj, bestScore = up->score ;
	    int b2 = up->a1 - 1 ;
	    for (i = ii - 1, vp = arrp (ds->exons, i, DSX) ; i >= 0 ; vp--, i--)
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
	    for (i = jj + 1, vp = arrp (ds->exons, i, DSX) ; i <= ii ; i++, vp++)
	      vp->score = bestScore ;
	  }
      }  
  
  /* overlapping elements are killed if below 1% */
  iiMax = ssHappyFew (ds->exons) ;
  if (1)
    {
      for (ii = 0, up = arrp (ds->exons, 0, DSX) ; ii < iiMax ; ii++, up++)
	if (up->score > 100)
	  {
	    int jj ;
	    DSX *vp ;
	    for (jj = ii - 1, vp = up -1 ; jj >= 0 && vp->a1 == up->a1 ; vp--, jj--)
	      ;
	    for (jj++, vp++ ; jj < iiMax && vp->a1 < up->a2 ; jj++, vp++)
	      if (100 * vp->score < up->score)
		vp->cover = vp->score = 0 ;
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

  iiMax = ssHappyFew (ds->exons) ;
  /* number the segments */
  for (ii = 0, up = arrp (ds->exons, 0, DSX) ; ii < iiMax ; ii++, up++)
    up->nn = ii ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignCleanExons\n") ;
      showDs (sc, ds->exons) ;
      /* showDs (sc, ds->introns) ; */
      fprintf (stderr, "mrnaDesignCleanExons done\n") ;
    }
  ac_free (h) ;
   return ;
} /* mrnaDesignCleanExons */

/**********************************************************************************/
/* extend the path Down via the best supported route */
static int mrnaDesignExtendDownRaw (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxScore)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;
  
  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;

  for (vp = up, ii = nn0 ; ii >= 0 ; ii--, vp--)
    {
      if (vp->path == path)
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
      int score2 = ((vp->type & gI) ? vp->score/3 :  vp->score) ;
      if (score < score2) { score = score2 ; nn = vp->nn ; wp = vp ; }
    }

  if (! score)
    return length ;

  if (nStop && 
      (wp->type & (gS | gS0 | gReal5p))
      )
    return length ;
  
  if (! (wp->type & (gA | gReal3p)))
    {
      if (up->type & gReal3p)
	{
	  if (up->score > 10 * score)
	    return length ;
	}
      if (up->score > 100 * score)
	return length ;

      if (100 * score < bestScore)
	return length ;
    }

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->path = path ;
  
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
  length += mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore) ;

  return length ;
} /* mrnaDesignExtendDownRaw */

/***********/

static int mrnaDesignExtendDownLoopCDS (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxScore)
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
	if (! (vp->type &  (gCompleteCDS)))
	  continue ;
	if (score < vp->cover) { score = vp->cover ; nn = vp->nn ; wp = vp ; }
      }
  if (! score)
    return   mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn0, nIntron, nStart, nStop, bestScore, maxScore) ;

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  wp->path = path ;
  if (score > bestScore)
    bestScore = score ;
  length += mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore) ;

  return length ;
} /* mrnaDesignExtendDownLoopCDS */

/***********/

static int mrnaDesignExtendDown (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn, int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxScore)
{
  if (useCDS)
    return mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore) ;

  return mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxScore) ;
} /* mrnaDesignExtendDown */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static int mrnaDesignExtendUp (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0
			       , int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxScore)
{ 
  int length = 0 ;
  int a1, ii, nn = 0, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  if (up->nn != nn0) messcrash ("nn != nn0 in  mrnaDesignExtendUp") ;
  a1 = up->a1 - 1 ;
  if (up->type & gS)
    return 0 ;
  if (useCDS)
    for (nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
      {
	if (vp->a2 != a1) continue ;
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
	  if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
	}
    }

  if (! score)
    return length ;

  if (nStart &&
      (wp->type & (gA | gReal3p))
      )
    return length ;
    
  if (! (wp->type & (gS | gS0 | gReal5p)))
    {
      if (up->type & gReal5p)
	{
	  if (up->score > 10 * score)
	    return length ;
	}
      if (up->score > 100 * score)
	return length ;

      if (100 * score < bestScore)
	return length ;
    }

  /* find local exon score */
  if (nStart && ! (wp->type & (gS | gReal5p | gCompleteCDS)))
    if (100 * score < bestScore)
      return length ;
    
  if (! useCDS)
    {
      if (nStart && (wp->type & ( gA | gReal3p) && ! (wp->type & gCompleteCDS)))
	return 0 ;
      if (0 && nIntron && (wp->type & (gA | gReal3p)))
	return 0 ;
    }
  
  length = wp->a2 - wp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->path = path ;
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
  length += mrnaDesignExtendUp (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, useCDS, bestScore, maxScore) ;

  return length ;
} /* mrnaDesignExtendUp */

/**********************************************************************************/

static int mrnaDesignExport (S2M *s2m, SC *sc, DS *ds, Array segs, KEYSET ks, int path, Array smrnas)
{
  int a1 = 0, ii,jj,  nn, maxScore = 0 ;
  DSX *up, *up2 ;
  SMRNA *smrna ;
  HIT *vp = 0 ;
  BOOL debug = TRUE ;
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
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      up = arrp (exions, ii, DSX) ;
      if (up->type & (gS | gS0) ||
	  100 * up->score > maxScore)
	break ;
      up->score = 0 ;
    }
  if (ii > 0 && (up->type & gI) && (up[-1].type & gX) && up[-1].score == 0)
    up[-1].score = up->score ;
  for (ii = keySetMax (ks) - 1 ; ii >= 0 ; ii--)
    {
      up = arrp (exions, ii, DSX) ;
      if (up->type & (gA) ||
	  100 * up->score > maxScore)
	break ;
      up->score = 0 ;
    }
    
   if (debug)
     {
       fprintf (stderr, "mrnaDesignExport### path %d exons/introns\n", path) ;
       showDs (sc, exions) ;
       fprintf (stderr, "mrnaDesignExport### exons/introns done\n") ;
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

static BOOL mrnaDesignIsNewPath (Array ksPaths, KEYSET ks0, AC_HANDLE h)
{
  BOOL ok = TRUE ;
  KEYSET ks = keySetHandleCopy (ks0, h) ;
  int iMax = keySetMax (ks) ;
  keySetSort (ks) ;

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
  if (ok)
    array (ksPaths, arrayMax(ksPaths), KEYSET) = ks ;
  return ok ;
} /* mrnaDesignIsNewpath */

/**********************************************************************************/
static int mrnaDesignFindPaths (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  BOOL debug = TRUE ;
  int ii2, path, nIntron, nStart, nStop, useCDS ;
  int eeMax = ssHappyFew (ds->exons) ;
  Array segs, segs2 ;
  DSX *vp, *vp2 ; ;
  KEYSET ks = keySetHandleCreate (ds->h) ;
  Array ksPaths = arrayHandleCreate (20, KEYSET, ds->h) ;

  arraySort (ksPaths, dsScoreOrder) ; /* for computer happiness */

  if (debug)
     {
       fprintf (stderr, "mrnaDesignFindPaths start\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignFindPaths start done\n") ;
     }
  
  /* sort a copy by cover */
  segs = ds->exons ;
  segs2 = arrayHandleCopy (segs, ds->h) ;
  arraySort (segs2, dsCoverOrder) ;

  /* construct all paths */
  for (path = 0, useCDS = 1 ; useCDS >=0 ; useCDS--)
    {
      ii2 = 0 ;
      vp2 = arrp (segs2, 0, DSX) ;
      int maxScore = vp2->cover ;
      for ( ; ii2 < eeMax ; ii2++, vp2++)
	{
	  /* find highest scoring seg not yet incorporated in a path */
	  vp = arrp (segs, vp2->nn, DSX) ;
	  if (vp->path || ! vp->score || (100 * vp->score < maxScore) || (vp->type & (gDF | gFF)))
	    continue ;
	  if (useCDS && ! (vp->type &  gCompleteCDS))
	    continue ;
	  /* do not start on a mini exon */
	  if (vp->a2 < vp->a1 + 10 && (vp->type & gX) && (!(vp->type & (gS | gA))))
	    continue ;
	  
	  path++ ;
	  vp->path = path ;
	  if (debug)
	    {
	      fprintf (stderr, "+++New path %d start on %d %d score=%d maxScore=%d\n", path, vp->a1, vp->a2, vp->score, maxScore) ;
	    }
	  ks = keySetReCreate (ks) ;
	  keySet (ks, 0) = vp->nn ;
	  if (1) /* ! (vp->type & gB) */
	    {
	      int bestScore = ((vp2->type & gI) ? vp2->score/3 :  vp2->score) ;
	      int nnn = (vp->type & gX) ? 3 : 0 ;
	      nIntron = (vp->type & gI) ? 3 : 0 ;
	      nStart = 1 ;
	      nStop = (vp->type & (gA | gReal3p)) ? 1 : 0 ;
	      mrnaDesignExtendDown (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, bestScore, maxScore) ;
	      nStart = (vp->type & gS) ? 1 : 0 ;
	      nStop = 1 ;
	      if (keySetMax (ks) > 1)
		nnn |= 1 ;
	      if (mrnaDesignIsNewPath (ksPaths, ks, ds->h))
		{
		  int nn2 = keySetMax (ks) ;
		  mrnaDesignExtendUp (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, maxScore, maxScore) ;
		  if (keySetMax (ks) > nn2)
		    nnn |= 2 ;
		  if (nnn >= 0)
		    mrnaDesignExport (s2m, sc, ds, segs, ks, path, smrnas) ;
		}
	    }
	  /*
	  else if (mrnaDesignIsNewPath (ksPaths, ks, ds->h))
	    mrnaDesignExport (s2m, sc, ds, segs, ks, path, smrnas) ;
	  */
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
  const char **errors = 0 ;
  AC_DB db = ac_open_db (0, errors) ; /* local database, cannot fail */

  if (! keySetMax (ds->reads)) 
    return ;
  ks = query (ds->reads, "(IS XY_* || IS XW_*) && (ct_ac || other || small_deletion) && (ct_ac > 1 || other > 1 || (ct_ac && other)) && ! composite > 50)") ;
  if (keySetMax (ks))
    {
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "Sequence %s\n-D Is_read\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), errors, 0, ds->h) ;
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
  mrnaDesignGetGraphElements (&ds, s2m, sc, gmrna, smrnas) ;
  mrnaDesignCleanExons (s2m, sc, &ds) ;
  smrnas = arrayReCreate (smrnas, 12, SMRNA) ;
  mrnaDesignFindPaths (s2m, sc, &ds, smrnas) ;
  if (0) mrnaDesignSetCompletenessFlags (s2m, sc, gmrna, smrnas) ;
  mrnaDesignCleanUp (s2m, sc, &ds, smrnas) ;
  
  ac_free (h) ;
  return ok ;
} /*  mrnaDesignUsingCompositeStrategy */

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
