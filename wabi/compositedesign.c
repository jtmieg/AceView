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
		if ((gX | gI) & up->type) /* exon */
		  {
		    if (gX & up->type) /* exon */
		      dsx = arrayp (exons, jje++, DSX) ;
		    if (gI & up->type) /* exon */
		      dsx = arrayp (introns, jji++, DSX) ;
		    dsx->a1 = up->a1 + smrna->a1 - 1 ;
		    dsx->a2 = up->a2 + smrna->a1 - 1 ;
		    dsx->type = up->type & (gX | gI) ;
		  }
	      }
	  }  
    }
  if (0)  /* stupid idea, i resuse all est from all genes in the current gene */
    {
      BOOL isDown =  (sc->a1 < sc->a2) ? TRUE : FALSE ;
      /* extract all exon fragments and all introns */
      for (j = 0, up = arrp(s2m->plainHits, 0, HIT) ; j < arrayMax(s2m->plainHits) ; up++, j++)
	{
	  if ((gX | gI) & up->type) /* exon */
	    {
	      if (gX & up->type) /* exon */
		dsx = arrayp (exons, jje++, DSX) ;
	      if (gI & up->type) /* exon */
		dsx = arrayp (introns, jji++, DSX) ;
	      if (isDown)
		{
		  dsx->a1 = up->a1 - sc->a1 + 1 ;
		  dsx->a2 = up->a2 - sc->a1 + 1 ;
		}
	      else
		{
		  dsx->a1 = sc->a1 - up->a2 + 1 ;
		  dsx->a2 = sc->a1 - up->a1 + 1 ;
		}
	      dsx->type = up->type & (gX | gI) ;
	    }
	}
    }  

   if (debug)
     {
       fprintf (stderr, "... getElements Z\n") ;
       showDs (sc, exons) ;
     }
   arraySort (exons, dsA1Order) ;
   arrayCompress (exons) ;
   if (debug)
     {
       fprintf (stderr, "... getElements A\n") ;
       showDs (sc, exons) ;
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
       fprintf (stderr, "... getElements B\n") ;
       showDs (sc, ds->exons) ;
     }

   arraySort (introns, dsA1Order) ;
   arrayCompress (introns) ;
   if (debug)
     showDs (sc, introns) ; 

   ac_free (h) ;
   return ;
} /* mrnaDesignGetElements */

/**********************************************************************************/
/* get the counts in introns included in XI XY XW */	
static void mrnaDesignGetIntronSupport (KEY intron, DSX *vp)
{
  int ir ;
  OBJ Intron = bsCreate (intron) ;
  static Array units = 0 ;
  BSunit *uu ;

  units = arrayReCreate (units, 12, BSunit) ;
  if (Intron)
    {
      if (bsFindTag (Intron, _Other) || bsFindTag (Intron, _ct_ac))
	vp->type |= gSuspect ;
      else
	vp->type &= ~gSuspect ;
      if (bsGetArray (Intron, str2tag ("Group_U"), units, 6))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 6)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[2].f > vp->cover) /* seqs */
		vp->cover = uu[2].f ;
	      if (uu[4].f > vp->cover) /* tags */
		vp->cover = uu[4].f ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("Validated_U"), units, 2))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[1].i > vp->cover)
		vp->cover = uu[1].i ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("RNA_seq"), units, 1))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[0].i > vp->cover)
		vp->cover = uu[0].i ;
	    }
	}
      if (bsGetArray (Intron, str2tag ("Magic_any_any"), units, 1))
	{
	  for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
	    {
	      uu = arrp (units, ir, BSunit) ;
	      if (uu[0].i > vp->cover)
		vp->cover = uu[0].i ;
	    }
	}
      bsDestroy (Intron) ;
    }

  return ;
} /* mrnaDesignGetIntronSupport */
		    
/**********************************************************************************/
/* get the XI or XD Composite supporting the introns and their counts */
static void  mrnaDesignGetSupport  (DS *ds, S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  DSX *vp, *ssp, *ssp2 ;
  int ii, jj, ir, iss, iss2, a1, a2, b1, b2, ssMax, u1 ;
  BOOL isDown ;
  BOOL debug = FALSE ;
  BOOL hasData = FALSE ;
  Array units = 0 ;
  Array ss = 0, ss2 = 0 ;
  BSunit *uu ;
  OBJ Est = 0 ;
  KEYSET xks, xws = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  ss = arrayHandleCreate (256, DSX, h) ; iss = 0 ;
  xks = keySetHandleCreate (h) ;
  xws = keySetHandleCreate (h) ;

  /* grab the introns in the plainHits EST alignments */
  isDown = (sc->a1 < sc->a2) ? TRUE : FALSE ;
  units = arrayHandleCreate (12, BSunit, h) ;
  for (ii = 0, up = arrp (s2m->plainHits, 0, HIT) ; ii < arrayMax (s2m->plainHits) ; ii++, up++)
    {
      if (! keyFindTag (up->est, _Composite))
	continue ;
      a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
      a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
      /* identify the introns in the intron table */
      if (arrayMax (ds->introns) &&  (up->type & gI))
	{
	  if (up->x1 < 8 || up->x2 > up->clipEnd - 8)
	    continue ;
	  for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < arrayMax (ds->introns) ; jj++, vp++)
	    {
	      if (vp->a1 != a1 || vp->a2 != a2)
		continue ;
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 5))
		    {
		      keySet (ds->reads, keySetMax (ds->reads)) = up->est ;
		      if (! strncmp (name(up->est), "XD_", 3))
			{
			  if (arrayMax (units) >= 5)
			    {
			      uu = arrp (units, 0, BSunit) ;
			      if (vp->donor < uu[1].i) vp->donor = uu[1].i ; 
			      if (vp->acceptor < uu[3].i) vp->acceptor = uu[3].i ; 
			    }
			}
		      else
			{
			  if (arrayMax (units) >= 1)
			    {
			      int k ;
			      uu = arrp (units, 0, BSunit) ;
			      k =  uu[0].i ; 
			      if (bsFindTag (Est, _Other))
				k /= 5 ;
			      if (vp->cover < k)
			      vp->cover  = k ;
			    }
			}
		    }
		  if (bsGetArray (Est, _Intron, units, 1) &&
		      arrayMax (units) == 1
		      )
		    {
		      KEY intron = 0 ;
		      if (bsGetKey (Est, _Intron, &intron))
			mrnaDesignGetIntronSupport (intron, vp) ;
		    }
		  bsDestroy (Est) ;
		}
	      if ((up->type & gSuspect) && ! (vp->type & gJ))
		{ 
		  vp->type |= gJ ; /* divide only once */
		  vp->cover /= 5 ;
		  vp->donor /= 5 ;
		  vp->acceptor /= 5 ;
		  vp->score /= 5 ;
		}
	    }
	}
      /* register the support of the exon support */
      if (up->type & gX)
	{
	  if ( ! strncmp (name(up->est), "XC_", 3) || ! strncmp (name(up->est), "XF_", 3))  /* coverons */
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 + 10 ; ssp->a2 = a2 - 10 ; ssp->cover = uu[0].i ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}	  
	      continue ;
	    }      
	  
	  if ( ! strncmp (name(up->est), "XG_", 3) || ! strncmp (name(up->est), "XH_", 3) || ! strncmp (name(up->est), "XE_", 3))
	    {
	      int iss0 = iss ;
	      if ((Est = bsCreate (up->est)))
		{
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 3) && arrayMax (units) >= 3)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 3)
		      {
			uu = arrp (units, ir, BSunit) ;
			if (ir) 
			  { 
			    /* oldb1 = b1 ;  */
			    /*  oldc = uu[2].i > uu[-1].i ? uu[-1].i : uu[2].i ; */
			  } 
			b1 = a1 + uu[0].i - 1 ; 	
			b2 = a1 + uu[1].i - 1 ; 
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
			ssp->clipable = 1 ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}
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
			  if (ssp->cover < ssq->cover)
			    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
			  if (ssp->cover > ssq->cover)
			    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
			}
		      if (ssp->a1 > ssp->a2 && ssp->a2 == ssq->a1 + 1)
			{
			  if (ssp->cover < ssq->cover)
			    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
			  if (ssp->cover > ssq->cover)
			    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
			}
		    }
		  ssp = arrayp (ss, iss0, DSX) ;
		  ssp->a1 += (ssp->a1 < ssp->a2 ? +10 : -10) ;
		  ssp = arrayp (ss, iss - 1, DSX) ;
		  ssp->a2 += (ssp->a1 < ssp->a2 ? -10 : +10) ;
		  if (ssp->a1 > ssp->a2) ssp->cover = 0 ;
		  if (debug)
		    fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		}
	      continue ;
	    }	  
      
	  if (! strncmp (name(up->est), "XI_", 3) ||  ! strncmp (name(up->est), "XK_", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  BOOL j ;
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 1)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			/* only count as supported the region close to the intron */
			j = keySetInsert (xks, up->est) ^ up->reverse ;
			if (j) /* first exonic apparition */
			  ssp->a1 = a2 - 30 ;
			else
			  ssp->a2 = a1 + 30 ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }	  
      
	  if (! strncmp (name(up->est), "XW_", 3) || ! strncmp (name(up->est), "XY_", 3))
	    {  /* double introns incoded in the name */
	      if (keySetInsert (xks, up->est) && /* first exonic apparition */
		  (Est = bsCreate (up->est)))
		{
		  int cover = 0 ;
		  const char *ccp = strstr (name(up->est), "__") ;
		  if (! bsFindTag (Est, _Other) &&
		      ccp && bsGetData (Est, _Composite, _Int, &cover) && cover > 0)
		    {
		      int k, u1 = 0, u2 = 0, v1 = 0, v2 = 0 ;
		      char c = 0 ;
		      hasData = TRUE ;
		      k = sscanf (ccp + 2, "%d_%d_%d_%d%c", &u1, &u2, &v1, &v2, &c) ;
		      if (k == 4) /* found 6 numbers, zero terminated */
			{
			  int isDown1 = isDown ? 1 : -1 ;
			  int isDown2 = (u1 < u2 ? 1 : -1) ;
			  int da1 = sc->a1 ;

			  if (! ds->c1)
			    {
			      OBJ Cosmid = bsCreate (s2m->cosmid) ;
			      int c1 = 0 ;
			      KEY map ;

			      if (Cosmid)
				{
				  if (bsGetKey (Cosmid, _IntMap, &map) &&
				      bsGetData (Cosmid, _bsRight, _Int, &c1)
				      )
				    ds->c1 = c1 ;
				  bsDestroy (Cosmid) ;
				}
			    }
			  da1 = sc->a1 + ds->c1 - 1 ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* first exon */
			  ssp->a1 = isDown1 * (u1 - 30 * isDown2 - da1) + 1  ; 
			  if (0) ssp->a1 = isDown1 * (up->a1 - sc->a1) + 1 ;
			  ssp->a2 = isDown1 * (u1 - isDown2   - da1 ) + 1 ;
			  ssp->cover = cover ;
			  if (debug) 
			      fprintf(stderr, "A::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* second exon */
			  ssp->a1 =  isDown1 * (u2 + isDown2  - da1) + 1  ; 
			  ssp->a2 =  isDown1 * (v1 -isDown2  - da1) + 1  ; 
			  ssp->cover = cover ;
			  if (debug)
			      fprintf(stderr, "B::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;

			  ssp = arrayp (ss, iss++, DSX) ;  /* third exon */
			  ssp->a1 =  isDown1 * (v2 + isDown2  - da1) + 1   ; 
			  ssp->a2 =  isDown1 * (v2 + 30 * isDown2  - da1) + 1  ; 
			  if (0) ssp->a2 = isDown1 * (up->a2 - sc->a1) + 1 ;
			  ssp->cover = cover ;
			  if (debug)
			      fprintf(stderr, "C::%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
			}
		    }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }	  

	  if (! strncmp (name(up->est), "XD_", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  BOOL j ;
		  hasData = TRUE ;
		  if (bsGetArray (Est, _Composite, units, 5) && arrayMax (units) >= 5)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 5)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			j = keySetInsert (xks, up->est) ^ up->reverse ;
			if (j) /* first exonic apparition */
			  { ssp->a1 = a2 - 30 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ; }
			else
			  {  
			    j = keySetInsert (xws, up->est) ;
			    if (j)  /* second exonic apparition = supported central exon */
			      { ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[4].i ; }
			    else
			      { ssp->a1 = a1 ; ssp->a2 = a1 + 30 ; ssp->cover = uu[4].i ; }
			  }
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }	    
		  bsDestroy (Est) ; 
		}
	      continue ;
	    }	  
	  
	  if (! strncmp (name(up->est), "XA_", 3)) 
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 2) && arrayMax (units) >= 2)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			ssp->type |= (gA | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d\n", name(up->est), a1, a2, uu[0].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  if (! strncmp (name(up->est), "Xcds_", 5)) 
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			ssp->cover = 12 ; /* a nominal value, uu[0].i is the length of the CDS, not very useful */
			ssp->type |= (gCompleteCDS | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d\n", name(up->est), a1, a2, uu[0].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  if (! strncmp (name(up->est), "Xends_", 6))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			 if (! strncmp (name(up->est), "Xends_ELF_", 9) && isDown) 
			   ssp->type |= (gS | gX) ;
			 if (! strncmp (name(up->est), "Xends_ERF_", 9) && isDown) 
			   ssp->type |= (gReal3p | gX) ;
			 if (! strncmp (name(up->est), "Xends_ERR_", 9) && ! isDown) 
			   ssp->type |= (gS | gX) ;
			 if (! strncmp (name(up->est), "Xends_ELR_", 9) && ! isDown) 
			   ssp->type |= (gReal3p | gX) ;
			if (debug)
			  fprintf(stderr, "....%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  
	  if (! strncmp (name(up->est), "XSL", 3))
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 2) && arrayMax (units) >= 2)
		    for (ir = 0 ; ir < arrayMax (units) ; ir += 2)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			ssp->type |= (gS | gX) ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d:%d \n", name(up->est), a1, a2, uu[0].i, uu[4].i) ;
		      }	    
		  bsDestroy (Est) ;
		}
	      continue ;
	    }
	  
	  if (1)
	    {
	      if ((Est = bsCreate (up->est)))
		{
		  if (bsGetArray (Est, _Composite, units, 1) && arrayMax (units) >= 1)
		    for (ir = 0 ; ir < arrayMax (units) ; ir ++)
		      {
			uu = arrp (units, ir, BSunit) ;
			ssp = arrayp (ss, iss++, DSX) ;
			ssp->a1 = a1 ; ssp->a2 = a2 ; ssp->cover = uu[0].i ;
			if (debug)
			  fprintf(stderr, "%s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
		      }
		  bsDestroy (Est) ;
		}	  
	      continue ;
	    }
	}
    }
  if (! hasData)
    arrayMax (ss) = 0 ;
  if (! arrayMax (ss))
    goto done ;
  /* at each position choose the maximal value */
  arraySort (ss, dsA1Order) ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetSupport phase 1 done\n") ;
      showDs (sc, ss) ;
    }
  if (arrayMax (ds->introns))
    {
      float intronBonus = 3.0 ;
      for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < arrayMax (ds->introns) ; jj++, vp++)
	{
	  if (! (vp->type & gSuspect))
	    vp->cover *= intronBonus ;
	  else
	    vp->cover /= 2 ;
	}
    }
   
  /*  register all boundaries */
  arrayDestroy (s2m->compositeDesignCovering) ;
  s2m->compositeDesignCovering = ds->covering = ss2 = arrayHandleCreate (2 * arrayMax (ss), DSX, s2m->h) ;
  for (jj = iss = 0, ssp = arrp (ss, 0, DSX) ; iss < arrayMax (ss) ; iss++, ssp++)
    {
      if(ssp->a1 > ssp->a2) { ssp->cover = 0 ; continue ; }
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a1 ; ssp2->type = ssp->type & gDebut ; 
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a2 ; ssp2->type = ssp->type & gFin ; 
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a2 + 1 ;
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a1 - 1 ;
    }
  /* some coords appear in the exons but not in ss */
  for (iss = 0, ssp = arrp (ds->exons, 0, DSX) ; iss < arrayMax (ds->exons) ; iss++, ssp++)
    {
      if(ssp->a1 > ssp->a2) { ssp->cover = 0 ; continue ; }
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a1 ; ssp2->type = ssp->type & gDebut ; 
      ssp2 = arrayp (ss2, jj++, DSX) ;
      ssp2->a1 = ssp2->a2 = ssp->a2 ; ssp2->type = ssp->type & gFin ; 
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a2 + 1 ;
      ssp2 = arrayp (ss2, jj++, DSX) ; /* this point is needed to go back to zero at the end of the intervals */
      ssp2->a1 = ssp2->a2 = ssp->a1 - 1 ;
    }
  arraySort (ss2, dsA1Order) ; 
  for (iss = 0, ssp = arrp (ss2, 0, DSX) ; iss < arrayMax (ss2) ; iss++, ssp++)
    { /* so that compress will eliminate coords duplicates created with different flags */
      for (iss2 = iss+1, ssp2 = ssp + 1 ;  iss2 < arrayMax (ss2) && ssp2->a1 == ssp->a1 && ssp2->a2 == ssp->a2; iss2++, ssp2++)
	{ ssp2->type |= ssp->type ; ssp->type = ssp2->type ; }
    }
  arrayCompress (ss2) ;
  /* register the max values in each zone */ 
  for (jj = iss = iss2 = 0, ssp = arrp (ss, 0, DSX), ssp2 = arrp (ss2, 0, DSX) ; iss < arrayMax (ss) ; iss++, ssp++)
    {
      /* position ssp2 on ssp->a1 */
      while (iss2 > 0 && ssp2->a1 > ssp->a1) { iss2-- ; ssp2-- ; }
      while (iss2 < arrayMax (ss2) && ssp2->a1 < ssp->a1) { iss2++ ; ssp2++ ; }
      /* possibly increment the coverage */
      while (iss2 < arrayMax (ss2) && ssp2->a1 <= ssp->a2) 
	{ 
	  if (ssp2->cover < ssp->cover) /* 2024_02 */
	    ssp2->cover = ssp->cover ;
	  iss2++ ; ssp2++ ; 
	}
    }
  if (debug)
    {
      fprintf (stderr, "max value in each zone\n") ;
      showDs (sc, ss2) ;
    }

  /* look for weak introns donor/acceptor */
   ii = 0 ; ssp = arrp (ss2, 0, DSX) ; ssMax = arrayMax (ss2) ;
   if (arrayMax(ds->introns))
     for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; ii < ssMax && jj < arrayMax (ds->introns) ; jj++, vp++)
       {
	 DSX* vp2 ;
	 int jj2 ;
	 long int cover = vp->cover ;

         cover *= 100 ;

	 a1 = vp->a1 ;
	 a2 = vp->a2 ;
	 
	 /* compare the cover to the donor exon (a1 - 1) and the retained intron (a1) */
	 while (ii > 0 && ssp->a1 > a1-5) { ii-- ; ssp-- ; }
	 for ( ; ssp->a1 <= a1+5 && ii < ssMax ; ssp++, ii++)
	   if (ssp->a1 >= a1 - 5 && ssp->a1 <= a1 + 5)
	     if (ssp->cover > cover)
	       { vp->donor = -1 ; break ; }
	 /* compare to other introns starting at same site */
	 for (jj2 = jj - 1, vp2 = vp - 1 ; jj2 >= 0  && vp2->a1 == a1 ; jj2--, vp2--)
	    if (vp2->cover > cover)
	       { vp->donor = -1 ; break ; }
	 for (jj2 = jj + 1, vp2 = vp + 1 ; jj2 <  arrayMax (ds->introns)  && vp2->a1 == a1 ; jj2++, vp2++)
	    if (vp2->cover > cover)
	       { vp->donor = -1 ; break ; }
	 /* compare the cover to the donor exon (a2 + 1) and the retained intron (a2) */
	 while (ii > 0 && ssp->a1 > a2 - 5) { ii-- ; ssp-- ; }
	 for ( ; ssp->a1 <= a2 + 5 && ii < ssMax ; ssp++, ii++)
	   if (ssp->a1 >= a2 - 5  && ssp->a1 <= a2 + 5)
	     if (ssp->cover > cover)
	       { vp->acceptor = -1 ; break ; }
	 for (jj2 = jj + 1, vp2 = vp + 1 ; jj2 <  arrayMax (ds->introns)  && vp2->a1 < a2 ; jj2++, vp2++)
	    if (vp2->a2 == a2 && vp2->cover > cover)
	       { vp->acceptor = -1 ; break ; }
	 if (vp->donor == -1 && vp->acceptor == -1)
	   vp->cover = vp->score = 0 ; 
	 if (ii) { ii-- ; ssp-- ; }
       }
  
   
   /* compute the average support of each exon fragment */
   ii = 0 ; ssp = arrp (ss2, 0, DSX) ; ssMax = arrayMax (ss2) ;
   if (arrayMax (ds->exons) )
     for (jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < arrayMax (ds->exons) ; jj++, vp++)
       {
	 long int nn, nnn ;
	 nnn = nn = 0 ;
	 a1 = vp->a1 ;
	 a2 = vp->a2 ;
	 while (ii > 0 && ssp->a1 > a1) { ii-- ; ssp-- ; }
	 u1 = a1 ; if (u1 < ssp->a1) u1 = ssp->a1 ; nn = ssp->cover ; 
	 
	 for ( ; ii < ssMax && ssp->a1 < a1 ; ii++, ssp++) 
	   {
	     nn = ssp->cover ; 
	   }
	 
	 for ( ; ii < ssMax && ssp->a1 < a2 ; ii++, ssp++)
	   {
	     nnn += nn * (ssp->a1 - u1) ;
	     if (ssp->a1 >= u1) u1 = ssp->a1 ; 
	     nn = ssp->cover ;
	   }
	 if (0 && ii < ssMax && a2 == ssp->a1)
	   vp->type |= (ssp->type & gFin) ; 
	 if (a2 >= u1) 
	   nnn += nn * (a2 - u1 + 1) ;
	 
	 vp->cover = nnn / (a2 - a1 + 1) ;
	 if (nnn >= 0 && vp->cover < 1) vp->cover = 1 ;
       }
   
   if (debug)
     {
       fprintf (stderr, "mrnaDesignGetSupport phase 2 done\n") ;
       showDs (sc, ss) ;
     }
   
   if (1) /* Xcds is used to mask Xends */
     for (jj = 0, vp = arrp (ss, 0, DSX) ; jj < arrayMax (ss) ; jj++, vp++)
       if (vp->type & gCompleteCDS)
	 for (iss = 0, ssp = arrp (ss2, 0, DSX) ; iss < arrayMax (ss2) ; iss++, ssp++)
	   {
	     if (
		 ssp->a1 >= vp->a1 &&
		 ssp->a1 <= vp->a2
		 )
	       {
		 ssp->type |=  gCompleteCDS ;
		 ssp->type &= ~gReal3p ;  /* Xcds is used to mask Xends */
		 ssp->type &= ~gA ;  /* Xcds is used to mask Xends */
	       }
	   }
   
   if (debug)
     {
       fprintf (stderr, "mrnaDesignGetSupport phase 3 done\n") ;
       showDs (sc, ss2) ;
	}
   
   if (1) /* reestablish the polyA flags */
     for (jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < arrayMax (ds->exons) ; jj++, vp++)
       {
	 for (iss = 0, ssp = arrp (ss2, 0, DSX) ; iss < arrayMax (ss2) ; iss++, ssp++)
	   {
	     if ((ssp->type & gCompleteCDS) &&
		 ssp->a2 <= vp->a2 &&
		 ssp->a2 >= vp->a1
		 )
	       vp->type |= gCompleteCDS ;
	     if ((ssp->type & gS) &&
		 ssp->a2 < vp->a2 &&
		 ssp->a2 >= vp->a1
		 )
	       vp->type |= gS ;
	     if ((ssp->type & (gA | gReal3p)) &&
		 ssp->a1 > vp->a1 &&
		 ssp->a1 <= vp->a2
		 )
	       {
		 vp->type |= (ssp->type & (gA | gCompleteCDS | gReal3p)) ;
		 if (vp->type & gCompleteCDS &  gReal3p) vp->type &= ~gReal3p ;
		 if (vp->type & gCompleteCDS &  gA) vp->type &= ~gA ;
	       }
	   }
       }
   if (debug)
     {
       fprintf (stderr, "mrnaDesignGetSupport  done\n") ;
       showDs (sc, ds->exons) ;
     }
 done:
   ac_free (h) ;
   return ;
} /* mrnaDesignGetSupport */

/**********************************************************************************/
/* copy the exon up0 into ds->exons, possibly splitting it */
static void mrnaDesignSplitOneExon (DS *ds, SC *sc, DSX *up0)
{
  DSX *up, *xp ;
  
  if ( 
      (up0->type & gX) &&  
      (	   
       (up0->type & (gReal5p | gS | gS0)) ||
       (up0->type & (gReal3p | gA))
	   )
       )
    {
      int ii ;
      int iMax = arrayMax (ds->covering) ;
      int a1 = up0->a1 ;
      int a2 = up0->a2 ;
      int cover = 0 ;
      unsigned int type = 0 ;
      
      for (ii = 0, xp = arrp (ds->covering, 0, DSX) ; ii < iMax ; xp++, ii++)
	{
	  if (xp->a1 >= a1 && xp->a1 <= a2)
	    {
	      if  (ii && 
		   (
		    ( 5 * xp->cover < (xp-1)->cover || 
		      xp->cover > 5 * (xp-1)->cover) ||
		    ( (type | xp->type) & (gS | gS0 | gReal3p | gA)) ||
		    ( (xp->type & gDebut) && (type & gFin)) ||
		    ( (type & gDebut) && (xp->type & gFin))
		    )
		   )
		{
		  if (cover)
		    {
		      up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
		      up->a1 = a1 ;
		      up->a2 = xp->a1 - 1 ;
		      up->type = type | gX ;
		      if (up->a2 >= up->a1)
			up->cover = cover/(up->a2 - up->a1 + 1) ;
		    }
		  a1 = xp->a1 ;
		  cover = 0 ;
		  type = 0 ;
		}
	      type |= xp->type ;
	      cover += xp->cover * ((xp+1)->a1 - xp->a1) ;
	    }
	}
      up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
      up->a1 = a1 ;
      up->a2 = up0->a2 ;
      if (up->a2 + 1 > up->a1)
	{
	  up->cover = cover/(up->a2 - up->a1 + 1) ;
	  up->type = type | gX ;
	}
      else
	up->a2 = up->a1 ;
    }
  else
    {
      up = arrayp (ds->exons, arrayMax (ds->exons), DSX) ;
      *up = *up0 ;
    }
  
  return ;
} /* mrnaDesignSplitOneExon */

/**********************************************************************************/

static void mrnaDesignSplitExons (DS *ds, SC *sc)
{
  BOOL debug =  FALSE ;

  if (debug)
    {
      fprintf (stderr, "mrnaDesignSpliExon starts\n") ;
      showDs (sc, ds->covering) ;
      fprintf (stderr, "   .. original exons\n") ;
      showDs (sc, ds->exons) ;
    }

  if (arrayMax (ds->exons))
    {
      int ii ;
      DSX *up ;
      Array aa = arrayCopy (ds->exons) ;

      arrayMax (ds->exons) = 0 ;
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	mrnaDesignSplitOneExon (ds, sc, up) ; 
      arrayDestroy (aa) ;
    }



  if (arrayMax (ds->exons) &&
      ds->covering &&
      arrayMax (ds->covering)
      ) /* reestablish the polyA flags */
    {
      int iss, jj ;
      DSX *vp, *ssp ;
      
      for (jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < arrayMax (ds->exons) ; jj++, vp++)
	{
	  for (iss = 0, ssp = arrp (ds->covering, 0, DSX) ; iss < arrayMax (ds->covering) ; iss++, ssp++)
	    {
	      if ((ssp->type & gCompleteCDS) &&
		  ssp->a2 <= vp->a2 &&
		  ssp->a2 >= vp->a1
		  )
		vp->type |= gCompleteCDS ;
	      if ((ssp->type & gS) &&
		  ssp->a1 < vp->a2 &&
		  ssp->a1 >= vp->a1
		  )
		vp->type |= gS ;
	      if ((ssp->type & (gA | gReal3p)) &&
		  ssp->a2 > vp->a1 &&
		  ssp->a2 <= vp->a2
		  )
		{
		  vp->type |= (ssp->type & (gA | gCompleteCDS | gReal3p)) ;
		  if (vp->type & gCompleteCDS &  gReal3p) vp->type &= ~gReal3p ;
		  if (vp->type & gCompleteCDS &  gA) vp->type &= ~gA ;
		}
	    }
	}
    }
  
  if (debug)
    {
      fprintf (stderr, "   .. final exons\n") ;
      showDs (sc, ds->exons) ;
    }
  
 return ;
} /* mrnaDesignSplitExons  */

/**********************************************************************************/
/* extend the path Down via the best supported route */
static int mrnaDesignExtendDownRaw (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int *nIntronp, int *nStartp, int *nStopp)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;

  for (nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
    {
      if (vp->a1 != a2) continue ;
      if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }

  if (! score)
    return 0 ;

  if (*nStopp &&  !(wp->type & gCompleteCDS) && 
      (wp->type & (gS | gS0 | gReal5p))
      )
    return length ;
  
  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  if (1 || ! wp->path) wp->path = path ;
  
  if (wp->type & gI)
    (*nIntronp)++ ;
  if (wp->type & gS)
    (*nStartp)++ ;
  if (wp->type & (gA | gReal3p))
    (*nStopp)++ ;
  
  if (0) fprintf (stderr, "Path=%d  %d:%d >> %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
  length += mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntronp, nStartp, nStopp) ;

  return length ;
} /* mrnaDesignExtendDownRaw */

/***********/
static int mrnaDesignExtendDownLoopCDS (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int *nIntronp, int *nStartp, int *nStopp)
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
	if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
      }
  if (! score)
    return   mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn0, nIntronp, nStartp, nStopp) ;

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  if (1 || ! wp->path) wp->path = path ;
  
  length += mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntronp, nStartp, nStopp) ;

  return length ;
} /* mrnaDesignExtendDownLoopCDS */

/***********/

static int mrnaDesignExtendDown (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn, int *nIntronp, int *nStartp, int *nStopp, int useCDS)
{
  if (useCDS)
    return mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntronp, nStartp, nStopp) ;

  return mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntronp, nStartp, nStopp) ;
} /* mrnaDesignExtendDown */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static int mrnaDesignExtendUp (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int *nIntronp, int *nStartp, int *nStopp, int useCDS)
{ 
  int length = 0 ;
  int a1, ii, nn, score = 0 ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;

  up = arrp (segs, nn0, DSX) ;
  a1 = up->a1 - 1 ;
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
  
  if (score)
    {
      /* kill alternative branches  with too low relative score */
      if (0)
	for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
	  {  
	    if (vp->a2 != a1) continue ;
	    if (score > 200 * vp->score) { vp->score = 0 ; }
	  }

      if (! useCDS)
	{
	  if (*nStartp && (wp->type & ( gA | gReal3p) && ! (wp->type & gCompleteCDS)))
	    return 0 ;
	  if (0 && *nIntronp && (wp->type & (gA | gReal3p)))
	    return 0 ;
	}

      length = vp->a2 - vp->a1 + 1 ;
      keySet (ks, arrayMax(ks)) = nn ;

      if (1 || ! wp->path) wp->path = path ;
      if (! useCDS)
	{
	  if (wp->type & gI)
	    (*nIntronp)++ ;
	  if (wp->type & gS)
	    (*nStartp)++ ;
	  if (wp->type & (gA | gReal3p))
	    (*nStopp)++ ;
	}

      if (0) fprintf (stderr, "Path=%d  %d:%d << %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
      length += mrnaDesignExtendUp (sc, ds, segs, ks, path, nn, nIntronp, nStartp, nStopp, useCDS) ;
    }
  return length ;
} /* mrnaDesignExtendUp */

/**********************************************************************************/

static int mrnaDesignExport (S2M *s2m, DS *ds, Array segs, KEYSET ks, int path, Array smrnas)
{
  int a1 = 0, ii,jj,  nn, nn2, nn3, maxScore = 0 ;
  DSX *up, *up2, *up3 ;
  SMRNA *smrna ;
  HIT *vp = 0 ;
  BOOL debug = FALSE ;

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
    }
  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      nn = keySet (ks, ii) ; 
      up = arrp (segs, nn, DSX) ;
      up2 = 0 ;
      if (ii + 1 < keySetMax (ks))
	{
	  nn2 = keySet (ks, ii + 1) ; 
	  up2 = arrp (segs, nn2, DSX) ;
 	}
      if (nn == 0 || 
	  100 * up->score > maxScore ||
	  (up->type & (gI | gS | gS0)) ||
	  (up2 && (up2->type & gI)) ||
	  (up2 && (up2->score < (up2->type & gS ? 3 : 10) * up->score)) 	  
	  )
	break ;
      up->score = 0 ;
    }
  for (ii = keySetMax (ks) - 1 ; ii >= 0 ; ii--)
    {
      nn = keySet (ks, ii) ; 
      up = arrp (segs, nn, DSX) ;
      up2 = up3 = 0 ;
      if (ii > 0)
	{
	  nn2 = keySet (ks, ii - 1) ; 
	  up2 = arrp (segs, nn2, DSX) ;
 	}
      if (ii > 1)
	{
	  nn3 = keySet (ks, ii - 2) ; 
	  up3 = arrp (segs, nn3, DSX) ;
 	}
     if ( up3 && 
	  (up2->type & (gA | gReal3p)) &&
	  ! (up->type & (gI | gA | gReal3p)) &&
	  up3->score > 100 * up->score
	  )
       up->score = 0 ;	  

      if (nn == arrayMax (segs) - 1 ||
	  100 * up->score > maxScore ||
	  (up->type & (gI | gA | gReal3p)) ||
	  (up2 && (up2->type & gI)) ||
	  (up2 && (up2->score <  (up2->type & (gA | gReal3p) ? 3 : 10) * up->score))
	  )
	break ;
      up->score = 0 ;
    }

  for (ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      nn = keySet (ks, ii) ; 
      up = arrp (segs, nn, DSX) ;
      if (!ii && (up->type & gI))
	continue ;
      if (! up->score)
	continue ;
      if (!vp) a1 = smrna->a1 = up->a1 ; 
      smrna->a2 = up->a2 ; 
      if (!vp || ((up->type & (gX | gI)) != (vp->type & (vp->type & (gX | gI)))))
	{	
	  vp = arrayp (smrna->hits, jj++, HIT) ;
	  vp->a1 = up->a1 - a1 + 1 ;
	}
      vp->a2 = up->a2 - a1 + 1 ;
      vp->type |= up->type ;
      if (debug) fprintf(stderr, "Path %d :: %d :: %d :: %s %d %d\t%d\td%d:a%d:c%d\t%d%s\ttvp %d:%d\n"
	      , path, ii, nn, up->type & gX ? "exon" : "intron", up->a1, up->a2, up->score, up->donor, up->acceptor, up->cover
	      , up->path, path == up->path ? "*****" : ""
	      , vp->a1, vp->a2
	      ) ;
    }
  return path ;
} /* mrnaDesignExport G_t_NC_003281.10_1_4561 */

/**********************************************************************************/

static int mrnaDesignFindPaths (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  int ii, jj, ii2, maxScore, path, a2Max, nIntron, nStart, nStop, useCDS ;
  int eeMax = arrayMax (ds->exons) ;
  int iiMax = arrayMax (ds->introns) ;
  Array segs, segs2 ;
  DSX *up, *vp, *vp2 ; ;
  KEYSET ks = keySetHandleCreate (ds->h) ;
  BOOL debug = FALSE ;

  /* group exons/introns */
  segs = arrayHandleCreate (eeMax + iiMax +1, DSX, ds->h) ;
  if (arrayMax (ds->exons))
    for (ii = jj = 0, vp = arrp (ds->exons, 0, DSX) ; jj < eeMax ; ii++, jj++, vp++)
      {	    
	up =  arrayp (segs,ii, DSX) ;
	*up = *vp ; up->score = vp->cover ; up->nn = ii ;
	}
  if (arrayMax (ds->introns))
    for (jj = 0, vp = arrp (ds->introns, 0, DSX) ; jj < iiMax ; ii++, jj++, vp++)
      {	    
	up =  arrayp (segs,ii, DSX) ;
	*up = *vp ; up->score = vp->cover ; up->nn = ii ;
	if (up->score < vp->donor) up->score = vp->donor ;
	if (up->score < vp->acceptor) up->score = vp->acceptor ;
	up->score <<= 1 ;
	if (0 && ! up->score) up->score = 1 ;
      }
  arraySort (segs, dsA1Order) ;
  eeMax = arrayMax (segs) ; 
  
  if (0)
    showDs (sc, segs) ;
  
      /* kill the non extremal exon extensions which are below 15 bp and do not hook onto an intron */
  if (arrayMax (segs))
    for (ii = a2Max = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
      if (a2Max < vp->a2) a2Max = vp->a2 ;
  
  if (arrayMax (segs))
    for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
      {
	if (! (vp->type & (gDebut | gFin)) && (vp->type & gX) && vp->a2 < vp->a1 + 15)
	  {
	    int ok1 = 0, ok2 = 0, ok3 = 0, ok4 = 0 ;
	    int a1 = vp->a2 + 1, a2 = vp->a1 - 1 ;
	    if (vp->a1 == 1)
	      ok1++ ;
	    else
	      for (jj = ii - 1, vp2 = vp - 1 ; ok1 == 0 && jj >= 0 ; jj--, vp2--)
		{
		  if (vp2->a2 == a2)
		    {
		      if (1 || (vp2->type & gI) ||
			  (vp->score > .8 * vp2->score && vp->score < 1.2 * vp2->score))
			ok1 = 1 ;
		    }
		  if (0 &&  (vp2->type & gI) && vp2->a2 == vp->a2 &&  vp->a2 < vp->a1 + 8)
		    {
		      ok3 = 1 ;
		    }
		}
	    if (vp->a2 == a2Max) 
	      ok2 = 1 ;
	    else
	      for (jj = ii + 1, vp2 = vp + 1 ; ok2 == 0 && jj < eeMax && vp2->a1 <= a1 ; jj++, vp2++)
		{
		  if (vp2->a1 == a1)
		    {
		      if (1 || (vp2->type & gI) ||   /* vp2 intron sticcking to vp exon */
			  (vp->score > .8 * vp2->score && vp->score < 1.2 * vp2->score)) /* vp2 exon with similar score directly sticking to vp exon */
			ok2 = 1 ;
		    }
		  if (0 &&  (vp2->type & gI) && vp2->a1 == vp->a1 &&  vp->a2 < vp->a1 + 8)
		    {
		      ok4 = 1 ;
		    }
		  
		}
	    if (ok1 + ok2 < 2) 
	      vp->score = -1 ;
	    else if (0 && ok3 + ok4) /* disfavor intron leaks below 10bp, the wiggle resolution */
	      vp->score /= 5 ; /* already done, grep for 'disfavor'  */
	  }
      }
  /* if we have   2 introns with same acceptor    XXXXXXXX=====----XXXX and a matching === bit of exon, the ===  should have the support of the short intron */
  if (arrayMax (segs))
    for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
      {
	if (vp->type & gI)
	  {
	    int ok2 = 0 ;
	    for (jj = ii + 1, vp2 = vp + 1 ; ok2 == 0 && jj < eeMax && vp2->a1 < vp->a2 ; jj++, vp2++)
	      if ((vp2->type & gI) && vp2->a2 <= vp->a2)
		{
		  /* search for the === exon */
		  int kk ;
		  DSX *vp3 ;
		  for (kk = ii - 1, vp3 = vp - 1 ; ok2 == 0 && kk >= 0 && vp3->a1 < vp->a2 ; kk--, vp3--)
		    if ((vp3->type & gX) && vp3->a1 == vp->a1 && vp3->a2 == vp2->a1 - 1)
		      {
			if (vp3->score < vp2->score)
			  vp3->score = vp2->score ;
			ok2 = 1 ;
		      }
		}
	  }
	/* same problem single donor 2 acceptors */
	if (vp->type & gI)
	  {
	    int ok2 = 0 ;
	    for (jj = ii + 1, vp2 = vp + 1 ; ok2 == 0 && jj < eeMax && vp2->a1 == vp->a1 ; jj++, vp2++)
	      if ((vp2->type & gI) && vp2->a2 > vp->a2)
		{
		  /* search for the === exon */
		  int kk ;
		  DSX *vp3 ;
		  for (kk = ii + 1, vp3 = vp + 1 ; ok2 == 0 && kk < eeMax && vp3->a1 <= vp->a2 + 1 ; kk++, vp3++)
		    if ((vp3->type & gX) && vp3->a1 == vp->a2 + 1 && vp3->a2 == vp2->a2 - 1)
		      {
			if (vp3->score < vp->score)
			  vp3->score = vp->score ;
			ok2 = 1 ;
		      }
		}
	  }
	
	/* a retained intron which is over 50% of the neighbouring exon
	 * should have a score higher than the corrsponding intron
	 * otherwise if it is below 10% it should be dropped
	 * otherwise it should be flagged and not extended
	 */
	if (vp->type & gX)
	  {
	    int ok2 = 0 ;
	    DSX *wp = 0 ;
	    for (jj = ii + 1, vp2 = vp + 1 ; ok2 < 2 && jj < eeMax && vp2->a1 == vp->a1 ; jj++, vp2++)
	      if ((vp2->type & gI) && vp2->a2 >= vp->a2)
		{ wp = vp2 ; ok2 =  (vp2->a2 == vp->a2 ? 2 : 1) ; }
	    for (jj = ii - 1, vp2 = vp - 1 ; ok2 < 2 && jj >= 0 ; jj--, vp2--)
	      if ((vp2->type & gI) && vp2->a2 == vp->a2)
		{ wp = vp2 ;  ok2 = (vp2->a1 == vp->a1 ? 2 : 1) ; }
	    if (ok2) 
	      {
		int score1 = 0, score2 = 0 ;
		/* search for highest sore of neighbouring exon */
		for (jj = ii + 1, vp2 = vp + 1 ; jj < eeMax && vp2->a1 <= vp->a2 + 1 ; jj++, vp2++)
		  if (vp2->a1 == vp->a2 + 1 && vp2->score > score1)
		    {
		      if (vp2->a1 >= wp->a1 + 1 || (vp2->type & gI))
			score1 = vp2->score ;
		      if ((vp2->a2 <= wp->a2 || vp2->a2 - vp2->a1 < 10) && jj + 1 < eeMax)
			{
			  DSX *vp3 = vp2 + 1 ;
			  if (vp3->a1 == vp2->a2 + 1 && vp3->score > score1 &&
			      vp3->a1 >= wp->a1 + 1)
			    score1 = vp3->score ;
			}
		    }
		for (jj = ii - 1, vp2 = vp - 1 ; jj >= 0 ; jj--, vp2--)
		  if ((vp2->type) && vp2->a2 == vp->a1 - 1 && vp2->score > score2)
		    {
		      if (vp2->a2 <= wp->a1 - 1 || (vp2->type & gI))
			score2 = vp2->score ; 
		      if ((vp2->a1 >= wp->a1 || vp2->a2 - vp2->a1 < 10) && jj > 0)
			{
			  DSX *vp3 = vp2 - 1 ;
			  if (vp3->a2 == vp2->a1 - 1 && vp3->score > score2 &&
			      vp3->a2 <= wp->a1 - 1)
			    score2 = vp3->score ;
			}
		    }
		if (vp->a2 < vp->a1 + 8 && score1 * score2 == 0)  /* disfavor intron leaks below 10bp, the wiggle resolution */
		  vp->score /= 5 ;
		if ((! score1 || 20 * vp->score < score1) && (!score2 || 20 * vp->score < score2))
		  vp->score = -1 ; /* drop very low retained introns */
		else if (
			 ((score1 && 2 * vp->score > score1) || (score2 && 2 * vp->score > score2)) &&
			 vp->score <= wp->score
			 )
		  { vp->score = wp->score + 1 ; vp->type &= (~gB) ; }
		else if (vp->score >= wp->score)
		  { }
		else
		  { vp->type |= gB ; vp->type &=(~ gCompleteCDS) ; }
	      }
	  }
      }
  
  /* keep happy few */
  if (arrayMax (segs))
    {
      for (ii = jj = 0, vp = vp2 = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
	{
	  if (vp->score <= 0)
	    continue ;
	  if (vp2 < vp)
	    *vp2 = *vp ;
	  vp2++ ; jj++ ;
	}
      arrayMax (segs) = eeMax = jj ;
    
      /* flag  gCompleteCDS introns: those linking gCompleteCDS exons */
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
      
      for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
	vp->nn = ii ;
    }

 /* sort a copy by score */
  segs2 = arrayHandleCopy (segs, ds->h) ;
  arraySort (segs2, dsScoreOrder) ;

  if (debug)
    showDs (sc, segs) ;

  for (path = 0, useCDS = 1 ; useCDS >=0 ; useCDS--)
    {
      for (ii2 = 0, vp2 = arrp (segs2, 0, DSX), maxScore = ((vp2->type & gI) ? vp2->score/6 :  vp2->score) ; ii2 < eeMax && (100 * vp2->score > maxScore || vp2->score >= 1) ; ii2++, vp2++)
	{
	  /* find highest scoring seg not yet incorporated in a path */
	  vp = arrp (segs, vp2->nn, DSX) ;
	  if (vp->path || ! vp->score)
	    continue ;
	  if (useCDS && ! (vp->type &  gCompleteCDS))
	    continue ;
	  /* do not start on a retained intron */
	  if (0 && vp->type & gX)
	    {
	      DSX *wp, *wp2 ;
	      int jj2 ;
	      BOOL ok2 = TRUE ;
	      
	      for (jj2 = 0, wp2 = arrp (segs2, 0, DSX) ; ok2 && jj2 < ii2 ; jj2++, wp2++)
		{
		  wp =  arrp (segs, wp2->nn, DSX) ;
		  if ((wp->type & gI) && wp->a1 <= vp->a1 && wp->a2 >= vp->a2 && wp->score > 3 * vp->score && 4*vp->score < maxScore)
		    ok2 = FALSE ;  /* wp is a true  intron better and including  my stupid  vp retained intron (i.e. artefactual exon) */
		}
	      if (! ok2)
		continue ;
	    }
	  
	  path++ ;
	  vp->path = path ;
	  ks = keySetReCreate (ks) ;
	  keySet (ks, 0) = vp2->nn ;
	  if (! (vp->type & gB))
	    {
	      nIntron = (vp->type & gI) ? 1 : 0 ;
	      nStart = 1 ;
	      nStop = (vp->type & (gA | gReal3p)) ? 1 : 0 ;
	      mrnaDesignExtendDown (sc, ds, segs, ks, path, vp->nn, &nIntron, &nStart, &nStop, useCDS) ;
	      nStart = (vp->type & gS) ? 1 : 0 ;
	      nStop = 1 ;
	      mrnaDesignExtendUp (sc, ds, segs, ks, path, vp->nn, &nIntron, &nStart, &nStop, useCDS) ;
	    }
	  mrnaDesignExport (s2m, ds, segs, ks, path, smrnas) ;
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
  ks = query (ds->reads, "IS XY_* && (ct_ac || other || small_deletion) && (ct_ac > 1 || other > 1 || (ct_ac && other) || ! composite > 50)") ;
  if (keySetMax (ks))
    {
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "Sequence %s\n-D Is_read\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), errors, 0, ds->h) ;
    } 
  keySetDestroy (ks) ;
  
#ifdef JUNK
  ks = query (ds->reads, "Intron && ! from_gene  && ! composite > 50 ; >cdna_clone") ;
  keySetDestroy (ks) ;
  ks = query (ds->reads, "! from_gene  ; >cdna_clone") ;
  if (0 && keySetMax (ks))
    { /* this does not work because for now none of the reads is associated to a tg */
      txt = vtxtHandleCreate (ds->h) ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
	vtxtPrintf (txt, "cDNA_clone %s\nIgnore_this_clone\n\n", name(keySet(ks, ii))) ;
      ac_parse (db, vtxtPtr (txt), errors, 0, ds->h) ;
    } 
  keySetDestroy (ks) ;
#endif

  /* shift the coordinates */
  ii = smrnas ? arrayMax (smrnas) : 0 ;
  if (0 && ii)
    {
      int iiMax = ii, a1 = 100000000 ;
      SMRNA *smrna =  arrayp(smrnas,0, SMRNA) ;
      for (ii = 0 ; ii < iiMax ; ii++)
	{
	  if (a1 > smrna[ii].a1)
	    a1 = smrna[ii].a1 ;
	}
      if (a1 > 1)
	{
	  a1-- ;
	  for (ii = 0 ; ii < iiMax ; ii++)
	    {
	      smrna[ii].a1 -= a1 ;
	      smrna[ii].a2 -= a1 ;
	    }
	  if (0) { sc->sh.a1 -= a1 ; sc->sh.a2 -= a1 ;  }
	  if (sc->isUp)
	    { sc->a1 -= a1 ; sc->a2 -= a1 ; }
	  else
	    { sc->a1 += a1 ; sc->a2 += a1 ; }
	}
    }

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
  int iss, issMax = arrayMax (covering) ;
  DSX *ssp ;
  unsigned int topFlags = gS | gS0 | gReal5p ;
  unsigned int bottomFlags = gReal3p | gA ;
  unsigned int cFlags = isTop ? topFlags : bottomFlags ;
  int a1 = up->a1, a2 = up->a2 ;
  if (!(up->type & gX) || !(up->type & cFlags))
    return ;
  up->type &= ~cFlags ;
 
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
      
      /* search the flags of the first exon */
      if (arrayMax (smrna->hits))
	for (j = 0, up = arrp (smrna->hits, 0, HIT) ; j < 1 && j < arrayMax(smrna->hits);  up++, j++)
	  mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, TRUE) ;

       /* search the polyA flags of the last exon */
      if (arrayMax (smrna->hits))
	for (j = arrayMax(smrna->hits) - 1 ; j >= 0 && (up = arrp (smrna->hits, j, HIT)) ; j = -1)
	  mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, FALSE) ;
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
  mrnaDesignGetSupport (&ds, s2m, sc, gmrna, smrnas) ;
  mrnaDesignSplitExons (&ds, sc) ;
  if (debug)
    {
      showDs (sc, ds.exons) ;
      showDs (sc, ds.introns) ; 
    }
  smrnas = arrayReCreate (smrnas, 12, SMRNA) ;
  mrnaDesignFindPaths (s2m, sc, &ds, smrnas) ;
  mrnaDesignSetCompletenessFlags (s2m, sc, gmrna, smrnas) ;
  mrnaDesignCleanUp (s2m, sc, &ds, smrnas) ;
  
  ac_free (h) ;
  return ok ;
} /*  mrnaDesignUsingCompositeStrategy */

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
