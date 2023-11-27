/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * Septembre 2021
 * altintrons.c
 *   Input: INTRON_DB database, with schema as in 
 *       waligner/metadata/wspec.aceview/models.wrm

 * Strategy:
 *   For each intron, the database countains raw support counts under
 *        de_uno ?Run Int
 *   We select pairs of overlapping hence incompatible introns
 *   there are four types; same donor, same accptor, citroen, pinching
 */
#include <acedb.h>
#include <ac.h>


typedef struct gxStruct {
  AC_HANDLE h ; 
  AC_DB db ; 
  const char *outFileName, *inFileName ;
  const char *project ; 
  BOOL gzo, gzi ;
  
  BOOL hasAltIntron ;
  Array runs ;
  Array histos ;
  Array introns ;
  Array pairs ;
  Array donors ;
  Array acceptors ;
  DICT *runDict ;
} GX ;

typedef struct itrStruct {
  KEY intron ;  /* a database KEY */
  KEY gene ;
  KEY map ;
  int a1, a2 ;
  int type ;
  int nSD, nSA;
  Array aa ;
} ITR ;

typedef struct pairStruct {
  int intron1, intron2 ;
  int gene ;
  Array aa ;
} PAIR ;

typedef struct runStruct {
  int run ; /* self, needed after we sort the table */
  KEYSET r2g ; /* list of groups g of which r is a member */
  KEYSET g2r ; /* list of runs in this group */
} RC ;

/*************************************************************************************/
/*************************************************************************************/

static int gxGetRuns (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *runDict ;
  AC_DB db = gx->db ;
  Array runs = gx->runs ;
  const char *errors = 0 ;
  char *qq = hprintf (h, "Find Project \"%s\" ; >run", gx->project) ;
  AC_TABLE tbl = ac_bql_table (db, qq, 0, 0, &errors, h) ;
  int ir, irMax = tbl ? tbl->rows : 0 ;

  gx->runDict = runDict = dictHandleCreate (200, gx->h) ;	
  for (ir = 0 ; ir < irMax ; ir++) 
    dictAdd (runDict, ac_table_printable (tbl, ir, 0, "NULL"), 0) ;

  if (! gx->runs)
    runs = gx->runs = arrayHandleCreate (256, RC, gx->h) ;
  for (ir = 0 ; ir < irMax ; ir++) 
    {
      AC_OBJ Run = ac_table_obj (tbl, ir, 0, 0) ;
      int run ;
      if (dictFind (runDict, ac_name(Run), &run))
	{
	  AC_TABLE tbl2 = ac_obj_bql_table (Run, ">>Group", 0, &errors, h) ;
	  int jrMax = tbl2 ? tbl2->rows : 0 ;
      
	  if (jrMax)
	    {
	      RC *rc = arrayp (runs, run, RC) ;
	      KEYSET r2g = rc->r2g ;
	  
	      if (! r2g)
		r2g = rc->r2g = keySetHandleCreate (gx->h) ;
	      for (int g = 0, k =0, jr = 0 ; jr < jrMax ; jr++)
		if (dictFind (runDict, ac_table_printable (tbl2, jr, 0, 0), &g))
		  keySet (r2g, k++) = g ;
	    }
	}
    }

  fprintf (stderr, "# Found %d runs in project %s\n"
	   , dictMax (runDict) 
	   , gx->project
	   ) ;

  ac_free (h) ;
  return dictMax (runDict) ;
} /* gxGetRuns */

/*************************************************************************************/
typedef struct ddStruct {int run, n, runCumul ; int intron ; } DD ;

static int ddRunOrder (const void *a, const void *b)
{
  const DD *up = a ;
  const DD *vp = b ;
  int n = up->run - vp->run ; if (n) return n ;
  n = up->intron - vp->intron ; if (n) return n ;
  n = up->n - vp->n ; if (n) return -n ;
  return 0 ;
} /* ddRunOrder */

static int ddIntronOrder (const void *a, const void *b)
{
  const DD *up = a ;
  const DD *vp = b ;
  int n = up->intron - vp->intron ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;
  n = up->n - vp->n ; if (n) return -n ;
  return 0 ;
} /* ddIntronOrder */

/********************************************/
/************************ EXPORT ********************************/

static void gxShowAcceptors (ACEOUT ao, AC_OBJ obj, AC_TABLE introns, Array aa, DICT *runDict)
{         

  int iaMax = arrayMax (aa) ;
  for (int ia = 0 ; ia < iaMax ; ia++)
    {
      DD *up = arrp (aa, ia, DD) ;
      if (up->run)
	{
	  const char *ccpI = ac_table_printable (introns, up->intron - 1, 0, 0) ;
	  if (ao)
	    aceOutf (ao, "%s___%s\t%s\tiit\t%d\t%d\t%.2f\n", ac_name(obj), ccpI, dictName (runDict, up->run), up->n, up->runCumul, 100.0 *up->n/up->runCumul) ;
	  else
	    fprintf (stderr, "%s___%s\t%s\t%d\t%d\t%.2f\n", ac_name(obj), ccpI, dictName (runDict, up->run), up->n, up->runCumul, 100.0 *up->n/up->runCumul) ;
	}
    }
  if (! ao)
    fprintf(stderr, "\n\n") ;
} /* gxShowAcceptors */

/********************************************/

static int gxGetAcceptors (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *runDict = gx->runDict ;
  const char *errors = 0 ;
  AC_ITER iter ;
  AC_OBJ obj = 0 ;
  AC_HANDLE h2 = 0;
  int nn = 0 ;
  int nnn = 0 ;
  int minIntronCount = 20 ;
  int minRunCount = 20 ;
  char *qq = "select run,n from intron in @, run in intron->de_duo, n in run[1]" ; 
  Array aa = arrayHandleCreate (100, DD, h) ;
  DD *dd ;
  KEYSET runsOk = keySetHandleCreate (h) ;
  KEYSET nrOk = keySetHandleCreate (h) ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".differential_acceptors.tsf", gx->gzo, h) ;

  /*  iter = ac_dbquery_iter (gx->db, "Find Intron de_duo && IS *37832378; >A ", h) ; */
  iter = ac_dbquery_iter (gx->db, "Find Intron de_duo; >A ", h) ;
  while (ac_free (h2), ac_free (obj),  nn < 300000 && (obj = ac_iter_obj (iter)))
    {
      h2 = ac_new_handle () ;
      AC_TABLE introns = ac_tag_table (obj , "Intron", h2) ;
      int iiMax = introns ? introns->rows : 0 ;
      int iaMax = 0 ;
      int nRunOk = 0 ;
      int nIntronOk = 0 ;

      if (iiMax < 2)
	continue ;
      aa = arrayReCreate (aa, 100, DD) ;
      for (int ii = 0 ; ii < iiMax ; ii++)
	{
	  AC_OBJ Intron = ac_table_obj (introns, ii, 0, h2) ;
	  AC_TABLE duo = ac_obj_bql_table (Intron , qq, 0, &errors, h2) ;
	  int irMax = duo ? duo->rows : 0 ;
	  for (int ir = 0 ; ir < irMax ; ir++)
	    {
	      int run = 0 ;
	      if (dictFind (runDict, ac_table_printable (duo, ir, 0, "toto"), &run))
		{
		  int n = ac_table_int (duo, ir, 1, 0) ;
		  if (n)
		    {
		      dd = arrayp (aa, iaMax++, DD) ;
		      dd->run = run ; dd->intron = ii + 1 ; dd->n = n ; dd->runCumul = 0 ; 
		    }
		}
	    }
	}

      arraySort (aa, ddRunOrder) ;
      nRunOk = 0 ;
      for (int ia = 0 ; ia < iaMax ; ia++)
	{
	  int ka = ia, ja, nr = 0 ;
	  DD *vp, *up = arrp (aa, ia, DD) ;
	  int run = up->run ;
	  if (run)
	    {
	      for (vp = up, ja = ia ; ja < iaMax && vp->run == run  ; ja++, vp++)
		nr += vp->n ;
	      ka = ja - 1 ;
	      for (vp = up, ja = ia ; ja < iaMax && vp->run == run  ; ja++, vp++)
		if (nr < minRunCount)
		  vp->run = vp->intron = vp->n = 0 ; /* eliminate */
		else
		  vp->runCumul = nr ;
	      if (nr >= minRunCount)
		{
		  keySet (runsOk, nRunOk) = run ;
		  keySet (nrOk, nRunOk) = nr ;
		  nRunOk++ ;
		}
	    }
	  ia = ka ;
	}


      /* add a zero count for all interesting values*/
      for (int ii = 0 ; ii < iiMax ; ii++)
	for (int k = 0 ; k < nRunOk ; k++)
	  {
	    int run = keySet (runsOk, k) ;
	    dd = arrayp (aa, iaMax++, DD) ;
	    dd->run = run ; dd->intron = ii + 1 ; dd->n = 0 ; dd->runCumul = keySet(nrOk, k) ; 
	  }
      /* remove these fake values in favor of actual values */
      arraySort (aa, ddRunOrder) ;
      for (int ia = 0 ; ia < iaMax ; ia++)
	{
	  int ka = ia, ja ;
	  DD *vp, *up = arrp (aa, ia, DD) ;
	  int run = up->run ;
	  int intron = up->intron ;
	  if (run)
	    {
	      for (vp = up+1, ja = ia+1 ; ja < iaMax && vp->run == run && vp->intron == intron ; ja++, vp++)
		vp->run = vp->intron = vp->n = 0 ; /* eliminate */
	      ka = ja - 1 ;
	    }
	  ia = ka ;
	}
      /* check for intron support */
      arraySort (aa, ddIntronOrder) ;
      for (int ia = 0 ; ia < iaMax ; ia++)
	{
	  int ka = ia, ja, nr = 0 ;
	  DD *vp, *up = arrp (aa, ia, DD) ;
	  int intron = up->intron ;
	  if (intron)
	    {
	      for (vp = up, ja = ia ; ja < iaMax && vp->intron == intron ; ja++, vp++)
		nr += vp->n ;
	      ka = ja - 1 ;
	      if (nr < minIntronCount)
		for (vp = up, ja = ia ; ja < iaMax && vp->intron == intron  ; ja++, vp++)
		  vp->run = vp->intron = vp->n = 0 ; /* eliminate */
	    }
	  ia = ka ;
	}

      if (0)
	{
	  arraySort (aa, ddIntronOrder) ;
	  gxShowAcceptors (0, obj, introns, aa, runDict) ;
	}
      /* check for ratio variability */
      arraySort (aa, ddIntronOrder) ;
      for (int ia = 0 ; ia < iaMax ; ia++)
	{
	  int ka = ia, ja ;
	  DD *vp, *up = arrp (aa, ia, DD) ;
	  int intron = up->intron ;
	  float rMin = 99999, rMax = -9999 ;
	  if (intron)
	    {
	      for (vp = up, ja = ia ; ja < iaMax && vp->intron == intron ; ja++, vp++)
		{
		  float ratio = vp->n ; ratio /= vp->runCumul ;
		  if (ratio > rMax) rMax = ratio ;
		  if (ratio < rMin) rMin = ratio ;
		}
	      ka = ja - 1 ;
	      if (rMax < rMin + .20)
		{ /* eliminate */
		  for (vp = up, ja = ia ; ja < iaMax && vp->intron == intron ; ja++, vp++)
		    vp->run = vp->intron = vp->n = 0 ; /* eliminate */
		}
	      else
		nIntronOk++ ;
	    }
	  ia = ka ;
	}
      if (nRunOk >= 4 && nIntronOk > 1)
	nn += nIntronOk ;

      arraySort (aa, ddIntronOrder) ;
      if (nIntronOk && nnn++ < 3)
	gxShowAcceptors (0, obj, introns, aa, runDict) ;
      if (nIntronOk)
	gxShowAcceptors (ao, obj, introns, aa, runDict) ;
    }

  fprintf (stderr, "# Found %d differential acceptors in project %s\n"
	   , nn
	   , gx->project
	   ) ;

  ac_free (h) ;
  return nn ;
} /* getAcceptors */

/*************************************************************************************/

static int gxGetIntrons (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *runDict = gx->runDict ;
  AC_ITER iter ;
  AC_OBJ obj = 0 ;
  int nn = 0, ir ;
  ITR itr ;
  int kMin = 100 ;

  iter = ac_dbquery_iter (gx->db, "Find Intron de_uno", h) ;
  while (ac_free (obj),  nn < 300000 && (obj = ac_iter_obj (iter)))
    {
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE map = ac_tag_table (obj, "IntMap", h) ;
      AC_TABLE uno = ac_tag_table (obj, "de_uno", h) ;
      Array aa = 0 ;
      int ok = 0 ;
      
      memset (&itr, 0, sizeof (itr)) ;
      itr.intron = ac_obj_key (obj) ;
      if (map && map->rows && map->cols >= 3)
	{
	  itr.map = ac_table_key (map, 0, 0, 0) ;
	  itr.a1 = ac_table_int (map, 0, 1, 0) ;
	  itr.a2 = ac_table_int (map, 0, 2, 0) ;
	}
      itr.gene = ac_tag_key (obj, "Gene", 0) ;

      if (ac_has_tag (obj, "In_mRNA"))
	itr.type |= 4 ;
      if (ac_has_tag (obj, "Known_donor"))
	itr.type |= 1 ;
      if (ac_has_tag (obj, "Known_acceptor"))
	itr.type |= 2 ;
      if (itr.type > 4)
	itr.type = 4 ;
      if (uno && uno->rows)
	for (ir = 0 ; ir  < uno->rows ; ir++)
	  {
	    int run = 0 ;
	    const char *rNam = ac_table_printable (uno, ir, 0, "toto") ;
	    if (dictFind (runDict, rNam, &run))
	      {
		const char *cr = ac_table_printable (uno, ir, 1, 0) ;
		int k = cr ? atoi (cr) : 0 ;		
		if (k>0)
		  {
		    int k2 = -1, k1 = k ;
		    KEYSET ks = array (gx->histos, itr.type+ 8*run, KEYSET) ;
		    if (!ks)
		      ks = array (gx->histos, itr.type + 8*run, KEYSET) = keySetHandleCreate (gx->h) ;
		    while (k1 > 0)
		      { k1/=2 ; k2++ ;}
		    for (k1 = 0 ; k1 <= k2 ; k1++)
		      keySet (ks, k1)++ ;
		  }


		if (k >= kMin)
		  ok++ ;
		if (! aa)
		  aa = itr.aa = arrayCreate (dictMax (runDict)+1, int) ;
		keySet (aa, run) = k ;
	      }
	  }

      if (1)
	{
	  int run, k = 0 ;
	  k = ac_tag_int (obj, "NB", 0) ;
	  if (k  && dictFind (runDict, "NB", &run))
	    {
	      if (! aa)
		aa = itr.aa = arrayCreate (dictMax (runDict)+1, int) ;
	      keySet (aa, run) = k ;
	    }
	}

      if (ok > 0)
	{
	  AC_TABLE tbl ;
	  ITR *itrp = arrayp (gx->introns, KEYKEY (itr.intron), ITR) ;
	  *itrp = itr ;		
	  nn++ ;
	  
	  tbl = ac_tag_table (obj, "Same_donor", h) ;
	  itrp->nSD = tbl ? tbl->rows : 0 ;
	  for (ir = 0 ; tbl && ir  < tbl->rows ; ir++)
	    {
	      int intron2 = ac_table_key (tbl, ir, 0, 0) ;
	      if (intron2 && KEYKEY(intron2) < arrayMax (gx->introns) && 
		  arrp(gx->introns,KEYKEY(intron2), ITR)->intron)
		{
		  PAIR *pp = arrayp (gx->pairs, arrayMax (gx->pairs), PAIR) ;
		  pp->intron1 = itr.intron ;
		  pp->intron2 = intron2 ;
		  pp->gene = itr.gene ;
		}
	    }
	  tbl = ac_tag_table (obj, "Same_acceptor", h) ;
	  itrp->nSA = tbl ? tbl->rows : 0 ;
	  for (ir = 0 ; tbl && ir  < tbl->rows ; ir++)
	    {
	      int intron2 = ac_table_key (tbl, ir, 0, 0) ;
	      if (intron2 && KEYKEY(intron2) < arrayMax (gx->introns) && 
		  arrp(gx->introns,KEYKEY(intron2), ITR)->intron)
		{
		  PAIR *pp = arrayp (gx->pairs, arrayMax (gx->pairs), PAIR) ;
		  pp->intron1 = itr.intron ;
		  pp->intron2 = intron2 ;
		  pp->gene = itr.gene ;
		}
	    }
	}
      
      ac_free (h) ;
    }
  fprintf (stderr, "# Found %d introns %d pairs\n", arrayMax (gx->introns), arrayMax (gx->pairs)) ;
  
  ac_free (h) ;
  
  return nn ;
} /* gxGetIntrons */

/*************************************************************************************/
/*************************************************************************************/

static void gxOneAltIntrons (GX *gx)
{
  return ;
} /* gxAltIntrons */

/*************************************************************************************/

static void gxAltIntrons (GX *gx)
{
  gxOneAltIntrons (gx) ;
  return ;
} /* gxAltIntrons */

/*************************************************************************************/

static void gxAltHistoExport (GX *gx)
{
  DICT *runDict = gx->runDict ;
  int type, run, runMax = dictMax (gx->runDict) ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".support.tsf", gx->gzo, h) ;
  const char *types[] = {"New", "NewKnownD", "NewKnownA", "NewKnownDA","Known","Known","Known","Known","Known"} ;

  aceOutDate (ao, "###", "Number of introns with 1,2,4,8,16... supports in selected groups of runs") ;
  
  for (run = 1 ; run <= runMax ; run++)
    for (type = 0 ; type < 8 ; type++)
      {
	int k ;
	KEYSET ks = array (gx->histos, type+ 8*run, KEYSET) ;
	
	if (ks && keySetMax (ks))
	  {
	    aceOutf (ao, "\n%s\t%s\tiiiiiiiiiiiiiiiiiiiiiiii", dictName(runDict, run), types[type]) ;
	    for (k = 0 ; k < 24 ; k++)
	      aceOutf (ao, "\t%d", keySet (ks, k) ) ;
	  }
      }

  aceOutf (ao, "\n") ;
  ac_free (h) ;
  return ;
} /* gxAltHistoExport */

/*************************************************************************************/

static void gxAltIntronsExport (GX *gx)
{
  ITR *itr ;  
  PAIR *pp ;
  DICT *runDict = gx->runDict ;
  Array aa = gx->introns ;
  Array pairs = gx->pairs ;
  int pMax = pairs ? arrayMax (pairs) : 0 ;
  int ii, iiMax = aa ? arrayMax (aa) : 0 ;
  int run, runMax = dictMax (gx->runDict) ;
  int nn = 0 ;
  const char *types[] = {"New", "NewKnownD", "NewKnownA", "NewKnownDA","Known","Known","Known","Known","Known"} ;


  if (! iiMax)
    return ;
  
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".altIntrons.txt", gx->gzo, h) ;
  aceOutDate (ao, "###", "Number of reads supporting the introns in selected groups of runs") ;
  
  ACEOUT aoP = aceOutCreate (gx->outFileName, ".altPairs.txt", gx->gzo, h) ;
  aceOutDate (aoP, "###", "Branching ratios in selected groups of runs") ;
  
  aceOut (ao, "#Intron\tType\tLength\tNbAlternatifs\tGene\tChromosome\tDonor\tAcceptor") ;
  for (run = 1 ; run <= runMax ; run++)
    aceOutf (ao, "\t%s", dictName (gx->runDict, run)) ;
  
  for (ii = 0, itr = arrp (aa, 0, ITR) ; -nn < 12 && ii < iiMax  ; ii++, itr++)
    {
      if (itr->intron)
	{
	  BOOL ok = FALSE ;
	  int run, run2 ;

	  for (run = 1 ; !ok && run <= runMax ; run++)
	    for (run2 = run + 1 ; !ok && run2 <= runMax ; run2++)
	      {
		int x = itr->aa ? arr (itr->aa, run, int) : 0 ;
		int y = itr->aa ? arr (itr->aa, run2, int) : 0 ;
	  
		if (
		    (x > 2 * y  && (x > 300 || y > 10)) ||
		    (y > 2 * x  && (y > 300 || x > 10)) 
		    )
		  ok = TRUE ;
	      }
	  if (ok)
	    {
	      int ln = itr->a2 - itr->a1 ;
	      if (ln < 0) ln = -ln ;
	      ln++ ;
	      nn++ ;
	      aceOutf (ao, "\n%s", ac_key_name (itr->intron)) ;
	      aceOutf (ao, "\t%s", types[itr->type]) ;
	      aceOutf (ao, "\t%d", ln) ;
	      aceOutf (ao, "\t%d", itr->nSA + itr->nSD) ;
	      aceOutf (ao, "\t%s", itr->gene ? ac_key_name (itr->gene) : "") ;
	      aceOutf (ao, "\t%s", itr->map ? ac_key_name (itr->map) : "") ;
	      aceOutf (ao, "\t%d\t%d", itr->a1, itr->a2) ;
	      for (run = 1 ; run <= runMax ; run++)
		aceOutf (ao, "\t%d", itr->aa ? arr (itr->aa, run, int) : 0) ;
	    }
	}
    }
  aceOut (ao, "\n") ;

  nn = 0 ;
  
  aceOut (aoP, "# Intron1\tIntron2\tType1\tType2\tGene1\tGene2\tLength 1\tLength 2") ;
  for (run = 1 ; run <= runMax ; run++)
      aceOutf (aoP, "\tIntron 1 support in %s", dictName (runDict, run)) ;
  for (run = 1 ; run <= runMax ; run++)
      aceOutf (aoP, "\tIntron 2 support in %s", dictName (runDict, run)) ;


  for (run = 1 ; run <= runMax ; run++)
    aceOutf (aoP, "\tBranching ratio in %s", dictName (runDict, run)) ;

  for (nn = ii = 0, pp = arrp (pairs, 0, PAIR) ; -nn < 12 && ii < pMax  ; ii++, pp++)
    {
      int run2 ;
      int ii1 = pp->intron1 ;
      int ii2 = pp->intron2 ;
      ITR *itr1 = arrayp (aa, KEYKEY(ii1), ITR) ;
      ITR *itr2 = arrayp (aa, KEYKEY(ii2), ITR) ;
      BOOL ok = TRUE ;
      int z[runMax+1] ;
      for (run = 1 ; ok && run <= runMax ; run++)
	{
	  int a1 = arr(itr1->aa, run, int) ;
	  int b1 = arr(itr2->aa, run, int) ;
	  
	  if (a1 + b1 < 100)
	    ok = FALSE ;
	  else
	    z[run] = (100.0 * a1) / (1.0 * (a1 + b1)) ;
	}
      if (! ok) continue ;

      ok = FALSE ;
      for (run = 1 ; !ok && run <= runMax ; run++)
	for (run2 = run + 1 ; !ok && run2 <= runMax ; run2++)
	  if ((z[run] > z[run2] + 20 || z[run2] > z[run] + 20))
	    ok = TRUE ;

      if (ok)
	{
	  int ln1 = itr1->a2 - itr1->a1 ;
	  int ln2 = itr2->a2 - itr2->a1 ;
	  if (ln1 < 0) ln1 = -ln1 ;
	  if (ln2 < 0) ln2 = -ln2 ;
	  ln1++ ; ln2++ ;

	  nn++ ;
	  aceOutf (aoP, "\n%s", ac_key_name (ii1)) ;
	  aceOutf (aoP, "\t%s", ac_key_name (ii2)) ;
	  aceOutf (aoP, "\t%s", types[itr1->type]) ;
	  aceOutf (aoP, "\t%s", types[itr2->type]) ;
	  aceOutf (aoP, "\t%s", itr1->gene ? ac_key_name (itr1->gene) : "") ;
	  aceOutf (aoP, "\t%s", itr2->gene  ? (itr1->gene != itr2->gene ? ac_key_name (itr2->gene) : "idem") :"") ;
	  aceOutf (aoP, "\t%d\t%d", ln1, ln2) ;

	  for (run = 1 ; ok && run <= runMax ; run++)
	    {
	      int a1 = arr(itr1->aa, run, int) ;
	      aceOutf (aoP, "\t%d", a1) ;
	    }
	  for (run = 1 ; ok && run <= runMax ; run++)
	    {
	      int b1 = arr(itr2->aa, run, int) ;
	      aceOutf (aoP, "\t%d", b1) ;
	    }
	  for (run = 1 ; ok && run <= runMax ; run++)
	    aceOutf (aoP, "\t%d", z[run]) ;
	}
    }
  aceOut (aoP, "\n") ;
  fprintf (stderr, "# exported %d pairs\n", nn) ;  
  ac_free (h) ;
} /* gxAltIntronsExport */

/*************************************************************************************/

static void gxInit (GX *gx)
{
  AC_HANDLE h = gx->h ;
  
  gx->runDict = dictHandleCreate (32, h) ;
  gx->introns = arrayHandleCreate (50000, ITR, h) ;
  gx->pairs = arrayHandleCreate (5000, PAIR, h) ;
  gx->histos = arrayHandleCreate (30, KEYSET, h) ;
  return ;
} /* gxInit */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  fprintf  (stderr,
	    "// bestali.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   to analyse the alternative introns:\n"
	    "// Database\n"
	    "//   -db ACEDB : acedb database holding the semantics of the experiments\n"
	    "// Input: the program expects a .hits file exported by clipalign or by the -exportBest option\n"
	    "//   -i file_name [-gzi] : default stdin\n"
	    "//      if the file is called .gz or if the option -gzi is specified, invokes gunzip\n"
	    "//   -o out_file_name [-gzo] : default stdout\n"
	    "//      redirect the output, each option adds its own suffix to the file name\n"
	    "//      if -gzo is specified, invokes gzip and add a .gz suffix to the file name\n"
	    "// Actions\n"
	    "//   -to be defined :\n"
	    "//      xxx\n"
	    "// Help\n"
	    "//   -h -help --help : export this on line help\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' in the input files are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  GX gx ;

  AC_HANDLE h = ac_new_handle () ;
  const char *dbName = 0 ;
  const char *errors = 0 ;
  char commandBuf [4000] ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  memset (&gx, 0, sizeof (GX)) ;
  memset (commandBuf, 0, sizeof (commandBuf)) ;
  gx.h = h ;

  if (argc == 1)
    usage (0) ;
  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;


  /* optional arguments */
  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "--db", &dbName) ;
  getCmdLineOption (&argc, argv, "-o", &(gx.outFileName)) ;
  getCmdLineOption (&argc, argv, "-p", &(gx.project)) ;
  getCmdLineOption (&argc, argv, "--project", &(gx.project)) ;

  if (argc != 1)
    {
      int ix ;
      char *cp ;
      for (ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 3900 ; cp += strlen (cp), ix++)
	sprintf(cp, "%s ", argv[ix]) ;
      fprintf (stderr, "unknown argument, sorry\n") ;
      usage (commandBuf) ;
    }

  fprintf (stderr, "Start %s\n", timeShowNow ()) ;

  if (dbName)
    {
      gx.db = ac_open_db (dbName, &errors);
      if (! gx.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }
  else
    usage ("-db dbName missing, dbName should eb an existing INTRON_DB database") ;
  
  if (! gx.project)
    usage ("Missing argument --project") ;

  gxInit (&gx) ;

  if (1)
    {
      gxGetRuns (&gx) ;
      gxGetAcceptors (&gx) ;
    }
  if (0)
    {
      gxGetIntrons (&gx) ;
      gxAltIntrons (&gx) ;
      gxAltIntronsExport (&gx) ;
      gxAltHistoExport (&gx) ;
    }


  if (1)
    {
      /* the existence of the .done file proves that the code worked correctly */
      int mx ;
      messAllocMaxStatus (&mx) ;   
      fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
    }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

