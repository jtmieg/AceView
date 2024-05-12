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
#include <wiggle.h>


typedef struct gxStruct {
  AC_HANDLE h ; 
  AC_DB db ; 
  const char *outFileName, *inFileName, *deMrna ;
  const char *wiggleDir ;
  const char *project ; 
  const char *chrom ;
  BOOL gzo, gzi ;
  BOOL setFeet ; 
  BOOL setDA ;
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
	    aceOutf (    ao, "%s___%s\t%s\tiit\t%d\t%d\t%.2f\n", ac_name(obj), ccpI, dictName (runDict, up->run), up->n, up->runCumul, 100.0 *up->n/up->runCumul) ;
	  else
	    fprintf (stderr, "%s___%s\t%s\tiit\t%d\t%d\t%.2f\n", ac_name(obj), ccpI, dictName (runDict, up->run), up->n, up->runCumul, 100.0 *up->n/up->runCumul) ;
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

static int gxSetIntronLnIntMapFeet (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  AC_ITER iter ;
  AC_OBJ obj = 0, chromObj = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  const char *errors = 0 ;
  const char *chromDna = 0 ;
  char chromName[128] ;
  char feet[6] ;
  const char *gooddies ="gt_ag,gc_ag,at_ac,ct_ac" ;

  memset (chromName, 0, 128) ;

  iter = ac_dbquery_iter (gx->db, "Find Intron type || ! Length || ! type", h) ;
  while (ac_free (obj),  obj = ac_iter_obj (iter))
    {
      AC_HANDLE h1 = ac_new_handle () ;
      char *chrom = strnew (ac_name(obj), h1), *cq ;
      int k = 0, a1 = 0, a2 = 0 ;
      int sliding = 0 ;
      vtxtClear (txt) ;

      cq = strstr (chrom, "__") ;
      if (cq)
	{ 
	  *cq = 0 ; cq += 2 ;
	  k = sscanf (cq,  "%d_%d", &a1, &a2) ;
	  if (a1 <= 0 || a2 <= 0)
	    k = 0 ;
	}
      if (k == 2)
	{
	  int ln = a2 - a1 ;
	  if (ln < 0) ln = - ln ;
	  ln++ ;
	  vtxtPrintf (txt, "Intron %s\nIntMap %s %d %d\nLength %d\n"
		      , ac_name (obj), chrom, a1, a2, ln
		      ) ;
	}
      else
	vtxtPrintf (txt, "-D Intron %s\n\n", ac_name(obj)) ;

      /* get intron feet  */
      if (k == 2 && ! ac_has_tag (obj, "Type"))
	{
	  if (! chromDna || strcmp (chrom, chromName))
	    {
	      if (strlen (chrom) > 127)
		messcrash ("ChromName too long %s", chromName) ;
	      strncpy (chromName, chrom, 127) ;
	      ac_free (chromDna) ;
	      chromObj = ac_get_obj (gx->db, "Sequence", chromName, h) ;
	      chromDna = ac_obj_dna (chromObj, h) ;
	    }
          if (chromDna)
	    {
	      /* get sliding */
	      for (int i = 0 ; chromDna[a1-1+i] == chromDna[a2+i] ; i++)
		sliding++ ;
	      for (int i = 0 ; chromDna[a1-2-i] == chromDna[a2-1-i] ; i++)
		sliding++ ;

	      /* get type */
	      if (a1 < a2)
		{
		  feet[0] = chromDna[a1-1] ;
		  feet[1] = chromDna[a1] ;
		  feet[2] = '_' ;
		  feet[3] = chromDna[a2-2] ;
		  feet[4] = chromDna[a2-1] ;
		  feet[5] = 0 ;
		}
	      else
		{
		  feet[0] = complementLetter(chromDna[a1-1]) ;
		  feet[1] = complementLetter(chromDna[a1-2]) ;
		  feet[2] = '_' ;
		  feet[3] = complementLetter(chromDna[a2]) ;
		  feet[4] = complementLetter(chromDna[a2-1]) ;
		  feet[5] = 0 ;
		}
	      if (strcasestr (gooddies, feet))
		vtxtPrintf (txt, "%s\n", feet) ;
	      else
		vtxtPrintf (txt, "Other %s\n", feet) ;
	      if (sliding)
		vtxtPrintf (txt, "Sliding %d\n", sliding) ;
	      else
		vtxtPrintf (txt, "-D Siding\n") ;
	    }
	}
	  /* set echo */
      if (k == 2)
	{
	  char *echo = hprintf (h1, "%s__%d_%d", chrom, a2, a1) ;
	  KEY key = ac_get_key (gx->db, "Intron", echo) ;
	  if (key)
	    {
	      vtxtPrintf (txt, "Has_echo %s\n", ac_key_name(key)) ;
	      if (! strcmp (feet, "ct_ac") || ! strcmp (feet, "ct_gc"))
		vtxtPrintf (txt, "Is_echo\n") ;
	    }
	}
    
      if (vtxtPtr (txt))
	ac_parse (gx->db, vtxtPtr (txt), &errors, 0, h1) ;
      ac_free (h1) ;
    }

  ac_free (h) ;
  return nn ;
} /* gxGetIntronLnIntMapFeet */

/*************************************************************************************/

static int gxSetDA (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  AC_ITER iter ;
  AC_OBJ Intron = 0, chromObj = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  const char *errors = 0 ;
  const char *chromDna = 0 ;
  char chromName[128] ;
  char feet[128], foot[3] ;
  const char *goodD ="gt,gc,at,ct" ;
  const char *goodA ="ag,ac" ;
  
  memset (chromName, 0, 128) ;
  
  iter = ac_dbquery_iter (gx->db, "Find Intron ! D || ! A", h) ;
  while (ac_free (Intron),  Intron = ac_iter_obj (iter))
    {
      AC_HANDLE h1 = ac_new_handle () ;
      char *chrom = strnew (ac_name(Intron), h1), *strand, *cq ;
      int k = 0, a1 = 0, a2 = 0, A1, A2, d1, d2, d3, d4, A3, A4 ;
      vtxtClear (txt) ;
      nn++ ;
      if (nn % 10000 == 1) fprintf(stderr, ".... setDA %d %s\n", nn, ac_name(Intron)) ;
      cq = strstr (chrom, "__") ;
      if (cq)
	{ 
	  *cq = 0 ; cq += 2 ;
	  k = sscanf (cq,  "%d_%d", &a1, &a2) ;
	  if (a1 <= 0 || a2 <= 0)
	    k = 0 ;
	}
      if (k != 2)
	continue ;
      
      /* get genome */
      if (! chromDna || strcmp (chrom, chromName))
	{
	  if (strlen (chrom) > 127)
	    messcrash ("ChromName too long %s", chromName) ;
	  strncpy (chromName, chrom, 127) ;
	  ac_free (chromDna) ;
	  chromObj = ac_get_obj (gx->db, "Sequence", chromName, h) ;
	  chromDna = ac_obj_dna (chromObj, h) ;
	}
      if (! chromDna)
	continue ;
      /* get donor/acceptor feet  */
      if (a1 < a2)
	{
	  strand = "f" ;
	  d1 = a1 - 50 ; d2 = a1 + 49 ; d3 = a1 - 1 ; d4 = a1 ;
	  A1 = a2 - 49 ; A2 = a2 + 50 ; A3 = a2 ; A4 = a2 + 1 ;
	}
      else
	{
	  strand = "r" ;
	  d1 = a1 + 50 ; d2 = a1 - 49 ; d3 = a1 + 1 ; d4 = a1 ;
	  A1 = a2 + 49 ; A2 = a2 - 50 ; A3 = a2 ; A4 = a2 + 1 ;
	}
      
      vtxtPrintf (txt, "Intron %s\n-D DA\n\n", ac_name (Intron)) ;
      
      /* donor */
      vtxtPrintf (txt, "Donor %s__%d_%s\n", chromName, d4, strand) ;
      vtxtPrintf (txt, "IntMap %s %d %d\n", chromName, d3, d4) ;
      vtxtPrintf (txt, "Intron %s\n", ac_name (Intron)) ;
      if (d1 < d2)
	{
	  int i, j ;
	  for (j = 0, i = d1 ; i <= d2 ; i++, j++)
	    {
	      if (j == 50)
		{ feet[50] = feet[51] = '-' ; j = 52 ; }
	      feet [j] = chromDna[i-1] ;
	    }
	  feet[j] = 0 ;
	}
      else
	{
	  int i, j ;
	  for (j = 0, i = d1 ; i >= d2 ; i--, j++)
	    {
	      if (j == 50)
		{ feet[50] = feet[51] = '-' ; j = 52 ; }
	      feet [j] = ace_lower(complementLetter(chromDna[i-1])) ;
	    }
	  feet[j] = 0 ;
	}
      foot[0] = feet[52] ; foot[1] = feet[53] ; foot[2] = 0 ;
      if (! strstr (goodD, foot))
	vtxtPrintf (txt, "Other") ;
      vtxtPrintf (txt, " %s\n", foot) ;
      vtxtPrintf (txt, "Motifs %s\n\n", feet) ;
      
      /* acceptor */
      vtxtPrintf (txt, "Acceptor %s__%d_%s\n", chromName, A3, strand) ;
      vtxtPrintf (txt, "IntMap %s__%d_%d\n", chromName, A3, A4) ;
      vtxtPrintf (txt, "Intron %s\n", ac_name (Intron)) ;
      if (A1 < A2)
	{
	  int i, j ;
	  for (j = 0, i = A1 ; i <= A2 ; i++, j++)
	    {
	      if (j == 50)
		{ feet[50] = feet[51] = '-' ; j = 52 ; }
	    feet [j] = chromDna[i-1] ;
	    }
	  feet[j] = 0 ;
	}
      else
	{
	  int i, j ;
	  for (j = 0, i = A1 ; i >= A2 ; i--, j++)
	    {
	      if (j == 50)
		{ feet[50] = feet[51] = '-' ; j = 52 ; }
	      feet [j] = ace_lower(complementLetter(chromDna[i-1])) ;
	    }
	  feet[j] = 0 ;
	}
      foot[0] = feet[48] ; foot[1] = feet[49] ; foot[2] = 0 ;
      if (! strstr (goodA, foot))
	vtxtPrintf (txt, "Other") ;
      vtxtPrintf (txt, " %s\n", foot) ;
      vtxtPrintf (txt, "Motifs %s\n\n", feet) ;
      
      if (vtxtPtr (txt))
	ac_parse (gx->db, vtxtPtr (txt), &errors, 0, h1) ;
      ac_free (h1) ;
    }
  
  ac_free (h) ;
  return nn ;
} /* gxSetDA */

/*************************************************************************************/

static int gxSetDAsupport (GX *gx) 
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  int pass ;
  AC_ITER iter ;
  AC_OBJ Donor = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  const char *errors = 0 ;

  for (pass = 0 ; pass < 2 ; pass++)
    {
      if (pass == 0)
	iter = ac_dbquery_iter (gx->db, "Find Donor", h) ;
      else
	iter = ac_dbquery_iter (gx->db, "Find Acceptor", h) ;
      
      while (ac_free (Donor),  Donor = ac_iter_obj (iter))
	{
	  AC_HANDLE h1 = ac_new_handle () ;
	  AC_TABLE introns = ac_tag_table (Donor, "Intron", h1) ;
	  const char *qq = "select r, n from d in @, ii in d->intron, r in ii->de_duo, n in r[1] where n>0" ;
	  AC_TABLE counts = ac_obj_bql_table (Donor, qq, 0, 0, h1) ; 
	  AC_TABLE genes = ac_obj_bql_table (Donor, "select g from d in @, ii in d->Intron, g in ii->gene", 0, 0, h1) ; 
	  AC_TABLE gfs = ac_obj_bql_table (Donor, "select g from d in @, ii in d->Intron, g in ii->From_genefinder", 0, 0, h1) ; 
	  AC_TABLE tgs = ac_obj_bql_table (Donor, "select g from d in @, ii in d->Intron, g in ii->From_gene", 0, 0, h1) ; 
	  AC_TABLE sTypes = ac_obj_bql_table (Donor, "select st from d in @, ii in d->Intron, st in ii->sType", 0, 0, h1) ; 
	  AC_TABLE mrnas = ac_obj_bql_table (Donor, "select m,m1 from d in @, ii in d->Intron, m in ii->In_mrna, m1 in m[1]", 0, 0, h1) ; 
	  int ir, jr ;

	  vtxtClear (txt) ;
	  nn++ ;
	  if (0*nn>1000)
	    { nn = 0 ; break ; }
	  if (nn % 10000 == 1) fprintf(stderr, ".... setSupport %d %s\n", nn, ac_name(Donor)) ;
	  
	  if (introns)
	    for (ir = 0 ; ir < introns->rows ; ir++)
	      for (jr = ir + 1 ; jr < introns->rows ; jr++)
		vtxtPrintf (txt , "Intron %s\n%s %s\n\n"
			    , ac_table_printable (introns, ir, 0, "xxx")
			    , pass == 0 ? "Same_donor" : "Same_acceptor"
			    , ac_table_printable (introns, jr, 0, "xxx")
			    ) 
		  ;
	  vtxtPrintf (txt , "\n%s %s\n"
		      , pass == 0 ? "Donor" : "Acceptor"
		      , ac_name (Donor) 
		      ) ;
	  if (genes)
	    for (ir = 0 ; ir < genes->rows ; ir++)
	      {
		vtxtPrintf (txt, "Gene %s\n"
			    , ac_table_printable (genes, ir, 0, "xxx") 
			    ) ;
	      }
	  if (gfs)
	    for (ir = 0 ; ir < gfs->rows ; ir++)
	      {
		vtxtPrintf (txt, "From_genefinder %s\n"
			    , ac_table_printable (gfs, ir, 0, "xxx") 
			    ) ;
	      }
	  if (tgs)
	    for (ir = 0 ; ir < tgs->rows ; ir++)
	      {
		vtxtPrintf (txt, "From_gene %s\n"
			    , ac_table_printable (tgs, ir, 0, "xxx") 
			    ) ;
	      }
	  if (mrnas)
	    for (ir = 0 ; ir < mrnas->rows ; ir++)
	      {
		vtxtPrintf (txt, "In_mRNA %s %d\n"
			    , ac_table_printable (mrnas, ir, 0, "xxx") 
			    , ac_table_int (mrnas, ir, 1, 0)
			    ) ;
	      }
	  if (sTypes)
	    for (ir = 0 ; ir < sTypes->rows ; ir++)
	      {
		vtxtPrintf (txt, "sType %s\n"
			    , ac_table_printable (sTypes, ir, 0, "xxx") 
			    ) ;
	    }
	  int nnn = 0 ;
	  if (counts)
	    for (ir = 0 ; ir < counts->rows ; ir++)
	      {
		int n = 0 ;
		KEY r1 ;
		r1 = ac_table_key (counts, ir, 0, 0) ;
		for (jr = ir ; jr < counts->rows && r1 == ac_table_key (counts, jr, 0, 0) ; jr++)
		n += ac_table_int (counts, jr, 1, 0) ;
		ir = jr - 1 ; /* reposition */
		nnn += n ;
		vtxtPrintf (txt, "De_uno %s %d\n"
			    , ac_table_printable (counts, ir, 0, "xxx") 
			    , n
			    ) ;
	      }
	  if (nnn)
	    vtxtPrintf (txt, "RNA_seq %d\n\n", nnn) ;
	  
	  if (vtxtPtr (txt))
	    ac_parse (gx->db, vtxtPtr (txt), &errors, 0, h1) ;
	  if (*errors)
	    {
	      fprintf(stderr, "%s\n", errors) ;
	      invokeDebugger () ;
	    }
	  ac_free (h1) ;
	}
    }
  
  ac_free (h) ;
  h = ac_new_handle () ;
  if (1)
    {
      char *qq = "Find Intron Gene && (same_donor || same_acceptor)" ;
      AC_ITER iter = ac_dbquery_iter (gx->db, qq, h) ;
      AC_OBJ Intron = 0 ;
      while (ac_free (Intron),  Intron = ac_iter_obj (iter))
	{
	  AC_HANDLE h1 = ac_new_handle () ;
	  AC_TABLE tbl1 = ac_obj_bql_table (Intron, "select g, i2 from ii in @, g in ii->gene, i2 in ii->same_donor", 0, 0, h1) ;
	  AC_TABLE tbl2 = ac_obj_bql_table (Intron, "select g, i2 from ii in @, g in ii->gene, i2 in ii->same_acceptor", 0, 0, h1) ;

	  vtxtClear (txt) ;
	  for (int ir = 0 ; ir < tbl1->rows ; ir++)
	    vtxtPrintf (txt, "Gene %s\nIntron %s\n\n"
			, ac_table_printable (tbl1, ir, 0, "xxx")
			, ac_table_printable (tbl1, ir, 1, "xxx")
			) ;
	  
	  for (int ir = 0 ; ir < tbl2->rows ; ir++)
	    vtxtPrintf (txt, "Gene %s\nIntron %s\n\n"
			, ac_table_printable (tbl1, ir, 0, "xxx")
			, ac_table_printable (tbl1, ir, 1, "xxx")
			) ;

	  ac_parse (gx->db, vtxtPtr (txt), &errors, 0, h1) ;
	}
    }
    
  ac_free (h) ;
  return nn ;
} /* gxSetDAsupport */

/*************************************************************************************/
typedef struct flatStruct { int da, run, n1, n2 ;} FLAT ;

static int gxFlatOrder (const void *a, const void *b)
{
  const FLAT *up = a ;
  const FLAT *vp = b ;
  int n ;
  n = up->da - vp->da ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;
  return 0 ;
} /* gxFlatOrder */

/*************************************************************************************/
/* parse  intron support count from mRNA alignemnts as found in tmp/INTRONRUNS
 * take the max between genomic and mrna counts
 */
static void gxDeMrna (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (gx->deMrna, 0, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  char *cp ;
  DICT *runDict = dictHandleCreate (1000, h) ;
  KEYSET ks = keySetHandleCreate (h) ;
  
  if (ai)
    {
      AC_HANDLE h1 = 0 ;
      AC_OBJ Intron = 0 ;
      AC_TABLE tbl ;
      char intronNam[256] ;
      int nModif = 0 ;

      memset (intronNam, 0, 256) ;
      aceInSpecial (ai, "\n\t") ;
      while (aceInCard (ai))
	{
	  int n1 = 0, n2 = 0 ;
	  int run ;
	  
	  cp = aceInWord (ai) ;
	  if (!cp || *cp == '#')
	    continue ;
	  if (strcasecmp (cp, intronNam))
	    {
	      ac_free (h1) ;
	      h1 = ac_new_handle () ;
	      Intron = 0 ; tbl = 0 ;
	      ks = keySetReCreate (ks) ;
	      if (nModif)
		ac_parse (gx->db, vtxtPtr (txt), 0, 0, h1) ;
	      nModif = 0 ;
	      vtxtClear (txt) ;
	      vtxtPrintf (txt, "\nIntron %s\n", cp) ;
	      strncpy (intronNam, cp,255) ;
	      Intron = ac_get_obj (gx->db, "Intron", cp, h1) ;

	      if (Intron)
		{
		  tbl = ac_tag_table (Intron, "de_duo", h1) ;
		  if (tbl)
		    for (int ir = 0 ; ir < tbl->rows ; ir++)
		      {
			dictAdd (runDict, ac_table_printable (tbl, ir, 0, "xxx"), &run) ;
			keySet (ks, run) = ac_table_int (tbl, ir, 1, 0) ;
		      }		      
		}
	    }
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  dictAdd (runDict, cp, &run) ;
	  n1 = 0 ;
	  if (run < keySetMax (ks))
	    n1 = keySet (ks, run) ;
	  
	  aceInStep (ai, '\t') ;
	  aceInInt (ai, &n2) ;
	  if (n2 > n1)
	    {
	      nModif++ ;
	      vtxtPrintf (txt, "de_duo %s %d\n"
			  , dictName (runDict, run)
			  , n2
			  ) ;
	    }
	}

      if (nModif)
	ac_parse (gx->db, vtxtPtr (txt), 0, 0, h1) ;
      ac_free (h1) ;
    }
  
  ac_free (h) ;
} /* gxDeMrna */

/*************************************************************************************/

static void gxFlatScan (BigArray aa, int run, AC_TABLE das, WIGGLE *sx)
{
  long int ii = bigArrayMax (aa) ;
  int da, daMax = das->rows ;
  Array www = arr (sx->aaa, 0, Array) ;
  WIGGLEPOINT *wp = arrp (www, 0, WIGGLEPOINT);
  int step = wp[1].x - wp[0].x ;
  int x0 = wp[0].x ;
  int xMax = arrayMax (www) ;

  if (xMax && step)
    for (da = 0 ; da < daMax && -da < 100 ; da++)
      {
	int x = ac_table_int (das, da, 0, 0) ;
	x = (x - x0)/step ;
	if (x > 1 && x < xMax - 2)
	  {
	    FLAT *up = bigArrayp (aa, ii++, FLAT) ;
	    up->da = da ; up->run = run ; 
	    wp = arrp (www, x - 2, WIGGLEPOINT) ;
	    up->n1 = wp->y ;
	    wp = arrp (www, x + 2, WIGGLEPOINT) ;
	    up->n2 = wp->y ;
	  }
      }
  return ;
} /* gxFlatScan */

/*************************************************************************************/

static int gxSetDAflatSupport (GX *gx) 
{
  const char *qq[4] ;
  qq[0] =  "select a1, da from da in  ?Donor where da ~ \"*_f\", chrom in da->intMap, a1 in chrom[1]" ;
  qq[1] =  "select a1, da from da in  ?Acceptor where da ~ \"*_f\", chrom in da->intMap, a1 in chrom[1]" ;
  qq[2] =  "select a1, da from da in  ?Donor where da ~ \"*_r\", chrom in da->intMap, a1 in chrom[1]" ;
  qq[3] =  "select a1, da from da in  ?Acceptor where da ~ \"*_r\", chrom in da->intMap, a1 in chrom[1]" ;

  for (int pass = 0 ; pass < 4 ; pass += 2)
    {
      
      AC_HANDLE h = ac_new_handle () ;
      int runMax = gx->runs ? dictMax (gx->runDict) + 1 : 0 ;
      BigArray aaD, aaA ;
      vTXT txt = vtxtHandleCreate (h) ;
      const char *errors = 0 ;
      AC_TABLE  dasD = ac_bql_table (gx->db, qq[pass+0], 0, 0, &errors, h) ;
      AC_TABLE  dasA = ac_bql_table (gx->db, qq[pass+1], 0, 0, &errors, h) ;
      
      aaD = bigArrayHandleCreate (100000, FLAT, h) ;
      aaA = bigArrayHandleCreate (100000, FLAT, h) ;
      for (int run = 1 ; run < runMax && -run < 3 ; run++)
	{
	  AC_HANDLE h1 = ac_new_handle () ;

	  char *fNam = hprintf(h1, "%s/%s/%s/R.chrom.u.%c.BF.gz", gx->wiggleDir, dictName(gx->runDict, run), gx->chrom, pass/2 == 0 ? 'f' : 'r') ;
	  WIGGLE sx ;
	  memset (&sx, 0, sizeof (WIGGLE)) ; 
	  sx.ai = aceInCreate (fNam, 0, h1) ;
	  sx.noRemap = TRUE ;
	  sx.out_step = 10 ;
	  sx.aaa = arrayHandleCreate (12, Array, h1) ;
	  array (sx.aaa, 0, Array) = arrayHandleCreate (100000, WIGGLEPOINT, h1) ;
	  if (sx.ai)
	    {
	      sxWiggleParse (&sx, 0, 0) ; /* whole wiggle, no zoning */
	      ac_free (sx.ai) ;

	      gxFlatScan (aaD, run, dasD, &sx) ;
	      gxFlatScan (aaA, run, dasA, &sx) ;
	    }
	  ac_free (h1) ;
	}
      for (int pass2 = 0 ; pass2 < 2 ; pass2++)
	{
	  BigArray aa = (pass2 == 0 ? aaD : aaA) ;
	  AC_TABLE das = (pass2 == 0 ? dasD : dasA) ;
	  if (bigArrayMax (aa)) 
	    {
	      long int ii, jj, iMax = bigArrayMax (aa) ;
	      bigArraySort (aa, gxFlatOrder) ;
	      for (ii = 0 ; ii < iMax ; ii++)
		{
		  FLAT *vp, *up = bigArrp (aa, ii, FLAT) ;
		  int da = up->da ;
		  
		  vtxtClear (txt) ;
		  vtxtPrintf (txt, "%s %s\n"
			      , pass2 == 0 ? "Donor" : "Acceptor"
			      , ac_table_printable (das, da, 1, "xxx")
			      ) ;
		  for (jj = ii, vp = up ; jj < iMax && vp->da == da ; jj++, vp++)
		    {
		      int n1, n2 ;
		      switch (pass + pass2)
			{
			case 0: n1 = vp->n1 ; n2 = vp->n2 ; break ;
			case 1: n1 = vp->n2 ; n2 = vp->n1 ; break ;
			case 2: n1 = vp->n2 ; n2 = vp->n1 ; break ;
			case 3: n1 = vp->n1 ; n2 = vp->n2 ; break ;
			}
		      if (n1 + n2 >= 10)
			vtxtPrintf (txt, "Sponge %s %d %d\n", dictName (gx->runDict, vp->run), n1, n2) ;
		      else
			vtxtPrintf (txt, "-D Sponge %s\n", dictName (gx->runDict, vp->run)) ;
		    }
		  ii = jj - 1 ;
		  vtxtPrintf (txt, "\n") ;
		  ac_parse (gx->db, vtxtPtr (txt), &errors, 0, h) ;
		  if (*errors)
		    {
		      fprintf(stderr, "%s\n", errors) ;
		      invokeDebugger () ;
		    }
		}
	    }
	}
      ac_free (h) ;
    }  /* loop on the 2 strands */
  return 0 ;
} /* gxSetDAflatSupport */

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
      AC_HANDLE h1 = ac_new_handle () ;
      AC_TABLE map = ac_tag_table (obj, "IntMap", h1) ;
      AC_TABLE uno = ac_tag_table (obj, "de_uno", h1) ;
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
	  
	  tbl = ac_tag_table (obj, "Same_donor", h1) ;
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
	  tbl = ac_tag_table (obj, "Same_acceptor", h1) ;
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
      
      ac_free (h1) ;
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
	    "// altintrons.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, May 2024, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "//   to analyse the alternative introns:\n"
	    "// Database\n"
	    "//   --db ACEDB : acedb database holding the semantics of the experiments\n"
	    "//   -p --project project_name : only use the runs belonging to the project\n"
	    "// Input: the program expects a .hits file exported by clipalign or by the -exportBest option\n"
	    "//   -i file_name [--gzi] : default stdin\n"
	    "//      if the file is called .gz or if the option -gzi is specified, invokes gunzip\n"
	    "//   -o out_file_name [--gzo] : default stdout\n"
	    "//      redirect the output, each option adds its own suffix to the file name\n"
	    "//      if -gzo is specified, invokes gzip and add a .gz suffix to the file name\n"
	    "// Actions\n"
	    "//   --setFeet :\n"
	    "//      for all introns, set IntMap, length, type gt_ag ...\n"
	    "//   --setDA :\n"
	    "//      for all introns, create the associated donors/acceptors\n"
	    "//      xxx\n"
	    "//   -w --wiggleDir dirName:\n"
	    "//      for each donor acceptor at position x, parse wiggle value around x+-13\n"
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

  gx.gzi = getCmdLineBool (&argc, argv, "--gzi") ;
  gx.gzo = getCmdLineBool (&argc, argv, "--gzo") ;
  gx.setFeet = getCmdLineBool (&argc, argv, "--setFeet") ;
  gx.setDA = getCmdLineBool (&argc, argv, "--setDA") ;

  /* optional arguments */
  getCmdLineOption (&argc, argv, "-db", &dbName) ;
  getCmdLineOption (&argc, argv, "--db", &dbName) ;
  getCmdLineOption (&argc, argv, "-o", &(gx.outFileName)) ;
  getCmdLineOption (&argc, argv, "-w", &(gx.wiggleDir)) ;
  getCmdLineOption (&argc, argv, "--wiggleDir", &(gx.wiggleDir)) ;
  getCmdLineOption (&argc, argv, "-p", &(gx.project)) ;
  getCmdLineOption (&argc, argv, "--project", &(gx.project)) ;
  getCmdLineOption (&argc, argv, "--chrom", &(gx.chrom)) ;
  getCmdLineOption (&argc, argv, "--deMrna", &(gx.deMrna)) ;

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

  if (gx.setDA && ! gx.project)
    usage ("Missing argument --project needed if --setDA") ;
  if (gx.wiggleDir && !gx.chrom)
    usage ("Missing argument --chrom needed if --wiggleDir") ;

  if (gx.deMrna)
    gxDeMrna (&gx) ;
  
  if (gx.setFeet)
      gxSetIntronLnIntMapFeet (&gx) ;
  if (gx.setDA)
    {
      gxSetDA (&gx) ;
      gxSetDAsupport (&gx) ;
    }
  if (gx.wiggleDir)
    {
      gxGetRuns (&gx) ;
      gxSetDAflatSupport (&gx) ;
    }
  if (0)
    {
      gxGetRuns (&gx) ;

    }
  if (0)
    {
      gxGetAcceptors (&gx) ;
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

