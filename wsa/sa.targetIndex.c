/*
 * sa.targetIndex.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to the creation, writing and reading of the target index
*/

#include "sa.h"

/**************************************************************/
/*************************************************************************************/
/* check the existence of the target files
 * identify their absolute file names
 * associate each file to a target class and to its optional parameters
 */
Array saTargetParseConfig (PP *pp)
{
  const char *tConfigFileName = pp->createIndex ? pp->tConfigFileName : hprintf (pp->h, "%s/tConfig", pp->indexName) ;
  Array tcs = arrayHandleCreate (64, TC, pp->h) ;
  TC *tc = 0 ;
  int nn = 0 ;
  
  if (pp->tFileName)
    {
      char *cr ;
      const char *cp = pp->tFileName ;
      cr = filName (cp, 0, "r") ;
      if (! cr)
	messcrash ("\nCannot open the target file -t %s\n", cp) ;
      tc = arrayp (tcs, nn++, TC) ;
      tc->fileName = strnew (cr, pp->h) ;
      tc->targetClass = 'G' ;
      
      tc->format = FASTA ; /* default */
      if (strstr (tc->fileName, ".fasta")) tc->format = FASTA ;
      if (strstr (tc->fileName, ".fna")) tc->format = FASTA ;
      if (strstr (tc->fileName, ".fa")) tc->format = FASTA ;

      /* user can override the defaults */
      if (pp->raw) tc->format = RAW ;
      if (pp->fasta) tc->format = FASTA ;
    }
  else if (tConfigFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *fDict = dictHandleCreate (64, h) ;
      ACEIN ai = aceInCreate (tConfigFileName, 0, h) ;
      int nn = 0, line = 0 ;

      while (aceInCard (ai))
	{
	  char cc, *cq, *cr, *cp = aceInWord (ai) ;
	  
	  line++ ;
	  if (! cp || ! *cp || *cp == '#')
	    continue ;
	  /* target class */
	  if (!cp || !*cp || *cp == '#')
	    continue ;
	  cc = *cp ;
	  if (!(cc >= 'A' && cc <= 'Z')) cc = 0 ;
	  if (cp[1]) cc = 0 ;
	  if (! cc)
	    messcrash ("\n\nThe target class must be specified as a single character (A-Z), not %s,  at line %d of -T target config file %s\n Please try sortalign --help\n"
		       , cp
		       , line
		       , pp->tConfigFileName
		       ) ;
	  if (! strchr ("GMCREATIBV", cc))
	    messcrash ("\n\nThe target class must be specified as a single character [GMCREATIBV], not %c,  at line %d of -T target config file %s\n Please try sortalign --help\n"
		       , cc
		       , line
		       , pp->tConfigFileName
		       ) ;
	    
	  tc = arrayp (tcs, nn++, TC) ;
	  tc->targetClass = cc ;
	  /* file names */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (!cp || !*cp || *cp == '#')
	    messcrash ("\nNo file name at line %d of file -T %s\n try sortalign --help\n"
		       , line
		       , pp->tConfigFileName
		       ) ;
	  cr = filName (cp, 0, "r") ;
	  if (! cr)
	    messcrash ("\nCannot open the target file -t %s\n", cp) ;
	  if (! dictAdd (fDict, cr, 0))
	    messcrash ("\nDuplicate target file name %s\n at line %d of file -T %s\n try sortalign --help\n"
		       , cr
		       , line
		       , pp->tConfigFileName
		       ) ;
	  tc->fileName = strnew (cr, pp->h) ;
	  tc->format = FASTA ; /* default */
	  if (cc == 'I')
	    {
	      if (strstr (cp, ".introns") || strstr (cp, ".gff") || strstr (cp, ".gtf"))
		tc->format = INTRONS ;
	      else
		messcrash ("\n\nThe Introns must be specified via a .gtf or a .gff file.\n try sortalign --help\n"
			   , cp
			   , line
			   , pp->tConfigFileName
		       ) ;
	    }
	  /* options */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  
	  while (cp)
	    {
	      if (*cp == '#')
		break ;
	      cq = strchr (cp, ',') ;
	      if (cq)
		*cq++ = 0 ;
	      if (! strcasecmp (cp, "fasta")) tc->format = FASTA ;
	      if (! strcasecmp (cp, "raw")) tc->format = RAW ;
	      cp = cq ;
	    }
	} 
      ac_free (h) ;
    }
  printf ("Found %d target files\n", arrayMax (tcs)) ;
  return tcs ;
} /* saTargetParseConfig */    
  
/**************************************************************/

static void storeTargetIndex (PP *pp, int tStep) 
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  long int nn = 0 ;
  int NN = pp->nIndex ;
  BB *bbG = &(pp->bbG) ;
  
  /* export the code words */
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s/cws.sortali.%d", pp->indexName, k) ;
      bigArrayMapWrite (bbG->cwsN[k], fNam) ;
      nn += bigArrayMax (bbG->cwsN[k]) ;
    }
  fprintf (stderr, "genomeCreateBinary exported %ld seed records\n", nn) ;

  /* export the global DNA */
  fNam = hprintf (h, "%s/dna.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->globalDna, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld target bases\n", bigArrayMax (bbG->globalDna)) ;

  /* export the complement of the global DNA */
  fNam = hprintf (h, "%s/dnaR.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->globalDnaR, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld complemented bases\n", bigArrayMax (bbG->globalDnaR)) ;

  /* memory map the coordinates */
  fNam = hprintf (h, "%s/coords.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->dnaCoords, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld coordinates\n", bigArrayMax (bbG->dnaCoords)) ;

  /* memory map the sequence identifiers (chromosome names) */
  /* transfer the dict to the seqids char array */
  int iMax = dictMax (bbG->dict) ;
  BigArray seqIds = bigArrayHandleCreate (256 * (iMax + 1), char, h) ;
  char *cp = bigArrayp (seqIds, 256 * (iMax + 1) - 1, char) ; /* make room */
  char buf[256] ;

  /* global title */
  cp = bigArrayp (seqIds, 0, char) ;
  memset (buf, 0, 256) ;
  char *signature = hprintf (h, "Sort Align version 1 step %d", tStep) ;
  strcpy (buf, signature) ;
  memcpy (cp, buf, 255) ;
  cp[255] = tStep ;
  /* list of names */
  for (int i = 1 ; i <= iMax ; i++)
    {
      cp += 256 ;
      memset (buf, 0, 256) ;
      strncpy (buf, dictName (bbG->dict, i), 255) ;
      memcpy (cp, buf, 256) ;
    }
  fNam = hprintf (h, "%s/ids.sortali", pp->indexName) ;
  bigArrayMapWrite (seqIds, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %d identifiers\n", dictMax (bbG->dict)) ;

  ac_free (h) ;
} /* storeTargetIndex */

/**************************************************************/

static BigArray GenomeAddSkips (const PP *pp, BigArray cws, BB *bb, int kk)
{
  long int iMax ; 
  long int jMax ; 
  long int i, j, k ;
  AC_HANDLE h = bb->h ;
  int maxRepeats = pp->maxTargetRepeats ;
    
  BigArray aa ;
  CW *up, *vp, *wp, *upMax ;
  unsigned int wordMax = 0xffffffff ;

  if (0)
    {
      int k = 0 ;
      vp = bigArrp (cws, 0, CW) ;
      for (j = 0 ; k < 6 && j < bigArrayMax(cws) ; j++, vp++)
	{
	  if (vp->nam == 44878)
	    {
	      k++ ;
	      fprintf (stderr, "==++ %d %d %u\n", vp->nam, vp->pos, vp->seed) ;
	    }
	}
    }
  /* remove highly repeated words and register number of repeats */
  if (1)
    {
      long int ks[21], cumul = 0 ;
      up = bigArrp (cws, 0, CW) ;
      vp = bigArrp (cws, 0, CW) ;

      iMax = bigArrayMax (cws) ;
      upMax = up + iMax ;
      memset (ks, 0, sizeof(ks)) ;
      for (i = 0, j = 0 ; i < iMax ; up++, i++)
	{
	  /* int tc = *dictName(pp->bbG.dict,up->nam >> 1) ; */
	  int m, n = 0, nI = 0, nX = 0 ;
	  wp = up ;

	  if (maxRepeats && pp->knownExons)
	    {
	      while (wp < upMax && wp->seed == up->seed)
		{
		  nI += ((wp->intron >> 31) & 0x1) ? 1 : 0 ;
		  nX += (((wp->intron >> 30) & 0x3) == 0x1) ? 1 : 0 ;
		  wp++ ;
		}
	    }
	  else
	    {
	      while (wp < upMax && wp->seed == up->seed)
		wp++ ;
	    }
	  n = wp - up ;
	  

	  if (!maxRepeats || n < maxRepeats)
	    {
	      for (wp = up, m = 0 ; m < n ; wp++, m++)
		{
		  if (((wp->intron >> 31) & 0x1) == 0x0)
		    wp->intron = n ;
		  if (vp < wp)
		    *vp = *wp ;
		  j++ ; vp++ ;
		}
	    }
	  else if (nI && nI < maxRepeats)  /* was nX < maxRepeats, but we stall on the RefEq XR */
	    {
	      for (wp = up, m = 0 ; m < n ; wp++, m++)
		{
		  nI = ((wp->intron >> 31) & 0x1) ? 1 : 0 ;
		  nX = (((wp->intron >> 30) & 0x3) == 0x1) ? 1 : 0 ;
		  nX = 0 ; /* do not enhance the exon seeds */
		  if (nI + nX)
		    {
		      if (((wp->intron >> 31) & 0x1) == 0x0)
			wp->intron = n ;
		      if (vp < wp)
			*vp = *wp ;
		      j++ ; vp++ ;
		    }
		}
	    }
	  up += n - 1 ; i += n - 1 ;
	  
	  if (n > 20) n = 20 ;
	  ks[n]++ ;
	}
      bigArrayMax (cws) = j ;

      if (0)
	{
	  AC_HANDLE h = ac_new_handle () ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, ".repeated_16_mers_in_target", FALSE, h) ;
	  aceOutf (ao, "#N\tWord\tInstances\tCumul");
	  for (i = 1 ; i <= 20 ; i++)
	    {
	      cumul += i * ks[i] ;
	      aceOutf (ao, "\n%ld\t%ld\t%ld\t%ld", i, ks[i], i * ks[i], cumul) ;
	    }
	  aceOut (ao, "\n") ;
	  ac_free (h) ;
	}
    }
  iMax = bigArrayMax (cws) ;
  if (! iMax) iMax = 1 ; /* insure non void */
  long int jMax0 = iMax + iMax/mstep1 + 1 ;
  aa = bigArrayHandleCreate (jMax0, CW, h) ;
  /* add skipping info */
  up = bigArrp (cws, 0, CW) ;
  vp = bigArrayp (aa, jMax0 - 1, CW) ; 
  vp = bigArrp (aa, 0, CW) ; jMax = 0 ;
  for (long int ii = 0 ; ii < iMax ; ii += mstep1)
    {
      vp->intron = ii + mstep4 < iMax ? (up + mstep4)->seed : wordMax ;
      vp->nam = ii + mstep3 < iMax ? (up + mstep3)->seed : wordMax ;
      vp->pos = ii + mstep2 < iMax ? (up + mstep2)->seed : wordMax ;
      vp->seed = ii + mstep1 < iMax ? (up + mstep1)->seed : wordMax ;

      k = iMax - ii ;
      if (k > mstep1)
	k = mstep1 ;
      vp++ ; 
      memcpy (vp, up, k * sizeof (CW)) ;
      vp += k ; up += k ;
      jMax += k + 1 ; /* k + 1_for_the_jumper */
      if (jMax > jMax0)
	messcrash ("add skipps error ") ;
    }
  bigArrayMax (aa) = jMax ;

  if (0)
    {
      int k = 0 ;
      vp = bigArrp (aa, 0, CW) ;
      for (j = 0 ; k < 6 &&j < jMax ; j++, vp++)
	{
	  if (vp->nam == 44878)
	    {
	      k++ ;
	      fprintf (stderr, "==== %d %d %u\n", vp->nam, vp->pos, vp->seed) ;
	    }
	}
    }
  return aa ;
} /* GenomeAddSkips */

/**************************************************************/
/* parse, code, sort the genome and create the index on disk
 * the human index takes around 18 GigaBytes
 */
static long int saTargetIndexCreateDo (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array tArray = pp->tArray ;
  BigArray cwsN[pp->nIndex] ;
  int nMax = arrayMax (tArray) ;
  TC *tc = 0 ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  long int nn = 0 ;
  int nTc = 0 ;
  BB *bbG = &(pp->bbG) ;
  
  memset (bbG, 0, sizeof (BB)) ;
  bbG->isGenome = TRUE ;
  t1 = clock () ;
  printf ("+++ %s: Parse the target files\n", timeBufShowNow (tBuf)) ;

  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  
  /* parse all targets into a single bbG->seqs array, with targetClass prefix in the sequence name */
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;
      if (tc->targetClass == 'I')
	continue ; /* we need to parse the genome before the introns */
      if (tc->targetClass == 'S')
	continue ; /* we need to parse the genome before the introns */
      nTc++ ;
    }

  for (int nn = 0, ntc = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;
      RC rc ;
      int step ;

      /* we need to parse the genome before the annotations */
      if (tc->targetClass == 'I')
	continue ; 
      if (tc->targetClass == 'S')
	continue ; 
      ntc++ ;
      
      memset (&rc, 0, sizeof (RC)) ;
      rc.fileName1 = tc->fileName ;
      rc.format = tc->format ;
      rc.run = nn + 1 ;
      saSequenceParse (pp, 0, tc, bbG, ntc == nTc ? 2 :1 ) ; /* 2 for last non-intron target */
      step = (bbG->length < 1<<20) ? 1 : 3 ;
      if (pp->tStep)
	step = pp->tStep ;
      if (step > pp->tStep)
	pp->tStep = step ;
      if (0)
	{
	  int mem = 0, mx = 0 ; /* megaBytes */
	  messAllocStatus (&mem) ;
	  messAllocMaxStatus (&mx) ;
	  fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
	}
    }
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;

      if (tc->targetClass == '*')
	saIntronParser (pp, tc) ;
    }
  
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;

      if (tc->targetClass == 'S')
	messcrash ("\n Unrecognized class S in traget configuration, please try sortalign --help") ;
    }
  
  t2 = clock () ;
  saCpuStatRegister ("1.Parse targets" , pp->agent, bbG->cpuStats, t1, t2, arrayMax (bbG->dnas) - 1) ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "===== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }

  t1 = clock () ;
  printf ("%s : extract the target seeds\n" , timeBufShowNow (tBuf)) ;
  saCodeSequenceSeeds (pp, bbG, pp->tStep, TRUE) ;
  if (pp->knownIntrons)
    saCodeIntronSeeds (pp, bbG) ;
  int NN = pp->nIndex ;
  for (int k = 0 ; k < NN ; k++)
    {
      long int n1 = bigArrayMax (bbG->cwsN[k]) ;
      cwsN[k] = bbG->cwsN[k] ;
      bbG->cwsN[k] = 0 ;
      nn += n1 ; 
    }
  for (int k = 0 ; k < NN ; k++)
    {
      long int n1 = bigArrayMax (cwsN[k]) ;
      if (1)
	{
	  int mem = 0, mx = 0 ; /* megaBytes */
	  messAllocStatus (&mem) ;
	  messAllocMaxStatus (&mx) ;
	  fprintf (stderr, "=== k=%d , %ld/%ld words %.1f %%,  Allocated %d Mb, max %d Mb\n", k, n1, nn, 100.0*n1/nn,  mem, mx) ;
	}
  }
  t2 = clock () ;
  saCpuStatRegister ("2.Extract target seeds" , pp->agent, bbG->cpuStats, t1, t2, bbG->nSeqs) ;
  t1 = clock () ;

  printf ("%s : sort the target seeds\n" , timeBufShowNow (tBuf)) ;
  for (int k = 0 ; k < NN ; k++)
    saSort (cwsN[k], 1) ; /* cwOrder */

  t2 = clock () ;
  saCpuStatRegister ("3.Sort seeds" , pp->agent, bbG->cpuStats, t1, t2, nn) ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  t1 = clock () ;

  printf ("%s : write the index to disk\n" , timeBufShowNow (tBuf)) ;
  bbG->cwsN = halloc (NN * sizeof(BigArray), bbG->h) ;
  for (int kk = 0 ; kk < NN ; kk++)
    {
      bbG->cwsN[kk] = GenomeAddSkips (pp, cwsN[kk], bbG, kk) ;
      bigArrayDestroy (cwsN[kk]) ;
    }

  storeTargetIndex (pp, pp->tStep) ;
  
  t2 = clock () ;
  saCpuStatRegister ("4.Write the index to disk" , pp->agent, bbG->cpuStats, t1, t2, nn) ;
  if (1)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }

  ac_free (h) ;
  return nn ;
} /* saTargetIndexCreateDo */

/**************************************************************/
/* The human genome index consumes around 18 Gigabytes of RAM */
void saTargetIndexCreate (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  /* check that input files were provided */

  pp->tMaxTargetRepeats = pp->maxTargetRepeats ;
  saTargetIndexCreateDo (pp) ;

  /* create short utility files in the IDX index directory */
  ACEOUT ao = aceOutCreate (filName (pp->indexName, "/seedLength", "w") , 0, 0, h) ;
  aceOutf (ao, "%d\t%d\t%d\n# SeedLength\ttStep\tmaxTargetRepeats\n"
	   , pp->seedLength, pp->tStep, pp->maxTargetRepeats
	   ) ;

  /* copy the actual config file used to create the index */
  if (pp->tConfigFileName)
    system (hprintf(h, "\\cp %s %s/tConfig\n", pp->tConfigFileName, pp->indexName)) ;
  else
    system (hprintf(h, "touch %s/tConfig\n", pp->indexName)) ;
  ac_free (h) ;
  return ;
}  /* saTargetIndexCreate */

/**************************************************************/
/**************************************************************/

static long int genomeParseBinary (const PP *pp, BB *bbG)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  BigArray seqIds = 0;
  DICT *dict = 0 ;
  long int ii, iMax, nn = 0 ;
  int NN = pp->nIndex ;
  
  clock_t t1, t2 ;

  t1 = clock () ;
  /* initialise the bbG (genome) block */
  bbG->h = ac_new_handle () ;
  bbG->cpuStats = arrayHandleCreate (128, CpuSTAT, bbG->h) ;
  bbG->dict = dict = dictHandleCreate (128, bbG->h) ;
  
  /* memory map the target DNA, seed Index, and Identifiers */
  BOOL READONLY = TRUE ;
  /* TRUE: memory mapping,
   * FALSE: read the data from disk into memory,
   *        100s slower in human
   */
  bbG->cwsN = halloc (NN * sizeof (BigArray), bbG->h) ;
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s.%d", pp->tFileBinaryCwsName, k) ;
      bbG->cwsN[k] = bigArrayMapRead (fNam, CW, READONLY, bbG->h) ; /* memory map the seed index */
      nn += bigArrayMax (bbG->cwsN[k]) ;
    }

  fNam = pp->tFileBinaryDnaName ;
  bbG->globalDna = bigArrayMapRead (fNam, unsigned char, READONLY, bbG->h) ; /* memory map the DNA */

  fNam = pp->tFileBinaryDnaRName ;
  bbG->globalDnaR = bigArrayMapRead (fNam, unsigned char, READONLY, bbG->h) ; /* memory map the reversed complemented DNA */

  fNam = pp->tFileBinaryCoordsName ;
  bbG->dnaCoords = bigArrayMapRead (fNam, unsigned int, READONLY, bbG->h) ; /* memory map the shared coordinates of the individual chromosomes in the globalDna/globalDnaR arrays */

  fNam = pp->tFileBinaryIdsName ;
  seqIds = bigArrayMapRead (fNam, char, FALSE, h) ; /* memory map the words */

  /* seqids is a char array, we need to transfer it to a dict */
  iMax = bigArrayMax (bbG->dnaCoords)/2 - 2 ;


  char *cp = bigArrp (seqIds, 0, char) ;
  int tStep = bigArr (seqIds, 255, char) ;
  /* check for common divisors */
  if (0 && /* no longer mandatory, since we select the best phasing */
      tStep > 1 && pp->iStep > 1)
    {
      int k ;
      for (int i = 2 ; i <= tStep ; i++)
	{
	  k = tStep/i ;
	  if (k * i == tStep) /* i divides tStep */
	    {
	      k = pp->iStep/i ;
	      if (k * i == pp->iStep) /* i divides iStep */
		messcrash ("\nThe target is indexed with step=%d, the requested read step is step=%d, these number are not relative primes, there will be systematic false negatives,\n please set the argument --istep %d (default 2) to a different value", tStep, pp->iStep, pp->iStep) ;
	    }
	}
    }

  bbG->step = tStep ;
  /* check the version */
  char *signature = hprintf (h, "Sort Align version 1 step %d", tStep) ;
  if (strcmp (cp, signature))
    messcrash ("\nCould not read the correct signature in the index files\n\texpected %s\n\treceived %s\nPlease destroy the index files %s.*.sortali and rerun sortalign --createIndex"
	       , signature
	       , cp
	       , fNam
	       ) ;
  /* create the dictionary of the chromosome identifiers */
  for (ii = 1 ; ii <= iMax ; ii++)
    {
      int n ;
      
      cp += 256 ;
      dictAdd (dict, cp, &n) ;
      if (n != ii)
	messcrash ("\nIndexing error, n=%d ii=%d : %s,  sorry, please rerun sortalign --createIndex\n", n, ii, cp) ;
    }

  /* create ancilary target dna arrays,
   * their memory is shared with the globalDna array
   * the coordinate of the seeds refer to the globalDna array
   */
  
  bbG->length = 0 ;
  bbG->nSeqs = iMax - 1 ;
  bbG->dnas = arrayHandleCreate (iMax + 1, Array, bbG->h) ;
  bbG->dnasR = arrayHandleCreate (iMax + 1, Array, bbG->h) ;
  /* entry zero is fake, because we index via a dictionary */
  array (bbG->dnas, 0, Array) = 0 ;
  array (bbG->dnasR, 0, Array) = 0 ;
  for (ii = 1 ; ii <= iMax ; ii++)
    {
      /* coords[iMax+1] is valid and initialised to bbG->length */
      unsigned int x1 = bigArr (bbG->dnaCoords, 2*ii, unsigned int) ;
      unsigned int x2 = bigArr (bbG->dnaCoords, 2*ii + 1, unsigned int) ;
      Array dna = arrayHandleCreate (8, unsigned char, bbG->h) ;
      Array dnaR = arrayHandleCreate (8, unsigned char, bbG->h) ;

      array (bbG->dnas, ii, Array) = dna ;
      array (bbG->dnasR, ii, Array) = dnaR ;
      /* manipulate the ancilary dna arrays */
      arrayLock (dna) ; /* protect dna->base. It must not be freed */
      messfree (dna->base) ;
      dna->base = bigArrp(bbG->globalDna, x1, char) ; 
      dna->max = dna->dim = x2 - x1 ;
      bbG->length += dna->max ;
      arrayLock (dnaR) ; /* protect dnaR->base */
      messfree (dnaR->base) ;
      dnaR->base = bigArrp(bbG->globalDnaR, x1, char) ; 
      dnaR->max = dnaR->dim = x2 - x1 ;
    }
  
#ifdef JUNK
  if (pp->wiggle)
    {
      bbG->wiggles = arrayHandleCreate (4*iMax, Array, bbG->h) ;
      for (ii = 1 ; ii <= iMax ; ii++)
	{
	  Array dna = array (bbG->dnas, ii, Array) ;
	  int ln = arrayMax (dna)/WIGGLE_STEP + 1 ;
	  array (bbG->wiggles, 4*ii, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 1, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 2, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 3, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	}
    }
#endif
  
  /*  Get thread CPU time at end */
  t2 = clock () ;
  saCpuStatRegister ("1.memMapTargets" , pp->agent, bbG->cpuStats, t1, t2, nn) ; 
  ac_free (h) ;
  return nn ; 
} /* genomeParseBinary */

/**************************************************************/

void saTargetIndexGenomeParser (const void *vp)
{
  const PP *pp = vp ;
  BB bbG = pp->bbG ;
  char tBuf[25] ;
  
  clock_t t2 =0,       t1 = clock () ;
  printf ("+++ %s: Start genome parser\n", timeBufShowNow (tBuf)) ;
  memset (&bbG, 0, sizeof (BB)) ;
  long int nn = genomeParseBinary (pp, &bbG) ;
  t2 = clock () ;


  saCpuStatRegister ("1.GParserDone" , pp->agent, bbG.cpuStats, t1, t2, nn) ;
  channelPut (pp->gmChan, &bbG, BB) ;
  channelClose (pp->gmChan) ;
  printf ("--- %s: Stop binary genome parser\n", timeBufShowNow (tBuf)) ;

  return ;
} /* genomeParser */

/**************************************************************/
/**************************************************************/
/**************************************************************/
