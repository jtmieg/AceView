/*
 * sa.config.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module parses the various configuration files
 * which define the sequence files to be aligned, the target genome
 * and its annotations (gens, exons, introns ...)
 */

#include "sa.h"


/*********************************************************************/

void saConfigIsIndexAccessible (PP *pp)
{
  if (! pp->indexName)
    saUsage ("Missing parameter --index", 0, 0) ;
  else
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *cp = pp->indexName ;
      char *cq ;
      if (pp->createIndex)
	{
	  cq = hprintf (h, "mkdir %s ; \\rm %s/* ; echo test > %s/sortalign.testFile", cp, cp, cp) ;
	  system (cq) ;
	}
      cq = hprintf (h, "%s/sortalign.testFile", pp->indexName) ;
      ACEIN ai = aceInCreate (cq, FALSE, h) ;
      if (! ai)
	saUsage (hprintf (h, "Cannot create and access the index directory: %s", pp->indexName), 0, 0) ;
      ac_free (h) ;
    }
} /* saConfigIsIndexAccessible */

/*********************************************************************/

void saConfigIsOutDirWritable (PP *pp)
{
  if (pp->outFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      
      char *cp = strrchr (pp->outFileName, '/') ;
      if (cp)
	{
	  *cp = 0 ;
	  char *cmd = hprintf (h, "mkdir -p %s", pp->outFileName) ;
	  system (cmd) ;
	  *cp = '/' ;
	}
      
      ACEOUT ao = aceOutCreate (pp->outFileName, ".test", 0, h) ;
      if (ao)
	{
	  int n = randint (), n2 = 0 ;
	  char * nam = strnew (aceOutFileName (ao), h) ;
	  
	  aceOutf (ao, "%d\n", n ) ;
	  ac_free (ao) ;
	  
	  ACEIN ai = aceInCreate (nam, 0, h) ;
	  if (ai)
	    {
	      if (aceInCard (ai))
		aceInInt (ai, &n2) ;
	      if (n != n2)
		saUsage (hprintf (h, "Unable to open a test file -o %s.test", pp->outFileName), 0, 0) ;
	      ac_free (ai) ;
	    }
	  char *cmd = hprintf (h, "rm %s", nam) ;
	  system (cmd) ;
	}
      ac_free (h) ;
    }
  
  return ;
} /* saConfigOutFiles */

/*********************************************************************/

/* check the existence of the sequence files to be aligned
 * identify their absolute file names
 * associate each file to a run name and to its optional parameters
 */
Array saConfigGetRuns (PP *pp, Array runStats)
{
  Array rcs = arrayHandleCreate (64, RC, pp->h) ;
  RC *rc = 0 ;
  int nn = 0 ;
  
  if (pp->inFileName)
    {   /* Split the individual file names, they are coma separated 
	 * Make the names absolute, and load them in the file-parser fp channel
	 *  all names must be loaded at once, otherwise we may have a deadlock
	 *  because the agents needing genome parser completion are not yet launched 
	 */
      char *buf = strnew (pp->inFileName, pp->h) ;
      char *cp, *cq, *cr, *cs ;
      int run = 0 ;
      char *filName2 = 0 ;
      if (! pp->runName)
	pp->runName = "runX" ;
      dictAdd (pp->runDict, pp->runName, &run) ;
      cp = buf ; *(cp + strlen (cp)) = ',' ;
      while ((cq = strchr (cp, ',')))
	{
	  cs = strchr (cp, ':') ;
	  if (cs)
	    { /* switch runName */
	      *cs = 0 ;
	      dictAdd (pp->runDict, cp, &run) ;
	      cp = cs + 1 ;
	    }
	  *cq = 0 ;
	  if (! *cp)
	    continue ;
	  cs = strchr (cp, '+') ;
	  if (cs)
	    {
	      *cs = 0 ;
	      cr = filName (cs+1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open input file %s\n", cs+1), 0, 0) ;
	      filName2 = strnew (cr, pp->h) ;
	    }
	  

	  rc = arrayp (rcs, nn++, RC) ; 
	  rc->fileName1 = strnew (cp, pp->h) ;
	  rc->fileName2 = filName2 ;
	  rc->run = run ;

	  rc->format = FASTA ; /* default */
	  int ln = strlen (rc->fileName1) ;
	  if (strstr (rc->fileName1, ".fasta.gz") == rc->fileName1 + ln - 9) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastq.gz") == rc->fileName1 + ln - 9) rc->format = FASTQ ;
	  else if (strstr (rc->fileName1, ".fasta") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastq") == rc->fileName1 + ln - 6) rc->format = FASTQ ;
	  else if (strstr (rc->fileName1, ".fna.gz") == rc->fileName1 + ln - 7) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fa.gz") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastc.gz") == rc->fileName1 + ln - 9) rc->format = FASTC ;
	  else if (strstr (rc->fileName1, ".fastc") == rc->fileName1 + ln - 6) rc->format = FASTC ;
	  else if (strstr (rc->fileName1, ".fna") == rc->fileName1 + ln - 4) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fa") == rc->fileName1 + ln - 3) rc->format = FASTA ;
	  else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	  
	  /* user can override the defaults */
	  if (pp->raw) rc->format = RAW ;
	  if (pp->fasta) rc->format = FASTA ;
	  if (pp->fastq) rc->format = FASTQ ;
	  if (pp->fastc) rc->format = FASTC ;
	  if (pp->sra) rc->format = SRA ;

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open input file %s\n", cp), 0, 0) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    {
	      rc->pairedEnd = TRUE ;
	      rs->nFiles++ ;
	    }
	  if (rc->format == FASTC)
	      rc->pairedEnd = TRUE ;
	  
	  cp = cq + 1 ;
	}
    }
  else if (pp->inConfigFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *fDict = dictHandleCreate (64, h) ;
      ACEIN ai = aceInCreate (pp->inConfigFileName, 0, h) ;
      int nRuns = 0, run = 0 ;
      int line = 0 ;

      if (! ai)
	saUsage (hprintf (h, "Cannot find file -I %s specified on the command line", pp->inConfigFileName), 0, 0) ;
      while (aceInCard (ai))
	{
	  char *cq, *cr, *cp = aceInWord (ai) ;
	  
	  line++ ;
	  if (! cp || ! *cp || *cp == '#')
	    continue ;
	  /* file names */
	  rc = arrayp (rcs, nn++, RC) ; 
	  nRuns++ ;
	      /* file pairs */
	  cq = strchr (cp, '+') ;
	  if (cq) 
	    {
	      *cq++ = 0 ;
	      if (*cq)
		{
		  cr = filName (cq, 0, "r") ;
		  if (! cr)
		    saUsage (hprintf (h, "Cannot open input file %s\n", cq), 0, 0) ;
		  rc->pairedEnd = TRUE ;
		  rc->fileName2 = strnew (cr, pp->h) ;
		}
	    }
	  rc->fileName1 = strnew (cp, pp->h) ;
	  if (! dictAdd (fDict, cp, 0))
	    saUsage (hprintf (h, "Duplicate file name %s\n at line %d of file -T %s"
		       , cr
		       , line
		       , pp->inConfigFileName
			      ), 0, 0) ;

	  /* run name */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (cp && *cp && *cp != '#')
	    dictAdd (pp->runDict, cp, &run) ;
	  else
	    dictAdd (pp->runDict, hprintf (h, "r.%d", nRuns), &run) ;
	  rc->run = run ;
	  rc->RNA = TRUE ; /* default */

	  /* options */
	  cp = aceInWord (ai) ;
	  while (cp)
	    {
	      if (*cp == '#')
		break ;
	      cq = strchr (cp, ',') ;
	      if (cq)
		*cq++ = 0 ;
	      if (! strcasecmp (cp, "fasta")) rc->format = FASTA ;
	      if (! strcasecmp (cp, "fastq")) rc->format = FASTQ ;
	      if (! strcasecmp (cp, "fastc")) rc->format = FASTC ;
	      if (! strcasecmp (cp, "raw")) rc->format = RAW ;
	      if (! strcasecmp (cp, "SRA")) rc->format = SRA ;

	      if (! strcasecmp (cp, "rna")) rc->RNA = TRUE ;
	      if (! strcasecmp (cp, "dna")) rc->RNA = FALSE ;
	      
	      cp = cq ;
	    }

	  if (! rc->format)
	    {
	      int ln = strlen (rc->fileName1) ;
	      rc->format = FASTA ; /* default */
	      if (strstr (rc->fileName1, ".fasta.gz") == rc->fileName1 + ln - 9) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastq.gz") == rc->fileName1 + ln - 9) rc->format = FASTQ ;
	      else if (strstr (rc->fileName1, ".fasta") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastq") == rc->fileName1 + ln - 6) rc->format = FASTQ ;
	      else if (strstr (rc->fileName1, ".fna.gz") == rc->fileName1 + ln - 7) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fa.gz") == rc->fileName1 + ln - 6) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastc.gz") == rc->fileName1 + ln - 9) rc->format = FASTC ;
	      else if (strstr (rc->fileName1, ".fastc") == rc->fileName1 + ln - 6) rc->format = FASTC ;
	      else if (strstr (rc->fileName1, ".fna") == rc->fileName1 + ln - 4) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fa") == rc->fileName1 + ln - 3) rc->format = FASTA ;
	      else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	    }

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		saUsage (hprintf (pp->h, "Cannot open sequence file %s listed in configuration file -I %s \tDid you mean -i %s\n"
				  , cp, pp->inConfigFileName, pp->inConfigFileName
				  ), 0, 0) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    rs->nFiles++ ;
	  if (rc->fileName2 || rc->format == FASTC)
	      rc->pairedEnd = TRUE ;
	} 
      
      ac_free (h) ;
    }
  pp->nFiles = arrayMax (rcs) ;
  printf ("Found %d sequence file%s\n", arrayMax (rcs), arrayMax (rcs) > 1 ? "s" : "") ;
  return rcs ;
} /* saConfigGetRuns */

/*********************************************************************/
/*********************************************************************/

/* check the completeness of the target index directory
 */
int saConfigCheckTargetIndex (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = TRUE ;
  int k = 0 ;
  char *cp ;
  int NN = 0 ;
  while (1)
    {
      char *fNam = hprintf (h, "/cws.sortali.%d", k) ;
      k++ ;
      cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  NN++;
	  if (NN == 1)
	    {
	      char *cq = cp + strlen(cp) - 2 ;
	      *cq = 0 ;
	      pp->tFileBinaryCwsName = strnew (cp, pp->h) ;
	    }
	}
      else
	break ;
    }
  if (NN != 1 && NN != 2 && NN != 4 && NN != 8 && NN != 16 && NN != 32 && NN != 64)
    messcrash ("The number of indexes NN=%d must be  a power of 2, say 1 2 4 8 ...", NN) ;
  
  ACEIN ai = aceInCreate (filName (pp->indexName, "/seedLength", "r") , 0, h) ;
  if (ai)
    {
      int x = 0, y = 0, z = 0 ;
      if (aceInCard (ai))
	{
	  if (aceInInt (ai, &x))
	    pp->seedLength = x ;
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &y))
	    pp->tStep = y ;
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &z))
	    pp->tMaxTargetRepeats = z ;
	}
      ac_free (ai) ;
    }
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/dna.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryDnaName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/dnaR.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryDnaRName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/ids.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryIdsName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  cp = filName (pp->indexName, "/coords.sortali", "rb") ;
  if (cp)
    pp->tFileBinaryCoordsName = strnew (cp, pp->h) ;
  else
    ok = FALSE ;
  if (!ok)
    saUsage ("Some of the target binary index files are missing, please run sortalign --createIndex <indexName> ", 0, 0) ;
  
  if (pp->iStep > pp->tStep)      /* impose phasing */
    pp->iStep -= (pp->iStep % pp->tStep) ;
      
  ac_free (h) ;
  return NN ;
} /* saConfigCheckTargetIndex */

/*********************************************************************/
/*********************************************************************/
