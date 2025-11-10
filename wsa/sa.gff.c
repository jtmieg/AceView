/*
 * sa.gff.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * the parsing the gff annotations 

 */

#include "sa.h"

/**************************************************************/
/**************************************************************/

static int exonOrder (const void *va, const void *vb)
{
  const CW *up = va ;
  const CW *vp = vb ;
  int n ;
  n = (up->seed > vp->seed) - (up->seed < vp->seed) ; if (n) return n ;
  n = up->nam - vp->nam ; if (n) return n ;
  n = (up->pos > vp->pos) - (up->pos < vp->pos) ; if (n) return n ;
  return 0 ;
} /* exonOrder */

/**************************************************************/

long int saGffParser (PP *pp, BB *bbG, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nnE = 0, nnI = 0 ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom= 0 ;
  int a1, a2 ;
  int line = 0 ;
  const int twoMb = (0x1 << 21) ;
  BOOL isGff = FALSE ;
  BOOL isGtf = FALSE ;
  BOOL isIntrons = FALSE ;
  
  if (1)
    {
      const char *cp = tc->fileName ;
      int i = strlen (cp) ;
      if (i >= 8 && !strcmp (cp + i - 8, ".introns")) isIntrons = TRUE ;
      else if (i >= 11 && !strcmp (cp + i - 11, ".introns.gz")) isIntrons = TRUE ;
      else if (i >= 4 && !strcmp (cp + i - 4, ".gff")) isGff = TRUE ;
      else if (i >= 5 && !strcmp (cp + i - 5, ".gff3")) isGff = TRUE ;
      else if (i >= 4 && !strcmp (cp + i - 4, ".gtf")) isGtf = TRUE ;
      else if (i >= 7 && !strcmp (cp + i - 7, ".gff.gz")) isGff = TRUE ;
      else if (i >= 8 && !strcmp (cp + i - 8, ".gff3.gz")) isGff = TRUE ;
      else if (i >= 7 && !strcmp (cp + i - 7, ".gtf.gz")) isGtf = TRUE ;

      if (isIntrons) fprintf (stderr, "Parsing gff file %s\n", cp) ;
      else if (isGff) fprintf (stderr, "Parsing gff file %s\n", cp) ;
      else if (isGtf) fprintf (stderr, "Parsing gtf file %s\n", cp) ;
      else messcrash ("\n\tOption I in -T tConfig_file_name expects a gff or gtf file\n\tCannot parse file %s found in tConfig file %s\n"
		      , cp, pp->tConfigFileName ) ;
    }
      
  if (pp->seedLength < 16)
    messcrash ("\nSorry, to study the introns defined in file %s,\n the seed length must be at least 16, not %d\n"
	       , aceInFileName (ai)
	       , pp->seedLength
	       ) ;
  if (! pp->intronSeeds)   pp->intronSeeds = bigArrayHandleCreate (100000, CW, pp->h) ;
  if (! pp->exonSeeds) pp->exonSeeds = bigArrayHandleCreate (100000, CW, pp->h) ;

  nnE = bigArrayMax (pp->exonSeeds) ;
  nnI = bigArrayMax (pp->intronSeeds) ;
    
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      int type ;
      char *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || *cp == '#' || *cp == '%' || *cp == '/')
	continue ;
      chrom = 0 ; a1 = a2 = 0 ;
      
      if (! dictFind (bbG->dict, hprintf (h, "G.%s", cp), &chrom))
	messcrash ("\nUnknown chromosome name %s, not matching the G targets, at line %d of file %s\n"
		   , cp
		   , line
		   , aceInFileName (ai)
		   ) ;

      if (! isIntrons)
	{
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* method, drop it */

	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* type */
	  type = 0 ;
	  if (! strcasecmp (cp, "exon")) type = 1 ;
	  else if (! strcasecmp (cp, "intron")) type = 2 ;
	  if (! type) continue ;
	}
      else
	type = 2 ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;

      if (! isIntrons)
	{
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* some kind of score or a dot, drop it */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* strand */
	  if (! cp)
	    continue ;
	  else if (*cp == '+')
	    ;
	  else if (*cp == '-')
	    {
	      int a0 = a1 ; a1 = a2 ; a2 = a0 ;
	    }
	  else
	    continue ;
	}

      if (type == 1)
	{
	  CW *up = bigArrayp (pp->exonSeeds, nnE++, CW) ;
	  up->nam = chrom ;
	  up->pos = a1 < a2 ? a1 : a2 ;
	  up->intron = a1 < a2 ? a2 : a1 ;
	}

      if (type == 2)
	{
	  int da = (a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1) ;
	  if (da >= twoMb) /* our format only allows 21 bits for the intron length */ 
	    continue ;
	  CW *up = bigArrayp (pp->intronSeeds, nnI++, CW) ;
	  up->nam = chrom ;
	  up->pos = a1 ;
	  up->intron = a2 ;
	}
    }
  
  if (nnE) /* fuse exons by projection on the top strand of the genome */
    {
      CW *up, *vp ;
      long int ii, jj ;
      bigArraySort (pp->exonSeeds, exonOrder) ;

      for (up = vp = bigArrp (pp->exonSeeds, 0, CW), ii = jj = 0 ; ii < nnE ; ii++, up++)
	{
	  if (up->nam == vp->nam && up->pos <= vp->intron)
	    {
	      if (up->intron > vp->intron)
		vp->intron = up->intron ;
	      continue ;
	    }
	  if (jj < ii) *vp = *up ;
	  jj++ ; vp++ ;
	}
      bigArrayMax (pp->exonSeeds) = jj ; 
    }

  fprintf (stderr, "+++++++ Found %ld exons %ld introns in file %s\n", nnE, nnI, aceInFileName (ai)) ;
  ac_free (h) ;

  return nnE + nnI ;
} /* gffParser */

/**************************************************************/

long int saIntronParser (PP *pp, BB *bbG, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nn = 0 ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom= 0 ;
  int a1, a2, da ;
  int line = 0 ;
  const int twoMb = (0x1 << 21) ;

  if (pp->seedLength < 16)
    messcrash ("\nSorry, to study the introns defined in file %s,\n the seed length must be at least 16, not %d\n"
	       , aceInFileName (ai)
	       , pp->seedLength
	       ) ;
  pp->intronSeeds = bigArrayHandleCreate (100000, CW, pp->h) ;
  CW *up ;

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || *cp == '#')
	continue ;
      chrom = 0 ; a1 = a2 = 0 ;

      if (! dictFind (bbG->dict, hprintf (h, "G.%s", cp), &chrom))
	messcrash ("\nUnknown chromosome name %s, not matching the G targets, at line %d of file %s\n"
		   , cp
		   , line
		   , aceInFileName (ai)
		   ) ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;
      da = (a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1) ;
      if (da >= twoMb) /* our format only allows 21 bits for the intron length */ 
	continue ;
      up = bigArrayp (pp->intronSeeds, nn++, CW) ;
      up->nam = chrom ;
      up->pos = a1 ;
      up->intron = a2 ;
    }
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, aceInFileName (ai)) ;
  ac_free (h) ;

  return nn ;
} /* saIntronParser */

