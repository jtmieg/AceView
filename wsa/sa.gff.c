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
/* #define ARRAY_CHECK */

#include "sa.h"

/**************************************************************/
/**************************************************************/

static int exonMrnaOrder (const void *va, const void *vb)
{
  const EXONINTRON *up = va ;
  const EXONINTRON *vp = vb ;
  int n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->mrna - vp->mrna ; if (n) return n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  return 0 ;
} /* exonMrnaOrder */

/**************************************************************/

static int exonA1Order (const void *va, const void *vb)
{
  const EXONINTRON *up = va ;
  const EXONINTRON *vp = vb ;
  int n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->mrna - vp->mrna ; if (n) return n ;
  return 0 ;
} /* exonMrnaOrder */

/**************************************************************/

static int a1GeneOrder (const void *va, const void *vb)
{
  const GENE *up = va ;
  const GENE *vp = vb ;
  int n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;

  return 0 ;
} /* a1GeneOrder */

/**************************************************************/
#ifdef JUNK
static int geneA1Order (const void *va, const void *vb)
{
  const GENE *up = va ;
  const GENE *vp = vb ;
  int n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* geneA1Order */
#endif
/**************************************************************/
/**************************************************************/

long int saGffParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nnE = 0, nnI = 0, nnEX = 0 ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom = 0, gene = 0, mrna = 0 ;
  int a1, a2 ;
  int line = 0 ;
  BB *bbG = &(pp->bbG) ;
  const int twoMb = (0x1 << 21) ;
  BOOL isGff = FALSE ;
  BOOL isGtf = FALSE ;
  BOOL isIntrons = FALSE ;
  DICT *chromDict = bbG->dict ;
  DICT *geneDict = pp->geneDict ;
  DICT *mrnaDict = pp->mrnaDict ;

  if (! geneDict) geneDict = pp->geneDict = dictHandleCreate (1 << 15, pp->h) ;
  if (! mrnaDict) mrnaDict = pp->mrnaDict = dictHandleCreate (1 << 15, pp->h) ;

  if (! pp->knownIntrons)   pp->knownIntrons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;
  if (! pp->knownExons) pp->knownExons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;

  if (! pp->geneBoxes) pp->geneBoxes = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;
  if (! pp->geneExons) pp->geneExons = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;
  
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

      if (isIntrons) fprintf (stderr, "Parsing .intron file %s\n", cp) ;
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

  if (isIntrons)
    return saIntronParser (pp, tc) ;
  
  if (! pp->knownIntrons)   pp->knownIntrons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;
  if (! pp->knownExons) pp->knownExons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;

  nnE = bigArrayMax (pp->knownExons) ;
  nnI = bigArrayMax (pp->knownIntrons) ;
    
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


      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* method, drop it */
      
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* type */
      type = 0 ;
      if (! strcasecmp (cp, "exon")) type = 1 ;
      if (! type) continue ;
  

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;

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
      
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* probably . (some kind of score ?) drop it */
      if (! cp)
	continue ;
      
      aceInStep (ai, '\t') ;
      char cutter = 0 ;
      cp = aceInWordCut (ai, "\t", &cutter) ;  /* the metadata */
      if (! cp)
	continue ;
      mrna = 0 ;
      if (type == 1) /* identify the transcript */
	{ /* expect : transcript_id "NM_001364814.1"; */
	  char *cq = strstr (cp, "transcript_id") ;
	  if (cq)
	    {
	      AC_HANDLE h1 = ac_new_handle () ;
	      ACEIN aj = aceInCreateFromText (cq + 13, 0, h1) ;
	      aceInCard (aj) ;
	      char *cr = aceInWord (aj) ;
	      if (cr)
		dictAdd (mrnaDict, cr, &mrna) ;
	      ac_free (h1) ;
	    }
	}
      gene = 0 ;
      if (type == 1) /* identify the gene */
	{ /* expect : transcript_id "NM_001364814.1"; */
	  char *cq = strstr (cp, "gene_id") ;
	  if (cq)
	    {
	      AC_HANDLE h1 = ac_new_handle () ;
	      ACEIN aj = aceInCreateFromText (cq + 7, 0, h1) ;
	      aceInCard (aj) ;
	      char *cr = aceInWord (aj) ;
	      if (cr)
		dictAdd (geneDict, cr, &gene) ;
	      ac_free (h1) ;
	    }
	}
      
      if (type == 1)
	{
	  EXONINTRON *up = bigArrayp (pp->knownExons, nnE++, EXONINTRON) ;
	  up->gene = gene ;
	  up->mrna = mrna ;
	  up->chrom = chrom ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}

    }

  if (nnE) /* create the corresponding introns */
    {
      EXONINTRON *up, *vp ;
      long int ii, jj ;

      bigArraySort (pp->knownExons, exonMrnaOrder) ;
      for (up = bigArrp (pp->knownExons, 0, EXONINTRON), ii = 0 ; ii < nnE ; ii++, up++)
	if (up->mrna)
	  {
	    for (vp = up + 1, jj = ii + 1 ; jj < nnE && vp->mrna == up->mrna ; vp++, jj++) 
	      {
		int chrom = up->chrom ;
		int a1 = up->a1 < up->a2 ? vp[-1].a2 + 1 : vp[0].a2 - 1 ;
		int a2 = up->a1 < up->a2 ? vp[0].a1 - 1 : vp[-1].a1 + 1 ;
		int da = a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1 ;
		if (0)
		  fprintf (stderr, "Intron %ld : %d %d\n", nnI, a1, a2) ;
		if (da >= twoMb) /* our format only allows 21 bits for the intron length */ 
		  continue ;

		EXONINTRON *wp = bigArrayp (pp->knownIntrons, nnI++, EXONINTRON) ;
		wp->mrna = 0 ;
		wp->gene = up->gene ;
		if (a1 < a2)
		  {
		    wp->chrom = chrom << 1 ; 
		    wp->a1 = a1 ;
		    wp->a2 = a2 ;
		  }
		else
		  {
		    wp->chrom = (chrom << 1) | 0x1 ;
		    wp->a1 = a2 ;
		    wp->a2 = a1 ;
		  }
		vp->mrna = 0 ;
	      }
	    up->mrna = 0 ;
	  }
      bigArraySort (pp->knownIntrons, exonA1Order) ;
      bigArrayCompress (pp->knownIntrons) ;
      nnI = bigArrayMax (pp->knownIntrons) ;
    }


  if (nnE) /* construct the gene expression shadows */
    {
      EXONINTRON *up ;
      GENE *gb, *ge ;
      long int ii ;
      
      fprintf (stderr, "max knownExons = %ld\n", nnE) ;
      bigArraySort (pp->knownExons, exonMrnaOrder) ;
      fprintf (stderr, "max knownExons = %ld\n", nnE) ;
      bigArrayCompress (pp->knownExons) ;
      nnE = bigArrayMax (pp->knownExons) ;
      fprintf (stderr, "max knownExons = %ld\n", nnE) ;
      
      for (up = bigArrp (pp->knownExons, 0, EXONINTRON), ii = 0 ; ii < nnE ; ii++, up++)
	{
	  gene = up->gene ;
	  chrom = up->chrom ;
	  a1 = up->a1 ;
	  a2 = up->a2 ;
	  if (a1 < a2)
	    {
	      chrom = chrom << 1 ; 
	    }
	  else
	    {
	      chrom = (chrom << 1) | 0x1 ;
	      int a0 = a1 ; a1 = a2 ; a2 = a0 ;
	    }

	  if  (gene)
	    {
	      Array geneBoxes = array (pp->geneBoxes, chrom, Array) ;
	      Array geneExons = array (pp->geneExons, chrom, Array) ;
	      if (! geneBoxes) geneBoxes = array (pp->geneBoxes, chrom, Array) = arrayHandleCreate (1000, GENE, pp->h) ;
	      if (! geneExons) geneExons = array (pp->geneExons, chrom, Array) = arrayHandleCreate (1000, GENE, pp->h) ;
	      
	      /* register a geneBox */
	      gb = arrayp (geneBoxes, gene, GENE) ;
	      gb->gene = gene ;
	      gb->a1 = (! gb->a1 || gb->a1 > a1) ? a1 : gb->a1 ;
	      gb->a2 = (! gb->a2 || gb->a2 < a2) ? a2 : gb->a2 ;
	      
	      /* register a geneExon */
	      ge = arrayp (geneExons, nnEX++, GENE) ;
	      ge->gene = gene ;
	      ge->a1 = a1 ;
	      ge->a2 = a2 ;
	    }
	}
    }

  /* sort the gene boxes */
  if (nnE)
    {
      for (int chrom = 0 ; chrom < arrayMax (pp->geneBoxes) ; chrom++)
	{
	  Array geneBoxes = array (pp->geneBoxes, chrom, Array) ;
	  if (geneBoxes)
	    arraySort (geneBoxes, a1GeneOrder) ;
	}
      for (int chrom = 0 ; chrom < arrayMax (pp->geneExons) ; chrom++)
	{
	  Array geneExons = array (pp->geneExons, chrom, Array) ;
	  if (geneExons)
	    arraySort (geneExons, a1GeneOrder) ;
	}
    }
  
  if (nnE) /* fuse exons by projection on the top strand of the genome */
    {
      EXONINTRON *up, *vp ;
      long int ii, jj ;

      for (up = bigArrp (pp->knownExons, 0, EXONINTRON), ii = 0 ; ii < nnE ; ii++, up++)
	if (up->a1 > up->a2)
	  { int a0 = up->a1 ; up->a1 = up->a2 ; up->a2 = a0 ; }
      bigArraySort (pp->knownExons, exonA1Order) ;
      for (up = vp = bigArrp (pp->knownExons, 0, EXONINTRON), ii = jj = 0 ; ii < nnE ; ii++, up++)
	{
	  if (up->chrom == vp->chrom && up->a1 <= vp->a2 + 1)
	    {
	      if (up->a2 > vp->a2)
		vp->a2 = up->a2 ;
	      continue ;
	    }
	  if (jj < ii) *vp = *up ;
	  jj++ ; vp++ ;
	}
      bigArrayMax (pp->knownExons) = jj ;
      if (0)
	for (up = bigArrp (pp->knownExons, 0, EXONINTRON), ii = 0 ; ii < nnE ; ii++, up++)
	  fprintf (stderr, "Exon %ld : %d %d\n", ii, up->a1, up->a2) ;
    }

  fprintf (stderr, "+++++++ Found %ld exons %ld introns in file %s\n", nnE, nnI, aceInFileName (ai)) ;
  ac_free (h) ;

  return nnE + nnI ;
} /* gffParser */

/**************************************************************/

long int saIntronParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nn = 0 ;
  BB *bbG = &(pp->bbG) ;
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
  pp->knownIntrons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;
  EXONINTRON *up ;

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
      up = bigArrayp (pp->knownIntrons, nn++, EXONINTRON) ;
      if (a1 < a2)
	{
	  up->chrom = chrom << 1 ; 
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}
      else
	{
	  up->chrom = (chrom << 1) | 0x1 ;
	  up->a1 = a2 ;
	  up->a2 = a1 ;
	}
	
    }
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, aceInFileName (ai)) ;
  ac_free (h) ;

  return nn ;
} /* saIntronParser */

/**************************************************************/
/**************************************************************/
