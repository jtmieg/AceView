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
} /* exonOrder */

/**************************************************************/

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
} /* exonOrder */

/**************************************************************/
/**************************************************************/

long int saGffParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nnE = 0, nnI = 0 ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom= 0 ;
  int a1, a2 ;
  int line = 0 ;
  BB *bbG = &(pp->bbG) ;
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

      if (isIntrons) fprintf (stderr, "Parsing .intron file %s\n", cp) ;
      else if (isGff) fprintf (stderr, "Parsing gff file %s\n", cp) ;
      else if (isGtf) fprintf (stderr, "Parsing gtf file %s\n", cp) ;
      else messcrash ("\n\tOption I in -T tConfig_file_name expects a gff or gtf file\n\tCannot parse file %s found in tConfig file %s\n"
		      , cp, pp->tConfigFileName ) ;
    }

  if (isIntrons)
    return saIntronParser (pp, tc) ;
  
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


      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* method, drop it */
      
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* type */
      type = 0 ;
      if (! strcasecmp (cp, "exon")) type = 1 ;
      else if (! strcasecmp (cp, "intron")) type = 2 ;
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
      if (a1 < a2)
	{
	  up->nam = chrom << 1 ; 
	  up->pos = a1 ;
	  up->intron = a2 ;
	}
      else
	{
	  up->nam = (chrom << 1) | 0x1 ;
	  up->pos = a2 ;
	  up->intron = a1 ;
	}
	
    }
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, aceInFileName (ai)) ;
  ac_free (h) ;

  return nn ;
} /* saIntronParser */

/**************************************************************/

static int saSpongeParserDo (PP *pp, const char *filNam)
{
  AC_HANDLE h = ac_new_handle () ;
  int nnnE = 0, nnnB = 0 ;
  ACEIN ai = aceInCreate (filNam, 0, h) ;
  int mrna, gene, exon, chrom ;
  int oldMrna, oldGene, oldExon, oldChrom ;
  int a1, a2, da, b2 ;
  int line = 0 ;
  int nnES, nnI, nnE ;
  CW *up ; ;
  GENE *gb, *ge ;
  BB *bbG = &(pp->bbG) ;
  Array geneBoxes = 0, geneExons = 0 ;
  DICT *chromDict = bbG->dict ;
  DICT *geneDict = pp->geneDict ;
  DICT *mrnaDict = pp->mrnaDict ;
  const int twoMb = (0x1 << 21) ;
  
  if (! geneDict) geneDict = pp->geneDict = dictHandleCreate (1 << 15, pp->h) ;
  if (! mrnaDict) mrnaDict = pp->mrnaDict = dictHandleCreate (1 << 15, pp->h) ;

  if (! pp->intronSeeds)   pp->intronSeeds = bigArrayHandleCreate (100000, CW, pp->h) ;
  if (! pp->exonSeeds) pp->exonSeeds = bigArrayHandleCreate (100000, CW, pp->h) ;

  if (! pp->geneBoxes) pp->geneBoxes = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;
  if (! pp->geneExons) pp->geneExons = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;

  
  nnES = bigArrayMax (pp->exonSeeds) ;
  nnI = bigArrayMax (pp->intronSeeds) ;
  nnE = bigArrayMax (pp->geneExons) ;

  
  
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || *cp == '#')
	continue ;
      mrna = chrom = gene = exon = 0 ; a1 = a2 = 0 ;

      dictAdd (mrnaDict, cp, &mrna) ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &exon) ;

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;
      
      if (! dictFind (bbG->dict, hprintf (h, "G.%s", cp), &chrom) &&
	  ! dictFind (bbG->dict, hprintf (h, "M.%s", cp), &chrom) &&
	  ! dictFind (bbG->dict, hprintf (h, "R.%s", cp), &chrom) &&
	  ! dictFind (bbG->dict, hprintf (h, "C.%s", cp), &chrom) 
	  )
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
      chrom <<= 1 ;
      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; chrom |= 0x1 ;}
      da = a2 - a1 + 1 ;
      geneBoxes = array (pp->geneBoxes, chrom, Array) ;
      geneExons = array (pp->geneExons, chrom, Array) ;
      if (! geneBoxes) geneBoxes = array (pp->geneBoxes, chrom, Array) = arrayHandleCreate (1000, GENE, pp->h) ;
      if (! geneExons) geneExons = array (pp->geneExons, chrom, Array) = arrayHandleCreate (1000, GENE, pp->h) ;
      

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;

      dictAdd (geneDict, cp, &gene) ;

      /* register intronSeeds */
      if (mrna == oldMrna && oldGene == gene && oldChrom == chrom && oldExon + 1 == exon)
	if (da < twoMb) /* our format only allows 21 bits for the intron length */ 
	  { /* register  intronSeeds */
	    int u1 = b2 + 1 ;
	    int u2 = a1 - 1 ;
	    int da =  u2 - u1 + 1 ;
	    
	    if (da > 3 && da < twoMb) /* our format only allows 21 bits for the intron length */ 
	      {
		up = bigArrayp (pp->intronSeeds, nnI++, CW) ;
		up->nam = chrom ;
		up->pos = u1 ;
		up->intron = u2 ;
	      }
	  }

      /* register exonSeeds */
      up = bigArrayp (pp->exonSeeds, nnES++, CW) ;
      up->nam = chrom ;
      up->pos = a1 ;
      up->intron = a2 ;

      /* register a geneBox */
      gb = arrayp (geneBoxes, gene, GENE) ;
      gb->gene = gene ;
      gb->a1 = (! gb->a1 || gb->a1 > a1) ? a1 : gb->a1 ;
      gb->a2 = (! gb->a2 || gb->a2 < a2) ? a2 : gb->a2 ;

      /* register a geneExon */
      ge = arrayp (geneExons, nnE++, GENE) ;
      ge->gene = gene ;
      ge->a1 = a1 ;
      ge->a2 = a2 ;

      oldMrna = mrna ;
      oldGene = gene ;
      oldChrom = chrom ;
      oldExon = exon ;
      b2 = a2 ;
    }

  for (chrom = 0 ; chrom < arrayMax (pp->geneBoxes) ; chrom++)
    {
      geneBoxes = array (pp->geneBoxes, chrom, Array) ;
      geneExons = array (pp->geneExons, chrom, Array) ;

      
      if (geneBoxes)
	{
	  int ii, jj, iMax = arrayMax (geneBoxes) ;
	  GENE *gb, *gb2 ;

	  arraySort (geneBoxes, geneA1Order) ;
	  
	  gb = gb2 = arrp (geneBoxes, 0, GENE) ;
	  for (ii = 0, jj = 0 ; ii < iMax ; ii++, gb++)
	    {
	      if (gb->gene)
		{
		  if (gb2->gene != gb->gene)
		    {
		      if (gb2->gene) { gb2++ ; jj++ ;}
		      if (gb2 != gb) *gb2 = *gb ;
		    }
		  else
		    {
		      gb2->a1 = (gb2->a1 < gb->a1 ? gb2->a1 : gb->a1) ;
		      gb2->a2 = (gb2->a2 > gb->a2 ? gb2->a2 : gb->a2) ;
		    }
		}
	    }
	  arrayMax (geneBoxes) = jj + 1 ;
	  arraySort (geneBoxes, a1GeneOrder) ;

	  nnnB += arrayMax (geneBoxes) ;
	}

      if (geneExons)
	{
	  int ii, jj, iMax = arrayMax (geneExons) ;
	  GENE *gb, *gb2 ;

	  arraySort (geneExons, geneA1Order) ;

	  gb = gb2 = arrp (geneExons, 0, GENE) ;
	  for (ii = 0, jj = 0 ; ii < iMax ; ii++, gb++)
	    {
	      if (gb->gene)
		{
		  if (gb2->gene != gb->gene || gb->a1 > gb2->a2 + 1)
		    {
		      if (gb2->gene) { gb2++ ; jj++ ;}
		      if (gb2 != gb) *gb2 = *gb ;
		    }
		  else
		    {
		      gb2->a1 = (gb2->a1 < gb->a1 ? gb2->a1 : gb->a1) ;
		      gb2->a2 = (gb2->a2 > gb->a2 ? gb2->a2 : gb->a2) ;
		    }
		}
	    }
	  arrayMax (geneExons) = jj + 1 ;

	  arraySort (geneExons, a1GeneOrder) ;
	  
	  nnnE += arrayMax (geneExons) ;
	}

    }


  fprintf (stderr, "+++++++ Found %d genes %d projected exons in file %s\n", nnnB, nnnE, filNam) ;
  ac_free (h) ;

  return nnnB ;
} /* saSpongeParserDo */

/**************************************************************/

int saSpongeParser (PP *pp, TC *tc)
{
  return saSpongeParserDo (pp, tc->fileName) ;
} /* saSpongeParser */

/**************************************************************/

int saSpongeParserDirect (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fNam = hprintf (h, "%s/tConfig", pp->indexName) ;
  ACEIN ai = aceInCreate (fNam, 0, h) ;

  aceInSpecial (ai, "\n") ;
  
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      if (! cp || *cp == '#')
	continue ;
      if (*cp == 'I')
	{

	}
    }
  

  if (pp->tFileSpongeFileNameF)
    saSpongeParserDo (pp, pp->tFileSpongeFileNameF) ;
  if (pp->tFileSpongeFileNameR)
    saSpongeParserDo (pp, pp->tFileSpongeFileNameR) ;

  ac_free (h) ;
  return 0 ;
} /* saSpongeParserDirect */

/**************************************************************/
/**************************************************************/
