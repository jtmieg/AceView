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

static void showEX (Array aa, int start)
{
  if (aa)
    {
      EXONINTRON *up = aa && start < arrayMax (aa) ? arrp (aa, start, EXONINTRON) : 0 ;
      if (up)
	for (int i = start ; i < arrayMax (aa) && i < start + 50 ; up++, i++)
	  fprintf (stderr, "%d :: gene:%d m:%d chr:%d a=%d__%d\n", i, up->gene, up->mrna, up->chrom, up->a1, up->a2) ;
    }
  return ;
} /* showEX */

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

static int exonGeneOrder (const void *va, const void *vb)
{
  const EXONINTRON *up = va ;
  const EXONINTRON *vp = vb ;
  int n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  return 0 ;
} /* exonGeneOrder */

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
/**************************************************************/
/* memory map the annotated introns and gene boxes */
void saGffBinaryParser (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL READONLY = TRUE ;
  DICT *geneDict = pp->geneDict = dictHandleCreate (10000, pp->h) ;
  DICT *mrnaDict = pp->mrnaDict = dictHandleCreate (10000, pp->h) ;

  if (0)
    {
      char *fNam = "/mrnaNames.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  saDictMapRead (mrnaDict, cp) ;
	  fprintf (stderr, "saGffParser found %d mrna identifiers in file %s\n", dictMax (mrnaDict), cp) ;
	}
    }
  
  if (1)
    { /* used in the introns and gene expression reports */
      char *fNam = "/geneNames.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  saDictMapRead (geneDict, cp) ;
	  fprintf (stderr, "saGffParser found %d gene identifiers in file %s\n", dictMax (geneDict), cp) ;
	}
    }

  if (1) /*annotated intron coordinates and gene name */
    {
      char *fNam = "/knownIntrons.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  pp->knownIntrons = bigArrayMapRead (cp, EXONINTRON, READONLY, pp->h) ;
	  long int nnI = bigArrayMax (pp->knownIntrons) ;
	  fprintf (stderr, "saGffParser found %ld introns coordinates in file %s\n", nnI, cp) ;
	}
    }

  if (pp->wiggle)
    { /* gene expression records */
      long int nnE = 0 ;
      BigArray geneBoxes = 0 ;
      char *fNam = "/geneBoxes.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  geneBoxes = bigArrayMapRead (cp, EXONINTRON, READONLY, h) ;
	  nnE = bigArrayMax (geneBoxes) ;
	  fprintf (stderr, "saGffParser found %ld gene fragments coordinates in file %s\n", nnE, cp) ;
	}

      if (nnE) /* split per chromosomes */
	{
	  EXONINTRON *up, *vp ;
	  long int ii ;
	  pp->geneBoxes = arrayHandleCreate (dictMax (pp->bbG.dict) + 1, Array, pp->h) ;
	  for (ii = 0, up = bigArrp (geneBoxes, 0, EXONINTRON) ; ii < nnE ; ii++, up++)
	    {
	      int chrom = up->chrom ;
	      Array chromBoxes = array (pp->geneBoxes, chrom, Array) ;
	      if (! chromBoxes)
		chromBoxes = array (pp->geneBoxes, chrom, Array) = arrayHandleCreate (10000, EXONINTRON, pp->h) ;
	      vp = arrayp (chromBoxes, arrayMax (chromBoxes), EXONINTRON) ;
	      *vp = *up ;	      
	    }
	}
    }
} /* saGffBinaryParser */

/**************************************************************/

long int saGffParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  long int nnE = 0, nnI = 0, nnCDS = 0, nnGE = 0, nnME = 0, nnPP = 0, nnPP2 = 0 ;
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
  Array geneExons = 0 ;
  Array geneCDS = 0 ;
  Array geneBoxes = 0 ;
  Array geneBoxes2 = 0 ;
  Array mrnaExons = 0 ;
  char strand = 0 ;
  
  if (! geneDict) geneDict = pp->geneDict = dictHandleCreate (1 << 15, pp->h) ;
  if (! mrnaDict) mrnaDict = pp->mrnaDict = dictHandleCreate (1 << 15, pp->h) ;

  if (! pp->knownIntrons)   pp->knownIntrons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;
  if (! pp->geneBoxes) pp->geneBoxes = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;

  geneExons = arrayHandleCreate (100000, EXONINTRON, h) ;
  geneCDS = arrayHandleCreate (100000, EXONINTRON, h) ;
  geneBoxes = arrayHandleCreate (100000, EXONINTRON, h) ;
  mrnaExons = arrayHandleCreate (100000, EXONINTRON, h) ;
  
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
      if (! strcasecmp (cp, "CDS")) type = 2 ;
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
	strand = 0x0 ;
      else if (*cp == '-')
	strand = 0x1 ;
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
      if (type >= 1) /* identify the transcript */
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
      if (type >= 1) /* identify the gene */
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
      
      if ((type == 1 || type == 2) && gene)  /* guarantee that all CDS are also annotated as exons */
	{
	  EXONINTRON *up = arrayp (geneExons, nnGE++, EXONINTRON) ;
	  up->gene = (gene << 1) | strand ;
	  up->mrna = 1 ;
	  up->chrom = chrom ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}

      if (type == 2 && gene)
	{
	  EXONINTRON *up = arrayp (geneCDS, nnCDS++, EXONINTRON) ;
	  up->gene = (gene << 1) | strand ;
	  up->mrna = type ;
	  up->chrom = chrom ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}

      if (type == 1 && gene && mrna)
	{
	  EXONINTRON *up = arrayp (mrnaExons, nnME++, EXONINTRON) ;
	  up->gene = gene ;
	  up->mrna = mrna ;
	  up->chrom = (chrom << 1) | strand ; 
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}
    }

  if (nnME) /* create the mRNA specific corresponding introns */
    {
      EXONINTRON *up, *vp ;
      int ii, jj ;
      
      arraySort (mrnaExons, exonMrnaOrder) ;
      for (up = arrp (mrnaExons, 0, EXONINTRON), ii = 0 ; ii < nnME ; ii++, up++)
	if (up->mrna)  
	  {
	    for (vp = up + 1, jj = ii + 1 ; jj < nnME && vp->mrna == up->mrna ; vp++, jj++) 
	      {
		int chrom = up->chrom ;
		int a1 = vp[-1].a2 + 1 ;
		int a2 = vp[0].a1 - 1 ;
		int da = a2 -= a1 + 1 ;
		if (0)
		  fprintf (stderr, "Intron %ld : %d %d\n", nnI, a1, a2) ;
		if (da >= twoMb) /* our format only allows 21 bits for the intron length */ 
		  continue ;

		EXONINTRON *wp = bigArrayp (pp->knownIntrons, nnI++, EXONINTRON) ;
		wp->chrom = chrom ;
		wp->mrna = 0 ;
		wp->gene = up->gene ;
		wp->a1 = a1 ;
		wp->a2 = a2 ;
		vp->mrna = 0 ;
	      }
	    up->mrna = 0 ;
	  }
      bigArraySort (pp->knownIntrons, exonA1Order) ;
      bigArrayCompress (pp->knownIntrons) ;
      nnI = bigArrayMax (pp->knownIntrons) ;
      
      /* memory map the known introns */
      if (nnI)
	{
	  char *fNam = hprintf (h, "%s/knownIntrons.sortali", pp->indexName) ;
	  bigArrayMapWrite (pp->knownIntrons, fNam) ;
	  fprintf (stderr, "saGffParser exported %ld knownIntrons\n", bigArrayMax (pp->knownIntrons)) ;
	}
      ac_free (mrnaExons) ;
    }
  
  
  if (nnGE) /* merge overlapping or continous exons of the same gene */
    {
      EXONINTRON *up, *vp ;
      int ii, jj ;
      
      arraySort (geneExons, exonGeneOrder) ;
      /* merge */
      for (up = arrp (geneExons, 0, EXONINTRON), ii = 0 ; ii < nnGE ; ii++, up++)
	if (up->gene)  
	  {
	    for (vp = up + 1, jj = ii + 1 ; jj < nnGE && vp->gene == up->gene ; vp++, jj++) 
	      {
		if (vp->a1 <= up->a2 + 1) /* merge */
		  {
		    if (vp->a2 > up->a2)
		      up->a2 = vp->a2 ;
		    vp->gene = 0 ;
		  }
		else
		  break ;
	      }
	  }
      /* compress */
      for (up = vp = arrp (geneExons, 0, EXONINTRON), ii = jj = 0 ; ii < nnGE ; ii++, up++)
	if (up->gene)
	  {
	    if (up > vp) *vp = *up ;
	    vp++ ; jj++ ;
	  }
      nnGE = arrayMax (geneExons) = jj ;
    }
  
  if (nnCDS) /* merge overlapping or continous CDS of the same gene */
    {
      EXONINTRON *up, *vp ;
      int ii, jj ;
      
      arraySort (geneCDS, exonGeneOrder) ;
      /* merge */
      for (up = arrp (geneCDS, 0, EXONINTRON), ii = 0 ; ii < nnCDS ; ii++, up++)
	if (up->gene)  
	  {
	    for (vp = up + 1, jj = ii + 1 ; jj < nnCDS && vp->gene == up->gene ; vp++, jj++) 
	      {
		if (vp->a1 <= up->a2 + 1) /* merge */
		  {
		    if (vp->a2 > up->a2)
		      up->a2 = vp->a2 ;
		    vp->gene = 0 ;
		  }
		else
		  break ;
	      }
	  }
      /* compress */
      for (up = vp = arrp (geneCDS, 0, EXONINTRON), ii = jj = 0 ; ii < nnCDS ; ii++, up++)
	if (up->gene)
	  {
	    if (up > vp) *vp = *up ;
	    vp++ ; jj++ ;
	  }
      nnCDS = arrayMax (geneCDS) = jj ;
    }
  
  if (nnCDS + nnGE) /* create a set of CDS/UTR/INTRON boxes of each gene */
    {
      EXONINTRON *up, *vp, *wp ;
      int ii = 0, jj = 0 ;
      
      arraySort (geneExons, exonGeneOrder) ;
      arraySort (geneCDS, exonGeneOrder) ;
      up = nnGE ? arrp (geneExons, 0, EXONINTRON) : 0 ;
      vp = nnCDS ? arrp (geneCDS, 0, EXONINTRON) : 0 ;
      for (ii = 0 ; ii < nnGE ; ii++, up++)
	{
	  if (ii > 0 && up[0].gene == up[-1].gene && up[0].a1 > up[-1].a2 + 1) 
	    { /* create a UTR segment */
	      wp = arrayp (geneBoxes, nnPP++, EXONINTRON) ;
	      wp->gene = up->gene ;
	      wp->mrna = 1 ; /* INTRON */
	      wp->chrom = up->chrom ;
	      wp->a1 = up[-1].a2 + 1 ;
	      wp->a2 = up[0].a1 - 1 ;
	    }

	  
	  if (vp)
	    {
	      for ( ; jj > 0 && vp->gene > up->gene ; jj--, vp--)
		if (vp->gene == up->gene && vp->a2 < up->a1)
		  break ;
	      for ( ; jj < nnCDS && vp->gene < up->gene ; jj++, vp++)
		;
	    }
	  if (vp && vp->gene == up->gene)
	    for ( ; jj < nnCDS && vp->gene == up->gene && vp->a1 <= up->a2 ; jj++, vp++)
	      {
		if (vp->a1 > up->a1) /* create a UTR segment */
		  {
		    wp = arrayp (geneBoxes, nnPP++, EXONINTRON) ;
		    wp->gene = up->gene ;
		    wp->mrna = 2 ; /* UTR */
		    wp->chrom = up->chrom ;
		    wp->a1 = up->a1 ;
		    wp->a2 = vp->a1 - 1 ;
		    up->a1 = vp->a1 ;
		  }
		if (vp->a1 >= up->a1 && vp->a2 <= up->a2) /* create a CDS segment */
		  {
		    wp = arrayp (geneBoxes, nnPP++, EXONINTRON) ;
		    wp->gene = up->gene ;
		    wp->mrna = 4 ; /* CDS */
		    wp->chrom = up->chrom ;
		    wp->a1 = vp->a1 ;
		    wp->a2 = vp->a2 ;
		    up->a1 = vp->a2 + 1 ;
		  }
	      }
	  if (up->a1 < up->a2)   /* create a UTR segment */
	    {
	      wp = arrayp (geneBoxes, nnPP++, EXONINTRON) ;
	      wp->gene = up->gene ;
	      wp->mrna = 2 ; /* UTR */
	      wp->chrom = up->chrom ;
	      wp->a1 = up->a1 ;
	      wp->a2 = up->a2 ;
	      up->a1 = up->a2 + 1 ;
	    }
	}
    }

  if (nnPP) /* identify the ambiguous segments, including on same of opposite strands */
    {
      EXONINTRON *up, *vp, *wp ;
      int ii = 0, jj = 0 ;

      geneBoxes2 = arrayHandleCreate (nnPP * 1.3, EXONINTRON, h) ;

      arraySort (geneBoxes, exonA1Order) ;
      up = nnPP ? arrp (geneBoxes, 0, EXONINTRON) : 0 ;

      for (ii = 0 ; ii < nnPP ; ii++, up++)
	{
	  for (vp = up + 1, jj = ii + 1 ; jj < nnPP && vp->chrom == up->chrom && vp->a1 <= up->a2 ; jj++, vp++)
	    {
	      if (up->a1 < vp->a1)   /* create a non ambiguous segment */
		{
		  wp = arrayp (geneBoxes2, nnPP2++, EXONINTRON) ;
		  wp->gene = up->gene ;
		  wp->mrna = up->mrna ; /* type */
		  wp->chrom = up->chrom ;
		  wp->a1 = up->a1 ;
		  wp->a2 = vp->a1 - 1 ;
		  up->a1 = vp->a1 ;
		}
	      if (1) /* create one or even two ambiguous segment */
		{
		  int a2 = up->a2 <= vp->a2 ? up->a2 : vp->a2 ;
		  wp = arrayp (geneBoxes2, nnPP2++, EXONINTRON) ;
		  wp->chrom = up->chrom ;
		  wp->a1 = up->a1 ;
		  wp->a2 = a2 ;
		  up->a1 = vp->a1 = a2 + 1 ;
		
		  if (up->type == 1 && vp->type == 1)
		    {   /* double introns create an intronic segment attributed to no gene */
		      wp->gene = 0 ;
		      wp->mrna = 1 ;
		    }
		  else if (up->type > 1 && vp->type == 1)
		    {   /* attribute to up gene */
		      wp->gene = up->gene ;
		      wp->mrna = up->mrna ;
		    }
		  else if (up->type == 1 && vp->type > 1)
		    {   /* attribute to vp gene */
		      wp->gene = vp->gene ;
		      wp->mrna = vp->mrna ;
		    }
		  else if (up->type > 1 && vp->type > 1)
		    {   /* create 2 ambiguous segments */
		      wp->gene = up->gene ;
		      wp->mrna = up->mrna + 8 ;
		      wp = arrayp (geneBoxes2, nnPP2++, EXONINTRON) ;
		      wp->chrom = up->chrom ;
		      wp->a1 = up->a1 ;
		      wp->a2 = a2 ;
		      wp->gene = vp->gene ;
		      wp->mrna = vp->mrna + 8 ;
		    }
		}
	    }
	  if (up->a1 < up->a2)   /* create a non ambiguous segment */
	    {
	      wp = arrayp (geneBoxes2, nnPP2++, EXONINTRON) ;
	      wp->gene = up->gene ;
	      wp->mrna = up->mrna ; /* type */
	      wp->chrom = up->chrom ;
	      wp->a1 = up->a1 ;
	      wp->a2 = up->a2 ;
	    }
	}
    }

  if (nnPP2) /* construct and memory map the final gene boxes */
    {
      EXONINTRON *up, *vp ;
      int ii ;
      BigArray geneBoxes3 = bigArrayHandleCreate (nnPP2, EXONINTRON, h) ;
      arraySort (geneBoxes2, exonA1Order) ;
      for (ii = 0, nnE = 0,  up = arrp (geneBoxes2, 0, EXONINTRON) ; ii < nnPP2 ; ii++, up++)
	{
	  vp = bigArrayp (geneBoxes3, nnE++, EXONINTRON) ;
	  *vp = *up ;
	}

      char *fNam = hprintf (h, "%s/geneBoxes.sortali", pp->indexName) ;
      bigArrayMapWrite (geneBoxes3, fNam) ;
      fprintf (stderr, "saGffParser exported %ld gene box segments\n", bigArrayMax (geneBoxes3)) ;
    }

  if (geneDict && dictMax (geneDict))
    {
      char *fNam = hprintf (h, "%s/geneNames.sortali", pp->indexName) ;
      saDictMapWrite (geneDict, fNam) ;
      fprintf (stderr, "saGffParser exported %d gene identifiers\n", dictMax (geneDict)) ;
    }
  if (0 && mrnaDict && dictMax (mrnaDict))
    {
      char *fNam = hprintf (h, "%s/mrnaNames.sortali", pp->indexName) ;
      saDictMapWrite (mrnaDict, fNam) ;
      fprintf (stderr, "saGffParser exported %d mrna identifiers\n", dictMax (mrnaDict)) ;
    }

  fprintf (stderr, "+++++++ Found %ld boxes %ld introns in file %s\n", nnE, nnI, tc->fileName) ;
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
  EXONINTRON *up ;
  
  if (pp->seedLength < 16)
    messcrash ("\nSorry, to study the introns defined in file %s,\n the seed length must be at least 16, not %d\n"
	       , aceInFileName (ai)
	       , pp->seedLength
	       ) ;
  if (! pp->knownIntrons)
    pp->knownIntrons = bigArrayHandleCreate (100000, EXONINTRON, pp->h) ;

  nn = bigArrayMax (pp->knownIntrons) ;
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
		   , tc->fileName
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

  bigArraySort (pp->knownIntrons, exonA1Order) ;
  bigArrayCompress (pp->knownIntrons) ;
  
  /* memory map the known introns */
  if (pp->knownIntrons && bigArrayMax (pp->knownIntrons))
    {
      char *fNam = hprintf (h, "%s/knownIntrons.sortali", pp->indexName) ;
      bigArrayMapWrite (pp->knownIntrons, fNam) ;
      fprintf (stderr, "genomeCreateBinary exported %ld known_introns\n", bigArrayMax (pp->knownIntrons)) ;
    }

  showEX (0, 0) ; /* for compiler happiness */
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, tc->fileName) ;
  ac_free (h) ;

  return nn ;
} /* saIntronParser */

/**************************************************************/
/**************************************************************/
