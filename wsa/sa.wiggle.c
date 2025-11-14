/*
 * sa.wiggle.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * to the "wiggle" coverage plots
 *   construction of the different types of wiggle
 *   evaluation of gene expression
 *   annotating transcripts starts and ends
 *
 */

#define WIGGLE_STEP 10
#define WIGGLETYPEMAX 2 /* strand */
#include "sa.h"

typedef struct wigPosStruct { unsigned int pos ; unsigned short ln ; unsigned short weight ; } WP ;

/**************************************************************/

static int wpOrder (const void *va, const void *vb)
{
  const WP *up = va ;
  const WP *vp = vb ;
  int n ;
  n = (up->pos > vp->pos ? 1 : (up->pos < vp->pos ? -1 : 0)) ; if (n) return n ;

  return 0 ;
} /* wiggleOrder */

/*************************************************************************************/

void wiggleCumulate (BigArray aaa, BigArray aa)
{
  WP *up, *vp ;
  long int iMax = bigArrayMax (aa), iiMax = bigArrayMax (aaa) ;
  if (iMax)
    {
      bigArraySort (aa, wpOrder) ;
      up = bigArrayp (aaa, iMax + iiMax - 1, WP) ; /* make room */
      up = arrp (aaa, 0, WP) ;
      vp = arrp (aa, 0, WP) ;
      memcpy (up, vp, iMax * sizeof (WP)) ;
    }
  return ;
} /* wiggleCumulate */
  
/**************************************************************/
  
static int wiggleCreate (const PP *pp, BB *bb)
{
  ALIGN *ap ;
  long int ii, iMax = bigArrayMax (bb->aligns) ;
  int nsw = 0 ;
  int chrom = 0 ;
  BigArray wig = 0 ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  Array wiggles = bb->wiggles = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;
  char tBuf[25] ;
  int step = WIGGLE_STEP ;  /* examples s=10, 5, 1 */
  int demiStep = step/2 ;
  if (2*demiStep == step) demiStep-- ; /* examples d=4, 2, 0 */

  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int mult = ap->nTargetRepeats ;
      int weight = 720/mult ;
      if (weight)
	{
	  int a1 = ap->a1 ;
	  int a2 = ap->a2 ;
	  int w1 = (a1 + demiStep)/step ;
	  int w2 = (a2 + demiStep)/step ;

	  if (chrom != ap->chrom)
	    {
	      nsw++ ;
	      wig = 0 ;
	      chrom = ap->chrom ;
	      if (*dictName (pp->bbG.dict, chrom >> 1) == 'G')
		{
		  wig = array (wiggles, chrom, BigArray) ;
		  if (! wig)
		    wig = array (wiggles,  chrom, BigArray) = bigArrayHandleCreate (100000, WP, bb->h) ;
		}
	    }

	  if (w1 > w2) { int w0 = w1 ; w1 = w2 ; w2 = w0 ; }
	  if (wig && w1 < w2)
	    {
	      WP *wp = bigArrayp (wig, bigArrayMax (wig), WP) ;
	      wp->pos = w1 ; wp->ln = w2 - w1 + 1 ; wp->weight = weight ;
	    }
	}
    }
  fprintf (stderr, "%s: lane %d: wiggleCreate created %d wiggles, %d switching\n", timeBufShowNow (tBuf), bb->lane, arrayMax (wiggles), nsw) ;
  return  arrayMax (wiggles) ;
} /* wiggleCreate */

/*************************************************************************************/

void saWiggleCumulate (const PP *pp, BB *bb)
{
  Array wiggles = pp->wiggles ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int iwMax = wiggleCreate (pp, bb) ;
  BigArray *ap1, *ap0 = arrp (wiggles, 0, BigArray) ;

  if (iwMax > 2 * chromMax) messcrash ("iwMax too large ?") ;
  for (int iw = 0 ; iw < iwMax ; iw++)
    {
      BigArray aa = array (bb->wiggles, iw, BigArray) ;

      if (aa)
	{
	  BigArray aaa = array (wiggles, 2 * bb->run * chromMax + iw, BigArray) ;
	  if (! aaa)
	    aaa = array (wiggles, 2 * bb->run * chromMax + iw, BigArray) = bigArrayHandleCreate (10000, WP, pp->h) ;
	  if (1) wiggleCumulate (aaa, aa) ;
	}
    }
  ap1 = arrp (wiggles, 0, BigArray) ;
  if (ap1 != ap0)
    messcrash ("pp->wiggles was relocalized which is not allowed here because of multiuthreading") ;
  return ;
} /* saWiggleCumulate */

/**************************************************************/

static void wiggleExportOne (const PP *pp, int nw)
{
  BigArray wig = array (pp->wiggles, nw, BigArray) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int run = nw / (2 * chromMax) ;
  int chrom = (nw % (2 * chromMax)) ;
  char strand = ( nw & 0x1) ? 'r' : 'f' ;
  long int ii, iMax = bigArrayMax (wig) ;
  unsigned int pos0 ;
  Array genes = pp->genes ? array (pp->genes, chrom, Array) : 0 ;
  Array geneC = array (pp->geneCounts, nw, Array) ;
  int step = WIGGLE_STEP ;
  int demiStep = step/2 ;
  if (2*demiStep == step) demiStep-- ; /* examples d=4, 2, 0 */

  bigArraySort (wig, wpOrder) ;
    
  if (wig && iMax)
    {
      AC_HANDLE h = ac_new_handle () ;
      WP *wp0, *wp = arrp (wig, iMax - 1, WP) ;
      unsigned int posMax = wp->pos ;
      Array a = arrayHandleCreate (posMax + 1000, int, h) ;
      unsigned int *xp = arrayp (a, posMax, unsigned int) ;

      wp = arrp (wig, 0, WP) ;
      wp0 = arrp (wig, 0, WP) ;
      pos0 = wp0->pos ;
      for (ii = 0 ; ii < iMax ; ii++, wp++)
	{
	  if (wp->weight)
	    {
	      pos0 = pos0 ? pos0 : wp->pos ;
	      xp = arrayp (a, wp->pos + wp->ln - 1 - pos0, unsigned int) ;
	      xp -= wp->ln ;
	      for (int i = 0 ; i < wp->ln ; i++)
		xp[i] += wp->weight ;
	    }
	}

      if (arrayMax(a))
	{
	  const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) ;
	  const char *runNam = dictName (pp->runDict, run) ;
	  char *fNam = hprintf (h, ".%s.%s.u.%c.BF", runNam, chromNam, strand) ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
	  aceOutDate (ao, "##", "wiggle") ;
	  aceOutf (ao, "track type=wiggle_0\n") ;

	  aceOutf (ao, "fixedStep chrom=%s start=%d step=%d\n", chromNam, pos0 * step, step) ;

      	  xp = arrayp (a, 0, unsigned int) ;
	  for (int j = 0, jMax = arrayMax(a) ; j < jMax ; j++)
	    aceOutf (ao, "%u\n", xp[j]/720) ;
	}

      if (arrayMax(a) && genes)
	{
	  int ig = 0, igMax = arrayMax (genes) ;
	  GENE *gp = arrayp (genes, 0, GENE) ;
      	  xp = arrayp (a, 0, unsigned int) ;

	  for (int j = 0, jMax = arrayMax(a), x = WIGGLE_STEP * (j + pos0) ; j < jMax ; j++, x += WIGGLE_STEP)
	    {
	      while (ig > 0 && gp->a1 > x) { gp-- ; ig-- ; }
	      while (ig < igMax && gp->a2 < x) { gp++ ; ig++ ; }
	      while (ig < igMax && gp->a1 <= x + demiStep && gp->a2 >= x - demiStep)
		{ array (geneC, gp->gene, int) += xp[j]/720 ; gp++ ; ig++ ; }
	    }
	}

	  
      ac_free (h) ;
    }

  return ;
} /* wiggleExportOne */

/**************************************************************/

void wiggleExportAgent (const void *vp)
{
  const PP *pp = vp ;
  int nw ;
  
  while (channelGet (pp->wwChan, &nw, int))
    {
      wiggleExportOne (pp, nw) ;
      channelPut (pp->wwDoneChan, &nw, int) ;
    }
  channelCloseSource (pp->wwDoneChan) ;
  return ;
} /* wiggleExportAgent */

/**************************************************************/
/**************************************************************/

typedef struct gcStruct {
  int gene, run, count, dummy ; 
} __attribute__((aligned(16))) GC ;

/**************************************************************/

static int gcOrder (const void *va, const void *vb)
{
  const GC *up = va ;
  const GC *vp = vb ;
  int n ;
  
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;

  return 0 ;
} /* wiggleOrder */

/**************************************************************/

static void wiggleExportGeneCounts (const PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nw, wMax = arrayMax (pp->wiggles) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  BigArray allGeneC ;
  long int igc = 0, igcMax = 0 ;
  GC *gc ; 
  char tBuf[25] ;
  
  fprintf (stderr, "%s: start geneCounts export\n", timeBufShowNow (tBuf)) ;
  
  allGeneC = bigArrayHandleCreate (1000000, GC, h) ;
  for (nw = 0 ; nw < wMax ; nw++)
    {
      Array geneC = array (pp->geneCounts, nw, Array) ;
      int run = nw / (2 * chromMax) ;
      int gMax = geneC ? arrayMax (geneC) : 0 ;
      int gene, *xp ;
      
      if (gMax)
	for (gene = 0, xp = arrp (geneC, 0, int) ; gene < gMax ; gene++, xp++)
	  if (*xp > 0)
	    {
	      gc = bigArrayp (allGeneC, igcMax++, GC) ;
	      gc->gene = gene ;
	      gc->run = run ;
	      gc->count = *xp ;
	    }
    }
  bigArraySort (allGeneC, gcOrder) ;
  
  ACEOUT ao = aceOutCreate (pp->outFileName, ".geneCounts.tsf", pp->gzo, h) ;
  aceOutDate (ao, "##", "wiggle") ;
  for (igc = 0, gc = bigArrp (allGeneC, 0, GC) ; igc < igcMax ; igc++, gc++)
    aceOutf (ao, "%s\t%s\ti\t%d\n"
	     , dictName (pp->geneDict, gc->gene)
	     , dictName (pp->runDict, gc->run)
	     , gc->count
	     ) ;
  
  fprintf (stderr, "%s: stop geneCounts export\n", timeBufShowNow (tBuf)) ;
  ac_free (h) ;
  return ;
} /* wiggleExportGeneCounts */

/**************************************************************/
/**************************************************************/

void saWiggleExport (PP *pp, int nAgents)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  BOOL debug = FALSE ;
  char tBuf[25] ;

  if (pp->genes)
    pp->geneCounts = arrayHandleCreate (wMax, Array, h) ;
  fprintf (stderr, "%s: start exportation of  %d wiggles\n", timeBufShowNow (tBuf), wMax) ;
  /* parallelize: open a channel and start agents */
  pp->wwChan = channelCreate (wMax + 1, int, h) ;
  channelDebug (pp->wwChan, debug, "wwChan") ;
  pp->wwDoneChan = channelCreate (wMax+1, int, h) ;
  channelDebug (pp->wwDoneChan, debug, "wwChan") ;

  for (int ii = 0 ; ii < nAgents ; ii++)
    {
      pp->agent = ii ;
      wego_go (wiggleExportAgent, pp, PP) ; channelAddSources (pp->wwDoneChan, 1) ;
    }
  
  /* load the channel to start execution */ 

  int k = 0, n = 0 ;
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      Array wig = arr (pp->wiggles, nw, Array) ;
      if (wig)
	{
	  array (pp->geneCounts, nw, Array) = arrayHandleCreate (2048, dictMax (pp->geneDict) + 1, h) ;
	  channelPut (pp->wwChan, &nw, int) ;
	}
    }
  channelClose (pp->wwChan) ;

  /* synchronize */
  while (channelGet (pp->wwDoneChan, &k, int))
    n++ ;

  fprintf (stderr, "%s: stop wiggle export\n", timeBufShowNow (tBuf)) ;

  if (pp->genes)
    wiggleExportGeneCounts (pp) ;
  
  ac_free (h) ;
  return ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
 
