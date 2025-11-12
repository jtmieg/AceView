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

  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int w1 = ap->w1, w2 = ap->w2 ;
      int mult = ap->nTargetRepeats ;
      int weight = 720/mult ;
      
      if (! weight)
	continue ;
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

static void wiggleExportOne (const PP *pp, int nn)
{
  BigArray wig = array (pp->wiggles, nn, BigArray) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int run = nn / (2 * chromMax) ;
  int chrom = (nn % (2 * chromMax)) >> 1 ;
  char strand = ( nn & 0x1) ? 'r' : 'f' ;
  long int ii, iMax = bigArrayMax (wig) ;
  unsigned int pos0 ;
  
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

      if (arrayMax (a))
	{
	  const char *chromNam = dictName (pp->bbG.dict, chrom) ;
	  const char *runNam = dictName (pp->runDict, run) ;
	  char *fNam = hprintf (h, ".%s.%s.u.%c.BF", runNam, chromNam, strand) ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
	  aceOutDate (ao, "##", "wiggle") ;
	  aceOutf (ao, "track type=wiggle_0\n") ;

	  aceOutf (ao, "fixedStep chrom=%s start=%d step=%d\n", chromNam, pos0 * WIGGLE_STEP , WIGGLE_STEP) ;
      

      	  xp = arrayp (a, 0, unsigned int) ;
	  for (ii = 0 ; ii < iMax ; ii++)
	    aceOutf (ao, "%u\n", xp[ii]/720) ;
	}
      ac_free (h) ;
    }

  return ;
} /* wiggleExportOne */

/**************************************************************/

void wiggleExportAgent (const void *vp)
{
  const PP *pp = vp ;
  int nn ;
  
  while (channelGet (pp->wwChan, &nn, int))
    {
      wiggleExportOne (pp, nn) ;
      channelPut (pp->wwDoneChan, &nn, int) ;
    }
  channelCloseSource (pp->wwDoneChan) ;
  return ;
} /* wiggleExportAgent */

/**************************************************************/

void saWiggleExport (PP *pp, int nAgents)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  BOOL debug = FALSE ;
  char tBuf[25] ;
  
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

  int k = 0, n = 0, nn = 0 ;
  for (int w = 0 ; w < wMax ; w++)
    {
      Array wig = arr (pp->wiggles, w, Array) ;
      if (wig)
	{
	  channelPut (pp->wwChan, &w, int) ;
	  nn++ ;
	}
    }
  channelClose (pp->wwChan) ;

  /* synchronize */
  while (channelGet (pp->wwDoneChan, &k, int))
    n++ ;

  fprintf (stderr, "%s: stop wiggle export\n", timeBufShowNow (tBuf)) ;
  ac_free (h) ;
  return ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
