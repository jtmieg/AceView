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

/*************************************************************************************/

void wiggleCumulate (BigArray aaa, BigArray aa)
{
  long int i, iMax = bigArrayMax (aa) ;
  if (iMax)
    {
      i = bigArray (aaa, iMax - 1, int) ; /* make room */
      for (i = 0 ; i < iMax ; i++)
	bigArray (aaa, i, int) += bigArr (aa, i, int) ;
    }
  return ;
} /* wiggleCumulate */
  
/**************************************************************/
  
static int wiggleCreate (const PP *pp, BB *bb)
{
  ALIGN *ap ;
  long int ii, iMax = bigArrayMax (bb->aligns) ;
  int chrom = 0 ;
  BigArray wig = 0 ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  Array wiggles = bb->wiggles = arrayHandleCreate (2 * chromMax, BigArray, bb->h) ;

  saSort (bb->aligns, 4) ;
  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int w1 = ap->w1, w2 = ap->w2 ;
      int mult = ap->nTargetRepeats ;
      int weight = 720/mult ;
      
      if (! weight)
	continue ;
      if (chrom != ap->chrom)
	{
	  wig = 0 ;
	  chrom = ap->chrom ;
	  if (*dictName (pp->bbG.dict, chrom >> 1) == 'G')
	    {
	      wig = array (wiggles, chrom, BigArray) ;
	      if (! wig)
		wig = array (wiggles,  chrom, BigArray) = bigArrayHandleCreate (100000, int, bb->h) ;
	    }
	}
      if (wig && w1 < w2)
	for (int i = w1 ; i <= w2 ; i++)
	  bigArray (wig, i, int) += weight ;
    }
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
	    aaa = array (wiggles, 2 * bb->run * chromMax + iw, BigArray) = bigArrayHandleCreate (100000, int, pp->h) ;
	  saSort (aa, 4) ;
	  wiggleCumulate (aaa, aa) ;
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
  int step = WIGGLE_STEP ;
  BigArray wig = array (pp->wiggles, nn, BigArray) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  int run = nn / (2 * chromMax) ;
  int chrom = (nn % (2 * chromMax)) >> 1 ;
  char strand = ( nn & 0x1) ? 'r' : 'f' ;
  
  if (wig && arrayMax (wig))
    {
      AC_HANDLE h = ac_new_handle () ;
      const char *chromNam = dictName (pp->bbG.dict, chrom) ;
      const char *runNam = dictName (pp->runDict, run) ;
      char *fNam = hprintf (h, ".%s.%s.u.%c.BF", runNam, chromNam, strand) ;
      ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
      aceOutDate (ao, "##", "wiggle") ;
      aceOutf (ao, "track type=wiggle_0\n") ;
      aceOutf (ao, "fixedStep chrom=%s start=%d step=%d\n", chromNam, step, step) ;
      
      for (long int jj = 1 ; jj < bigArrayMax (wig) ; jj++)
	aceOutf (ao, "%d\n", bigArr (wig, jj, int)) ;
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
      wego_go (wiggleExportAgent, pp, PP) ;
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
  while (n < nn)
    {
      channelGet (pp->wwDoneChan, &k, int) ;
      n++ ;
    }
  fprintf (stderr, "%s: stop wiggle export\n", timeBufShowNow (tBuf)) ;
  ac_free (h) ;
  return ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
