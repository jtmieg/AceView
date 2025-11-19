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
  const int step = pp->wiggle_step ;  /* examples s=10, 5, 1 */
  const int demiStep = step/2 - (step % 1);

  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int mult = ap->nTargetRepeats ? ap->nTargetRepeats : 1 ;
      int weight = 720/mult ;
      if (weight)
	{
	  int a1 = ap->a1 ;
	  int a2 = ap->a2 ;
	  int w1 = (a1 + demiStep)/step ;
	  int w2 = (a2 + demiStep)/step ;
	  int aChrom = ap->chrom ^ (ap->read & 0x1) ; 

	  if (1) aChrom ^= 0x1 ; /* negative run */

	  /*  chr 0  read 0 faux,
	   *  chr 0 read 1 ok
	   *  chr 1 read 1 faux
	   *  chr 1 read 0 ok
	   */
	  if (0 && !(ap->chrom & 0x1)) continue ;
	  if (0 && (ap->read & 0x1)) continue ;
	  if ( ((ap->chrom ^ ap->read) & 0x1)) continue ; 
	  if (chrom != aChrom)
	    {
	      nsw++ ;
	      wig = 0 ;
	      chrom = aChrom ;
	      if (*dictName (pp->bbG.dict, chrom >> 1) == 'G' ||
		  *dictName (pp->bbG.dict, chrom >> 1) == 'M' ||
		  *dictName (pp->bbG.dict, chrom >> 1) == 'C'
		  )
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
  if (0) fprintf (stderr, "%s: lane %d: wiggleCreate created %d wiggles, %d switching\n", timeBufShowNow (tBuf), bb->lane, arrayMax (wiggles), nsw) ;
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
  long int cumul = 0 ;
  long int exonic = 0, intronic = 0, intergenic = 0 ;
  unsigned int pos0 ;
  Array geneExons = pp->geneExons ? array (pp->geneExons, chrom, Array) : 0 ;
  Array geneBoxes = pp->geneBoxes ? array (pp->geneBoxes, chrom, Array) : 0 ;
  Array geneC = pp->geneExons ? array (pp->geneExonCounts, nw, Array) : 0 ;
  Array geneB = pp->geneBoxes ? array (pp->geneBoxCounts, nw, Array) : 0 ;
  const int step = pp->wiggle_step ;
  const int demiStep = step/2 - (step % 1);

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
	    {
	      unsigned int w = xp[j] ;
	      aceOutf (ao, "%u\n", w / 720) ;
	      cumul += w * step ;
	    }
	}
      
      if (arrayMax(a) && geneBoxes)
	{
	  int ie = 0, ie1 = 0, ieOld = 0, ieMax = arrayMax (geneExons) ;
	  int ib = 0, ib1 = 0, ibOld = 0, ibMax = arrayMax (geneBoxes) ;
	  GENE *ge = arrayp (geneExons, 0, GENE) ;
	  GENE *gb = arrayp (geneBoxes, 0, GENE) ;
      	  xp = arrayp (a, 0, unsigned int) ;
	  
	  /* the candidate genes that may cover position x are
	   * not earlier than the first gene igOld  covering the previous position
	   * not later than the fist gene ig1 starting after x
	   */
	  for (int j = 0, jMax = arrayMax(a), x = step * (j + pos0) ; j < jMax ; j++, x += step) 
	    {
	      int isGeneTr = 0 ;
	      int weight = step * xp[j] ; 

	      if (! weight)
		continue ;
	      /* find the first gene starting rtight of x + demi step */
	      for (gb = arrp (geneBoxes, ib1, GENE) ; ib1 < ibMax && gb->a1 < x + demiStep ; ib1++, gb++)
		;
	      /* find the genes intersecting x and reset ibOld */
	      for (ib = ibOld, gb = arrp (geneBoxes, ib, GENE), ibOld = 0  ; ib < ib1 ; ib++, gb++)

		{
		  if (gb->a1 <= x + demiStep && gb->a2 >= x - demiStep)
		    {
		      ibOld = ibOld ? ibOld : ib ;
		      array (geneB, gb->gene, int) += weight ;
		      isGeneTr = 1 ;
		    }
		}
	      ibOld = ibOld ? ibOld : ib1 ;
	      
	      if (isGeneTr)  /* we are ina gene, are we in exons */
		{
		  /* find the first exon starting rtight of x + demi step */
		  for (ge = arrp (geneExons, ie1, GENE) ; ie1 < ieMax && ge->a1 < x + demiStep ; ie1++, ge++)
		    ;
		  /* find the exons intersecting x and reset ibOld */
		  for (ie = ieOld, ge = arrp (geneExons, ie, GENE), ieOld = 0  ; ie < ie1 ; ie++, ge++)
		    {
		      if (ge->a1 <= x + demiStep && ge->a2 >= x - demiStep)
			{
			  ieOld = ieOld ? ieOld : ie ;
			  array (geneC, ge->gene, int) += weight ;
			  isGeneTr = 2 ;
			}
		    }
		}
	      switch (isGeneTr)
		{
		case 2: exonic += weight ; break ;
		case 1: intronic += weight ; break ;
		default: intergenic += weight ; break ; 
		}
	  }
	}

      ac_free (h) ;
    }
  array (pp->wiggleCumuls, nw, long int) = cumul ;
  array (pp->exonics, nw, long int) = exonic ;
  array (pp->intronics, nw, long int) = intronic ;
  array (pp->intergenics, nw, long int) = intergenic ;
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
  int gene, run, boxCount, exonCount ; 
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

static long int wiggleExportGeneCounts (const PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nw, wMax = arrayMax (pp->wiggles) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  BigArray allGeneC ;
  long int igc = 0, jgc, igcMax = 0, nnn = 0 ;
  GC *gc, *gc2 ; 
  char tBuf[25] ;
  
  fprintf (stderr, "%s: start geneCounts export\n", timeBufShowNow (tBuf)) ;
  
  allGeneC = bigArrayHandleCreate (1000000, GC, h) ;
  for (nw = 0 ; nw < wMax ; nw++)
    {
      Array geneC = array (pp->geneExonCounts, nw, Array) ;
      Array geneB = array (pp->geneBoxCounts, nw, Array) ;
      int run = nw / (2 * chromMax) ;
      int gCMax = geneC ? arrayMax (geneC) : 0 ;
      int gBMax = geneB ? arrayMax (geneB) : 0 ;
      int gene, *xp ;
      if (gCMax)
	for (gene = 0, xp = arrp (geneC, 0, int) ; gene < gCMax ; gene++, xp++)
	  if (*xp > 0)
	    {
	      gc = bigArrayp (allGeneC, igcMax++, GC) ;
	      gc->gene = gene ;
	      gc->run = run ;
	      gc->exonCount = *xp ;
	    }
      if (gBMax)
	for (gene = 0, xp = arrp (geneC, 0, int) ; gene < gBMax ; gene++, xp++)
	  if (*xp > 0)
	    {
	      gc = bigArrayp (allGeneC, igcMax++, GC) ;
	      gc->gene = gene ;
	      gc->run = run ;
	      gc->boxCount = *xp ;
	      nnn += *xp ;
	    }
    }
  bigArraySort (allGeneC, gcOrder) ;
  for (igc = 0, gc = bigArrp (allGeneC, 0, GC), gc2 = gc ; igc < igcMax ; igc++, gc++)
    {
      if (gc2->gene != gc->gene || gc2->run != gc->run)
	{
	  gc2++ ; jgc++ ;
	  if (gc2 < gc) *gc2 = *gc ;
	}
      else
	{
	  gc2->boxCount += gc->boxCount ;
	  gc2->exonCount += gc->exonCount ;
	}
    }
  bigArrayMax (allGeneC) = jgc ;
  
  ACEOUT ao = aceOutCreate (pp->outFileName, ".geneCounts.tsf", pp->gzo, h) ;
  aceOutDate (ao, "##", "wiggle") ;
  for (igc = jgc = 0, gc = bigArrp (allGeneC, 0, GC), gc2 = gc ; igc < igcMax ; igc++, gc++)
    if (gc->boxCount)
      aceOutf (ao, "%s\t%s\tii\t%d\t%d\n"
	       , dictName (pp->geneDict, gc->gene)
	       , dictName (pp->runDict, gc->run)
	       , gc->boxCount/720, gc->exonCount/720
	       ) ;
  
  fprintf (stderr, "%s: stop geneCounts export total count %ld\n", timeBufShowNow (tBuf), nnn) ;
  ac_free (h) ;
  return nnn ;
} /* wiggleExportGeneCounts */

/**************************************************************/

static void wiggleExportWiggleStats (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  int chromMax = dictMax (pp->bbG.dict) + 1 ;
  long int nnn = 0 ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".wiggleCumuls.tsf", pp->gzo, h) ;
  aceOutDate (ao, "##", "wiggleCumuls in Million Bases") ;
  aceOutDate (ao, "#", "Target\tRun\tFormat\tCumul\tExonic\tIntronic\tIntergenic in Million Bases") ;

  /* export the wiggle cumuls per run. chromosome, strand */
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      int run = nw / (2 * chromMax) ;
      int chrom = (nw % (2 * chromMax)) ;
      char strand = ( nw & 0x1) ? 'r' : 'f' ;
      long int cumul = array (pp->wiggleCumuls, nw, long int) ;
      long int exonic = array (pp->exonics, nw, long int) ;
      long int intronic = array (pp->intronics, nw, long int) ;
      long int intergenic = array (pp->intergenics, nw, long int) ;

      if (cumul)
	{
	  const char *chromNam = dictName (pp->bbG.dict, chrom >> 1) ;
	  const char *runNam = dictName (pp->runDict, run) ;

	  RunSTAT *rc = arrayp (pp->runStats, run, RunSTAT) ;
	  rc->wiggleCumul += cumul ;
	  rc->exonic += exonic ;
	  rc->intronic += intronic ;
	  rc->intergenic += intergenic ;
	  nnn += cumul ;
	  aceOutf (ao, "%s.%c\t%s\tiiii\t%ld\t%ld\t%ld\t%ld\n"
		   , chromNam, strand
		   , runNam
		   , cumul / 720 
		   , exonic / 720 , intronic / 720 , intergenic / 720 
		   ) ;	       
	}
    }

  /* export the cumul per run and globally */
  pp->wiggleCumul = 0 ;
  pp->exonic = 0 ;
  pp->intronic = 0 ;
  pp->intergenic = 0 ;
  for (int run = 1 ; run <= dictMax (pp->runDict) ; run++)
    {
      RunSTAT *rc = arrayp (pp->runStats, run, RunSTAT) ;

      rc->wiggleCumul /= 720 ;
      rc->exonic /= 720 ;
      rc->intronic /= 720 ;
      rc->intergenic /= 720 ;
      
      pp->wiggleCumul += rc->wiggleCumul ;
      pp->exonic += rc->exonic ;
      pp->intronic += rc->intronic ;
      pp->intergenic += rc->intergenic ;

      const char *runNam = dictName (pp->runDict, run) ;
      
      aceOutf (ao, "%s\t%s\tiiii\t%ld\t%ld\t%ld\t%ld\n"
	       , "Any"
	       , runNam
	       , rc->wiggleCumul
	       , rc->exonic
	       , rc->intronic
	       , rc->intergenic
	       ) ;
    }
  
  ac_free (h) ;
  return ;
} /* wiggleExportWiggleStats */

/**************************************************************/
/**************************************************************/

void saWiggleExport (PP *pp, int nAgents)
{
  AC_HANDLE h = ac_new_handle () ;
  int wMax = arrayMax (pp->wiggles) ;
  BOOL debug = FALSE ;
  char tBuf[25] ;

  pp->wiggleCumuls = arrayHandleCreate (wMax, long int, h) ;
  pp->exonics = arrayHandleCreate (wMax, long int, h) ;
  pp->intronics = arrayHandleCreate (wMax, long int, h) ;
  pp->intergenics = arrayHandleCreate (wMax, long int, h) ;
      
  if (pp->geneBoxes)
    {
      pp->geneExonCounts = arrayHandleCreate (wMax, Array, h) ;
      pp->geneBoxCounts = arrayHandleCreate (wMax, Array, h) ;
    }

  int k = 0, n = 0 ;
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      Array wig = arr (pp->wiggles, nw, Array) ;
      if (wig)
	n++ ;
    }
  fprintf (stderr, "%s: start exportation of  %d wiggles\n", timeBufShowNow (tBuf), n) ;
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

  k = n = 0 ;
  for (int nw = 0 ; nw < wMax ; nw++)
    {
      Array wig = arr (pp->wiggles, nw, Array) ;
      if (wig)
	{
	  if (pp->geneBoxCounts)
	    {
	      array (pp->geneBoxCounts, nw, Array) = arrayHandleCreate (2048, dictMax (pp->geneDict) + 1, h) ;
	      array (pp->geneExonCounts, nw, Array) = arrayHandleCreate (2048, dictMax (pp->geneDict) + 1, h) ;
	    }
	  channelPut (pp->wwChan, &nw, int) ;
	}
    }
  channelClose (pp->wwChan) ;

  /* synchronize */
  while (channelGet (pp->wwDoneChan, &k, int))
    n++ ;

  fprintf (stderr, "%s: stop wiggle export\n", timeBufShowNow (tBuf)) ;

  if (pp->geneBoxes)
    wiggleExportGeneCounts (pp) ;

  wiggleExportWiggleStats (pp) ;
  
  ac_free (h) ;
  return ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
 
