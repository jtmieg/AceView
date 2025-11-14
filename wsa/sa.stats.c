/*
 * sa.stats.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module registers and exports detailled run statistics in tsf format
 */

#include "sa.h"

/**************************************************************/

void saCpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, long int n)
{
  int s = arrayMax (cpuStats) ;
  CpuSTAT *bs = arrayp (cpuStats, s, CpuSTAT) ;
  bs->nB++ ;
  strncpy (bs->nam, nam, 30) ;
  bs->agent = agent ;
  bs->n += n ;
  bs->tA += t2 - t1 ;  /* active time */
  return ;
} /* saCpuStatRegister */

/**************************************************************/

void saCpuStatCumulate (Array aa, Array a)
{
  int i = arrayMax (aa) ;
  int j = arrayMax (a) ;
  for (int k = 0 ; k < j ; k++)
    array (aa, i + k, CpuSTAT) = array (a, k, CpuSTAT) ;
  return ;
} /* saCpuStatCumulate */

/**************************************************************/

static int cpuStatsOrder (const void *va, const void *vb)
{
  const CpuSTAT *up = va, *vp = vb ;
  int n = 0 ;

  n = strcmp (up->nam, vp->nam) ; if (n) return n ;
  n = up->agent - vp->agent ; if (n) return n ; 
  return 0 ;
} /* cpuStatsOrder */

/**************************************************************/

void saCpuStatExport (const PP *pp, Array stats)
{
  int i, j, iMax = arrayMax (stats) ;
  CpuSTAT *bt, *bs = arrayp (stats, 0, CpuSTAT) ;

  arraySort (stats, cpuStatsOrder) ;
  printf ("\n# Action\tAgent\tnB\tn\ttime (s)") ;
  if (0)
    {  /* details per agent, only useful for debugging and optimizing */
      for (i = 0 ; i < iMax ; i++, bs++)
	{
	  if (! bs->nB)
	    continue ;
	  for (j = i + 1, bt = bs + 1 ; j < iMax ; j++, bt++)
	    {
	      if (strcmp (bs->nam, bt->nam) || bs->agent != bt->agent)
		break ;
	      bs->nB += bt->nB ; bs->n += bt->n, bs->tA += bt->tA ; bt->nB = 0 ;
	    }
	  printf ("\n%s\t%d\t%d\t%ld\t%.0f", bs->nam, bs->agent, bs->nB, bs->n, bs->tA * 1.0 / CLOCKS_PER_SEC) ;
	}
      printf ("\n") ;
    }
  bs = arrayp (stats, 0, CpuSTAT) ;
  for (i = 0 ; i < iMax ; i++, bs++)
    {
      if (! bs->nB)
	continue ;
      for (j = i + 1, bt = bs + 1 ; j < iMax ; j++, bt++)
	{
	  if (strcmp (bs->nam, bt->nam))
	    break ;
	  if (bt->nB)
	    { bs->nB += bt->nB ; bs->n += bt->n, bs->tA += bt->tA ; bt->nB = 0 ; }
	}
      printf ("\n%s\t%d\t%d\t%ld\t%.0f", bs->nam, bs->agent, bs->nB, bs->n, bs->tA  * 1.0 / CLOCKS_PER_SEC) ;
    }
  printf ("\n") ;
  return ;
} /* saCpuStatExport */

/**************************************************************/
/* cumulate int global runStats the content of bb->runStats */
void saRunStatsCumulate (int run, Array aa, RunSTAT *vp)
{
  RunSTAT *up = arrayp (aa, run, RunSTAT) ;
    
  up->run = run ;
  up->nPairs += vp->nPairs ;
  up->nCompatiblePairs += vp->nCompatiblePairs ;
  up->nCirclePairs += vp->nCirclePairs ;
  up->nOrphans += vp->nOrphans ;
  up->n2ChromsPairs += vp->n2ChromsPairs ;
  up->nAlignedPairs += vp->nAlignedPairs ;
  up->nReads += vp->nReads ;
  up->nBase1 += vp->nBase1 ;
  up->nBase2 += vp->nBase2 ;
  up->nIntronSupports += vp->nIntronSupports ;
  up->nIntronSupportPlus += vp->nIntronSupportPlus ;
  up->nIntronSupportMinus += vp->nIntronSupportMinus ;
  for (int i = 0 ; i < 5 ; i++)
    up->NATGC[i] += vp->NATGC[i] ;
  for (int i = 0 ; i < 11 ; i++)
    up->nMultiAligned[i] += vp->nMultiAligned[i] ;
  for (int i = 0 ; i < 256 ; i++)
    up->nAlignedPerTargetClass[i] += vp->nAlignedPerTargetClass[i] ;
  up->nPerfectReads += vp->nPerfectReads ;
  up->nAlignments += vp->nAlignments ;
  up->nBaseAligned1 += vp->nBaseAligned1 ;
  up->nBaseAligned2 += vp->nBaseAligned2 ;

  up->nErr += vp->nErr ;
  if (vp->errors)
    for (int i = 0 ; i < arrayMax (vp->errors) ; i++)
      array (up->errors, i, int) += array (vp->errors, i, int) ;
  for (int i = 0 ; i < 256 ; i++)
    {
      up->GF[i] += vp->GF[i] ;
      up->GR[i] += vp->GR[i] ;
    }
  return ;
} /* saRunStatsCumulate */

/**************************************************************/

void saRunStatExport (const PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".runStats.tsf", 0, h) ;
  int run, runMax = arrayMax (runStats) ;
  RunSTAT *up ;
  
  aceOutDate (ao, "##", "Run statistics") ;
  aceOut (ao, "#") ;
  for (run = 1 ; run < runMax ; run++)
    aceOutf (ao, "%c%s", run == 1 ? ' ' : '\t', dictName (pp->runDict, run)) ;
  aceOutf (ao, "\nPairs") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nPairs) ;
  aceOutf (ao, "\nAlignedPairs") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nAlignedPairs) ;
  aceOutf (ao, "\nCompatiblePairs") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nCompatiblePairs) ;
  aceOutf (ao, "\nCirclePairs") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nCirclePairs) ;
  aceOutf (ao, "\nReads") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nReads) ;
  aceOutf (ao, "\nAligned_Reads") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nMultiAligned[0]) ;
  aceOutf (ao, "\nMissmatches") ;
  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
    aceOutf (ao, "\t%d", up->nErr) ;
  for (int ii = 1 ; ii < 5 ; ii++)
    {
      aceOutf (ao, "\n%d Alis", ii) ;
      for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
	aceOutf (ao, "\t%d", up->nMultiAligned[ii]) ;
    }
  for (int ii = 1 ; ii < 256 ; ii++)
    {
      up = arrp (runStats, 0, RunSTAT) ;
      if (up->nAlignedPerTargetClass[ii])
	{
	  aceOutf (ao, "\nAligned in class %c", ii) ;
	  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
	    aceOutf (ao, "\t%d", up->nAlignedPerTargetClass[ii]) ;
	}
    }
  for (int ii = 1 ; ii < 256 ; ii++)
    {
      up = arrp (runStats, 0, RunSTAT) ;
      if (up->nAlignedPerTargetClass[ii])
	{
	  aceOutf (ao, "\nStranding in class %c", ii) ;
	  for (run = 1,  up = arrp (runStats, run, RunSTAT) ; run < runMax ; up++, run++)
	    {
	      int f = up->GF[ii] ;
	      int r = up->GR[ii] ;
	      int t = f + r ;
	      if (t) aceOutf (ao, "\t%.3f"
			      , 100.0 * f/t
			      ) ;
	    }
	}
    }
  aceOut (ao, "\n") ;
  
  ac_free (h) ;
} /* saRunStatExport */

