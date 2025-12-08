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

/*************************************************************************************/

static void samStatsExports (const PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *METHOD = pp->method ? pp->method : "01_SortAlign" ;
  Array aa = runStats ;
  RunSTAT *s0 = arrayp (aa, 0, RunSTAT) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".s2g.samStats", 0, h) ;	
  const char *run = pp->runName ? pp->runName : "xxx" ;
  
  long int nRawReads = pp->nRawReads ? pp->nRawReads : s0->nReads ;
  long int nRawBases = pp->nRawBases ? pp->nRawBases : s0->nBase1 + s0->nBase2 ;
  
  aceOutf (ao, "%s\t%s\tnPairs\t%ld\n", run, METHOD, s0->nPairs) ;
  aceOutf (ao, "%s\t%s\tnPairsAligned\t%ld\t%.2f%%\n", run, METHOD, s0->nPairsAligned, (100.0 * s0->nPairsAligned)/(nRawReads/2 + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCompatiblePairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCompatiblePairs, (100.0 * s0->nCompatiblePairs)/(s0->nPairsAligned + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCirclePairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCirclePairs, (100.0 * s0->nCirclePairs)/(s0->nPairsAligned + .000001)) ;
  aceOutf (ao, "%s\t%s\tnIncompatiblePairs\t%ld\t%.2f%%\n", run, METHOD, s0->nIncompatiblePairs, (100.0 * s0->nIncompatiblePairs)/(s0->nPairsAligned + .000001)) ;
  
  
  aceOutf (ao, "\n%s\t%s\tnRawReads\t%ld\n", run, METHOD, nRawReads) ;
  aceOutf (ao, "%s\t%s\tnReads\t%ld\n", run, METHOD, s0->nReads) ;  
  aceOutf (ao, "%s\t%s\tnPerfectReads\t%ld\t%.2f%%\n", run, METHOD, s0->nPerfectReads, (100.0 * s0->nPerfectReads)/(nRawReads + .000001)) ;
  aceOutf (ao, "%s\t%s\tnAlignedReads\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0], (100.0 * s0->nMultiAligned[0])/(nRawReads + .000001)) ;
  
  aceOutf (ao, "\n%s\t%s\tnAlignments\t%ld\t%.2f per aligned read\n", run, METHOD, s0->nAlignments, (1.0 * s0->nAlignments)/(s0->nMultiAligned[0] + .000001)) ;
  aceOutf (ao, "%s\t%s\tnMultiAligned\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0] - s0->nMultiAligned[1], (100.0 * (s0->nMultiAligned[0] - s0->nMultiAligned[1]))/(s0->nReads + .000001)) ;

  aceOutf (ao, "\n%s\t%s\tnRawBases\t%ld\n", run, METHOD, nRawBases) ;
  aceOutf (ao, "%s\t%s\tnBases\t%ld\n", run, METHOD, s0->nBase1 + s0->nBase2) ;
  aceOutf (ao, "%s\t%s\tnAlignedBases\t%ld\t%.2f%%\n", run, METHOD, s0->nBaseAligned1 + s0->nBaseAligned2, 100.0 * (s0->nBaseAligned1 + s0->nBaseAligned2) / (nRawBases + .000001)) ;
  aceOutf (ao, "%s\t%s\tnErrors\t%ld\t%.6f%%\n", run, METHOD, s0->nErr, (100.0 * s0->nErr)/(s0->nBaseAligned1 + s0->nBaseAligned2 + 0.00000001)) ;

  aceOutf (ao, "\n%s\t%s\tnSupportedIntrons\t%ld\n", run, METHOD, arrayMax (pp->confirmedIntrons)) ;
  aceOutf (ao, "%s\t%s\tnIntronSupports\t%ld\t%ld\t%ld\t%.3f\n", run, METHOD
	   , s0->nIntronSupportPlus + s0->nIntronSupportMinus
	   , s0->nIntronSupportPlus
	   , s0->nIntronSupportMinus
	   , s0->intronStranding
	   ) ;
  
  long int nUnaligned = nRawReads - s0->nMultiAligned[0] ;
  aceOutf (ao, "\n%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD, 0
	   , nUnaligned
	   , 100.0 * nUnaligned / (nRawReads + .000001)
	   ) ;
  for (int j = 1 ; j < 11 ; j++)
    if (s0->nMultiAligned[j])
      {
	aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD
		 , j, s0->nMultiAligned[j]
		 , 100.0 * s0->nMultiAligned[j] / (nRawReads + .000001)
		 ) ;
      }
  long int verif = 0 ;
  for (int j = 1 ; j < 11 ; j++)
    verif += s0->nMultiAligned[j] ;
  aceOutf (ao, "nReads = %ld , sum of multiAli = %ld, verif = %ld\n", nRawReads, verif, nRawReads - verif) ; 

  ac_free (h) ;
  return ;
} /* saSamStatsExport */

/**************************************************************/
/**************************************************************/
/* cumulate int global runStats the content of bb->runStats */
void saRunStatsCumulate (int run, PP *pp, BB *bb)
{
  Array aa = pp->runStats ;
  RunSTAT *up = arrayp (aa, run, RunSTAT) ;
  RunSTAT *vp = &(bb->runStat) ;
      
  up->run = run ;
  up->nPairs += vp->nPairs ;
  up->nCompatiblePairs += vp->nCompatiblePairs ;
  up->nIncompatiblePairs += vp->nIncompatiblePairs ;
  up->nCirclePairs += vp->nCirclePairs ;
  up->nOrphans += vp->nOrphans ;
  up->n2ChromsPairs += vp->n2ChromsPairs ;
  up->nPairsAligned += vp->nPairsAligned ;
  up->nReads += vp->nReads ;
  up->nBase1 += vp->nBase1 ;
  up->nBase2 += vp->nBase2 ;
  up->nSupportedIntrons += vp->nSupportedIntrons ;
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
  up->nSupportedIntrons = saSupportedIntrons (pp, run) ;
  
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
  
  aceOutDate (ao, "##", "Run statistics") ;
  for (run = 1 ; run < runMax ; run++)
    {
      RunSTAT *up = arrp (runStats, run, RunSTAT) ;

      if (up->nReads)
	{
	  const char *runNam = dictName (pp->runDict, run) ;
	  
	  aceOutf (ao, "%s\tPairs\ti\t%ld\n", runNam, up->nPairs) ;
	  aceOutf (ao, "%s\tAlignedPairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nPairsAligned
		   , 100.0 * up->nPairsAligned / (.000001 + up->nPairs)
		   ) ;
	  aceOutf (ao, "%s\tCompatiblePairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nCompatiblePairs
		   , 100.0 * up->nCompatiblePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tCirclePairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nCirclePairs
		   , 100.0 * up->nCirclePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tIncompatiblePairs\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nIncompatiblePairs
		   , 100.0 * up->nIncompatiblePairs / (.000001 + up->nPairsAligned)
		   ) ;
	  aceOutf (ao, "%s\tReads\ti\t%ld\n", runNam, up->nReads) ;
	  aceOutf (ao, "%s\tAligned_Reads\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nMultiAligned[0]
		   , 100.0 * up->nMultiAligned[0] / (.000001 + up->nReads)
		   ) ;
	  aceOutf (ao, "%s\tPerfect_Reads\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nPerfectReads
		   , 100.0 * up->nPerfectReads / (.000001 + up->nMultiAligned[0])
		   ) ;

	  aceOutf (ao, "%s\tBases\tiii\t%ld\t%ld\t%ld\n"
		   , runNam
		   , up->nBase1 + up->nBase2
		   , up->nBase1
		   , up->nBase2
		   ) ;
	  aceOutf (ao, "%s\tnAligned_Bases\tififif\t%ld\t%.3f\t%ld\t%.3f\t%ld\t%.3f\n"
		   , runNam
		   , up->nBaseAligned1 + up->nBaseAligned2
		   , 100.0*(up->nBaseAligned1 + up->nBaseAligned2)/(.0001 + up->nBase1 + up->nBase2)
		   , up->nBaseAligned1
		   , 100.0 * up->nBaseAligned1/(.000001 + up->nBase1)
		   , up->nBaseAligned2
		   , 100.0 * up->nBaseAligned2/(.000001 + up->nBase2)
		   ) ;

	  up->intergenic = up->wiggleCumul - up->exonic - up->intronic ;
	  aceOutf (ao, "%s\tnExonic_intronic_intergenic_Bases\tififif\t%ld\t%.3f\t%ld\t%.3f\t%ld\t%.3f\n"
		   , runNam
		   , up->exonic
		   , 100.0*(up->exonic)/(.000001 + up->exonic + up->intronic + up->intergenic)
		   , up->intronic
		   , 100.0*(up->intronic)/(.000001 + up->exonic + up->intronic + up->intergenic)
		   , up->intergenic
		   , 100.0*(up->intergenic)/(.000001 + up->exonic + up->intronic + up->intergenic)
		   ) ;

	  aceOutf (ao, "%s\tATGCN\tiiiii\t%ld\t%ld\t%ld\t%ld\n"
		   , runNam
		   , up->NATGC[0]
		   , up->NATGC[1]
		   , up->NATGC[2]
		   , up->NATGC[3]
		   , up->NATGC[4]
		   ) ;

	  aceOutf (ao, "%s\tRead_length\tiii\t%.0f\t%0.f\t%.0f\n"
		   , runNam
		   , (up->nBase1 + up->nBase2) / (- .00000001 + up->nReads) 
		   , (up->nBase1) / (- .00000001 + (up->nPairs ? up->nPairs : up->nReads) )
		   , (up->nBase2) / (- .00000001 + (up->nPairs ? up->nPairs : up->nReads) )
		   ) ;

	  aceOutf (ao, "%s\twiggleCumul\ti\t%ld\n", runNam, up->wiggleCumul) ;
	  aceOutf (ao, "%s\tMissmatches\tif\t%ld\t%.3f\n"
		   , runNam
		   , up->nErr
		   , 100.0 * up->nErr /(.0001 + up->nBaseAligned1 + up->nBaseAligned2)
		   ) ;
	  aceOutf (ao, "%s\tIntronSupport\tiiif\t%ld\t%ld\t%ld\t%.2f%%\n"
		   , runNam
		   , up->nIntronSupportPlus + up->nIntronSupportMinus
		   , up->nIntronSupportPlus
		   , up->nIntronSupportMinus
		   , up->intronStranding
		   ) ;
	  aceOutf (ao, "%s\tSupportedIntrons\ti\t%ld\n"
		   , runNam
		   , up->nSupportedIntrons
		   ) ;

	  for (int ii = 1 ; ii < 2 ; ii++)
	    aceOutf (ao, "%s\tReads_Aligned_once\tif\t%ld\t%.3f\n"
		     , runNam
		     , up->nMultiAligned[ii]
		     , 100.0 * up->nMultiAligned[ii]/(.000001 + up->nMultiAligned[0]) 
		     ) ;
	  for (int ii = 2 ; ii < 11 ; ii++)
	    if (up->nMultiAligned[ii]) 
	      aceOutf (ao, "%s\tReads_Multi_aligned__%d\tif\t%ld\t%.3f\n"
		       , runNam, ii
		       , up->nMultiAligned[ii]
		       , 100.0 * up->nMultiAligned[ii]/(.000001 + up->nMultiAligned[0])
		       ) ;
	  
	  for (int ii = 1 ; ii < 256 ; ii++)
	    {
	      if (up->nAlignedPerTargetClass[ii])
		aceOutf (ao, "%s\tAligned_in_class_%c\tif\t%ld\t%.3f\n"
			 , runNam, ii
			 , up->nAlignedPerTargetClass[ii]
			 , 100.0 * up->nAlignedPerTargetClass[ii]/(.000001 + up->nMultiAligned[0]) 
			 ) ;
	    }
	  for (int ii = 1 ; ii < 256 ; ii++)
	    {
	      if (up->nAlignedPerTargetClass[ii])
		{
		  int f = up->GF[ii] ;
		  int r = up->GR[ii] ;
		  int t = f + r ;
		  if (t)		  
		    aceOutf (ao, "%s\tStranding_in_class_%c\tf\t%.3f\n", runNam, ii
			     , 100.0 * f/t
			     ) ;
		}
	    }
	}
    }
  ac_free (h) ;

  samStatsExports (pp, runStats) ;

  return ;    
} /* saRunStatExport */

/**************************************************************/
