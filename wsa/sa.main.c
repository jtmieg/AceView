/*
 * sa.main.c : sortalign RNA aligner

 * Created April 18, 2025
 * In collaboration with Greg Boratyn, NCBI

 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality

 * The algorithm extract words from the union of the targets and from the reads creating 2 long tables
 * Both tables are sorted
 *  sorting the targets is a preprocessing, so the time does not count
 *  sorting the reads is done datablock by datablock
 * The two sorted tables are 'merged' which is local in memory
 *  exporting a sequential table.
 * That table is sorted per read
 *  alignments are extended using the magicblast jumper algorithm
 * Finally auxiliary tables, like introns, or wiggles, are exported on demand

 * The code is heavilly parallelized using the channel paradigm borrowed from the Go language
 */

/*
  #define ARRAY_CHECK
  #define MALLOC_CHECK
*/
/* Examples: try   sortalign -h */


#include "sa.h"

static int NN = 1 ;

static void showHitsDo (HIT *hit, long int iMax) ;


/**************************************************************/
/************** utilities *************************************/
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

static void cpuStatCumulate (Array aa, Array a)
{
  int i = arrayMax (aa) ;
  int j = arrayMax (a) ;
  for (int k = 0 ; k < j ; k++)
    array (aa, i + k, CpuSTAT) = array (a, k, CpuSTAT) ;
  return ;
} /* cpuStatCumulate */

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

static void cpuStatExport (const PP *pp, Array stats)
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
} /* cpuStatExport */

/**************************************************************/
/* cumulate int global runStats the content of bb->runStats */
static void runStatsCumulate (int run, Array aa, RunSTAT *vp)
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
} /* runStatsCumulate */

/**************************************************************/

static void runStatExport (const PP *pp, Array runStats)
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
} /* runStatExport */

/**************************************************************/

static void showCountChroms (Array countChroms)
{
  int ii, iMax = countChroms ? arrayMax (countChroms) : 0 ;
  for (ii = 0 ; ii < iMax && ii < 3000 ; ii++)
    {
      COUNTCHROM *zp = arrp (countChroms, ii, COUNTCHROM) ;
      if (zp->chrom)
	printf (".. %d:  weight=%.2f\tseeds=%d\t1:%d/2:%d/4:%d/8:%d/16:%d/32:%d\tindex %d %d\tchrom=%d\tpos=%u %u\tda=%u\tx= %u %u\n"
		, ii, zp->weight, zp->seeds
		, zp->seed1
		, zp->seed2
		, zp->seed4
		, zp->seed8
		, zp->seed16
		, zp->seed32
		, zp->i1, zp->i2
		, zp->chrom
		, zp->a1, zp->a2
		, zp->a2 > zp->a1 ? zp->a2 - zp->a1 : zp->a1 - zp->a2
		, zp->x1, zp->x2
		) ;
    }
} /* showCountChroms */

/**************************************************************/

static void showHits (BigArray hits)
{
  long int ii, iMax = hits ? bigArrayMax (hits) : 0 ;
  int n = NTARGETREPEATBITS ;
  int mask = (1 << n) - 1 ;
  for (ii = 0 ; ii < iMax && ii < 2700 ; ii++)
    {
      HIT *hit = bigArrp (hits, ii, HIT) ;
      printf (".. %ld:  r=%u\t%u\t%u\t%u\tmult %d %s\n"
	      , ii, hit->read, hit->chrom, hit->a1, hit->x1 >> (n+3), hit->x1 & mask
	      , (hit->x1 >> n) & 0x7 ? "intron" : ""
	      ) ;
    }
  printf ("......... max %ld\n", iMax) ;
} /* showHits */

/**************************************************************/

static void showHitsDo (HIT *hit, long int iMax)
{
  if (hit)
    for (long int ii = 0 ; ii < iMax && ii < 30 ; ii++, hit++)
      {
	printf (".. r=%u\t%u\t%u\t%u\n"
		, hit->read, hit->chrom, hit->a1, hit->x1
		) ;
      }
  printf (".........\n") ;
} /* showHitsDo */

/**************************************************************/

static void showAli (Array aligns)
{
  if (aligns)
    {
      long int ii, iMax = arrayMax (aligns) ;
      
      for (ii = 0 ; ii < iMax && ii < 100 ; ii++)
	{
	  ALIGN *ali = arrp (aligns, ii, ALIGN) ;
	  if (ali->score)
	    printf (".. chain %d\tchainScore %d\tchainErr %d chainAli %d\tscore %d\tr=%u : %d:%d/%d\tchr=%u : %d:%d/%d\tnErr %d\t %d::p %d n %d\tdonor %d acc %d  a1-x1=%d\n"
		    , ali->chain, ali->chainScore
		    , ali->chainErr, ali->chainAli
		    , ali->score
		    , ali->read, ali->x0, ali->x1, ali->x2
		    , ali->chrom, ali->a0, ali->a1, ali->a2
		    , ali->nErr
		    , ali->id, ali->previous, ali->next
		    , ali->donor, ali->acceptor
		    , ali->a1 - ali->x1
		    ) ;
	}
      printf (".........\n") ;
    }
  return ;
} /* showAli */

/**************************************************************/
/**************************************************************/
/* cumulate the number of BB to be analyzed */  
static void npCounter (const void *vp)
{
  const PP *pp = vp ;
  int nPuts = 0, nn = 0, n = 0 ;
  char tBuf[25] ;

  int nf = pp->nFiles ;
  
  while (n < nf)
    {
      channelGet (pp->npChan, &nPuts, int) ;
      nn += nPuts ;
      printf ("--- %s: npCounter Parsed %d/%d files\t%d blocs\tcumul %d blocks to be analysed\n",
	      timeBufShowNow (tBuf), ++n, nf, nPuts, nn) ;
    }
  printf ("--- %s: npCounter Processed %d files, closing pcChan at %d blocks\n", timeBufShowNow (tBuf), nf, nn) ;
  channelCloseSource  (pp->plChan) ;
  return ;
} /* npCounter */

/**************************************************************/
/**************************************************************/

static void loadRegulator (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  int nn = 0, nMax = pp->nBlocks ;
  
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->plChan, &bb, BB))
    {
      nn++ ;
      while  (nn > nMax)
	{
	  int nd = channelCount (pp->doneChan) ;
	  nMax = pp->nBlocks - nd ;
	  if (nn <= nMax)
	    break ;
	  sleep (1) ;
	}
      channelPut (pp->lcChan, &bb, BB) ;
    }
  channelClose (pp->lcChan) ; /* there is a single instance of the load regulator */
  return ;
} /* load regulator */

/**************************************************************/
/**************************************************************/
/**************************************************************/

static void readParser (const void *vp)
{
  const PP *pp = vp ;
  RC rc ;
  
  while (channelGet (pp->fpChan, &rc, RC))
    saSequenceParse (pp, &rc, 0, 0, 0) ;
  return ;
} /* Readparser */

/**************************************************************/

#ifdef JUNK
include <sys/resource.h>

void* f(void* arg) {
    struct rusage usage_start, usage_end;
    // Get resource usage at start
    getrusage(RUSAGE_THREAD, &usage_start);
    
    // Simulate memory-intensive work
    int size = 1000000;
    char* buffer = malloc(size);
    for (int i = 0; i < size; i += 4096) { // Touch pages
        buffer[i] = 1;
    }
    
    // Get resource usage at end
    getrusage(RUSAGE_THREAD, &usage_end);
    
    // Calculate page faults
    long minor_faults = usage_end.ru_minflt - usage_start.ru_minflt;
    long major_faults = usage_end.ru_majflt - usage_start.ru_majflt;
    
    printf("Thread %ld: Minor faults: %ld, Major faults: %ld\n", 
           (long)arg, minor_faults, major_faults);
    free(buffer);
    return NULL;
}
#endif

/**************************************************************/

#ifndef YANN
static void codeWords (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
  

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->lcChan, &bb, BB))
    {
      long int nn= 0 ;

      if (pp->debug) printf ("+++ %s: Start code words\n", timeBufShowNow (tBuf)) ;

      t1 = clock () ;
      saSequenceParseGzBuffer (pp, &bb) ;
      saCodeSequenceSeeds (pp, &bb, pp->iStep, FALSE) ;
      t2 = clock () ;

      for (int i = 0 ; i < NN ; i++)
	nn += bigArrayMax (bb.cwsN[i]) ;
      nnn += nn ;
      saCpuStatRegister ("3.CodeWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
      if (pp->debug) printf ("--- %s: Stop code words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;

      channelPut (pp->csChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close csChan at %d,  coded %ld words\n", n, nnn) ;
  channelCloseSource (pp->csChan) ;

  return ;
} /* codeWords */
#endif

/**************************************************************/
/**************************************************************/

#ifndef YANN
static void sortWords (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
  int k ;
  
  memset (&bb, 0, sizeof (BB)) ;
  while ((k = channelGet (pp->csChan, &bb, BB)))
    {
      if (k < 0)
	sleep (1) ;
      else
	{
	  long int nn = 0 ;
	  if (pp->debug) printf ("+++ %s: Start sort words\n", timeBufShowNow (tBuf)) ;
	  
	  t1 = clock () ;
	  for (int k = 0 ; k < NN ; k++)
	    if (bb.cwsN[k])
	      {
		saSort (bb.cwsN[k], 1) ; /* cwOrder */
		nn += bigArrayMax (bb.cwsN[k]) ;
	      }
	  t2 = clock () ;
	
	  saCpuStatRegister ("4.SortWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
	  if (pp->debug) printf ("--- %s: Stop sort words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;
	  channelPut (pp->smChan, &bb, BB) ;
	}
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close smChan at %d,  sorted %ld words\n", n, nnn) ;
  
  channelCloseSource (pp->smChan) ;

  return ;
} /* sortWords */
#endif

/**************************************************************/

static long int  matchHitsDo (const PP *pp, BB *bbG, BB *bb)
{
  BOOL debug = FALSE ;
  BOOL useIntronSeeds = ! pp->ignoreIntronSeeds ;
  BigArray hitsArray = bigArrayHandleCreate(64, BigArray, bb->h);
  long int nn = 0, k = 0, kkk = 0 ;
  int nHA = 0, hMax = 100000 ;
  const long unsigned int mask26 = (1L << 26) - 1 ;
  BigArray hits = bigArrayHandleCreate(hMax, HIT, bb->h);
  bigArray (hitsArray, nHA++, BigArray) = hits ;
  const int seedLength = pp->seedLength ;
  const int intronBonus = 1 ;
  int absoluteMax = 0x1 << NTARGETREPEATBITS ;
  int absoluteX1Max = 0x1 << ( 31 - 3 - NTARGETREPEATBITS) ;
  int maxTargetRepeats = pp->maxTargetRepeats  ;
  
  for (int kk = 0; kk < NN ; kk++)
    {
      long int i = 0, iMax = bigArrayMax (bbG->cwsN[kk]);
      long int j = 0, jMax = bigArrayMax (bb->cwsN[kk]);


      if (!iMax || !jMax)
	continue ;
      
      unsigned int mask = step1 - 1 ;
      const CW *restrict cw = bigArrp(bbG->cwsN[kk], 0, CW) ;
      const CW *restrict rw = bigArrp(bb->cwsN[kk], 0, CW) ;
      const CW *restrict cw1;
      const CW *restrict cwMax = cw + iMax ;
      HIT *restrict hit;
        
      while  (i < iMax && j < jMax)
	{
	  if (0 && kk == 1 && rw->seed == 185667857)
	    printf("(rw->seed == 185667857)\n") ;
	  if ((i & mask) == 0)
	    {
	      if (rw->seed <= (unsigned int) cw->seed)
		{
		  cw++ ;
		  i++ ;
		  bb->skips0++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->pos)
		{
		  cw += step1 ;
		  i += step1 ;
		  bb->skips1++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->nam)
		{
		  cw += step2 ;
		  i += step2 ;
		  bb->skips2++;
		  continue;
		}
	      else if (rw->seed <= (unsigned int) cw->intron)
		{
		  cw += step3 ;
		  i += step3 ;
		  bb->skips3++;
		  continue;
		}
	      else
		{
		  cw += step4 ;
		  i += step4 ;
		  bb->skips4++;
		  continue;
		}
	    }
	  if (cw->seed < rw->seed)
	    { i++; cw++; bb->skipsNotFound++ ; }
	  else if (cw->seed > rw->seed)
	    { j++ ; rw++ ; }
	  else
	    {
	      long int a1, x1, i1 = i ;
	      bb->skipsFound++ ;
	      int nTargetRepeats = cw->intron ;
	      
	      if (0 &&   /* avoid, this kills the intron seeds */
		  nTargetRepeats > maxTargetRepeats)
		{ j++ ; rw++ ; continue ; }
	      if (nTargetRepeats >= absoluteMax)
		nTargetRepeats = absoluteMax - 1 ; /* we will report absoluteMax (i.e. 31) even is value is higher */

	      for (cw1 = cw ; cw1 < cwMax ; i1++, cw1++)
		{
		  if ((i1 & mask) == 0)
		    continue ;
		  if (cw1->seed != rw->seed)
		    break ;
		  /* success, non intron case */
		  if (((cw1->intron >> 31) & 0x1) == 0x0)
		    {
		      BOOL readUp = rw->nam & 0x1 ;
		      BOOL chromUp = cw1->nam & 0x1 ;
		      int nTargetRepeats = cw1->intron ;
		      
		      if (1 && nTargetRepeats > maxTargetRepeats)
			continue ; 
		      if (0 && useIntronSeeds) continue ;
		      nn++ ;
		      hit = bigArrayp (hits, k++, HIT) ;
		      hit->read = rw->nam >> 1 ;
		      a1 = cw1->pos ;
		      x1 = rw->pos ;
		      chromUp ^= readUp ; /* we want to be on strand plus of the read */
		      readUp = 0 ;
		    
		      if (! chromUp)  /* plus strand of the genome */
			{
			  hit->a1 =
			    a1            /* position of the first base of the seed */
			    - x1          /* locate the virtual position of base 0 of the read */
			    + 1 ;         /* avoid zero */
			  hit->x1 = x1 << 3 ;  /* reserve 3 bits for the intron seeds */
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | nTargetRepeats ;  /* all intron seeds are valuable */			
			}
		      else    /* minus strand of the genome */
			{
			  hit->a1 =
			    a1            /* position of the first base of the seed */
			    + (seedLength - 1)  /* go the last base of the seed */
			    + x1   /* locate the virtual position of base 0 of the read */
			    + 1 ;  /* avoid zero */
			  hit->x1 = x1 << 3 ;  /* reserve 3 bits for the intron seeds */
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | nTargetRepeats ;  /* all intron seeds are valuable */			
			}
		    }
		  else  if (useIntronSeeds)  /* INTRON */
		    {
		      BOOL readUp = rw->nam & 0x1 ;
		      BOOL chromUp = cw1->nam & 0x1 ; 

		      unsigned int z = cw1->intron ;
		      unsigned int isIntronDown = (z >> 28) & 0x4 ;
		      int da1 =  z & 0xf ; /* nb of letters in first exon : 4....11 */
		      int da  =  ((z >> 4) & mask26) ;  /* intron length < 32Mb */

		      chromUp ^= readUp ; /* we want to be on strand plus of the read */
		      readUp = 0 ;
		      hit = 0 ;
		      
		      if (! chromUp)  /* plus strand on the read and on the genome */
			{
			  a1 = cw1->pos ;       /* first base of intron in the genome, in bio coordinates */
		          x1 = rw->pos + da1 ;  /* matching base on the read */
			  
			  /* Create a hit to the last two bases of the donor exon (x1-2 / a1-2) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (0 && isIntronDown == 0) hit->chrom |= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 - 2) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives the intron strand */
			    | 0x1             /* donor site */
			    ;
			  hit->a1 =
			    (a1 - 2)            /* position of the first base of the seed */
			    - (x1 - 2)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			  
			  /* Create a hit to the first two bases of the acceptor exon (x1/ a1 + da) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam & 0xfffffffe ; /* to select plus strand, kill the last bit */
			  if (0 && isIntronDown == 0) hit->chrom |= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives the intron strand */
			    | 0x2             /* acceptor site */
			    ;
			  hit->a1 =
			    (a1 + da)            /* position of the first base of the seed */
			    - (x1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  if (hit->x1 > absoluteX1Max)
			    messcrash ("read coordinate x1=%d too large, please edit the source code", x1) ;
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			}

		     else  /* plus strand on the read and minus strand  on the genome */
			{
			  a1 = cw1->pos ;       /* first base of intron in the genome, in bio coordinates */
		          x1 = rw->pos + (seedLength - da1) - 1 ;  /* matching base on the read */
			  
			  /* Create a hit to the first two bases of the acceptor exon (x1+1,x1+2 / a1-1,a1-2) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (0 && isIntronDown == 0) hit->chrom ^= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 + 1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives th intron strand */
			    | 0x2             /* acceptor site */
			    ;
			  hit->a1 =
			    (a1 - 1)       /* position of the first base of the seed */
			    + (x1 + 1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			  
			  /* Create a hit to the last two bases of the donor exon (x1-1,x1/ a1 + da+1,a1+da) */
			  nn++ ;
			  hit = bigArrayp (hits, k++, HIT) ;
			  hit->read = rw->nam >> 1 ;
			  hit->chrom = cw1->nam | 0x1 ; /* to select minus strand, set the last bit */
			  if (0 && isIntronDown == 0) hit->chrom ^= 0x1 ; /* to cluster on the correct strand of the chroms */
			  hit->x1 =
			    ((x1 - 1) << 3)   /* reserve 3 bits for the intron seeds */
			    | isIntronDown    /* bit 3 gives th intron strand */
			    | 0x1             /* donor site */
			    ;
			  hit->a1 =
			    (a1 + da + 1)            /* position of the first base of the seed */
			    + (x1 - 1)          /* locate the virtual position of base 0 of the read */
			    - intronBonus       /* prefer intron match to exon match */
			    + 1 ;               /* avoid zero */
			  hit->x1 = ( hit->x1 << NTARGETREPEATBITS) | 0x1 ;  /* all intron seeds are valuable */
			}
		    }
		  else
		    continue ;
		  
		  hit->chrom ^= (hit->read & 0x1) ; /* flip chrom for the second read of a pair */
		  if (0 && rw->seed == 168430082)
		    printf("(rw->seed == 168430082)\n") ;
		
		  if (k >= hMax)
		    {
		      bigArrayMax (hits) = k  ;
		      hits = bigArrayHandleCreate(hMax, HIT, bb->h);
		      bigArray (hitsArray, nHA++, BigArray) = hits ;
		      kkk += k ;
		      k = 0 ;
		    }
		}
	      j++ ; rw++ ;
	    }
	}
      bigArrayDestroy (bb->cwsN[kk]) ; /* the read words are no longer needed */
      bigArrayMax (hits) = k ;
      kkk += k ;
    }
  bb->hits = hitsArray ;

  if (debug)
    fprintf (stderr, "..MatchHitsDo found %ld matches\n", kkk) ;

  return nn ;
}  /* matchHitsDo */

/**************************************************************/
/**************************************************************/

#ifndef YANN
static void matchHits (const void *vp)
{
  const PP *pp = vp ;
  BB bb ;
  BB bbG = pp->bbG;
  char tBuf[25] ;
  long int nnn = 0 ;
  
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  /* grab and match a block of reads */
  while (channelGet (pp->smChan, &bb, BB))
    {
      if (bb.length)
	{
	  if (pp->debug) printf ("+++ %s: Start match %ld bases againt %ld target bases\n", timeBufShowNow (tBuf), bb.length, bbG.length) ;


	  t1 = clock () ;
	  long int nn = matchHitsDo (pp, &bbG, &bb) ;
	  if (pp->debug) printf ("--- %s: Stop match hits constructed %ld arrays\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	  t2 = clock () ;

	  nnn += nn ;
	  saCpuStatRegister ("5.MatchHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
	}
      channelPut (pp->moChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close moChan at %d,  found %ld hits\n", n, nnn) ;
  channelCloseSource (pp->moChan) ;

  return ;
} /* matchHits */
#endif

/**************************************************************/
/**************************************************************/
/**************************************************************/

static int alignOrder (const void *va, const void *vb)
{
  const ALIGN *up = va ;
  const ALIGN *vp = vb ;
  int n ;
  n = up->read - vp->read ; if (n) return n ;
  n = up->chainScore - vp->chainScore ; if (n) return -n ;
  n = up->chain - vp->chain ; if (n) return n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->x1 - vp->x1  ; if (n) return n ;
  n = up->x2 - vp->x2  ; if (n) return -n ;
  n = up->nErr - vp->nErr ; if (n) return n ;
  n = up->a1 - vp->a1  ; if (n) return n ;
  n = up->a2 - vp->a2  ; if (n) return n ;
  n = up->nTargetRepeats - vp->nTargetRepeats ; if (n) return 1 ;
  return 0 ;
} /* alignOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static int hitReadPosOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n, n1, n2 ;

  n = ((up->read > vp->read) - (up->read < vp->read)) ; if (n) return n ;
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; if (n) return n ;
  n1 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n2 = vp->a1 + (vp->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n = n1 - n2 ;
  return n ;
} /* hitReadPosOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static int countChromOrder (const void *va, const void *vb)
{
  const COUNTCHROM *up = va ;
  const COUNTCHROM *vp = vb ;
  int n ;
  n = up->weight - vp->weight ; if (n) return -n ;
  n= up->seeds - vp->seeds ; if (n) return -n ;
  n= up->chrom - vp->chrom ; if (n) return n ;
  n= up->a1 < up->a2 ? 1 : -1 ;
  
  return n ;
} /* countChromOrder */

/**************************************************************/

static void sortHitsFuse (const PP *pp, BB *bb)
{
  BigArray aa = bb->hits ;
  long int iMax = aa ? bigArrayMax (aa) : 0 ;
  char tBuf[25] ;
  
  if (pp->debug) printf ("+++ %s: Start sort hits merging %ld arrays\n", timeBufShowNow (tBuf), iMax) ;


  if (iMax == 0)
    bb->hits = 0 ;
  else if (iMax == 1)
    bb->hits = bigArray (aa, 0, BigArray) ;
  else
    {
      BigArray a ;
      long int n = 0 ;
      HIT *up ;
      for (long int i = 0 ; i < iMax ; i++)
	{
	  a = bigArray (aa, i, BigArray) ;
	  n += bigArrayMax (a) ;
	}
      bb->hits = bigArrayHandleCreate (n, HIT, bb->h) ;
      up = bigArrayp (bb->hits, n - 1, HIT) ; /* make room */
      up = bigArrayp (bb->hits, 0, HIT) ;
      for (int i = 0 ; i < iMax ; i++)
	{
	  a = bigArray (aa, i, BigArray) ;
	  n = bigArrayMax (a) ;
	  if (n > 0)
	    {
	      memcpy (up, bigArrp (a, 0, HIT), n * sizeof (HIT)) ;
	      up += n ;
	    }
	  bigArrayDestroy (a) ;
	}
      bigArrayDestroy (aa) ;
    }
  return ;
} /* sortHitsFuse */

/**************************************************************/
/* sort the hits */
#ifndef YANN
static void sortHits (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->moChan, &bb, BB))
    {
      if (pp->align && bb.hits)
	{
	  t1 = clock () ;
	  sortHitsFuse (pp, &bb) ;
	  if (bb.hits)
	    {
	      saSort (bb.hits, 3) ; /* hitPairOrder */
	      t2 = clock () ;

	      long int nn = bigArrayMax (bb.hits) ;
	      nnn += nn ;
	      saCpuStatRegister ("6.SortHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
	      if (pp->debug) printf ("--- %s: Stop sort hits %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	    }
	}
      channelPut (pp->oaChan, &bb, BB) ;
    }

  if (1)
    {
      memset (&bb, 0, sizeof (BB)) ;
      bb.isGenome = TRUE ;
      channelPut (pp->doneChan, &bb, BB) ; /* destroy bbG.cws, all Matches are already computed */ 
    }
  channelCloseSource (pp->oaChan) ;
  channelCloseSource (pp->doneChan) ;
  
  return ;
} /* sortHits */
#endif

/**************************************************************/
/**************************************************************/
#include <stdint.h>

/**************************************************************/

static BOOL alignExtendHit (Array dna, Array dnaG, Array dnaGR, Array err
			    , BOOL isDown, int chromLength
			    , int *a1p, int *a2p, int *x1p, int *x2p
			    , int errCost
			    , BOOL isIntron
			    , int errMax, int minAli
			    )
{
  int nN = 0, dx = 0, nerr = 0 ;
  int x1 = *x1p, x2 = *x2p ;
  int a1 = *a1p, a2 = *a2p ;
  arrayMax (err) = 0 ;

  errMax = isIntron ? 0 : errMax ;
  
  

  if (! isDown)
    {
      Array dummy = dnaG ; dnaG = dnaGR ; dnaGR = dummy ;
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
      *a1p = a1 ; *a2p = a2 ;
    }

  arrayMax (err) = 0 ;
  aceDnaDoubleTrackErrors (dna, &x1, &x2, TRUE /* bio coordinates, extend = TRUE */
			   , dnaG, dnaGR, &a1, &a2
			   , &nN, err, MAXJUMP, errMax, TRUE, 0) ;
  if (x1 > *x1p || x2 < *x2p)
    return FALSE ;
  dx = x2 - x1 + 1 ;
  if (*a1p < *a2p && (a1 > a2 || a1 > *a1p || a2 < *a2p))
    return FALSE ;
  if (*a1p > *a2p && (a1 < a2 || a1 < *a1p || a2 > *a2p))
    return FALSE ;
  /* reclip */
  /* clip left errors */
  nerr = arrayMax (err) ;
  if (nerr)
    {
      int i, j, y1 = x1 - 1, y2 = x2 - 1 ;  /* natural coordinates */
      A_ERR *up = arrp (err, 0, A_ERR) ;
      A_ERR *vp = up ;
      BOOL wonder = TRUE ;
      
      for (i = j = 0 ; i < nerr ; i++, up++)
	{
	  if (wonder && errCost > (up->iShort - y1)) /* clip this error */
	    {
	      x1 = up->iShort + 1 ; /* +1 for bio coordinates */
	      a1 = up->iLong + 1 ;  /* +1 for bio coordinates */
	      switch (up->type)
		{
		case AMBIGUE: break ;
		case INSERTION_TRIPLE: x1 += 3 ; break ;
		case TROU_TRIPLE: a1 += 3 ;break ;
		case INSERTION_DOUBLE: x1 += 2 ; break ;
		case TROU_DOUBLE: a1 += 2 ; break ;
		case INSERTION: x1++ ; break ;
		case  TROU: a1++ ; break ;
		default : a1++ ; x1++ ; break ;
		}
	      y1 = x1 - 1 ; /* natural coordinates */
	      continue ;
	    }
	  else  /* register all remaining errors */
	    {	      
	      wonder = FALSE ;
	      if (j < i) *vp = *up ;
	      j++ ; vp++ ;
	    }
	}
      nerr =  arrayMax (err) = j ; /* some 5' errors are dropped */
	    
      wonder = TRUE ;
      if (nerr)
	for (i = nerr - 1, up = arrp (err, i, A_ERR) ; i >= 0 ; i--, up--)
	  {
	    if (wonder && errCost > y2 - up->iShort) /* clip this error */
	      {
		nerr-- ;
		x2 = up->iShort - 1 + 1 ; /* +1 for bio coordinates */
		a2 = up->iLong  - 1 + 1 ; /* +1 for bio coordinates */
		y2 = x2 - 1 ; /* natural coordinates */
		continue ;
	      }
	  }
      arrayMax(err) = nerr ; /* some 3' errors are dropped */
    }
  
  if (x1 > *x1p || x2 < *x2p)
    return FALSE ;
  if (*a1p < *a2p && (a1 > a2 || a1 > *a1p || a2 < *a2p))
    return FALSE ;
  if (*a1p > *a2p && (a1 < a2 || a1 < *a1p || a2 > *a2p))
    return FALSE ;

  if (! isDown)
    {
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
    }
  
  *x1p = x1 ; *x2p = x2 ;
  *a1p = a1 ; *a2p = a2 ;
  dx = x2 - x1 + 1 ;

  minAli = isIntron ? 8 : 22 ;
  if (dx < minAli)
    return FALSE ;
  if (errMax >= 0 &&  nerr > errMax)
    return FALSE ;
  if (dx < nerr * errCost)
    return FALSE ;

  return TRUE ;
} /* alignExtendHit */

/**************************************************************/
/* merge align->errors into bb->errors */
/* sorry, we use this extremelly obfuscated way to encode the errors
 * rather than a simple clear utilisation of vp->type + actual bases
 * just to maintain a small max value for the variable type.
 * The intention is to help the CPU cache by maintaining a very small aa array
 *
 * There is also a complex way to choose the strand
 * If we use paired end sequencing, we want to switch the strand of read2
 * In addition, if we sequence RNA (as opposed to DNA) and
 * we map antistrand to a known gene, we flip read1 and again read2
 */
static void mergeErrors (Array aa, Array bb, unsigned int flip)
{
  A_ERR *vp = arrp (bb, 0, A_ERR) ;
  int i, iMax = arrayMax (bb) ;

  for (i = 0 ; i < iMax ; i++, vp++)
    {
      int type = 0 ;
      switch (vp->type)
	{
	case ERREUR : /* use 4 right bits type < 15 */
	  switch (vp->baseLong)
	    {
	    case A_: break ;
	    case T_: type = 0xc ;break ;
	    case G_: type = 0x4 ; break ;
	    case C_: type = 0x8 ; break ;
	    }
	  switch (vp->baseShort)
	    {
	    case A_: break ;
	    case T_: type |= 0x3 ;break ;
	    case G_: type |= 0x1 ; break ;
	    case C_: type |= 0x2 ; break ;
	    }
	  type = 0xf & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	  break ;
	case INSERTION:
	  switch (vp->baseShort)
	    {
	    case A_: type = 0x0 ; break ;
	    case T_: type = 0x3 ; break ;
	    case G_: type = 0x1 ; break ;
	    case C_: type = 0x2 ; break ;
	    }
	  type = 0x3 & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	  type |= 0x10 ;
	  break ;
	case TROU:
	  switch (vp->baseLong)
	    {
	    case A_: type = 0x0 ; break ;
	    case T_: type = 0xc ;break ;
	    case G_: type = 0x4 ; break ;
	    case C_: type = 0x8 ; break ;
	    }
	  type = 0xc & (type ^ flip) ; /* flip the last 4 bits : A<->T   G<->C */
	  type |= 0x20 ;
	  break ;
	case INSERTION_DOUBLE:
	  type = 0x32 ;
	  break ;
	case INSERTION_TRIPLE:
	  type = 0x33 ;
	  break ;
	case TROU_DOUBLE:
	  type = 0x38 ;
	  break ;
	case TROU_TRIPLE:
	  type = 0x3c ;
	  break ;
	default:
	  break ;
	}
      if (type)
	(array (aa, type, int)) += 1 ;
    }
  return ;
} /* mergeErrors */

/**************************************************************/

static void runErrorsCumulate (int run, Array runErrors, Array bb)
{
  Array aa = array (runErrors, run, Array) ;
  int iMax = arrayMax (bb) ;
  if (aa && iMax)
    for (int i = 0 ; i < iMax ; i++)
      array (aa, i, int) += array (bb, i, int) ;
} /* runErrorsCumulate */

/**************************************************************/

static void reportRunErrors (const PP *pp, Array runStats, Array runErrors)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".runErrors.tsf", 0, h) ;
  int type, run, runMax = arrayMax (runErrors) ;

  for (run = 1 ; run < runMax ; run++)
    {
      char buf[8] ;
      const char *AGCT = "AGCT" ;
      const char *rNam = dictName (pp->runDict, run) ;      
      RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
      Array errors = array (runErrors, run, Array) ;

      aceOutDate (ao, "##", "Mismatches") ;
      aceOutf (ao, "%s\tAny\ti\t%d\n"
	       , rNam, rs->nErr
	       ) ;

      if (! errors) continue ;
      /* substitutions */
      for (int longBase = 0 ; longBase < 4 ; longBase++)
	{
	  buf[0] = AGCT[longBase] ;
	  buf[1] = '>' ;
	  for (int shortBase = 0 ; shortBase < 4 ; shortBase++)
	    {
	      type = (longBase << 2) | shortBase ;
	      buf[2] = AGCT[shortBase] ;
	      buf[3] = 0 ;

	      aceOutf (ao, "%s\t%s\ti\t%d\n"
		       , rNam, buf, array (errors, type, int)
		       ) ;
	    }
	}
      /* insertions */
      memcpy (buf, "Ins", 3) ;
      for (int shortBase = 0 ; shortBase < 4 ; shortBase++)
	{
	  type = 0x10 + shortBase ;
	  buf[3] = AGCT[shortBase] ;
	  buf[4] = 0 ;
	  aceOutf (ao, "%s\t%s\ti\t%d\n"
		   , rNam, buf, array (errors, type, int)
		   ) ;
	}
      /* deletions */
      memcpy (buf, "Del", 3) ;
      for (int longBase = 0 ; longBase < 4 ; longBase++)
	{
	  type = 0x20  | (longBase << 2) ;
	  buf[3] = AGCT[longBase] ;
	  buf[4] = 0 ;
	  aceOutf (ao, "%s\t%s\ti\t%d\n"
		   , rNam, buf, array (errors, type, int)
		   ) ;
	}
    }
  ac_free (h) ;
} /* runErrorsCumulate */

/**************************************************************/
/**************************************************************/

static void alignFormatLeftOverhang (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR)
{
  int x1 = up->x1 ;

  if (x1 > 1)
    {
      int dx = x1 - 1 ;
      char buf[dx + 1] ;
      if (dx > 30) dx = 30 ;
      memcpy (buf, arrp (dna, x1 - 1 - dx, char), dx) ; buf[dx] = 0 ;
      for (int i = 0 ; i < dx ; i++)
	buf[i] = dnaDecodeChar[(int)buf[i]] ;
      dictAdd (bb->dict, buf, &up->leftOverhang) ;
    }
  return ;
} /* alignFormatLeftOverhang */

/**************************************************************/

static void alignFormatRightOverhang (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR)
{
  int x2 = up->x2 ;
  int ln = arrayMax (dna) ; /* length to align */

  if (x2 < ln)
    {
      int dx = ln - x2 ;
      char buf[dx + 1] ;
      if (dx > 30) dx = 30 ;
      memcpy (buf, arrp (dna, x2, char), dx) ; buf[dx] = 0 ;
      for (int i = 0 ; i < dx ; i++)
	buf[i] = dnaDecodeChar[(int)buf[i]] ;
      dictAdd (bb->dict, buf, &up->rightOverhang) ;
    }
  return ;
} /* alignFormatRightOverhang */

/**************************************************************/

static void alignFormatErrors (const PP *pp, BB *bb, ALIGN *up, Array dna, Array dnaG, Array dnaGR)
{
  A_ERR *ep = arrp (up->errors, 0, A_ERR) ;
  int ii, nerr = arrayMax (up->errors) ;
  vTXT txt1 = bb->txt1 ;
  vTXT txt2 = bb->txt2 ;
  char *sep ;
  int xShort, xLong  ;
  BOOL isUp = (up->a1 > up->a2) ;
  const BOOL debug = FALSE ;
  if (debug)
    aceDnaShowErr (up->errors) ;
  vtxtClear (txt1) ;
  vtxtClear (txt2) ;
  
  for (ii = 0, sep = "" ; ii < nerr ; ii++, ep++, sep = ",")
    {
      xShort = ep->iShort + 1 ;
      xLong = ep->iLong + 1 ;

      if (xShort <= 0 || xLong <= 0 || xShort >= arrayMax (dna))
	continue ;
      switch (ep->type)
	{
	  case TYPE80:
	    {
	      char cc1a, cc2a ;
	      
	      cc1a = arr (dnaG, (isUp ? xLong - 0 : xLong - 2), unsigned char) ;
	      cc2a = arr (dnaG, (isUp ? xLong - 1 : xLong - 1), unsigned char) ;
	      
	      vtxtPrintf (txt1,"%s%d:%c%c>oo"
			  , sep
			  , xShort - 1
			  , isUp ? dnaDecodeChar[(int)cc2a] : dnaDecodeChar[(int)cc1a]
			  , isUp ? dnaDecodeChar[(int)cc1a] : dnaDecodeChar[(int)cc2a]
			  ) ;
	      vtxtPrintf (txt2,"%s%d:%c%c>oo"
			  , sep
			  , isUp ? xLong - 1 : xLong - 1
			  , isUp ? dnaDecodeChar[(int)complementBase[(int)cc1a]] : dnaDecodeChar[(int)cc1a]
			  , isUp ? dnaDecodeChar[(int)complementBase[(int)cc2a]] : dnaDecodeChar[(int)cc2a]
			  ) ;
	      
	    }
	    break ;
	    
	case AMBIGUE:
	case ERREUR:
	  {
	    char ccS, ccL, ccSR, ccLR ;
	    /* int xLongR = arrayMax (dnaG) - xLong + 1 ; */
	    
	    ccS = ep->baseShort ;
	    ccL = ep->baseLong ;
	    
	    ccSR = isUp ? complementBase[(int)ccS] : ccS ; 
	    ccLR = isUp ? complementBase[(int)ccL] : ccL ; 

	    ccS = dnaDecodeChar[(int)ccS] ;
	    ccL = dnaDecodeChar[(int)ccL] ;
	    ccSR = dnaDecodeChar[(int)ccSR] ;
	    ccLR = dnaDecodeChar[(int)ccLR] ;
	    
	    vtxtPrintf(txt1, "%s%d:%c>%c"
		       , sep
		       , xShort 
		       , ccL, ccSR
		       ) ;

	    vtxtPrintf(txt2, "%s%d:%c>%c"
		       , sep
		       , xLong 
		       , ccLR, ccS
		       ) ;
	  }  
	  break ;
	case TROU: 
	  {
	    char *ss = "-", cc1a, cc2a, cc1ac, cc2ac ;
	    
	    cc1a = arr(dnaG, xLong - 1, unsigned char) ;
	    cc2a = isUp ? complementBase[(int)cc1a] : cc1a ; 
	    
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;

	    vtxtPrintf (txt2,"%s%d:%s%c"
			, sep
			, xLong 
			, ss 
			, cc1ac 
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c"
			, sep
			, xShort
			, ss
			, cc2ac
			) ;
	  }
	  break ;

	case TROU_DOUBLE:
	  {
	    char *ss = "--", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;
	    cc1a = arr(dnaG, xLong - 1 , unsigned char) ;
	    cc1b = arr(dnaG, xLong - 1 + 1, unsigned char) ;
	    
	    cc2a = isUp ? complementBase[(int)cc1b] : cc1a ; 
	    cc2b = isUp ? complementBase[(int)cc1a] : cc1b ; 
	    
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc1bc = dnaDecodeChar[(int)cc1b] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    cc2bc = dnaDecodeChar[(int)cc2b] ;
	    
	    vtxtPrintf (txt2, "%s%d:%s%c%c"
			, sep
			, xLong 
			, ss 
			, cc1ac, cc1bc
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c%c"
			, sep
			, xShort + (isUp ? 1 : 0)
			, ss
			, cc2ac, cc2bc
			) ;
	  }
	  break ;
	case TROU_TRIPLE:
	  {
	    char *ss = "---", cc1a, cc1b, cc1c, cc2a, cc2b, cc2c, cc1ac, cc1bc, cc1cc, cc2ac, cc2bc, cc2cc ;
	    cc1a = arr(dnaG, xLong - 1 , unsigned char) ;
	    cc1b = arr(dnaG, xLong - 1 + 1 , unsigned char) ;
	    cc1c = arr(dnaG, xLong - 1 + 2 , unsigned char) ;
	    
	    cc2a = isUp ? complementBase[(int)cc1c] : cc1a ; 
	    cc2b = isUp ? complementBase[(int)cc1b] : cc1b ; 
	    cc2c = isUp ? complementBase[(int)cc1a] : cc1c ; 
	    
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc1bc = dnaDecodeChar[(int)cc1b] ;
	    cc1cc = dnaDecodeChar[(int)cc1c] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    cc2bc = dnaDecodeChar[(int)cc2b] ;
	    cc2cc = dnaDecodeChar[(int)cc2c] ;
	    
	    vtxtPrintf (txt2, "%s%d:%s%c%c%c"
			, sep
			, xLong 
			, ss 
			, cc1ac, cc1bc, cc1cc
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c%c%c"
			, sep
			, xShort + (isUp ? 1 : 0)
			, ss
			, cc2ac, cc2bc, cc2cc
			) ;
	  }
	  break ;
	case INSERTION: 
	  {
	    char *ss = "+", cc1a, cc2a, cc1ac, cc2ac ;

	    cc1a = arr (dna, xShort - 1, char) ;
	    cc2a = isUp ? complementBase[(int)cc1a] : cc1a ; 
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    
	    vtxtPrintf (txt2,"%s%d:%s%c"
			, sep
			, xLong 
			, ss 
			, cc2ac
			  ) ;
	    
	    vtxtPrintf (txt1,"%s%d:%s%c"
			, sep
			, xShort 
			, ss
			, cc1ac
			) ;
	  }
	  break ;

	case INSERTION_DOUBLE:
	  {
	    char *ss = "++", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;

	    cc1a = arr (dna, xShort - 1 + (isUp ? -1 : 0), char) ;
	    cc1b = arr (dna, xShort - 1 + (isUp ? 0 : 1), char) ;

	    cc2a = isUp ? complementBase[(int)cc1b] : cc1a ; 
	    cc2b = isUp ? complementBase[(int)cc1a] : cc1b ; 
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc1bc = dnaDecodeChar[(int)cc1b] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    cc2bc = dnaDecodeChar[(int)cc2b] ;
	    
	    vtxtPrintf (txt2, "%s%d:%s%c%c"
			, sep
			, xLong
			, ss 
			, cc2ac, cc2bc
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c%c"
			, sep
			, xShort + (isUp ? -1 : 0) 
			, ss
			, cc1ac, cc1bc
			) ;
	  }
	  break ;
	  
	case INSERTION_TRIPLE:
	  {
	    char *ss = "+++", cc1a, cc1b, cc1c, cc2a, cc2b, cc2c, cc1ac, cc1bc, cc1cc, cc2ac, cc2bc, cc2cc ;

	    cc1a = arr (dna, xShort - 1 + (isUp ? -2 : 0), char) ;
	    cc1b = arr (dna, xShort - 1 + (isUp ? -1 : 1), char) ;
	    cc1c = arr (dna, xShort - 1 + (isUp ? 0 : 2), char) ;
	    
	    cc2a = isUp ? complementBase[(int)cc1c] : cc1a ; 
	    cc2b = isUp ? complementBase[(int)cc1b] : cc1b ; 
	    cc2c = isUp ? complementBase[(int)cc1a] : cc1c ; 
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc1bc = dnaDecodeChar[(int)cc1b] ;
	    cc1cc = dnaDecodeChar[(int)cc1c] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    cc2bc = dnaDecodeChar[(int)cc2b] ;
	    cc2cc = dnaDecodeChar[(int)cc2c] ;
	    
	    vtxtPrintf (txt2, "%s%d:%s%c%c%c"
			, sep
			, xLong
			, ss 
			, cc2ac, cc2bc, cc2cc
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c%c%c"
			, sep
			, xShort + (isUp ? -2 : 0)
			, ss
			, cc1ac, cc1bc, cc1cc
			) ;
	    
	  }
	  break ;
      
	default:
	  vtxtPrintf(txt2, "%s%d:%c>%c"
		     , sep
		     , xLong 
		     , 'Z', 'Z'
		     ) ;
	  
	  vtxtPrintf(txt1, "%s%d:%c>%c"
		     , sep
		     , xShort 
		     , 'Z', 'Z'
		     ) ;
	}
    }
  dictAdd (bb->errDict, vtxtPtr (txt1), &up->errShort) ;
  dictAdd (bb->errDict, vtxtPtr (txt2), &up->errLong) ;
  
  return ;
} /* alignFormatErrors */

/**************************************************************/
#ifdef JUNK

static void alignClipErrorLeft (ALIGN *vp, int errCost)
{
  A_ERR *ep ;
  int i, iMax = arrayMax (vp->errors) ;
  int bestMax = -1 ;
  int x1 = vp->x1 ;
  
  if (iMax)
    {
      for (i = 0, ep = arrp (vp->errors, i, A_ERR) ; i < iMax ; i++, ep++)
	{
	  int dx = ep->iShort - x1 ;
	  if (dx < (i+1) * errCost)
	    { bestMax = i ; x1 = ep->iShort ; }
	}
      if (bestMax > -1)
	{
	  ep = arrp (vp->errors, bestMax, A_ERR) ;
	  if (0 && ep->iShort > vp->x2 - 5)
	    return ;
	  vp->x1 = ep->iShort + 1 ; /* last exact base */
	  if (vp->a1 < vp->a2)
	    vp->a1 = ep->iLong + 1 ;
	  else
	    {
	      vp->a1 = ep->iLong - 1 ;
	    }
	  switch (ep->type)
	    {   
	    case INSERTION:
	      vp->a1 += ((vp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_DOUBLE:
	      vp->x1 += 1 ;
	      vp->a1 += ((vp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_TRIPLE:
	      vp->x1 += 2 ;
	      vp->a1 += ((vp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case TROU:
	      vp->x1 -- ;
	      break ;
	    case TROU_DOUBLE:
	      vp->x1 -- ;
	      vp->a1 += ((vp->chrom & 0x1) ? -1: 1) ;
	      break ;
	    case TROU_TRIPLE:
	      vp->x1 -- ;
	      vp->a1 += ((vp->chrom & 0x1) ? -2 : 2) ;
	      break ;
	    default:
	      break ;
	    }
	  int k = bestMax + 1 ;
	  for (i = 0, ep = arrp (vp->errors, i, A_ERR) ; i < iMax - k ; i++, ep++)
	    *ep = *(ep + k) ;
	  arrayMax (vp->errors) -= k ;	    
	}
    }
  
  
  return ;
} /* alignClipErrorLeft */
#endif

/**************************************************************/

static void alignClipErrorRight (ALIGN *vp, int errCost)
{
  A_ERR *ep ;
  int i, iMax = arrayMax (vp->errors) ;
  int bestMax = iMax ;
  int x2 = vp->x2 ;
  
  if (iMax)
    {
      for (i = iMax - 1, ep = arrp (vp->errors, i, A_ERR) ; i >= 0 ; i--, ep--)
	{
	  int dx = x2 - ep->iShort ;
	  if (dx < (bestMax - i) * errCost)
	    { bestMax = i ; x2 = ep->iShort; }
	}
      if (bestMax < iMax)
	{
	  ep = arrp (vp->errors, bestMax, A_ERR) ;
	  vp->x2 = ep->iShort ; /* last exact base */
	  if (vp->a1 < vp->a2)
	    vp->a2 = ep->iLong ;
	  else
	    vp->a2 = ep->iLong ;
	  arrayMax (vp->errors) = bestMax ;
	}
    }

  return ;
} /* alignClipErrorRight */

/**************************************************************/
/* locate again the chains to eliminate the chain == -1 killed exons */
static int alignLocateChains (Array bestUp, Array aa, int myRead)
{
  int i1, i2, jj, iMax = arrayMax (aa) ;
  ALIGN *up, *vp ;
  
  bestUp = arrayReCreate (bestUp, keySetMax (bestUp), KEY) ;
  if (iMax)
    {
      for (i1 = i2 = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain > 0 && up->a1 != up->a2)
	    {
	      if (vp < up) *vp = *up ;
	      if (vp->read == myRead && vp->chain > jj)
		{
		  jj = vp->chain ;
		array (bestUp, jj, int) = i2 ;
		}
	      vp++ ; i2++ ;
	    }
	}
      iMax = arrayMax (aa) = i2 ;
    }
      
  up = arrayp (aa, iMax, ALIGN) ;
  memset (up, 0, sizeof (ALIGN)) ; /* force a null record */
  arrayMax (aa) = iMax ;

  return iMax ;
} /* alignLocateChains */

/**************************************************************/

static void alignAdjustIntrons (const PP *pp, BB *bb, Array bestUp, Array aa)
{
  ALIGN *up, *vp, *wp ;
  int chromA = 0 ;
  Array dnaG = 0 ;
    
  for (int ic = 1 ; ic < arrayMax (bestUp) ; ic++)
    {
      /* adjust introns */
      up = vp = arrp (aa, array (bestUp, ic, int), ALIGN) ; 
      int chain = up->chain ;
      if (up->chrom != chromA)
	{
	  chromA = up->chrom ;
	  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
	}

      
      wp = (vp[1].chain == chain ? vp + 1 : 0) ;
      while (vp && wp)
	{
	  if (vp > up && vp->x1 > wp->x1)
	    {
	      memset (vp, 0, sizeof (ALIGN)) ; vp->chain = -1 ;
	      saIntronsOptimize (bb, up, wp, dnaG) ;
	    }
	  else
	    saIntronsOptimize (bb, vp, wp, dnaG) ;
	  if (wp->chain != -1)   /* happens if the exons were merged */
	    vp = wp ;
	  wp = (vp[0].chain && wp[1].chain == vp[0].chain ? wp + 1 : 0) ;
	}
    }

} /* alignAdjustIntrons */

/**************************************************************/

static void alignAdjustExons (const PP *pp, BB *bb, Array bestUp, Array aa, Array dna)
{
  ALIGN *up, *vp ;
  int chromA = 0, chain = 0 ;
  AC_HANDLE h1 = ac_new_handle () ;
  Array dnaR = 0 ;
  Array dnaG = 0 ;


  for (int ic = 1 ; ic < arrayMax (bestUp) ; ic++)
    {
      /* adjust introns */
      up = vp = arrp (aa, array (bestUp, ic, int), ALIGN) ; 
      chain = up->chain ;
      if (up->chrom != chromA)
	{
	  chromA = up->chrom ;
	  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
	}

      
      /* count all errors */
      int nErr = 0 ;
      int da = 0, nAli = 0 ;
      int lnShort = arrayMax (dna) ;
      for (vp = up ; vp->chain == chain || vp->chain == -1 ; vp++)
	{
	  da = vp->a2 - vp->a1 ;
	  nErr += vp->nErr ;
	  nAli += da > 0 ? da + 1 : -da + 1 ;
	}

      if (nErr || nAli < lnShort)
	{
	  KEYSET ks = keySetHandleCreate (h1) ;
	  Array dnaI = 0, dnaShort = dna, errors = arrayHandleCreate (20, A_ERR, bb->h) ;
	  int ie, ia = 0, da ;
	  int ln = arrayMax (dnaG) ;
	  char *cp, *cq ;
	  A_ERR *ep ;
	  ALIGN zp ;
	  BOOL isDown ;
	  memset (&zp, 0, sizeof (A_ERR)) ;
		  
	  /* reconstruct the image of the transcript */

	  isDown = (up->chrom & 0x1) ? FALSE : TRUE ;
      	  nAli = nAli + 200 ;

	  dnaI = arrayHandleCreate (nAli + 1, char, h1) ;
	  array (dnaI, nAli, char) = 0 ; /* add a terminal zero */
	  arrayMax (dnaI) = nAli ;
	  if (! isDown)
	    {
	      int nvp = 0 ;
	      if (! dnaR)
		{
		  dnaR = dnaHandleCopy (dna, h1) ;
		  reverseComplement (dnaR) ;
		}
	      dnaShort = dnaR ;
	      for (vp = up ; vp->chain == chain ; vp++)
		{
		  int dummy = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = dummy ;
		  dummy = vp->x1 ; vp->x1 = lnShort - vp->x2 + 1 ; vp->x2 = lnShort - dummy + 1 ;
		  nvp++ ;
		}
	      for (int i = 0 ; i < nvp/2 ; i++)
		{
		  ALIGN wp = up[i] ;
		  up[i] = up[nvp - i - 1] ;
		  up[nvp -i - 1] = wp ;
		}
	    }
	  if (1)
	    {
	      int jj = up->a1 > 100 ? 100 : up->a1 - 1 ;
	      keySet (ks, ia++) = up->a1 - jj ;
	      cp = arrp (dnaI, 0, char) ;
	      if (jj)
		{
		  cq = arrp (dnaG, up->a1 - 1 - jj, char) ;
		  memcpy (cp, cq, jj) ;
		  cp += jj ;
		}
	      zp.x1 = up->x1 ;
	      zp.a1 = jj + 1 ;
	      for (vp = up ; vp->chain == chain || vp->chain == -1 ; vp++)
		{
		  da = vp->a2 - vp->a1 + 1 ;
		  if (vp->chain == -1 || da < 1) continue ;
		  keySet (ks, ia++) = vp->a1 - jj - 1 ;
		  cq = arrp (dnaG, vp->a1 - 1 , char) ;
		  memcpy (cp, cq, da) ;
		  jj += da ;
		  cp += da ;
		  if (vp[1].chain == chain)
		    {
		      int du = vp[1].x1 - vp->x2 - 1 ;
		      BOOL isDonor = TRUE ;
		      if (du >= 0)
			{ /* in the overlap OR the bases form the donor and acceptor exon */
			  cq = arrp (dnaG, vp->a2 - 1, char) ;
			  char *cr = arrp (dnaG, vp[1].a1 - 1, char) ;
			  if (cr[-2] == A_ && cr[-1] == G_ && cq[du + 1] == G_ && cq[du+2] == T_)
			    { /* copy the extension of the donor exon */
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -1] == G_ && cq[du+2-1] == T_)
			    { /* copy the extension of the donor exon with deletion */
			      du-- ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 2 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -2] == G_ && cq[du+2-2] == T_)
			    { /* copy the extension of the donor exon with double deletion */
			      du-=2 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >=3 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 1 -3] == G_ && cq[du+2-3] == T_)
			    { /* copy the extension of the donor exon with triple deletion */
			      du-=3 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 2] == G_ && cq[du+3] == T_)
			    { /* copy the extension of the donor exon with insertion */
			      du++ ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >= 1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 3] == G_ && cq[du+4] == T_)
			    { /* copy the extension of the donor exon with double insertion */
			      du+=2 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (du >=1 && cr[-2] == A_ && cr[-1] == G_ && cq[du + 4] == G_ && cq[du+5] == T_)
			    { /* copy the extension of the donor exon with triple insertion */
			      du+=3 ;
			      memcpy (cp, cq+1, du) ;
			    }
			  else if (cq[1] == G_ && cq[2] == T_ && cr[- du - 2] == A_ && cr[-du - 1] == G_)
			    { /* copy the extension of the acceptor exon */
			      isDonor = FALSE ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 2] == A_ && cr[-du - 1] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 3] == A_ && cr[-du - 2] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      du++ ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 4] == A_ && cr[-du - 3] == G_)
			    { /* copy the extension of the acceptor exon with insertion */
			      isDonor = FALSE ;
			      du+=2 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 1] == A_ && cr[-du - 0] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-- ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du - 0] == A_ && cr[-du + 1] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-=2 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else if (du >= 0 && cq[1] == G_ && cq[2] == T_ && cr[- du + 1] == A_ && cr[-du + 2] == G_)
			    { /* copy the extension of the acceptor exon with deletion */
			      isDonor = FALSE ;
			      du-=3 ;
			      /* do not copy it will be automatic since we shift vp[1] below */
			    }
			  else
			    memset (cp, N_, du) ;
			  if (isDonor)
			    {
			      jj += du ;
			      cp += du ;
			      vp->a2 += du ; vp->x2 += du ;
			    }
			  else
			    {
			      vp[1].x1 -= du ; vp[1].a1 -= du ;
			    }
			}
		    }
		}
	      vp-- ; ia-- ;
	      zp.x2 = vp->x2 ;
	      zp.a2 = jj ;

	      if (vp->x2 < lnShort)
		{
		  cp = arrp (dnaI, jj, char) ;
		  da = ln > vp->a2 + 100 ? 100 : ln - vp->a2 ;
		  if (da)
		    {
		      cp = arrayp (dnaI, jj + da, char) ; /* may reallocate */
		      cp = arrp (dnaI, jj, char) ;
		      cq = arrp (dnaG, vp->a2, char) ;
		      memcpy (cp, cq, da ) ;
		      jj += da ;
		    }
		  arrayMax (dnaI) = jj ;
		}
	      
	      /* align the read on the genomic image of the transcript */
	      zp.errors = errors ;
	      arrayMax (zp.errors) = 0 ;
	      if (1)
	      {
		int x2 = zp.x2,  a2 = zp.a2 ;
		aceDnaDoubleTrackErrors (dnaShort, &(zp.x1), &(zp.x2), TRUE   /* isDown = TRUE */    /* MAXJUMP */
					 , dnaI, 0, &(zp.a1), &(zp.a2)
					 , 0, zp.errors, MAXJUMP, -3, TRUE, 0) ; /* bio coordinates, extend right */
		if (zp.x2 < x2 - 50)
		  {
		    zp.x2 = x2 ; zp.a2 = a2 ;
		    arrayMax (zp.errors) = 0 ;
		    aceDnaDoubleTrackErrors (dnaShort, &(zp.x1), &(zp.x2), TRUE   /* isDown = TRUE */    /* MAXJUMP */
					     , dnaI, 0, &(zp.a1), &(zp.a2)
					     , 0, zp.errors, MAXJUMP2, -1, FALSE, 0) ; /* bio coordinates, jump 8 but do not extend */
		  }
	      }
	      alignClipErrorRight (&zp, pp->errCost) ;
	      /* remap */
	      int ja ;
	      for (vp = up, ja = 1 ; vp->chain == chain || vp->chain == -1 ; vp++)
		{
		  int dz = keySet (ks, ja)  ;
		  int j = 0 ;
		  da = vp->a2 - vp->a1 + 1 ;
		  if (vp->chain == -1 || da < 1) continue ;
		  if (vp->a1 > zp.a2 + dz)
		    { vp->chain = -1 ; continue ; }
		  if (vp->a2 > zp.a2 + dz) { vp->a2 = zp.a2 + dz ; vp->x2 = zp.x2 ; }
		  
		  if (vp->errors)
		    arrayMax (vp->errors) = 0 ; 
		  for (ie = 0 ; ie < arrayMax (zp.errors) ; ie++)
		    {
		      ep = arrp (zp.errors, ie, A_ERR) ;
		      int es = ep->iShort + 1 ;
		      if (es >= vp->x1 && es <= vp->x2)
			{
			  if (! vp->errors)
			    vp->errors = arrayHandleCreate (8, A_ERR, bb->h) ;
			  A_ERR *eq = arrayp (vp->errors, j++, A_ERR) ;
			  *eq = *ep ;
			  eq->iLong += dz ;
			}
		    }
		  ja++ ;
		  vp->nErr = vp->errors ? arrayMax (vp->errors) : 0 ;
		}
	    }
	  if (! isDown)
	    {
	      int nvp = 0 ;
	      for (vp = up ; vp->chain == chain ; vp++)
		{
		  int dummy = vp->a1 ; vp->a1 = vp->a2 ; vp->a2 = dummy ;
		  dummy = vp->x1 ; vp->x1 = lnShort - vp->x2 + 1 ; vp->x2 = lnShort - dummy + 1 ;
		  nvp++ ;

		  /* flip the rror positions */
		  if (vp->nErr)
		    {
		      int i ;
		      for (i = 0, ep = arrp (vp->errors, 0, A_ERR) ; i < vp->nErr ; i++, ep++)
			{
			  ep->iShort = lnShort - ep->iShort - 1 ;
			  ep->baseShort = complementBase [(int)ep->baseShort] ;
			  ep->sens = -1 ;
			}
		    }
		}
	      for (int i = 0 ; i < nvp/2 ; i++)
		{
		  ALIGN wp = up[i] ;
		  up[i] = up[nvp - i - 1] ;
		  up[nvp-i - 1] = wp ;
		}
	    }
	}
    }

  ac_free(h1) ;
  
} /* alignAdjustExons */

/**************************************************************/
/* Dynamic programming of path score */

static void alignSelectBestDynamicPath (const PP *pp, BB *bb, Array aaa, Array aa, Array dna, int chromA, Array dnaG, Array dnaGR, Array bestUp) 
{
  AC_HANDLE h = 0 ;
  int ii, jj, i1, i2, iMax ;
  ALIGN *up, *vp, *wp ;
  int i02 = 0 ;
  int chainAli = 0 ;
  int chainScore = 0 ;
  int bestChainScore = 0 ;
  int bestI1 = 0 ;
  int chain = 0 ;
  int maxIntron = pp->maxIntron ;
  int minAli = pp->minAli ;
  Array aaNew = 0 ;
  int errCost = pp->errCost ;
  int bigErrCost = 8 ; /* errCost ; */
  int myRead = arrp (aa, 0, ALIGN)->read ;
  iMax = arrayMax (aa) ;
  arraySort (aa, alignOrder) ;

  /* create scores */
  if (iMax)
    {
      for (ii = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; ii < iMax ; ii++, up++)
	{
	  int ali = up->x2 - up->x1 + 1 ;
	  int score = ali - up->nErr * errCost ;
	  
	  if (score > 0)
	    {
	      up->chain = 0 ;
	      up->chainScore = 0 ;
	      up->ali = ali ;
	      up->id = jj + 1 ;
	      up->score = up->ali - up->nErr * bigErrCost ;
	      if (vp < up) *vp = *up ;
	      if (1)
		{
		  wp = vp - 1 ;
		  if (jj == 0 || vp->a1 != wp->a1 || vp->a2 != wp->a2 || vp->x1 != wp->x1 || vp->x2 != wp->x2)
		    {
		      vp++ ; jj++ ;
		    }
		  else if (jj)
		    {
		      if (vp->donor && ! wp->donor)
			wp->donor = vp->donor ;
		      if (vp->acceptor && ! wp->acceptor)
			wp->acceptor = vp->acceptor ;
		      }
		}
	    }
	}
    }

  iMax = arrayMax (aa) = jj ;
  if (! iMax) return ;
  
  /* for each exon moving left to right, construct the chains */
  i2 = 0 ; vp = arrp (aa, 0, ALIGN) ; /* preposition */
  
  for (i1 = 0, up = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
    {
      int chrom = up->chrom ;
      int x1 = up->x1 ;
      int x2 = up->x2 ;
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int bestPrevious = 0 ;
      int bestPreviousScore = 0 ;
      int bestPreviousAli = 0 ;
      BOOL isDown = ! (chrom & 0x1) ;
      BOOL foundI2 = FALSE ;
      
      i2 = i02 ; vp = arrp (aa, i2, ALIGN) ; /* preposition */
      for (foundI2 = FALSE ; pp->splice && i2 < iMax ; i2++, vp++)
	{
	  if (i1 == i2)
	    continue ;
	  if (vp->chrom > chrom || (vp->chrom == chrom && vp->x1 > x2 - MAXJUMP2))
	    break ;
	  if (vp->chrom < chrom || vp->x2 < x1 - MAXJUMP2)
	    {
	      if (!foundI2) i02 = i2 + 1 ;
	      continue ;
	    }  
	  foundI2 = TRUE ;
	  if (vp->chrom == chrom
	      && vp->x2 >= x1 - 8 && vp->x2 < x2 && vp->x1 < x2
	      &&
	      (
	       ( isDown && vp->a1 < a2 && vp->a2 + maxIntron > a1 && vp->a2 - vp->x2 < a1 - x1 + MAXJUMP2) ||
	       ( ! isDown && vp->a1 > a2 && vp->a2 - maxIntron < a1 && vp->a2 + vp->x2 > a1 + x1  - MAXJUMP2)
	       )
	      )
	    {
	      int dx = vp->x2 - x1 + 1 > 0 ? vp->x2 - x1 + 1 : 0 ;
	      if (vp->chainScore - dx > bestPreviousScore)
		{
		  bestPreviousScore = vp->chainScore - dx ;
		  bestPreviousAli = vp->chainAli - dx ;
		  bestPrevious = i2 + 1 ;
		}
	      }
	}
      
      up->chainScore = up->score + bestPreviousScore ;
      up->chainAli = up->ali + bestPreviousAli ;
      up->previous = bestPrevious ;
      if (up->chainScore > bestChainScore)
	{
	  bestChainScore = up->chainScore ;
	  bestI1 = i1 ;
	}
    }
  
  /* attribute a chain number recursivelly to the best paths */
  chain = 0 ;
  i1 = bestI1 ;
  h = ac_new_handle () ;
  aaNew = arrayHandleCreate (iMax, ALIGN, h) ;

  while (bestChainScore > 0)
    {
      /* find the top of the chain */
      up = arrp (aa, bestI1, ALIGN) ;
      chainScore = up->score ;
      chainAli = up->chainAli ;
      while (up->previous)
	{
	  vp = arrp (aa, up->previous - 1, ALIGN) ;
	  vp->next = up->id ;
	  int dx = vp->x2 - up->x1 + 1 ;
	  chainScore += vp->score + (dx> 0 ? -dx : 0) ;
	  chainAli += vp->ali + (dx> 0 ? -dx : 0) ;

	  i1 = up->previous ;
	  up = vp ;
	}
      /* register this chain in aaNew, kill it in aa */
      chain++ ;
      int j0 = 0, jj = arrayMax (aaNew) ;
      array (bestUp, chain, int) = jj ;  /* best chain */
      wp = 0 ;
      while (up)
	{
	  up->chainScore =  chainScore ;
	  up->chainAli =  chainAli ;
	  vp = arrayp (aaNew, jj++, ALIGN) ;
	  *vp = *up ;
	  
	  vp->id = jj ;
	  if (j0)
	    {
	      vp[-1].next = jj ;
	      vp->previous = j0 ;
	    }
	  else
	    {
	      wp = vp ;
	      wp->chainX1 = up->x1 ;
	    }
	  j0 = jj ;
	  up->chain = vp->chain = chain ;
	  wp->chainX2 = vp->x2 ;
	  up = (up->next ? arrayp (aa, up->next - 1, ALIGN) :0) ;
	}
      /* filter on chainX1/chainX2 */
      BOOL ok = TRUE ;
      if (wp->chainX2 - wp->chainX1 + 1 < minAli)
	ok = FALSE ;
      for (int ic = 1 ; ok && ic < chain ; ic++)
	{
	  int iw = arr (bestUp, ic, int) ;
	  up = arrp (aaNew, iw, ALIGN) ;
	  int z1 = (up->chainX1 > wp->chainX1 ? up->chainX1 : wp->chainX1) ;
	  int z2 = (up->chainX2 < wp->chainX2 ? up->chainX2 : wp->chainX2) ;
	  int dz = z2 - z1 ;
	  int du = up->chainX2 - up->chainX1 ;
	  int dw = wp->chainX2 - wp->chainX1 ;
	  int tc = *dictName(pp->bbG.dict,wp->chrom >> 1) ;
	  if (
	      (2 * dz > du || 2 * dz > dw) /* significant overlap */
	      && (up->chainScore > 2 * wp->chainScore   || up->chainScore > wp->chainAli + pp->bonus[tc]) /* 10: allow for class bonus */
	      )
	    ok = FALSE ;
	}
      if (! ok)
	{  /* destroy this chain in aaNew */
	  int iw = arr (bestUp, chain, int) ;
	  arrayMax (aaNew) = iw ;
	  arrayMax (bestUp) = chain ;
	  chain-- ;

	  /* flag the bad chain */
	  up = arrp (aa, bestI1, ALIGN) ;
	  up->chain = -1 ;
	  while (up->previous)
	    {
	      up = arrp (aa, up->previous - 1, ALIGN) ;
	      up->chain = -1 ;
	    }	  
	}
      
      /* edit the score of the other exons
       * and look for new best score
       * this applies even for destroyed chains 
       */
      bestI1 = 0 ;
      bestChainScore = 0 ;
      for (i1 = 0, up = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain) /* allready used */
	    continue ;
	  if (up->previous && (vp = arrp (aa, up->previous - 1, ALIGN)) && !vp->chain)
	    up->chainScore = up->score + vp->chainScore ;
	  else
	    { up->chainScore = up->score ; up->previous = 0 ; } /* disconnect */
	  if (up->chainScore > bestChainScore)
	    { bestChainScore = up->chainScore ; bestI1 = i1 ; }
	}
    }

  /* transfer the sorted chains back in aa */
  iMax = arrayMax (aaNew) ;
  if (iMax)
    for (i1 = i2 = 0, up = arrp (aaNew, 0, ALIGN) ; i1 < iMax ; i1++, up++)
      {
	if (up->chain > 0)
	  {
	    vp = arrp (aa, i2++, ALIGN) ;
	    *vp = *up ;
	  }
      }
  iMax =   arrayMax (aa) = i2 ;
  arrayMax (bestUp) = 0 ;
  if (iMax)
    for (i1 = i2 = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
      {
	if (1)
	  {
	    if (vp < up) *vp = *up ;
	    if (vp->chain > jj)
	      {
		jj = vp->chain ;
		array (bestUp, jj, int) = i2 ;
	      }
	    vp++ ; i2++ ;
	  }
      }
  iMax = arrayMax (aa) = i2 ;      /* readjust */
  
  /* adjust introns and scores
   * both operations need to know the genome 
   */
  
  /* locate the chains */
  iMax = alignLocateChains (bestUp, aa, myRead) ;

  if (iMax)
    {
      /* adjust introns */
      alignAdjustIntrons (pp, bb, bestUp, aa) ;
      iMax = alignLocateChains (bestUp, aa, myRead) ;
      /* adjust exons */
      alignAdjustExons (pp, bb, bestUp, aa, dna) ;
      iMax = alignLocateChains (bestUp, aa, myRead) ;
      
      /* Compute the clean chain score */
      bestChainScore = 0 ;
      for (int ic = 1 ; ic < arrayMax (bestUp) ; ic++)
	{
	  ii =array (bestUp, ic, int) ;
	  up = arrp (aa, ii, ALIGN) ; 
	  
	  int tc = *dictName(pp->bbG.dict,up->chrom >> 1) ;
	  int chain = up->chain ;
	  int chainX1 = up->chainX1 ;
	  int chainX2 = up->chainX2 ;
	  
	  int chainAli = up->chainAli = 0 ;
	  int chainErr = up->chainErr = 0 ;
	  int chainScore =	pp->bonus[tc] ;
	  int chainA1 = up->a1 ;
	  int chainA2 = up->a2 ;
	  
	  for (vp = up, jj = ii ; jj < iMax && vp->chain == chain ; jj++, vp++)
	    {
	      vp->targetClass = tc ;
	      vp->ali = vp->x2 - vp->x1 + 1 ;
	      vp->score = vp->ali - errCost * vp->nErr ;
	      if (vp->score < -10)
		continue ;
	      chainX2 = vp->x2 ;
	      chainAli += vp->ali ;
	      chainErr += vp->nErr ;
	      chainScore += vp->score ;
	      
	      if (chainA1 > vp->a1) chainA1 = vp->a1 ;
	      if (chainA1 > vp->a2) chainA1 = vp->a2 ;
	      if (chainA2 < vp->a1) chainA2 = vp->a1 ;
	      if (chainA2 < vp->a2) chainA2 = vp->a2 ;
	      
	    }
	  if (chainAli > arrayMax (dna))
	    chainAli = arrayMax (dna) ;
	  /* filter */
	  if (chainScore < pp->minScore ||
	      chainAli < pp->minAli ||
	      100 * chainAli < pp->minAliPerCent * arrayMax (dna) ||
	      100 * chainErr > pp->errRateMax * chainAli
	      )
	    chainScore = chainAli = chainErr = 0 ;
	  
	  if (bestChainScore < chainScore)
	    bestChainScore = chainScore ;
	  
	  /* set the chain values in all exons */
	  if (up->a1 > up->a2)
	    { int dummy = chainA1 ; chainA1 = chainA2 ; chainA2 = dummy ; }
	  for (vp = up, jj = ii ; jj < iMax && vp->chain == chain ; jj++, vp++)
	    {
	      vp->chainScore = chainScore ;
	      vp->chainAli = chainAli ;
	      vp->chainErr = chainErr ;
	      vp->chainX1 = chainX1 ;
	      vp->chainX2 = chainX2 ;
	      vp->chainA1 = chainA1 ;
	      vp->chainA2 = chainA2 ;
	    }
	  up = vp - 1 ; ii = jj - 1 ;
	}
    }
  
  /* clean up the destroyed chains and adjust the chain numbers */
  iMax = arrayMax (aa) ;

  if (iMax)
    {
      arraySort (aa, alignOrder) ;
      
      int newChain = 0 ;
      for (i1 = i2 = chain = 0,  up = vp = arrp (aa, 0, ALIGN) ; i1 < iMax ; i1++, up++)
	{
	  if (up->chain && up->chainScore > 0 && up->score > -10)
	    {
	      int k = 0 ;
	      if (1)
		{
		  for (int i3 = 1 ; up->score && i3 <= newChain ; i3++)
		    {
		      wp = arrp (aa, arr (bestUp, i3, int), ALIGN) ;
		      if (wp->chainScore > up->chainScore)
			{
			  int z1 = (up->chainX1 > wp->chainX1 ? up->chainX1 : wp->chainX1) ;
			  int z2 = (up->chainX2 < wp->chainX2 ? up->chainX2 : wp->chainX2) ;
			  int dz = z2 - z1 ;
			  int du = up->chainX2 - up->chainX1 ;
			  int dw = wp->chainX2 - wp->chainX1 ;
			  
			  if  (2 * dz > du || 2 * dz > dw) /* significant overlap */
			    {
			      int c = up->chain ;
			      for (wp = up ; i1 < iMax && wp->chain == c ; wp++)
				{ k++ ; wp->chainScore = wp->chain = wp->score = 0 ; }
			      i1 += k - 1 ; up += k - 1 ;
			    }
			}
		      
		    }
		  if (up->score)
		    {
		      int c = up->chain ;
		      newChain++ ;
		      array (bestUp, newChain, int) = i2 ;
		      for (wp = up ; i1 < iMax && wp->chain == c ; wp++)
			if (wp->score)
			  {
			    wp->chain = newChain ;
			    if (vp < wp) *vp = *wp ;
			    
			    i2++ ; vp++ ;
			    i1++ ; up++ ;
			  }
		      i1-- ; up-- ;
		    }		
		}
	    }
	}
      arrayMax (bestUp) = newChain ? newChain + 1 : 0 ;
      arrayMax (aa) = i2 ;
    }
  
  /* register the alignments */
  int kMax = arrayMax (aaa) ;
  iMax = arrayMax (aa) ;
  if (iMax)
    {
      up = arrp (aa, 0, ALIGN) ;
      /*   bitSet (bb->isAligned, up->read) ; */
      for (ii = 0 ; ii < iMax ; ii++, up++)
	{
	  if (up->score > -10)
	    {
	      vp = arrayp (aaa, kMax++, ALIGN) ;
	      *vp = *up ;
	      vp->nChains = arrayMax (bestUp) ;
	    }
	}
    }
  iMax = alignLocateChains (bestUp, aaa, myRead) ;  
  ac_free (h) ;
  
  return ;
} /* alignSelectBestDynamicPath */

/**************************************************************/
/* Establish chain scores, select best */
static void  alignDoRegisterOnePair (const PP *pp, BB *bb, BigArray aaa, Array aa, int read, Array bestUp)
				
{
  ALIGN *up, *vp ;

  int ii ;
  int iMax = alignLocateChains (bestUp, aa, read) ;  
  int nChains = 0 ;
  char allTc[256] ;
  Array dna = 0, dna1 = 0, dna2 = 0 ;
  int chromA = 0 ;
  Array dnaG = 0, dnaGR = 0 ;
  int read1 = read & (~0x1) ;
  int read2 = read | 0x1 ;

  dna1 = arr (bb->dnas, read1, Array) ;
  dna2 = arr (bb->dnas, read2, Array) ;
  
  /* create chains */

  /* register the wiggle boundaries */
  if (pp->wiggle)
    {
      int wiggleStep = WIGGLE_STEP ;  /* examples s=10, 5, 1 */
      int demiStep = wiggleStep/2 ;
      if (2*demiStep == wiggleStep) demiStep-- ; /* examples d=4, 2, 0 */
      
      iMax = arrayMax (aa) ;
      if (iMax)
	for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
	  if (read == up->read)
	    {
	      int a1 = up->a1 ;
	      int a2 = up->a2 ;
	      if (up->read & 0x1) { int a0 = a1 ; a1 = a2 ; a2 = a0 ;}
	      up->w1 = (a1 + demiStep)/wiggleStep ;
	      up->w2 = (a2 + demiStep)/wiggleStep ;
	    }
    }

  /* overhangs */
  iMax = arrayMax (aa) ;
  if (iMax)
    for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
      if (read == up->read)
	{
	  dna = up->read & 0x1 ? dna2 : dna1 ;
	  if (up->chrom != chromA)
	    {
	      chromA = up->chrom ;
	      dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
	      dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
	    }
	  
	  if (up->x1 && up->x1 == up->chainX1)	
	    alignFormatLeftOverhang (pp, bb, up, dna, dnaG, dnaGR) ;
	  if (up->x2 == up->chainX2)	
	    alignFormatRightOverhang (pp, bb, up, dna, dnaG, dnaGR) ;
	}
  
  /* format the errors */
  iMax = arrayMax (aa) ;
  if (iMax)
    for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
      if (read == up->read)
	{
	  if (arrayExists (up->errors))
	    {
	      unsigned int flip = 0 ;
	      
	      dna = up->read & 0x1 ? dna2 : dna1 ;
	      if (up->chrom != chromA)
		{
		  chromA = up->chrom ;
		  dnaG = arr (pp->bbG.dnas, chromA >> 1, Array) ;
		  dnaGR = arr (pp->bbG.dnasR, chromA >> 1, Array) ;
		}
	      
	      if (bb->runStat.nPairs && (up->read & 0x1))
		flip = 0x0f ; /* will flip last 4 bits */
	      if (arrayMax (up->errors))
		mergeErrors (bb->errors, up->errors, flip) ;
	      if (arrayMax (up->errors))
		alignFormatErrors (pp, bb, up, dna, dnaG, dnaGR) ;	  
	      if (! pp->sam)
		arrayDestroy (up->errors) ;
	    }
	}
  
  /* global statistics */
  /* stranding : once per target class */
  memset (allTc, 0, sizeof (allTc)) ;
  
  for (int ic = 0 ; ic < arrayMax (bestUp) ; ic++)
    {
      up = arrp (aa, array (bestUp, ic, int), ALIGN) ; 
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int tc = up->targetClass ;
      BOOL s = a1 < a2 ;
      
      if (! tc || allTc[tc])
	continue ;
      allTc[tc] = 1 ;
      if (up->read & 0x1)
	s = !s ;
      if (s)
	bb->runStat.GF[tc]++ ;
      else
	bb->runStat.GR[tc]++ ;
      
      bb->runStat.nAlignedPerTargetClass[tc]++ ;
    }
  
  /* increase the block stats */
  nChains = 0 ;
  if (arrayMax (bestUp))
    {
      memset (allTc, 0, sizeof (allTc)) ;
      
      up = arrp (aa, array (bestUp, 1, int), ALIGN) ; 
      int tc0 = up->targetClass ;
      bb->nAli++ ;
      bb->runStat.nAlignedPerTargetClass[0]++ ;
      bb->runStat.nMultiAligned[0]++ ;

      nChains = 1 ;

      bb->runStat.nErr += up->chainErr ;
      bb->aliDx += up->chainAli ;
      bb->aliDa += up->chainAli ;
      bb->runStat.nAlignments++ ;
      
      dna = up->read & 0x1 ? dna2 : dna1 ;
      if (up->chainErr == 0 && up->chainAli == arrayMax (dna))
	bb->runStat.nPerfectReads++ ;
      
      if (up->read & 0x1)
	bb->runStat.nBaseAligned2 += up->chainAli ;
      else
	bb->runStat.nBaseAligned1 += up->chainAli ;

      for (int ic = 2 ; ic < arrayMax (bestUp) ; ic++)
	{
	  vp = arrp (aa, array (bestUp, ic, int), ALIGN) ;
	  int tc = vp->targetClass ;
	  if (read == vp->read)
	    if (tc == tc0)
	      { /* count multiali only in main class */
		int z1 = (up->chainX1 > vp->chainX1 ? up->chainX1 : vp->chainX1) ;
		int z2 = (up->chainX2 < vp->chainX2 ? up->chainX2 : vp->chainX2) ;
		int dz = z2 - z1 ;
		int du = up->chainX2 - up->chainX1 ;
		int dv = vp->chainX2 - vp->chainX1 ;

		if (2 * dz > du || 2 * dz > dv) /* significant overlap */
		  {
		    bb->runStat.nAlignments++ ;
		    nChains++ ;
		  }
		else
		  {
		    bb->runStat.nErr += vp->chainErr ;
		    bb->aliDx += vp->chainAli ;
		    bb->aliDa += vp->chainAli ;
		    if (vp->read & 0x1)
		      bb->runStat.nBaseAligned2 += vp->chainAli ;
		    else
		      bb->runStat.nBaseAligned1 += vp->chainAli ;
		  }
	      }
	}
      if (nChains)
	bb->runStat.nMultiAligned[nChains > 10 ? 10 : nChains]++ ;
    }
  /* register the alignments */
  long int kMax = bigArrayMax (aaa) ;

  iMax = arrayMax (aa) ;
  if (iMax)
    {
      up = arrp (aa, 0, ALIGN) ;
      /*   bitSet (bb->isAligned, up->read) ; */
      for (ii = 0 ; ii < iMax ; ii++, up++)
	if (read == up->read)
	  {
	    vp = bigArrayp (aaa, kMax++, ALIGN) ;
	    *vp = *up ;
	    vp->nChains = nChains ;
	  }
    }
  return ;
} /* alignSelectBestChain */

/**************************************************************/
static void alignDoOneRead (const PP *pp, BB *bb
			    , Array aaa, BigArray hits
			    , Array aa, Array err, Array bestUp)
{   
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  HIT * restrict hit ;
  ALIGN *ap ;
  long int ii, iMax = bigArrayMax (hits), kMax = 0 ;
  int a1, a2, x1, x2 ;
  int b1, b2, y1, y2, ha1, readOld = 0, chromOld = 0, readA = 0, chromA = 0, read1 = 0, iiGood = 0 ;
  BOOL isDownOld = TRUE ;
  Array dna = 0, dnaG = 0, dnaGR = 0 ;
  int errMax = pp->errMax ; /* 999999 ; */
  int chromLength = 0 ;
  /*
    HIT * restrict h1 ;
    int r1 = 0, nh1 = 0, chrom1 = 0 ;
  */
  int errCost = pp->errCost ;
  /*   unsigned int uu = 0 ; */
  int donor = 0, acceptor = 0 ;
  const int intronBonus = 1 ;
  int nTargetRepeats  = 1 ;
  int nTargetRepeatsOld = 0 ;
  const int nTRmask = (0x1 << NTARGETREPEATBITS) - 1 ;

  for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; ii++, hit++)
    {
      int read = hit->read ;
      int chrom = hit->chrom ;
      BOOL isDown = TRUE ;
      BOOL isIntron = ((hit->x1  >> NTARGETREPEATBITS )  & 0x7) ? TRUE : FALSE ;

      if (! read || ! chrom)
	continue ;
      if (ii < iMax  && ! memcmp (hit, hit + 1, sizeof (HIT)))
	continue ;
      if (read != read1)
	{
	  read1 = read ;
	  if (arrayMax (aa))
	    { /* create chain scores */
	      alignSelectBestDynamicPath (pp, bb, aaa, aa, dna, chromA, dnaG, dnaGR, bestUp) ;
	    }
	  arrayMax (aa) = kMax = 0 ;
	}

      if (read != readA)
	{ readA = read ; dna = arr (bb->dnas, read, Array) ; }
      if (chrom != chromA)
	{
	  chromA = chrom ;
	  dnaG = arr (pp->bbG.dnas, chrom >> 1, Array) ;
	  dnaGR = arr (pp->bbG.dnasR, chrom >> 1, Array) ;
	  chromLength = arrayMax (dnaG) ;
	}

      x1 = hit->x1 ;
      nTargetRepeats = (x1 & nTRmask) ; 
      x1 = x1 >> NTARGETREPEATBITS ;
      BOOL isIntronDown = (x1 >> 2) & 0x1 ;
      isDown = (chrom & 0x1)  ? FALSE : TRUE ;
      donor = x1 & 0x1 ;
      acceptor = x1 & 0x2 ;
      x1 = x1 >> 3 ; x2 = x1 + 1 ;    /* bio coordinates */   
      if (isDown)   /* plus strand of the genome */
	{
	  a1 =
	    hit->a1
	    + x1
	    + (isIntron ? intronBonus : 0) 
	    - 1 ;        /* compensate avoid zero */
	  a2 = a1 + 1 ;
	  if (donor) donor = a1 + 2 ; /* first base of intron */
	  if (acceptor) acceptor = a1 - 1 ; /* last base of intron */
	}
      else   /* minus strand of the genome */
	{
	  a1 =
	    hit->a1
	    - x1
	    + (isIntron ? intronBonus : 0) 
	    - 1 ;        /* compensate avoid zero */
	  a2 = a1 - 1 ;
	  if (donor) donor = a1 - 2 ; /* first base of intron */
	  if (acceptor) acceptor = a1 + 1 ; /* last base of intron */
	}

      if (! isIntronDown)
	{ donor = - donor ; acceptor = - acceptor ; }
      if (1 && read == readOld && chrom == chromOld && isDown == isDownOld &&
	  x1 >= y1 && x2 <= y2 && hit->a1 < ha1 + 3 &&   /* MAXJUMP */
	  (
	   (isDown && a1 >= b1 && a2 <= b2) ||
	   (! isDown && a1 <= b1 && a2 >= b2)
	   )
	  &&
	  nTargetRepeats >= nTargetRepeatsOld >> 1
	  )
	{
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tDoublet of %d\t%s\t%d\n", ii, read, x1, x2, chrom, a1, a2, iiGood, dictName (pp->bbG.dict, chrom >> 1), hit->a1) ;
	  hit->read = 0 ;  /* remove doublets */
	  if (kMax)
	    {
	      if (donor)
		ap->donor = donor ;
	      if (acceptor)
		ap->acceptor = acceptor ;
	    }
	}
      else 
	{
	  chromOld = 0 ;
	  int a0 = a1, x0 = x1 ;
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tbefore align\t%s\t%u\n"
			      , ii, read, x1, x2, chrom, a1, a2, dictName (pp->bbG.dict, chrom >> 1), hit->a1) ;
	  if (alignExtendHit (dna, dnaG, dnaGR, err, isDown, chromLength, &a1, &a2, &x1, &x2, errCost, isIntron, errMax, 22))
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tAccepted\t%s, u=%u, nErr=%d\n"
				  , ii, read, x1, x2, chrom, a1, a2
				  , dictName (pp->bbG.dict, chrom >> 1)
				  , hit->a1
				  , arrayMax (err)
				  ) ;
	      ap = arrayp (aa, kMax++, ALIGN) ;
	      memset (ap, 0, sizeof (ALIGN)) ;
	      ap->read = read ;
	      ap->chrom = chrom ;
	      ap->a0 = a0 ;
	      ap->a1 = a1 ;
	      ap->a2 = a2 ;
	      ap->x0 = x0 ;
	      ap->x1 = x1 ;
	      ap->x2 = x2 ;
	      ap->nTargetRepeats = nTargetRepeats ;
	      ap->donor = donor ;
	      ap->acceptor = acceptor ;
	      ap->readLength = arrayMax (dna) ;
	      ap->nErr = arrayMax (err) ;
	      if (ap->nErr)
		ap->errors = arrayHandleCopy (err, bb->h) ;
	      readOld = read ;
	      chromOld = chrom ;
	      isDownOld = isDown ;
	      nTargetRepeatsOld = nTargetRepeats ;
	      b1 = a1 ; b2 = a2 ; y1 = x1 ; y2 = x2 ; ha1 = hit->a1 ; iiGood = ii ;
	    }
	  else
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tRejected\tu=%d\n", ii, read, x1, x2, chrom, a1, a2, hit->a1) ;
	      hit->read = 0 ; /* remove false positive */
	    }
	}
    }
  if (arrayMax (aa))
    { /* create chain scores */
      alignSelectBestDynamicPath (pp, bb, aaa, aa, dna, chromA, dnaG, dnaGR, bestUp) ;
    }
  ac_free (h) ;
  return ;
} /* alignDoOneRead */

/**************************************************************/
/**************************************************************/
typedef struct pairStruct { ALIGN *up, *vp ; int score, chrom, a1, a2 ; } PAIR ;
static int pairOrder (const void *va, const void *vb)
{
  const PAIR *up = va ;
  const PAIR *vp = vb ;
  int n ;
  n = up->score - vp->score ; if (n) return -n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;

  return 0 ;
} /* pairOrder */

/**************************************************************/
static void alignDoOnePair (const PP *pp, BB *bb
			    , BigArray aaaa, BigArray hits
			    , Array aa, Array err)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT *hit ;
  int ii, iMax = arrayMax (hits) ;
  int read1, read2 ;
  Array aaa = arrayHandleCreate (128, ALIGN, h) ;
  Array bestUp1 = arrayHandleCreate (8, int, h) ;
  Array bestUp2 = 0 ;

  read1 = bigArr (hits, 0, HIT).read ;
  read2 = bigArr (hits, iMax -1, HIT).read ;

  if (read1 != read2)
    {
      int iMax1 ;
      for (ii = 0, hit = bigArrp (hits, 0, HIT) ; ii < iMax ; hit++, ii++)
	if (hit->read == read2)
	  break ;
      arrayMax (hits)  = iMax1 = ii ;
      alignDoOneRead (pp, bb, aaa, hits, aa, err, bestUp1) ;
      arrayMax (hits) = iMax ;
      if (iMax > iMax1)
	for (ii = iMax1, hit = bigArrp (hits, ii, HIT) ; ii < iMax ; hit++, ii++)
	  *(hit - iMax1) = *hit ;
      arrayMax (hits) = iMax - iMax1 ;
      arrayMax (aa) = arrayMax (err) = 0 ;
      bestUp2 = arrayHandleCreate (8, int, h) ;
      iMax1 = arrayMax (aaa) ;
      alignDoOneRead (pp, bb, aaa, hits, aa, err, bestUp2) ;
    }
  else
    alignDoOneRead (pp, bb, aaa, hits, aa, err, bestUp1) ;

  if (bestUp2) /* we have a pair, good  example polyA_B_1 read 144/145*/
    {
      ALIGN *up, *vp ;
      int iMax1 = arrayMax (bestUp1) ;
      int iMax2 = arrayMax (bestUp2) ;

      if (iMax1 && iMax2)
	{
	  Array pairs = arrayHandleCreate (iMax1 * iMax2 , PAIR, h) ;
	  PAIR *px ;
	  int i1, i2 ;
	  int jj = 0 ;
	  for (i1 = 1 ; i1 < iMax1 ; i1++)
	    for (i2 = 1 ; i2 < iMax2 ; i2++)
	      {
		int j1 = arr (bestUp1, 1, int) ;
		int j2 = arr (bestUp2, 1, int) ;

		up = arrp (aaa, j1, ALIGN) ;
		vp = arrp (aaa, j2, ALIGN) ;

		if ((up->chrom ^ vp->chrom) == 0x1)
		  {
		    if ((up->chrom & 0x1) == 1) { ALIGN *zp = up ; up = vp ; vp = zp ;}
		    int da = vp->chainA1 - up->chainA1 ;
		    int db = vp->chainA2 - up->chainA2 ;

		    if (da > 0 && db < 1000000) /* true pair */
		      {
			px = arrayp (pairs, jj++, PAIR) ;
			px->up = up ;
			px->vp = vp ;
			px->chrom = up->chrom ;
			px->a1 = up->a1 ;
			px->a2 = vp->a1 ;
			px->score = up->chainScore + vp->chainScore ;
		      }
		  }
	      }
	  int jMax = jj ;
	  if (jMax)
	    {
	      ALIGN *wp, *zp ;
	      int m ;
	      arraySort (pairs, pairOrder) ;
	      PAIR *px0 = arrayp (pairs, 0, PAIR), *qx = 0, *rx ;
	      int bestScore = px0->score, kk = 0 ;
	      Array aaa1 = arrayHandleCreate (arrayMax (aaa), ALIGN, h) ;

	      for (jj = 0, px = px0 ; jj < jMax && px->score == bestScore ; jj++, px++)
		{
		  if (qx && px->up == qx->up && px->a2 > qx->a2)
		    continue ; /* eliminate vp2 in   ---> <--- <--- config */
		  BOOL ok = TRUE ;
		  for (m = jj + 1, rx = px+1 ; /* eliminate up1 in   ---> ---> <--- config */
		       ok && m < jMax && rx->a1 < px->a2 && rx->chrom == px->chrom && rx->score == bestScore ;
		       m++, rx++ )
		    ok = FALSE ;
		  if (!ok)
		    continue ;
		  qx = px ;

		  bb->runStat.nAlignedPairs++ ;
		  bb->runStat.nCompatiblePairs++ ;
		  
		  up = px->up ;
		  for (vp = up ; vp->chain == up->chain && vp->read == up->read ; vp++)
		    {
		      wp = px->vp ;
		      vp->pairScore = px->score ;
		      vp->mateChrom = wp->chrom ;
		      vp->mateA1 = wp->chainA1 ;
		      vp->mateA2 = wp->chainA2 ;

		      zp = arrayp (aaa1, kk++, ALIGN) ;
		      *zp = *up ;
		    }	  


		  up = px->vp ;
		  for (vp = up ; vp->chain == up->chain && vp->read == up->read ; vp++)
		    {
		      wp = px->up ;
		      vp->pairScore = px->score ;
		      vp->mateChrom = wp->chrom ;
		      vp->mateA1 = wp->chainA1 ;
		      vp->mateA2 = wp->chainA2 ;

		      zp = arrayp (aaa1, kk++, ALIGN) ;
		      *zp = *up ;
		    }	  
		}
	      arrayMax (aaa) = 0 ;
	      if (kk)
		{
		  zp = arrayp (aaa, kk, ALIGN) ;
		  memcpy (arrp (aaa, 0, ALIGN), arrp (aaa1, 0, ALIGN), kk * sizeof (ALIGN)) ;
		}
	    }
	}
    }
  /* alignDoSelectBestPair */ 
  alignDoRegisterOnePair (pp, bb, aaaa, aaa, read1, bestUp1) ;
  if (bestUp2) alignDoRegisterOnePair (pp, bb, aaaa, aaa, read2, bestUp2) ;
  ac_free (h) ;
} /* alignDoOnePair */

/*
  run -x Aligners/011_SortAlignG5R5/IDX.hs2013.18.81 --maxTargetRepeats 81 -I Fasta/RNA_PolyA_B_1/RNA_PolyA_B_1.Config --align --method 011_SortAlignG5R5  --run RNA_PolyA_B_1  -o RESULTS/011_SortAlignG5R5/RNA_PolyA_B_1/RNA_PolyA_B_1 --step 5 --nRawReads 313433604 --nRawBases 47328474204 --numactl --nB 1
read1 = 144 line=5884
*/
  
/**************************************************************/

static void alignDo (const PP *pp, BB *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT * restrict hit ;
  HIT *h1, *h2 ;
  long int ii, jj, iMax = bigArrayMax (bb->hits) ;
  BigArray hits = bigArrayHandleCreate (256, HIT, h) ;
  BigArray hits2 = bigArrayHandleCreate (256, HIT, h) ;  
  Array aa = arrayHandleCreate (128, ALIGN, h) ;
  Array err = arrayHandleCreate (256, A_ERR, h) ;
  BigArray aaa = bigArrayHandleCreate (iMax, ALIGN, h) ;
  Array countChroms = arrayHandleCreate (256, COUNTCHROM, h) ;
  int n = NTARGETREPEATBITS ;
  int mask = (1 << n) - 1 ;
  
  bb->confirmedIntrons = arrayHandleCreate (64000, INTRON, bb->h) ;
  /*
    bb->isAligned = bitSetHandleCreate (bb->nSeqs, bb->h) ;
  */

  for (ii = 0, hit = iMax ? bigArrp (bb->hits, 0, HIT) : 0 ; ii < iMax ; ii++, hit++)
    {
      int nn = 1, read = hit->read, pair = read >> 1 ;
      for (jj = ii + 1, h1 = hit + 1 ; jj < iMax && (h1->read >> 1) == pair ; jj++, h1++)
	nn++ ;
      if (nn >= 1) /* this read has n+1 hit */
	{ /* create  a copy of the hits of that read */
	  h2 = bigArrayp (hits, nn - 1, HIT) ; /* make room */
	  bigArrayMax (hits) = nn ;
	  h2 = bigArrayp (hits, 0, HIT) ; /* make room */
	  memcpy (h2, hit, nn * sizeof(HIT)) ;

	  /* bb->hits, hence its slice hits, are in hitPairOrder : pair, chrom, position */

	  if (1)
	    {
	      int chrom = 0, mult ;
	      int k, kk = 0, a1 = 0 ;
	      HIT *up  ;
	      COUNTCHROM *zp, *zp0 ;
	      arrayMax (countChroms) = 0 ;

	      /* establish zones */
	      for (k = kk = 0, up = bigArrayp (hits, 0, HIT) ; k < nn ; up++, k++)
		{
		  BOOL isIntron = ((up->x1  >> NTARGETREPEATBITS )  & 0x7) ? TRUE : FALSE ;
		  if (0 && isIntron) continue ;
		  if (up->chrom != chrom || up->a1 > a1 + 1000000)
		    {
		      a1 = up->a1 ;

		      zp = arrayp (countChroms, kk++, COUNTCHROM) ;
		      zp->i1 = zp->i2 = k ;
		      zp->a1 = zp->a2 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
		      zp->chrom = chrom = up->chrom ;
		    }
		  else
		    {  /* hits to this chrom are between [zp->a1,zp->a2] */
		      zp->i2 = k ;
		      zp->a2 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
		      a1 = up->a1 ;
		    }
		}

	      /* compute scores as number of seeds in the read x coordinates */
	      if (kk > 0)
		{
		  for (k = 0, zp = arrayp (countChroms, 0, COUNTCHROM)  ; k < kk ; k++, zp++)
		    {
		      int ii, i1 = zp->i1, i2 = zp->i2 + 1 ;
		      int x0 = -999999999 ;

		      zp->seeds = 0 ;
		      zp->seed1 = 0 ;
		      zp->seed2 = 0 ;
		      zp->seed4 = 0 ;
		      zp->seed8 = 0 ;
		      zp->seed16 = 0 ;
		      zp->seed32 = 0 ;
		      zp->weight = 0 ;
		      
		      bigArraySortSlice (hits, i1, i2, hitReadPosOrder) ;
		      for (ii = i1, up = bigArrp (hits, ii, HIT) ; ii < i2 ; ii++, up++)
			{
			  int x1 = up->x1 >> NSHIFTEDTARGETREPEATBITS ;
			  if (x1 > x0)
			    {
			      x0 = x1 ;
			      mult = up->x1 & mask ;
			      /* zp->weight += mult < 4 ? 720/mult : 720/(4 * mult * mult) ; */
			      if (!zp->seeds) zp->x1 = zp->x2 = x1 ;
			      else zp->x2 = x1 ;
			      zp->seeds++ ;
			      if (mult == 1)       { zp->seed1++ ;  zp->weight += 4 ; }
			      else if (mult == 2)  { zp->seed2++ ;  zp->weight += 3 ; }
			      else if (mult <= 4)  { zp->seed4++ ;  zp->weight += 2 ; }
			      else if (mult <= 8)  { zp->seed8++ ;  zp->weight += 1 ; }
			      else if (mult <= 16) { zp->seed16++ ; zp->weight += 1 ; }
			      else if (mult <= 32) { zp->seed32++ ; zp->weight += 1 ; }
			      if (0) zp->weight = zp->seeds ;
			    }
			}
		    }
		  arraySort (countChroms, countChromOrder) ;
		  /* copy the relevant hits */
		  bigArrayMax (hits2) = 0 ;

		  int mm = 0, k8 = 1 ;
		  zp = zp0 = arrayp (countChroms, 0, COUNTCHROM) ;
		  if (zp->seed1 < zp->seed2) k8 = 2 ;
		  if (zp->seed1 < zp->seed4) k8 = 4 ;
		  if (zp->seed1 < zp->seed8) k8 = 8 ;
		  if (zp->seed1 < zp->seed16) k8 = 8 ;
		  if (zp->seed1 < zp->seed32) k8 = 8 ;
		  if (0 && zp->seed1 < 2 * zp->seed8) k8 = 8 ;
		  for (k = 0 ; k < kk ; k++, zp++)
		    {
		      /* keep at most 2 chromosomes */
		      if (
			  ( k > k8 && zp->weight < zp0->weight && zp->seeds < 3) ||
			  (k >= k8 && 3 * zp->weight < zp0->weight) ||
			  (k < k8 && 4 * zp->weight < zp0->weight)
			  )
			break ;
		      for (int i = zp->i1 ; i <= zp->i2 ; i++)
			{
			  HIT *up = bigArrp (hits, i, HIT) ;
			  HIT *vp = bigArrayp (hits2, mm++, HIT) ;
			  *vp = *up ;
			  vp->chrom ^= (vp->read & 0x1) ;
			}
		    }
		  if (k + 4 < kk) arrayMax(countChroms) = k + 4 ;
		  if (0) arrayMax(countChroms) = 1 ;
		  arrayMax (hits2) = mm ;
		  arrayMax (aa) = arrayMax (err) = 0 ;
		  /* switch chroms and reorder */
		  if (0) hits2 = bb->hits ;
		  saSort (hits2, 2) ; /* hitReadOrder */
		  if (0)  showCountChroms (countChroms) ;

		  alignDoOnePair (pp, bb, aaa, hits2, aa, err) ;
		}
	    }
	}
      hit += nn - 1 ; ii += nn - 1 ;
    }

  bb->aligns = bigArrayHandleCopy (aaa, bb->h) ; /* resize */

  ac_free (h) ;
  return ;
} /* alignDo */

/**************************************************************/

#ifndef YANN
static void align (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->oaChan, &bb, BB))
    {
      if (pp->align && bb.hits)
	{
	  if (pp->debug) printf ("--- %s: Start align %lu seeds\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;

	  t1 = clock () ;
	  alignDo (pp, &bb) ;
	  nnn += bb.nAli ;
	  t2 = clock () ;
	  
	  saCpuStatRegister ("7.Align_r", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  if (pp->debug) printf ("--- %s: Stop align %lu ali, %lu mismatches\n", timeBufShowNow (tBuf), bb.nAli, bb.nerr) ;
	}
      channelPut (pp->aeChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  if (pp->debug) printf ("..... close aeChan at %d,  found %ld ali\n", n, nnn) ;

  channelCloseSource (pp->aeChan) ;

  return ;
} /* align */
#endif

/**************************************************************/
/********************************************************************/

static void sortAlignTableCaption (const PP *pp, ACEOUT ao)
{
  char tBuf[25] ;
  aceOutf (ao, "## %s Magic aligner\n", timeBufShowNow(tBuf)) ;
  if (1)
    aceOutf (ao, "## Author: Danielle et Jean Thierry-Mieg, Greg Boratyn, NCBI, mieg@ncbi.nlm.nih.gov\n") ;
  if (0)
    aceOutf (ao, "## Latest version and documentation are  available from https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/Software\n") ;
  aceOutf (ao, "## Sequences are 1-based\n") ;
  if(0)
    {
      aceOutf (ao, "## 5\'All DNA sequences are read on the strand going away from the edges of the alignment\n") ;
      aceOutf (ao, "## So for donors and acceptors the match part (identical on read and target) and the overhang part are read on the opposite strands,\n"
	       "## The advantage of this unusual convention is that the beginning of all these words is independent of the length of the read, \n"
	       "## independent of the strand of the read or of the target, and the overhang at the donor site is identical to the \n"
	       "## match at the acceptor site and vice versa, so the donor acceptor pairs are easy to recognize.\n"
	       ) ;
    }
  aceOutf (ao, "##Read\tAlignment score\tRead multiplicity") ;
  aceOutf (ao, "\tRead length to be aligned, i.e. removing eventual adaptors, leading Ts or trailing As not matching the target") ;
  aceOutf (ao, "\tRead length aligned") ;
  aceOutf (ao, "\tfrom base") ;
  aceOutf (ao, "\tto base, coordinates of the alignment in the read") ;
  
  
  aceOutf (ao, "\tTarget_class (A_mito ... Z_genome) defines a Bayesian hierarchy of most desirable targets, used by some of the post-processing steps. Theses names are prioritized alphabetically") ;
  aceOutf (ao, "\tGene, as given in the target fasta file as : >target_name|Gene|gene_name, to avoid confusions, RefSeq gene names are prefixed by X__") ;
  aceOutf (ao, "\tUniqueness, number of genes or genomic sites on which the read aligns at best score, if the strand of the best alignements varies, the number is negative, hence -2 likely indicates a single alignment on 2 genes antisense to each other") ;
  aceOutf (ao, "\tTarget, the identifier given in the -t fasta file, for example a chromosome or a trancscritpt") ;
  
  aceOutf (ao, "\tfrom base") ;
  aceOutf (ao, "\tto base, coordinates of the alignment in the target") ;
  
  /*
    if (pp->solid)
    aceOutf (ao, "\tNumber of SOLiD corrected bases") ;
    else
  */
  aceOutf (ao, "\tNumber of N bases") ;
  aceOutf (ao, "\tNumber of mismatches") ;
  
  
  
  aceOutf (ao, "\tMissmatches: read position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the read") ;
  aceOutf (ao, "\tMissmatches: target position:type,... single base variation  (base in target > base in read); single base insertion (+ base inserted in read); single base deletion (- target base missing from read), given in the coordinates and the orientation of the target") ;
  
  aceOutf (ao, "\t5' overhang read backwards, i.e complementary to the unaligned 5' part of the read, max %d bp", pp->OVLN) ;
  aceOutf (ao, "\t3' overhang read forward, i.e. identical to the unaligned 3' part of the read, max %d bp", pp->OVLN) ;
  
  aceOutf (ao, "\tTarget sequence immediately upstream of the alignment") ;
  aceOutf (ao, "\tTarget sequence immediately downstream of the alignment") ;
  aceOutf (ao, "\tFragment length as measured on target in case of paired end sequencing or 0 if not available at this stage of the program") ;
  aceOutf (ao, "\n") ;
  
  aceOutf (ao, "#Run/Read\tScore\tRead multiplicity\tLength to align\tAligned length\tfrom base\tto base in read") ;
  aceOutf (ao, "\tTarget class\tGene (X__name for RefSeq)\tTarget multiplicity\tTarget name (_name for RefSeq)\tfrom base\tto base on target") ;
  
  aceOutf (ao, "\tNumber of %s bases\tNumber of mismatches\tPosition:type in read\tPosition:type in target\t5'overhang (reverse-complemented)\t3'Overhang\tTarget prefix\tTarget suffix\tPair length", pp->solid ? "corrected" : "N") ;
  
  aceOutf (ao, "\n") ;
} /* sortAlignTableCaption */

/**************************************************************/

static void exportDo (const PP *pp, BB *bb)
{
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ii, aMax = bigArrayMax (bb->aligns) ;
  DICT *dict = bb->dict ;
  DICT *errDict = bb->errDict ;
  DICT *dictG = pp->bbG.dict ;
  BOOL pairedEnd = bb->rc.pairedEnd ;
  const char *run = dictName (pp->runDict, bb->run) ;
  AC_HANDLE h = ac_new_handle () ;
  char *runNam = hprintf (h, ".%s.%d.hits", run, bb->lane) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, runNam, pp->gzo, h) ;
  aceOutDate (ao, "###", "sortaling hits") ;
  
  sortAlignTableCaption (pp, ao) ;

  for (ii = 0 ; ii < aMax ; ii++, ap++)
    {
      int read = ap->read ;
      int chrom = ap->chrom ;
      int x1 = ap->x1 ;
      int x2 = ap->x2 ;
      int a1 = ap->a1 ;
      int a2 = ap->a2 ;
      int nerr = ap->nErr ;
      int nN = ap->nN ;

      if (read)
	{
	  aceOutf (ao, "%s/%s%s", run, dictName (dict, read >> 1), (pairedEnd ? ((read & 0x1 )? "<" : ">") : "")) ; 
	  aceOutf (ao, "\t%d", ap->chainScore) ;
	  aceOutf (ao, "\t%d", 1) ; /* ap->multiplicity */
	  aceOutf (ao, "\t%d", ap->readLength) ;
	  aceOutf (ao, "\t%d\t%d\t%d", ap->chainAli, x1, x2) ;
	  
	  aceOutf (ao, "\t%c\t-", ap->targetClass) ;
	  aceOutf (ao, "\t%d", ap->nChains) ; /* ap->targetMultiplicity */
	  aceOutf (ao, "\t%s", dictName (dictG, chrom >> 1)) ; 
	  aceOutf (ao, "\t%d\t%d", a1, a2) ;
	  aceOutf (ao, "\t%d\t%d", nN, nerr) ;
	  aceOutf (ao, "\t%s", ap->errShort ? dictName (errDict, ap->errShort) : "-") ;
	  aceOutf (ao, "\t%s", ap->errLong ? dictName (errDict, ap->errLong) : "-") ;
	  aceOut  (ao, "\t-\t-") ;  /* prefix suffix in genome */
	  aceOutf (ao, "\t%s", ap->leftOverhang ? dictName (dict, ap->leftOverhang) : "-") ;
	  aceOutf (ao, "\t%s", ap->rightOverhang ? dictName (dict, ap->rightOverhang) : "-") ;
	  aceOutf (ao, "\tchain %d 1", ap->chain) ;
	  aceOutf (ao, "\t%d", ap->chainScore) ;
	  aceOut  (ao, "\n") ;
	}
    }

  ac_free (h) ;
  return ;
} /* exportDo */

/**************************************************************/
/**************************************************************/

static void export (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  clock_t  t1, t2, dt = 50 * CLOCKS_PER_SEC ;
  
  char tBuf[25] ;
  int n, nn = 0 ;
  float nG = 0 ;


  AC_HANDLE h = ac_new_handle () ;
  int runMax = dictMax (pp->runDict) ;
  ACEOUT aos[runMax + 1]  ;

  memset (aos, 0, sizeof (aos)) ;
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->aeChan, &bb, BB))
    {
      t1 = clock () ;
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  bigArraySort (bb.aligns, alignOrder) ;
	  if (! pp->sam)
	    exportDo (pp, &bb) ;
	  else
	    {
	      ACEOUT ao = aos[bb.run] ;
	      if (! ao)
		ao = aos[bb.run] = saSamCreateFile (pp, &bb, h) ;
	      saSamExport (ao, pp, &bb) ;
	    }
	  t2 = clock () ;
	  saCpuStatRegister ("8.Export_ali", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  n = channelCount (pp->plChan) ;
	  nG += bb.length / 1.0e+9 ;
	  ++nn ;
	  if (t2 - t1 > dt) 
	    {
	      printf ("--- %s: Exported %d/%d blocks, run %s, %lu reads, %lu ali, %.2f Gbases\n", timeBufShowNow (tBuf), nn, n, dictName (pp->runDict, bb.run), bb.nSeqs, bb.nAli, nG) ;
	      t1 = clock () ;
	    }
	}

      channelPut (pp->doneChan, &bb, BB) ;
    }

  n = channelCount (pp->plChan) ;
  if (0) printf ("--- %s: Export closes doneChan at %d\n", timeBufShowNow (tBuf), n) ;
  channelCloseSource (pp->doneChan) ;
  
  ac_free (h) ;
  return ;
} /* export */

/**************************************************************/

#ifdef YANN
static void wholeWork (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  BB bbG = pp->bbG;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->lcChan, &bb, BB))
    {
      long int nn = 0 ;

      t1 = clock () ;
      /* code words */
      saSequenceParseGzBuffer (pp, &bb) ;
      saCodeSequenceSeeds (pp, &bb, pp->iStep, FALSE) ;

      if (pp->debug) printf ("+++ %s: Start wholeWork %ld bases againt %ld target bases\n", timeBufShowNow (tBuf), bb.length, bbG.length) ;

      /* sort words */
      for (int k = 0 ; k < NN ; k++)
	if (bb.cwsN[k])
	  saSort (bb.cwsN[k], 1) ; /* cwOrder */
      
      /* match hits */
      if (bb.length)
	nn = matchHitsDo (pp, &bbG, &bb) ;
      nnn += nn ;

      /* sorthits */
      if (pp->align && bb.hits)
	{
	  sortHitsFuse (pp, &bb) ;
	  saSort (bb.hits, 3) ; /* hitPairOrder */
	  alignDo (pp, &bb) ;
	}

#ifdef JUNK      
      /* export */
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  bigArraySort (bb.aligns, alignOrder) ;
	  if (! pp->sam)
	    { if (1) exportDo (pp, &bb) ; }
	  else
	    {
	      AC_HANDLE h = ac_new_handle () ;
	      ACEOUT ao = bb.rc.aoSam  ;
	      char *VERSION = "0.1.1" ;
	      DICT *dictG = pp->bbG.dict ;
	      ao = aceOutCreate (pp->outFileName, hprintf (h, ".%s.%d.sam", dictName (pp->runDict, bb.run), bb.lane), pp->gzo, h) ;
	      aceOutf (ao, "@HD VN:1.5\tSO:queryname\n") ;
	      aceOutf (ao, "@PG ID:1\tPN:Magic\tVN:%s\n", VERSION) ;
	      
	      for (int chrom = 1 ; chrom <= dictMax (dictG) ; chrom++)
		{
		  Array dna = arr (pp->bbG.dnas, chrom, Array) ;
		  int ln = dna ? arrayMax (dna) : 0 ;
		  aceOutf (ao, "@SQ\tSN:%s\tLN:%d\n", dictName (dictG, chrom), ln) ;
		}
	      /* aceOutf (ao, "\tCL:%s", commandBuf) ; */
	      exportSamDo (ao, pp, &bb) ;
	      ac_free (h) ;
	    }
	}
#endif
      
      t2 = clock () ;
      saCpuStatRegister ("5.WholeWork", pp->agent, bb.cpuStats, t1, t2, nnn) ;
      channelPut (pp->aeChan, &bb, BB) ;
    }

  channelCloseSource (pp->aeChan) ;

  return ;
} /* wholeWork */
#endif

/*************************************************************************************/

static void reportRunStats (PP *pp, Array runStats)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *METHOD = pp->method ? pp->method : "01_SortAlign" ;
  Array aa = runStats ;
  int ii, iMax = arrayMax (aa) ;
  RunSTAT *s0 = arrayp (aa, 0, RunSTAT) ;
  ACEOUT ao = aceOutCreate (pp->outFileName, ".s2g.samStats", 0, h) ;	
  const char *run = pp->runName ? pp->runName : "xxx" ;

  printf ("\n####### Run Statistics") ;
  fprintf (stderr, "\n####### Run Statistics") ;
  
  printf ("\n# Run") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%s", s->run ? dictName (pp->runDict, s->run) : "Any") ;
    }
  
  printf ("\nFiles") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%d", s->nFiles) ;
    }
  
  printf ("\nPairs") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nPairs) ;
    }
  printf ("\nAlignedPairs") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nAlignedPairs) ;
    }
  
  fprintf (stderr, "AlignedPairs: %ld\n", s0->nAlignedPairs) ;
  
  printf ("\nCompatiblePairs") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nCompatiblePairs) ;
    }
  
  fprintf (stderr, "CompatiblePairs: %ld\n", s0->nCompatiblePairs) ;
  
  printf ("\nCirclePairs") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nCirclePairs) ;
    }
  
  fprintf (stderr, "CirclePairs: %ld\n", s0->nCirclePairs) ;
  
  printf ("\nReads") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nReads) ;
    }
    
  
  printf ("\nBase1") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nBase1) ;
    }
  
  printf ("\nBase2") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nBase2) ;
    }

  printf ("\nA") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->NATGC[1]) ;
    }

  printf ("\nT") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->NATGC[2]) ;
    }

  printf ("\nG") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->NATGC[3]) ;
    }

  printf ("\nC") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->NATGC[4]) ;
    }

  printf ("\nN") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->NATGC[0]) ;
    }

  printf ("\nLn1") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      float nR1 = s->nPairs ? s->nPairs : s->nReads ;
      printf ("\t%.1f", s->nBase1/nR1) ;
    }
  
  printf ("\nLn2") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      float nR2 = s->nPairs ? s->nPairs : 1 ;
      printf ("\t%.1f", s->nBase2/nR2) ;
    }
    
  printf ("\nAlignedReads\t%ld\t%.2f%%", s0->nMultiAligned[0], (100.0 * s0->nMultiAligned[0])/(s0->nReads + .000001)) ;
  printf ("\nMultiAligned\t%ld\t%.2f%%", s0->nMultiAligned[0] - s0->nMultiAligned[1], (100.0 * (s0->nMultiAligned[0] - s0->nMultiAligned[1]))/(s0->nReads + .000001)) ;
  printf ("\nPerfectReads\t%ld\t%.2f%%", s0->nPerfectReads, (100.0 * s0->nPerfectReads)/(s0->nReads + .000001)) ;
  printf ("\nAlignments\t%ld\t%.2f per aligned read", s0->nAlignments, (1.0 * s0->nAlignments)/(s0->nMultiAligned[0] + .000001)) ;
  printf ("\nCompatiblePairs\t%ld\t%.2f%%", s0->nCompatiblePairs, (100.0 * s0->nCompatiblePairs)/(s0->nPairs + .000001)) ;
  printf ("\nCirclePairs\t%ld\t%.2f%%", s0->nCirclePairs, (100.0 * s0->nCirclePairs)/(s0->nPairs + .000001)) ;
  printf ("\nReads supporting introns") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nIntronSupports) ;
    }
  fprintf (stderr, "Reads supporting introns : %ld, plus %ld, minus %ld, stranding %.2f%%\n"
	   , s0->nIntronSupports
	   , s0->nIntronSupportPlus
	   , s0->nIntronSupportMinus
	   , (100.0 * s0->nIntronSupportPlus + 0.00001)/(s0->nIntronSupportPlus+s0->nIntronSupportMinus - 0.00001)
	   ) ;
  
  printf ("\nSupportedIntrons") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	  continue ;
      printf ("\t%ld", s->nSupportedIntrons) ;
    }
  
  
  printf ("\nBase aligned") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nBaseAligned1) ;
    }
  
  printf ("\nBase2 aligned") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nBaseAligned2) ;
    }
  
  
   
  printf ("\nMismatches") ;
  for (ii = 0 ; ii < iMax ; ii++) 
    {	
      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
      if (! s->nReads)
	continue ;
      printf ("\t%ld", s->nErr) ;
    }
  
  for (int j = 0 ; j < 256 ; j++)
    if (s0->nAlignedPerTargetClass[j])
      {
	printf ("\nReads Aligned in %c", j ? j : '-' ) ;
	for (ii = 0 ; ii < iMax ; ii++) 
	  {	
	    RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	    if (! s->nReads)
	      continue ;
	    printf ("\t%ld", s->nAlignedPerTargetClass[j]) ;
	  }
      }
  
  for (int j = 0 ; j < 256 ; j++)
    {
      if (s0->nAlignedPerTargetClass[j])
	{
	  printf ("\nStranding in  %c", j ? j : '-' ) ;
	  for (ii = 0 ; ii < iMax ; ii++) 
	    {	
	      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	      int t = s->GF[j] + s->GR[j]  ;
	      printf ("\t") ;
	      if (t)
		printf ("\t%.3f", 100.0 * s->GF[j]/t) ;
	    }
      }
    }
  
    for (int j = 1 ; j < 11 ; j++)
      if (s0->nMultiAligned[j])
	{
	  printf ("\nReads with %d alignments", j) ;
	  for (ii = 0 ; ii < iMax ; ii++) 
	    {	
	      RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	      if (! s->nReads)
		continue ;
	      printf ("\t%ld", s->nMultiAligned[j]) ;
	    }
	}

  aceOutf (ao, "%s\t%s\tnRawReads\t%ld\n", run, METHOD, pp->nRawReads) ;
  aceOutf (ao, "%s\t%s\tnReads\t%ld\n", run, METHOD, s0->nReads) ;  
  aceOutf (ao, "%s\t%s\tnRawBases\t%ld\n", run, METHOD, pp->nRawBases) ;
  aceOutf (ao, "%s\t%s\tnBases\t%ld\n", run, METHOD, s0->nBase1 + s0->nBase2) ;
  aceOutf (ao, "%s\t%s\tnAlignedPairs\t%ld\n", run, METHOD, s0->nAlignedPairs) ;
  aceOutf (ao, "%s\t%s\tnCompatiblePairs\t%ld\n", run, METHOD, s0->nCompatiblePairs) ;
    aceOutf (ao, "%s\t%s\tnCirclePairs\t%ld\n", run, METHOD, s0->nCirclePairs) ;
  if (! pp->nRawReads) pp->nRawReads = s0->nReads ;
  if (! pp->nRawBases) pp->nRawBases = s0->nBase1 + s0->nBase2 ;
  
  aceOutf (ao, "%s\t%s\tnPerfectReads\t%ld\t%.2f%%\n", run, METHOD, s0->nPerfectReads, (100.0 * s0->nPerfectReads)/(pp->nRawReads + .000001)) ;
  aceOutf (ao, "%s\t%s\tnAlignedReads\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0], (100.0 * s0->nMultiAligned[0])/(pp->nRawReads + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCompatiblePairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCompatiblePairs, (100.0 * s0->nCompatiblePairs)/(pp->nRawReads/2 + .000001)) ;
  aceOutf (ao, "%s\t%s\tnCirclePairs\t%ld\t%.2f%%\n", run, METHOD, s0->nCirclePairs, (100.0 * s0->nCirclePairs)/(pp->nRawReads/2 + .000001)) ;
  aceOutf (ao, "%s\t%s\tnAlignments\t%ld\t%.2f per aligned read\n", run, METHOD, s0->nAlignments, (1.0 * s0->nAlignments)/(s0->nMultiAligned[0] + .000001)) ;
  aceOutf (ao, "%s\t%s\tnMultiAligned\t%ld\t%.2f%%\n", run, METHOD, s0->nMultiAligned[0] - s0->nMultiAligned[1], (100.0 * (s0->nMultiAligned[0] - s0->nMultiAligned[1]))/(s0->nReads + .000001)) ;
  aceOutf (ao, "%s\t%s\tnAlignedBases\t%ld\t%.2f%%\n", run, METHOD, s0->nBaseAligned1 + s0->nBaseAligned2, 100.0 * (s0->nBaseAligned1 + s0->nBaseAligned2) / (pp->nRawBases + .000001)) ;
  aceOutf (ao, "%s\t%s\tnErrors\t%ld\t%.6f%%\n", run, METHOD, s0->nErr, (100.0 * s0->nErr)/(s0->nBaseAligned1 + s0->nBaseAligned2 + 0.00000001)) ;
  long int nUnaligned = pp->nRawReads - s0->nMultiAligned[0] ;
  aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD, 0
	   , nUnaligned
	   , 100.0 * nUnaligned / (pp->nRawReads + .000001)
	   ) ;
  for (int j = 1 ; j < 11 ; j++)
    if (s0->nMultiAligned[j])
      {
	aceOutf (ao, "%s\t%s\tnMultiAligned %d times\t%ld\t%.2f%%\n", run, METHOD
		 , j, s0->nMultiAligned[j]
		 , 100.0 * s0->nMultiAligned[j] / (pp->nRawReads + .000001)
		 ) ;
      }
  long int verif = 0 ;
  for (int j = 1 ; j < 11 ; j++)
    verif += s0->nMultiAligned[j] ;
  aceOutf (ao, "nReads = %ld , sum of multiAli = %ld, verif = %ld\n", pp->nRawReads, verif, pp->nRawReads - verif) ; 
  printf ("\n") ;
  ac_free (h) ;
} /* reportRunStats */

/*************************************************************************************/
/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

void saUsage (char *message, int argc, const char **argv)
{
  int i ;

  if (message)
    {
      fprintf (stderr, "########## ERROR : %s\n", message) ;
      fprintf (stderr, "########## Please try : sortalign -h\n\n") ;
    }
  else if (argc == 0)
    {

      fprintf (stderr,
	       "// Usage: sortalign  <parameters> \n"
	       "//      try: -h --help --version\n"
	       "//\n"
	       "// EXAMPLES:\n"
	       "//      sortalign --createIndex XYZ -t target.fasta (needed once)\n"
	       "//      sortalign --index XYZ -i f_1.fastq.gz+f_2.fastq.gz --wiggle -o results/xxx \n"
	       "//      sortalign --index XYZ -i SRR35876976  -o results/xxx \n"
	       "//\n"
	       "// OBJECTIVE:\n"
	       "//    Sort align maps deep-sequencing files and can export at once alignments, coverage plots, introns.\n"
	       "//      On the first pass, sortalign analyses the (-t or -T) target(s) and creates a binary index.\n"
	       "//      On subsequent calls, sortalign memory maps the index and analyses the (-i or -I) sequence files.\n"
	       "//      All output files use the -o parameter (we recommend to give the name of a directory) as prefix.\n"
	       "//    The amount of data to be analysed can be large (say 100 Gbases), as it is processed in batches.\n"
	       "//      Depending on the  hardware, using 10 to 20 CPUs, it may process up to 1 Gigabases of human RNA-seq per minute\n"
	       "//      On human or mouse, the program requires around (30 + nB) gigabytes of RAM (see --nB --nA options).\n"
	       "//    Sequence files can be provided locall, or automatically downloaded from SRA as in the EXAMPLES above.\n"
	       "//\n"
	       "// TARGET CREATION:\n"
	       "// --createIndex <directory_name>\n"
	       "// --seedLength <int> : default 18, min 10, max 19\n"
	       "// --maxTargetRepeats <int> : default[81], do not index highly repeated target words\n"
	       "// -t <fasta_FileName>\n"
	       "//    Single target, for example a fasta file containing all the chromosomes of a given organism\n"
	       "// -T <target_configuration_fileName>\n"
	       "//    A tab delimited file describing a more complex protocol\n"
	       "//    Each line describes a target class, in a class specific format:\n"
	       "//      Column 1: one of G, M, R, C, A, I, T, B, V : target class.\n"
	       "//         G: Genome fasta files, RNA will be aligned on the genome jumping introns\n"
	       "//         M: Mitochondial fasta file\n"
	       "//         R: Ribosomal RNA fasta file\n"
	       "//         C: Control sequences (i.e. ERCC)\n"
	       "//         A: Adaptor fasta file \n"


	       "//         I: Known Introns and Transcripts,\n"
	       "//            The I file must be of type .gtf, .gff, .gff3 [.gz], or .introns.\n"
	       "//            In each case the the chromosome name (column 1) must match the G fasta file(s).\n"
	       "//            In gtf/gff format the code expects (at least) 7  columns chrom/method/tag/a1/a2/./[+|-]\n"
	       "//               tag : [intron|exon] other lines are ignored,\n"
	       "//               a1 a2 : for introns, the positions of the 2 G of the Gt_aG motif,\n"
	       "//               strand : [+|-] indicates the strand, a1 < a2, and the coordinates are 1-based.\n" 
	       "//            In .introns format the code expects 3 columns chrom/a1/a2\n"
	       "//               a1 a2 : give the positions of the 2 G of the Gt_aG motif,\n"
	       "//               on the top strand a1 < a2, on the bottom strand a1 > a2.\n"
	       "//            Please create the .introns file if the gff/gtf file does not contain intron lines\n"

	       "//         B: Frequent bacteria contaminants fasta file.\n"
	       "//         V: Frequent virus contaminants fasta file.\n"
	       "//      Column 2: f1[,f2] : one or several coma separated file names.\n"
	       "//                For a given class (G, ...) several files can also be specified on several lines.\n"
	       "//      Column 3: optional, the file format if not implied by the file name\n"
	       "//           fasta[default], fastq, fastc, raw.\n"
	       "//\n"
	       "// SEQUENCE FILES  TO BE ANALYZED:\n"
	       "// -x, --index <directory_name> : a directory of index files created  previously using --createIndex\n" 
	       "// -r, --run <runName>   [default -]  : global run name\n"
	       "// -i <sequence_fileName[s]> \n"
	       "//       example 1:   -i f0.fasta,f1.fasta.gz,f2.fastq,f3.fastq.gz,f4.raw,f5.fasta\n"
	       "//        Align files f0 to f5, possibly using different formats.\n"
	       "//       example 2:   -i f1_1.fastq+f1_2.fastq,f2_1.fastq.gz+f2_2.fastq.gz\n"
	       "//        Align read pairs, again possibly using different formats.\n"
	       "//     The file format is implied by the file name, or may be provided explicitly using\n"
	       "//       --raw   | --fasta  | --fastq | --fastc | -sra\n"
	       "//     for example when using  '-i - ' to pipe sequence files into the pipeline\n"
	       "//       example 3: zcat fx.*.fastq.gz | sortalign -x XYZ -i --fastq -o results/fx\n"
	       "// -I <config_fileName>\n"
	       "//     A more precise definition of a set of sequencing files to be analysed\n"
	       "//     Eech line contains 1 to 3 tab separated columns\n"
	       "//          FileName[s]  RunName Descriptors\n"
	       "//     Example:\n"
	       "//          f1.fasta run1 RNA,Nanopore\n"
	       "//          f.R1.gz+f.R2.gz  run2  DNA,fastq,Illumina\n"
	       "//        A table with a single column of sequence files is acceptable, the format will be deduced from the file names\n"
	       "//     1: file name, mandatory\n"
	       "//        For paired end sequencing, provide, as in the second example (f1.F+f1.R), two file names separated by a plus\n"
	       "//        A file named (.fa, .fasta, .fastq, .fastc) implies that format, the file may be gzipped (.gz)\n"
	       "//        A file named (SRR[0-9]*) implies direct SRA download, with optional caching in the ./SRA directory\n"
	       "//     2: Run name, optional,the name must not contain spaces or special characters, as it will be used to name subdirectories\n"
	       "//     3: Descriptors, optional, the options are coma separated, possibilities are\n"
	       "//        DNA/RNA : [default RNA], if DNA sortalign will not clip polyA or jump introns\n"
	       "//        Machine : [default Illumina], one of Illumina, PacBio, Nanopore,..\n"
	       "//        Adaptors=atagg,cctg   run specific adaptor sequences\n"
	       "//        Format : [default fasta], one of sra, raw, fasta, fastq, fastc, only needed if not implied by the file name\n"
	       "//  --sraCaching\n"
	       "//    When the requested sequences are is SRA format they are automatically downloaded from NCBI SRA archive\n"
	       "//    If --sraCaching is set, the (large) downloaded fasta.gz files are copied and saved in the ./SRA local directory\n"
	       "//    This is only useful if you intend to align these sequences several times or reuse them in other programs\n" 
	       "//  --sraDownload SRR123,SRR456,SRR999\n"
	       "//    Just download a set of SRR entries into the cache directory ./SRA/SRR123.fasta.gz ... \n"
	       "//    This command is optional, since sortalign provides direct SRA access\n"
	       "//  --gzi\n"
	       "//    Forces decompression of the input files, useful when piping into sortalign\n"
	       "//    All files named *.gz are automatically decompressed\n"
	       "// OUTPUT:\n"
	       "// -o <outFileNamePrefix> [default stdout]\n"
	       "//   All output files will share this prefix, large outputs are split\n"
	       "//   Examples : -o x/y/z  all output files will be named x/y/z.something\n"
       	       "//              -o x/y/   all output files will be named x/y/something\n"
	       "// --gzo : the output files will be gziped\n"      
	       "// --sam : export the alignment in sam format (default is self documented .hits format)\n"
	       "//\n"
	       "// REQUEST\n"
	       "// --align : [default] Extend the seed alignments hopefully to full read length\n"
	       "// --do_not_align : do not align, just test the validity of the index directory\n"
	       "//     The score s is computed as s = #aligned bases - errCost * #errors\n"
	       "//   --minAli <int> : [default 30] minimal length of the alignment\n"
	       "//   --minScore <int> : [default 30] minimal score  of the alignment\n"
	       "//   --errCost <int> : [default 8] cost of substitition, or short indel up to 3 bases\n"
	       "//   --errMax <int> : [default NA] maximal number of mismatches in any (partial) alignment\n"
	       "//   --errRateMax <int> : [default 10] maximal percentage of mismatches in any (partial) alignment\n"
	       "//   --no_splice : only accept comtinuous alignments, [by default search also spliced alignments]\n"
	       "//   --maxIntron [default 1000000] : max intron size\n"
	       "//   --ignoreIntronSeeds [default FALSE] : do not use the known intron provided in class I in the -T config file\n"
	       "// --wiggle  : Report target coverage wiggles in UCSC BF (fixed) format\n"
	       "// --intron  : Report intron support\n"
	       "// (--snp : not yet ready) Report candidate SNP counts (substitutions and short indels)\n"
	       "// STEPPING\n"
	       "// --step <int>, take a seed every <int> base\n"
	       "//   while createIndex, the default is 1 for targets < 1Mb, 3 for larger targets, 5 would accelerate the code\n"
	       "//   while analysing the reads, the default is 2. --step 1 would increase sensitivity\n"
	       "//\n"
	       "// OPTIONAL TECHNICAL PARAMETERS\n"
	       "//   Unless the program refuses to run, it is probably best not to specify these parameters\n"
	       "//   The program sortalign is parallelized and synchronized by GO-language like channels\n"	
	       "//     This is very different from the way most C programs are parallelized\n"
	       "//     All code layers execute at the same time, the work-load among layers is self balancing\n"
	       "//     and each agent and communication channel runs in its own virtual thread.\n"
	       "//   There is no limit to the size of the input, as data continuously flow in and out\n"
	       "// --nA or --nAgents <int> : [default 10] number of agent in each code layer\n"
	       "// --nB or --nBlocks <int> : [default 20] number of simultaneous data blocks circulating in the pipeline\n"
	       "// --bMax <int> : [default 20] (range 1--1024) max number of mega-bases in a data block.\n"
	       "// --NN <int> : [default 16] split the seed files in NN parts, allowed values (1,2,4,8,16,32,64)\n"
	       "//              seedLength > 16 imply NN >= 4^(sedLength - 16)\n"
	       "// --max_threads <int>  : [default 128] maximal number of simultaneous UNIX threads\n"
	       "//    Do not impose a low value, like --max_threads 8 or 16, the code will not run correctly\n"
	       "//    Possible values are \n"
	       "//        --nBlocks 1 : for a small test\n"
	       "//        --nAgents 10 --nBlocks 10 --max_threads 128 : default, uses less than 32 Gbytes of RAM\n"
	       "//        --nAgents 20 --nBlocks 30 --max_threads 512 : on a large machine\n"
	       "//    These parameters do not limit the amount of data to be processed, they affect the speed\n"
	       "//    If no message is emited within the first minute, the computer is overloaded,\n"
	       "//    please try to increase --max_threads to 256 or 512. alternativelly lower --nAgents and --nBlocks\n"
	       "//      The code is using virtual threads, at least one by agent and by communication channels.\n"

	       "// --verbose : all kinds of details, mostly usefull for debugging, are reported\n"
	       "//\n\n"
	       ) ;
    }
  else
    {
      fprintf (stderr,
	       "########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n##########\t It may have been misspelled or repeated\n") ; 
      fprintf (stderr, "##########\t Please try : sortalign -h\n\n") ;
    }
  exit (1) ;
} /* saUsage */

/*************************************************************************************/
/************ Public interface, get params, set channels and agents, launch **********/
/*************************************************************************************/

int main (int argc, const char *argv[])
{
  PP p ;
  BB bb ;
  Array cpuStats = 0 ;
  Array runStats = 0 ;
  Array runErrors = 0 ;
  Array inArray = 0 ;
  int maxThreads = 128 ;
  long unsigned int nHits = 0, nSeqs = 0, nBaseAligned = 0 ;
  long unsigned int nerr = 0 ;   /* cumulated number of errors */
  long unsigned int aliDx = 0 ; /* cumulated aligned read length */
  long unsigned int aliDa = 0 ;  /* cumulated genome coverage */
  /*   long unsigned int intronSupports = 0 ;  // cumulated intronSupports */ 
  long unsigned int nAli = 0 ;
  BOOL debug = FALSE ;
  int n = 0 ;
  long int skips0 = 0, skips1 = 0, skips2 = 0, skips3 = 0, skips4 = 0, skipsFound = 0, skipsNotFound = 0 ;
  char tBuf[25] ;
  char tBuf0[25] ;
  AC_HANDLE h ;
  int nAgents = 10 ;
  int channelDepth = 10 ;
  mytime_t t0, t1 ;
  
  freeinit () ; 
  messErrorInit (argv[0]) ;

  memset (&p, 0, sizeof (PP)) ;
  h = ac_new_handle () ;
  p.h = h ;

  /**************************  trivial parameters *********************************/
  
  if (argc < 2)
    saUsage (0, 0, argv) ;
  if (getCmdLineBool (&argc, argv, "-help")  ||
      getCmdLineBool (&argc, argv, "--help") ||
      getCmdLineBool (&argc, argv, "-h")
      )
    saUsage (0, 0, argv) ;

  if (getCmdLineBool (&argc, argv, "--version"))
    {
      fprintf (stderr, "sortalign version 0.0.42, novembre 2025") ;
      exit (0) ;
    }     

  p.debug  = getCmdLineText (h, &argc, argv, "--debug", 0) ;
  p.debug |= getCmdLineText (h, &argc, argv, "--verbose", 0) ;
  
  p.gzi = getCmdLineBool (&argc, argv, "--gzi") ;   /* decompress input files (implicit for files named .gz) */
  p.gzo = getCmdLineBool (&argc, argv, "--gzo") ;   /* compress most output files */

  /**************************  SRA downloader *********************************/
  
  {{  /* --sraDownloar SRR1,SRR2,,...
       * Just download SRA files from NCBI/SRA into the local SRA caching directory
       * This is not needed by the aligner, but is provided as a service to the user
       * since it costs us nothing as we have a real-time downloader available
       */
      
      const char *sraID = 0 ;
      if (getCmdLineText (h, &argc, argv, "--sraDownload", &sraID))
	{
	  char *buf = strnew (sraID, h) ;
	  char *cp = buf ;

	  while (cp)
	    {
	      char *cq = strchr (cp, ',') ;
	      if (cq) *cq = 0 ;
	      saSequenceParseSraDownload (cp) ;
	      cp = cq ? cq + 1 : 0 ;
	    }
	  exit (0) ;
	}
    }}

  /******************** NUMA harware  optimizer ******************************************/
  /******************** pin the threads to the least used node ***************************/
  /* numa (non unifirm mmory access) harware are composed of several nodes, each with many CPUs
   * moving data across nodes is slow and costly
   * This module finds the node with the least number of running threads
   * and relaunches the program pinned to that node
   * adding the command line parameter --numactl to prevent recursion
   *
   * This system could be useful in other C programs using multithreading and large memory
   */   

  if (! p.debug &&
      ! getCmdLineBool (&argc, argv, "--numactl")  &&
      !  (getenv("INVOCATION_NOTIFICATIONS") &&
	  strstr(getenv("INVOCATION_NOTIFICATIONS"), "numactl")) &&
      isExecutableInPath ("numactl")   /* see w1/utils.c */
      )
    {
      vTXT txt = vtxtHandleCreate (h) ;
      
      int nodes = -1 ;

      if (0)
	{       /* largest node */
	  FILE *f = fopen("/sys/devices/system/node/online","r");
	  if (f)
	    { /* it seems better to bind to the largest (less used) node */
	      unsigned a,b; 
	      if (fscanf(f,"%u-%u",&a,&b) == 2) nodes = b ; 
	      fclose(f);
	    }
	}
      
      if (1)
	{ /* node with least running threads */

	  int best_node = 0;
	  long long min_load = -1;

	  
	  DIR *d = opendir("/sys/devices/system/node"); 
	  if (d)
	    {	  
	      struct dirent *e;
	      while ((e = readdir(d)))
		{
		  int node;
		  if (sscanf(e->d_name, "node%d", &node) != 1) continue;
		  
		  char path[64];
		  snprintf(path, sizeof(path),
			   "/sys/devices/system/node/node%d/cpumap", node);
		  
		  FILE *f = fopen(path, "r");
		  if (!f) continue;
		  
		  /* count bits = number of online CPUs on this node */
		  unsigned long long map = 0;
		  int n = fscanf(f, "%llx", &map);
		  fclose(f);
		  if (n != 1) continue;
		  
		  long long load = 0;
		  for (unsigned long long m = map; m; m &= m-1) load++;
		  
		  /* count how many tasks are actually scheduled here */
		  load = 0;
		  snprintf (path, sizeof(path), "/proc/schedstat");
		  f = fopen(path, "r");
		  if (f)
		    {
		      char buf[256];
		      while (fgets(buf, sizeof(buf), f))
			{
			  int cpu;
			  if (sscanf(buf, "cpu%d %*s %*d %*d %*d %*d %*d %*d %*d %lld", &cpu, &load) == 2)
			    {
			      /* check if this cpu belongs to our node */
			      char cpu_path[64];
			      snprintf(cpu_path, sizeof(cpu_path),
				       "/sys/devices/system/node/node%d/cpu%d",
				       node, cpu);
			      if (! access(cpu_path, F_OK))
				load += 1;   // one more running thread
			    }
			}
		      fclose(f);
		    }
		  
		  if (min_load == -1 || load < min_load)
		    {
		      min_load = load;
		      best_node = node;
		    }
		}
	      closedir(d);
	    }
	  nodes = best_node;
	}


      vtxtPrintf (txt, "/usr/bin/numactl  --cpunodebind=%d --membind=%d ", nodes, nodes) ;
      for (int i = 0 ; i < argc ; i++)
	vtxtPrintf (txt, " %s " , argv[i]) ;
      vtxtPrintf (txt, " --numactl ") ;

      fprintf (stderr, "%s\n", vtxtPtr (txt)) ;
      return system (vtxtPtr (txt)) ;
    }
  
  /**************************  debugging modules, ignore *********************************/

  if ( getCmdLineInt (&argc, argv, "--checkSorter", &n))
    {                                           /* check the sorting algorithm */
      BigArray b, a = bigArrayHandleCreate (n, HIT, 0) ;

      
      for (int i = 0 ; i < n ; i++)
	bigArray (a, i, HIT).read = ( !(n &0x1) ? (n - i)%7 :  (unsigned int)randint()) ;  

      b = bigArrayHandleCopy (a, 0) ;
      saSort (b, 2) ; /* hitReadOrder */
      
      for (int i = 0 ; i < 50 && i < n ; i++)
	printf("%d %d\n"
	       , bigArr (a,i, HIT).read
	       , bigArr (b,i, HIT).read
	       ) ;
      for (int i = 0 ; i < n-1 ; i++)
	if (bigArr (b,i, HIT).read > bigArr (b,i+1, HIT).read)
	  messcrash ("\nerror line %d\n", i) ;
      showHits (0) ;
      showHitsDo (0, 0) ;
      exit (0) ;
    }

  /**************************  another debugging module, ignore ***************************/
  
  if (getCmdLineInt (&argc, argv, "--makeTest", &n))
    {
      ACEOUT ao = 0 ;
      BigArray aa2 =0, aa = bigArrayHandleCreate (6, int, h) ;
      for (int i = 0 ; i < 6 ; i++)
	bigArray (aa, i, int) = i+1 ;
      bigArrayMapWrite (aa, "mapTest1") ; 
      aa2 = bigArrayMapRead ("mapTest1", int, TRUE, h) ;
      for (int i = 0 ; i < 6 ; i++)
	printf ("\t%d", bigArray (aa2, i, int)) ;
      printf ("\n") ;
      
      /* create a test fasta file with a known 16 mer every 1000 bases */
      ao = aceOutCreate (hprintf (h, "testR.%d", n),  ".fasta", 0, h) ;
      char dna[1000 *n + n] ;
      char w[27] ;
      int i, j, k ;
      char *bb = "atgc" ;
      for (i = k = 0 ; i < n ; i++)
	{
	  for (j = 0 ; j < 26 ; j++)
	    dna[k++] = w[j] = bb[randint () % 4] ;
	  w[26] = 0 ; w[15] = dna[k-11] = 'c' ;
	  aceOutf (ao, ">s.%d\n%s\n", i, w) ;
	  for (j = 0 ; j < 974 ; j++)
	    dna[k++] = bb[randint () % 3] ;
	}
      dna[k++] = 0 ;
      if (n == 1) dna[27] = 0 ; /* expect 10 genome words */
      ao = aceOutCreate (hprintf (h, "testG.%d", n),  ".fasta", 0, h) ;
      aceOutf (ao, ">g.%d\n%s\n", n, dna) ;
      ac_free (h) ;
      exit (0) ;
    }

  /***************** creation of a new target index  *************/
  /* The target index must exist or be created
   *
   * The seed length is only used when creating a new target index
   * otherwise it is inherited from th index
   * The maximal number of repeats of a seed in the targets is used at creation
   * but one is allowed  provide a smaller number when aligning
   * which would drop the index seeds with more repetitions 
   */
  
  p.seedLength = 18 ; /* default */
  getCmdLineInt (&argc, argv, "--seedLength", &(p.seedLength));
  
  p.maxTargetRepeats = 31 ;  /* was 81  31 12 */
  getCmdLineInt (&argc, argv, "--maxTargetRepeats", &p.maxTargetRepeats) ;


  getCmdLineText (h, &argc, argv, "-x", &(p.indexName)) ;
  getCmdLineText (h, &argc, argv, "--index", &(p.indexName)) ;
  getCmdLineText (h, &argc, argv, "-t", &(p.tFileName)) ;
  getCmdLineText (h, &argc, argv, "-T", &(p.tConfigFileName)) ;
  p.createIndex = getCmdLineText (h, &argc, argv, "--createIndex", &(p.indexName)) ;
  
  if (p.createIndex)
    {
      if (! p.tFileName && ! p.tConfigFileName)
	saUsage ("--createIndex requires providing a target parameter -t or -T", argc, argv) ;
      if ( p.tFileName &&  p.tConfigFileName)
	saUsage ("conflicting parameters -t and -T, both define the targets", argc, argv) ;
    }
  
  if (! p.indexName)
    saUsage ("missing parameter -x or --createIndex", argc, argv) ;

  /******************* sequence files to be analysed , not used if --createIndex *************/

  getCmdLineText (h, &argc, argv, "-i", &(p.inFileName)) ;
  getCmdLineText (h, &argc, argv, "-I", &(p.inConfigFileName)) ;

  /* check that input files were provided */
  if (! p.createIndex &&
      ! p.inFileName && ! p.inConfigFileName)
    saUsage ("missing parameter -i or -I inputFileName(s)", argc, argv) ;
  if (p.inFileName && p.inConfigFileName)
    saUsage ("conflicting parameter -i and -I, both define the input files", argc, argv) ;

  /***************** requested analyses, not used if --createIndex *************/

  p.align = getCmdLineBool (&argc, argv, "--align") ;          /* left for back compatibility, over-ridden by --do_not_align */
  p.align = ! getCmdLineBool (&argc, argv, "--do_not_align") ; /* default is to align */
  p.ignoreIntronSeeds = getCmdLineBool (&argc, argv, "--ignoreIntronSeeds") ;

  p.sam = getCmdLineBool (&argc, argv, "--sam") ;
  if (p.sam)
    {
      p.exportSamSequence = TRUE ;
      p.exportSamQuality = TRUE ;
    }
  
  p.sraCaching = getCmdLineBool (&argc, argv, "--sraCaching");   /* cache the files downlaoded from NCBI/SRA */

  /* future options ***/
  p.wiggle = getCmdLineBool (&argc, argv, "--wiggle") ;
  p.snps = getCmdLineBool (&argc, argv, "--snp") ;
  p.introns = getCmdLineBool (&argc, argv, "--intron") ;
  
  /*****************  sequence file names and their formats  ************************/
  
  getCmdLineText (h, &argc, argv, "-o", &(p.outFileName)) ;
  p.runName = 0  ; /* default */

  /* force the sequence file formats,
   * in -I case it is better to provide the format run file by file in the third column of the iConfig file
   * but the best way is not to provide this parameter and let the the format be implied by the file names
          .fasta[.gz]
	  .fastq[.gz]
	  -fna[.gz]
	  .fa[.gz]
	  .fastc[.gz]
	  SRR*
   */
  p.fasta = getCmdLineBool (&argc, argv, "--fasta");
  p.fastq = getCmdLineBool (&argc, argv, "--fastq");
  p.fastc = getCmdLineBool (&argc, argv, "--fastc");
  p.raw = getCmdLineBool (&argc, argv, "--raw");
  p.sra = getCmdLineBool (&argc, argv, "--sra");

  /***************** run name  ***********/

  /* force the run name
   * in -I case it is better to provide the run name file by file in the second column of the tConfig file
   * several files can be associated to the same run
   */
  p.run = 0 ;
  getCmdLineText (h, &argc, argv, "-r", &(p.runName)) ;
  getCmdLineText (h, &argc, argv, "--run", &(p.runName)) ;

  /***************** method name, only used in some output files, convenient when optimizing parameters */
  
  getCmdLineText (h, &argc, argv, "--method", &(p.method)) ;
  
  /*****************  Optional technical parameters ************************/
  /*****************  number of parallel target indexes, only used when createIndex  ************************/

  NN = 16 ; /* default */
  getCmdLineInt (&argc, argv, "--NN", &(NN));
  
  if (NN != 1 && NN != 2 && NN != 4 && NN != 8 && NN != 16 && NN != 32 && NN != 64)
    messcrash ("The number of indexes NN=%d must be  a power of 2, say 1 2 4 8 ...", NN) ;

  if (p.seedLength > 19)
    messcrash ("\n parameter --seedlength %d cannot exceed 19, sorry\n", p.seedLength) ;
  else if (p.seedLength == 19 && NN < 64)
    NN = 64 ;
  else if (p.seedLength == 18 && NN < 16)
    NN = 16 ;
  else if (p.seedLength == 17 && NN < 4)
    NN = 4 ;
  else  if (p.seedLength < 12)
    messcrash ("\n parameter --seedlength %d must be at least 12, sorry\n", p.seedLength) ;

  /* harware check : the code can oly run on 64 bits machnes */
  if (sizeof (long unsigned int) != 8)
    messcrash ("The source code assumes that long unsigned ints use 64 bits not %d, sorry", 8 * sizeof (long unsigned int)) ;

  p.nIndex = NN ;

  /***************** amount or parallelization **************************/
  nAgents = 40 ;
  if (! getCmdLineInt (&argc, argv, "--nAgents", &(nAgents)))
    getCmdLineInt (&argc, argv, "--nA", &(nAgents)) ;

  p.nBlocks = 40 ;  /* max number of BB blocks processed in parallel */
  if (p.debug)      /* under debugger, it is more convenient to run with a single agent and a single block */
    p.nBlocks = 1 ; /* but we wish to be able to reset nBlocks even in debug mode */
  if (!getCmdLineInt (&argc, argv, "--nBlocks", &(p.nBlocks)))
    getCmdLineInt (&argc, argv, "--nB", &(p.nBlocks));
  if (p.nBlocks == 1)
    { nAgents = 1 ; }
  channelDepth = p.nBlocks ;
  maxThreads = 128 ;  /* UNIX  max on lmem12 machine */
  getCmdLineInt (&argc, argv, "--max_threads", &maxThreads) ;
  if (maxThreads < 24)
    maxThreads = 24 ;
  if (p.nBlocks == 1)
    maxThreads = 16 ;  /* if maxThreads == 8 and the single fasta file is split in 3 parts, the code stalls */
  if (p.createIndex)
    { nAgents = 1 ; maxThreads = 1 ; p.nBlocks = 1 ; }


  /****************** stepping when constructing sequence seeds ********************/
  p.iStep = 1 ;   /* read default */
  p.tStep = 0 ;   /* default 1 or 3 set in createIndex or read in existing index */
  if (p.createIndex)
    getCmdLineInt (&argc, argv, "--step", &(p.tStep)) ;
  else
    getCmdLineInt (&argc, argv, "--step", &(p.iStep)) ;

  getCmdLineLong (&argc, argv, "--nRawReads", &(p.nRawReads)) ;
  getCmdLineLong (&argc, argv, "--nRawBases", &(p.nRawBases)) ;

  /*****************  Aligner filters  ************************/
  p.minAli = p.minAliPerCent = p.minScore = -1 ;
  p.errMax = 1000000 ; /* on negative values, extendHits stops on first error */
  p.errRateMax = 10 ;
  p.OVLN = 30 ;
  p.splice = TRUE ;
  if (getCmdLineBool (&argc, argv, "--no_splice"))
    p.splice = FALSE ;
  p.errCost = 8 ;
  getCmdLineInt (&argc, argv, "--errCost", &(p.errCost)) ;
  getCmdLineInt (&argc, argv, "--errMax", &(p.errMax)) ;
  getCmdLineInt (&argc, argv, "--minScore", &(p.minScore)) ;
  getCmdLineInt (&argc, argv, "--minAli", &(p.minAli)) ;
  getCmdLineInt (&argc, argv, "--minAliPerCent", &(p.minAliPerCent)) ;
  getCmdLineInt (&argc, argv, "--errRatMax", &(p.errRateMax)) ;
  p.maxIntron = 1000000 ;
  getCmdLineInt (&argc, argv, "--maxIntron", &(p.maxIntron)) ;

  if (p.minScore < 0) p.minScore = 0 ;
  if (p.minAli < 0) p.minAli = 30 ;
  if (p.minAliPerCent < 0) p.minAliPerCent = 0 ;
  if (p.minAli < p.minScore) p.minAli = p.minScore ;

  p.BMAX = 20 ;
  getCmdLineInt (&argc, argv, "--bMax", &(p.BMAX)) ;
  if (p.BMAX < 1) p.BMAX = 1 ;
  if (p.BMAX > 1024) p.BMAX = 1024 ;
  
  /*****************  Check tat all parameters have been parsed *******************/
  if (argc > 1)
    saUsage (0, argc, argv) ;

  showAli (0) ; /* for compiler happiness */

  /****************** Check the existence of all file names ****************************************/ 

  /* check that the index directory is accessible */
  saConfigIsIndexAccessible (&p) ;

  /* check that the output directoy is writable */
  saConfigIsOutDirWritable (&p) ;   

  /*****************  Start working ***********************************************/
  t0 = timeNow () ;
  printf ("%s: Start\n", timeBufShowNow (tBuf0)) ;

  cpuStats = arrayHandleCreate (1024, CpuSTAT, h) ;
  runStats = arrayHandleCreate (1024, RunSTAT, h) ;
  runErrors = arrayHandleCreate (1024, Array, h) ;
  
  p.runDict = dictHandleCreate (16, p.h) ;
  p.targetClassDict = dictHandleCreate (16, p.h) ;

  dictAdd (p.targetClassDict, "rRNA", 0) ;
  dictAdd (p.targetClassDict, "Mito", 0) ;
  dictAdd (p.targetClassDict, "Genome", 0) ;
  dictAdd (p.targetClassDict, "Bacteria", 0) ;
  dictAdd (p.targetClassDict, "Virus", 0) ;

  /*******************  create the index ********************************************/

  if (p.createIndex)
    { /* The human genome index consumes around 18 Gigabytes of RAM */
      saTargetIndexCreate (&p) ;
      goto done ;
    }

  /*******************  otherwise verify the existence of the indexes ********************/
  
  NN = saConfigCheckTargetIndex (&p) ;

  /* check the existence of the input sequence files */
  inArray = saConfigGetRuns (&p, runStats) ;
  n = dictMax (p.runDict) + 1 ;
  p.runLanes = arrayHandleCreate (n, atomic_int, p.h) ;
  p.runLanesDone = arrayHandleCreate (n, atomic_int, p.h) ;
  array (p.runLanes, n - 1, atomic_int) = 0 ;
  array (p.runLanesDone, n - 1, atomic_int) = 0 ;
  
  /******************** launch the multiprocessing ***************************************/

  wego_max_threads (maxThreads) ;

  if (1)
    {
      /* Create the communication channels */
      p.fpChan = channelCreate (4096, RC, p.h) ; /* this chan nust be deep enough to accept all file names at once */
      channelDebug (p.fpChan, debug, "fpChan") ;
      p.npChan = channelCreate (16, int, p.h) ; /* count BB per inFile */
      channelDebug (p.npChan, debug, "npChan") ;
      p.gmChan = channelCreate (1, BB, p.h) ;
      channelDebug (p.gmChan, debug, "gmChan") ;
      p.plChan = channelCreate (1, BB, p.h) ;
      channelDebug (p.plChan, debug, "plChan") ;
      p.lcChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.lcChan, debug, "lcChan") ;
#ifndef YANN
      p.csChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.csChan, debug, "csChan") ;
      p.smChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.smChan, debug, "smChan") ;
      p.moChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.moChan, debug, "moChan") ;
      p.oaChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.oaChan, debug, "oaChan") ;
#endif
      p.aeChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.aeChan, debug, "aeChan") ;
      p.doneChan = channelCreate (channelDepth, BB, p.h) ;
      channelDebug (p.doneChan, debug, "doneChan") ;
      
      
      /* create the agents using  "wego_go()"
       *   Their execution is  triggered by their in-channels
       *   They export processed BB data-blocks on their out-channel
       *   They exit when their in-channel closes
       *
       * a non-writable COPY of p is passed to each agent
       * this allows to pass agent specific parameters
       * like p.agent : the identifier of an instance of the agent
       */
      
      /* The genome parsers start immediately */
      p.agent = nAgents ;
      wego_go (saTargetIndexGenomeParser, &p, PP) ;
      /* The load regulator maintains --nB data-blocks in the pipeline */
      wego_go (loadRegulator, &p, PP) ;
      
      /* Start the processing of the sequence files */
      p.nFiles = arrayMax (inArray) ; /* number of files to be processed */
      /* n Files must be processed
       * The readParsers will emit for each file, in parallel, n BB blocks
       * The number is passed to the npChan which feeds the npCounter
       * When n numbers have been received, the npCounter
       * knows the cumulated number N of BB blocks to be analyzed
       * and tells the plChan, hence the load regulator,
       * and recurssibvely all program layers
       * to close after having processed N BB blocks
       */
      wego_go (npCounter, &p, PP) ; channelAddSources (p.plChan, 1) ;
      /* Read preprocessing agents, they do not require the genome */
      for (int i = 0 ; i < p.nFiles && i < nAgents && i < 10 ; i++)
	{
	  fprintf (stderr, "Launch readParser %d\n", i) ;
	  p.agent = i ;
	  wego_go (readParser, &p, PP) ; channelAddSources (p.npChan, 1) ;
	}
      for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
	{
	  p.agent = i ;

#ifndef YANN
	  wego_go (codeWords, &p, PP) ; channelAddSources (p.csChan, 1) ;
	  wego_go (sortWords, &p, PP) ; channelAddSources (p.smChan, 1) ;
#endif
	}

      
      /* Load the fpChan, triggering the readParsers
       * and recursively the whole pipeline
       */
      for (int k = 0 ; k < p.nFiles ; k++)
	{
	  RC *rc = arrayp (inArray, k, RC) ;
	  array (runErrors, rc->run, Array) = arrayHandleCreate (32, int, p.h) ;
	  
	  channelPut (p.fpChan, rc, RC) ;
	}
      channelClose (p.fpChan) ;
      
      /* wait untill the genome is ready */
      channelGet (p.gmChan, &p.bbG, BB) ;
      if (! p.bbG.cwsN[0])
	messcrash ("matchHits received no target words") ;
      
      /* map the reads to the genome in parallel */
      for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
	{
	  p.agent = i ;
#ifdef YANN
	  wego_go (wholeWork, &p, PP) ;
	  channelAddSources (p.aeChan, 1) ;
#else

	  wego_go (matchHits, &p, PP) ;
	  channelAddSources (p.moChan, 1) ;
		  
	  wego_go (sortHits, &p, PP) ;
	  channelAddSources (p.oaChan, 1) ;
	  channelAddSources (p.doneChan, 1) ;
	      
	  if (!i || p.align) /* at least 1 agent */
	    {
	      p.agent = 3*i ;
	      wego_go (align, &p, PP) ;
	      p.agent = 3*i + 1 ;
	      wego_go (align, &p, PP) ;
	      p.agent = 3*i + 2 ;
	      wego_go (align, &p, PP) ;
	      p.agent = i ;
	      channelAddSources (p.aeChan, 3) ;
	    }
#endif
	  if (! p.sam || !i) /* only 1 export agent in sam case */
	    {
	      wego_go (export, &p, PP) ;
	      channelAddSources (p.doneChan, 1) ;
	    }
	}
    }

  BigArray confirmedIntrons  = bigArrayHandleCreate (1000, INTRON, h) ;
  int chromMax = dictMax (p.bbG.dict) + 1 ;
  int runMax = dictMax (p.runDict) + 2 ;
  p.wiggles = 0 ;
  if (p.wiggle)
    {
      p.wiggles = arrayHandleCreate (2 * chromMax * runMax, Array, h) ;
      array (p.wiggles, 2 * chromMax * runMax - 1, Array) = 0 ; /* initialize */
    }
  
  while (channelGet (p.doneChan, &bb, BB))
    {
      long int n = (bb.hits ? bigArrayMax (bb.hits) : 0) ;
      int laneToDo ;
      int laneDone ;
      if (bb.isGenome && p.bbG.cwsN)
	{
	  cpuStatCumulate (cpuStats, p.bbG.cpuStats) ;
	  for (int k = 0 ; k < NN ; k++)
	    {
	      if (p.bbG.cwsN) bigArrayDestroy (p.bbG.cwsN[k]) ;
	    }
	  ac_free (p.bbG.cwsN) ;
	  continue ;
	}
      if (p.debug) printf ("%s:Block done\n", timeBufShowNow (tBuf)) ;
      if (p.debug) printf ("Found %ld hits\n", n) ; 

      skips0 += bb.skips0 ;
      skips1 += bb.skips1 ;
      skips2 += bb.skips2 ;
      skips3 += bb.skips3 ;
      skips4 += bb.skips4 ;
      skipsFound += bb.skipsFound ;
      skipsNotFound += bb.skipsNotFound ;

      nHits += n ;
      nSeqs += bb.nSeqs ;
      nBaseAligned += bb.length ;
      nerr += bb.nerr ;
      nAli += bb.nAli ;
      aliDa += bb.aliDa ;
      aliDx += bb.aliDx ;
      if (bb.cpuStats)
	cpuStatCumulate (cpuStats, bb.cpuStats) ;
      bb.runStat.nIntronSupports = bb.confirmedIntrons ?
	arrayMax (bb.confirmedIntrons) : 0 ;
      bb.runStat.nIntronSupportPlus = bb.nIntronSupportPlus ;
      bb.runStat.nIntronSupportMinus = bb.nIntronSupportMinus ;
	
      if (bb.run)
	{
	  runStatsCumulate (0, runStats, &(bb.runStat)) ;
	  runStatsCumulate (bb.run, runStats, &(bb.runStat)) ;

	  runErrorsCumulate (0, runErrors, bb.errors) ;
	  runErrorsCumulate (bb.run, runErrors, bb.errors) ;

	  saIntronsCumulate (confirmedIntrons, bb.confirmedIntrons) ;
	  if (p.wiggle) saWiggleCumulate (&p, &bb) ;
	}
      
      laneToDo = atomic_fetch_add (arrp (p.runLanes, bb.run, atomic_int), 0) + 0 ;
      laneDone = atomic_fetch_add (arrp (p.runLanesDone, bb.run, atomic_int), 1) + 1 ;

      if (laneToDo == laneDone)
	{
	  fprintf (stderr, "run %s done\n", dictName (p.runDict, bb.run)) ;
	}

      /* release block  memory */
      if (bb.dnas)
	{
	  Array dna = 0 ;
	  int iMax = arrayMax (bb.dnas) ;
	  for (int i = 1 ; i < iMax ; i++)
	    {  /* the dna array points into the memory globalDna array */
	      dna = arr (bb.dnas, i, Array) ;
	      if (dna && dna->lock)
		{
		  dna->base = 0 ;
		  arrayUnlock (dna) ;
		  arr (bb.dnas, i, Array) = 0 ;
		}
	    }
	}
      if (bb.quals)
	{
	  Array qual = 0 ;
	  int iMax = arrayMax (bb.quals) ;
	  for (int i = 1 ; i < iMax ; i++)
	    { 
	      qual = arr (bb.quals, i, Array) ;
	      if (arrayExists (qual))
		ac_free (qual) ;
	      arr (bb.quals, i, Array) = 0 ;
	    }
	}
      if (bb.run)
	{
	  int ns = 0 ;
	  char tBuf[25], tBuf2[25] ;
	  bb.stop = timeNow () ;
	  timeDiffSecs (bb.start, bb.stop, &ns) ;
	  printf ("%s: run %d / slice %d done start %s elapsed %d s, nSeqs %ld nBases %.1g\n",  timeBufShowNow (tBuf), bb.run, bb.lane, timeShow (bb.start, tBuf2, 25), ns, bb.nSeqs, (double)bb.length) ; 
	}
      ac_free (bb.h) ;
    }
  cpuStatExport (&p, cpuStats) ;
  runStatExport (&p, runStats) ;
  saIntronsExport (&p, confirmedIntrons) ;
  if (p.wiggle)
    saWiggleExport (&p, nAgents) ;
  
  wego_log ("Done") ;
  wego_flush () ; /* flush the wego logs to stderr */

 done:
  if (1)
    {
      int ns = 0 ;
      printf ("%s: Program start\n",  tBuf0) ;
      t1 = timeNow () ;

      timeDiffSecs (t0, t1, &ns) ;
      printf ("%s: Program end elapsed %d s\n",  timeBufShowNow (tBuf), ns) ; 
    }

  printf ("\tTarget %d sequences %ld bases\n", p.bbG.dnas ? arrayMax (p.bbG.dnas) - 1 : 0, p.bbG.length) ;
  if (1 || p.debug) printf ("Skips: 0=%ld, %d=%ld, %d=%ld, %d=%ld, %d=%ld, found=%ld, notFound=%ld\n",
			    skips0, step1, skips1, step2, skips2, step3, skips3, step4, skips4, skipsFound, skipsNotFound);
  if (1 || p.debug) printf ("SeedLength %d, tStep=%d, iStep=%d, maxTargetRepeats read/target=%d/%d, nAgents=%d nBlocks=%d NN=%d\n"
			    , p.seedLength, p.tStep, p.iStep
			    , p.maxTargetRepeats, p.tMaxTargetRepeats
			    , nAgents, p.nBlocks, NN
			    ) ;
  if (arrayMax (runStats))
    reportRunStats (&p, runStats) ;
  if (p.align)
    reportRunErrors (&p, runStats, runErrors) ;
  
  /* release memory */
  if (p.bbG.dnas)
    {
      Array dna = 0 ;
      int iMax = arrayMax (p.bbG.dnas) ;
      for (int i = 1 ; i < iMax ; i++)
	{  /* the dna array points into the memory mapped globalTargetDna array */
	  dna = arr (p.bbG.dnas, i, Array) ;
	  if (dna->lock)
	    {
	      dna->base = 0 ;
	      arrayUnlock (dna) ;
	    }
	} 
      if (1) /* use valgrind if double freed errors occur */
	ac_free (p.bbG.h) ;
    }
  /* wego_log is the thread-safe way to pass messages to stderr */
  
  if (0)   ac_free (h) ; /* blocks on channel cond destroy */
  return 0 ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
 

