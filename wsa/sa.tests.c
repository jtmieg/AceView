/*
 * sa.seeds.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements several tests
*/

#include "sa.h"



/* create a random genome of so many Megabases
 * add a gt_ag introns at locations n*10kb + 101-102/198-199 // exon 200 300 // ... until 1000
 * export th corresponding gtf and corresponding mrna.fasta
 */

void saCreateRandomGenome (PP *pp, int nMb)
{
  AC_HANDLE h = ac_new_handle() ;
  int iMax = 10000* nMb ;
  int ii, jj ;
  char *atgc = "atgc" ;
  ACEOUT aog = aceOutCreate (pp->outFileName, "/random.genome.fasta", 0, h) ;
  ACEOUT aom = aceOutCreate (pp->outFileName, "/random.mrna.fasta", 0, h) ;
  ACEOUT aof = aceOutCreate (pp->outFileName, "/random.gtf", 0, h) ;
  ACEOUT aoc = aceOutCreate (pp->outFileName, "/random.tConfig", 0, h) ;
  
  if (iMax < 0) messcrash ("nMb = %d too large in saCreateRandomIndex", nMb) ;

  char buf[iMax+10] ;
  char mrna[1200] ;

  aceOutf (aoc, "G\t%s/random.genome.fasta\n", pp->outFileName) ;
  aceOutf (aoc, "I\t%s/random.gtf\n", pp->outFileName) ;
  
  for (ii = 0 ; ii < iMax ; ii++)
    buf[ii] = atgc[randint() % 4] ;
  buf[iMax] = 0 ;
  aceOutf (aog, ">chr1\n%s\n", buf) ;
  for (ii = 0 ; ii < nMb  ; ii++)
    {
      int a1, a2 ;
      for (jj = 0 ; jj < 5 ; jj++)
	{
	  a1 = 10000 * ii + 200 * jj ;
	  a2 = a1 + 100 ;
	  memcpy (mrna+100*jj, buf + a1, 100) ;
	  if (a1 > 2) buf[a1-2] = 'a' ;
	  if (a1 > 2) buf[a1-1] = 'g' ;
	  buf[a2] = 'g' ;
	  buf[a2+1] = 't' ;
	  aceOutf (aof, "chr1\trandom\texon\t%d\t%d\t.\t+\t0\tgene_id \"gene_%d ; transcript_id mrna_%d;\"\n", a1, a2, ii, ii) ;
	  if (jj > 0 && jj <4)
	    aceOutf (aof, "chr1\trandom\tCDS\t%d\t%d\t.\t+\t0\tgene_id \"gene_%d ; transcript_id mrna_%d;\"\n", a1, a2, ii, ii) ;
	}
      for (jj = 0 ; jj < 12 ; jj++)
	mrna[500 + jj] = 'a' ;
      mrna[500 + jj] = 0 ;
      aceOutf (aom, ">s%d\n%s\n", ii, mrna) ;
    }

  ac_free (h) ;
  exit (0) ;
}
