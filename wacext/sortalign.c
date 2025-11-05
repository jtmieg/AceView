/*
 * RNA aligner

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


  #define ARRAY_CHECK
  #define MALLOC_CHECK

/* Example: (for more details try: sortalign -h)
  cd worm_2024_RSMagic1

\rm -rf IDX.WG IDX.WGR
sortalign -t _a1.GR.fasta --createIndex IDX.WGR
sortalign -t _a1.G.fasta --createIndex IDX.WG

run -i _a1.fasta -x IDX.WGR --minAli 30 -o _a3.R.hits --align --max_threads 1

\rm -rf IDX.WX.1
bin/sortalign -T tConfig --createIndex IDX.WX.1 --NN 1

time sortalign -i _unc32.1.fastc -x IDX.WX.1 --minAli 30 -o _unc32.sort.hits --align
time sortalign -i _unc32.fastc -x IDX.WX.1 --minAli 30 -o _unc32.sort.hits --align

 */


#include "ac.h"
#include "channel.h"
#include <dirent.h>
#include <zlib.h>
#include <stdatomic.h>
#include "../wsra/sra_read.h"

#ifdef __SSE2__
#define VECTORIZED_MEM_CPY
#include <emmintrin.h> // SSE2
#endif

#define YANN 1

#define WIGGLE_STEP 10
#define MAXJUMP 3
#define MAXJUMP2 8

typedef struct nodeStruct { double x ; CHAN *cx, *cy, *cu, *cv, *done ; int k ; } NODE ;

typedef enum {FASTA=1, FASTQ, FASTC, RAW, SRA, SRACACHE, SRACACHE1, SRACACHE2, INTRONS} DnaFormat ;
typedef struct targetClassStruct {
  char targetClass ; /* single char a-z, A-Z */
  int bonus ;
  int priority ;
  DnaFormat format ;
  const char *fileName ;
  Array cws ;
} TC ;
		  
typedef struct runClassStruct {
  int run ;   /* index in pp->runDict */
  BOOL pairedEnd ;
  BOOL RNA ; /* 1: is RNA, 0 : is DNA (no introns) */
  int bonus ;
  DnaFormat format ;
  const char *fileName1 ;
  const char *fileName2 ;
  atomic_int lane ;
  ACEOUT aoSam ;
} RC ;
		  
typedef struct runStatStruct {
  int run ;   /* index in p.runDict */
  int nFiles ;
  long int nPairs, nReads ;
  long int nAlignedPairs, nCompatiblePairs, n2ChromsPairs, nOrphans, nCirclePairs ;
  long int nBase1, nBase2 ;
  long int nMultiAligned[11] ;
  long int nAlignedPerTargetClass[256] ;
  long int nPerfectReads ;
  long int nAlignments ;
  long int nPairsAligned, nBaseAligned1, nBaseAligned2 ;
  long int nClippedPolyA ;
  long int nClippedVectors ;
  long int nSupportedIntrons ;
  long int nIntronSupports ;
  long int nIntronSupportPlus ;
  long int nIntronSupportMinus ;    
  int minReadLength, maxReadLength ;
  char *adaptor1, *adaptor2 ;
  long int nN, nErr ;
  long int NATGC[5] ;
  int GF[256], GR[256] ; /* number of reads aligned per target_class on Forward/Reverse strand */
  Array errors ;  /* substitutions, insertions, deletions counts */
  /* coverage of long transcripts ? */
} RunSTAT ;
		  
static int NN = 1 ;
typedef struct bStruct {
  AC_HANDLE h ;
  int run, lane, readerAgent ;
  RC rc ;
  DICT *dict, *errDict ;
  mytime_t start, stop ;
  BigArray dnaCoords ;   /* offSets of the dna in the globalDna array */
  Array dnas ;           /* Array of const char Arrays */
  Array dnasR ;          /* Their reverse complement, only computed for the genome */
  Array quals ;
  BigArray globalDna ;   /* concatenation of all sequence DNAs separated by blocks of nnn */
  BigArray globalDnaR ;  /* concatenation of all reverse sequences, in the same order */
  BigArray hits ;        /* BigArray of read<->genome hits */
  BigArray *cwsN ;         /* BigArray of codeWords */
  Array confirmedIntrons ;
  long unsigned int nSeqs ;  /* number of sequences in bloc */
  long unsigned int length ; /* cumulated number of bases */
  long unsigned int nerr ;   /* cumulated number of errors */
  long unsigned int nAli ;   /* cumulated number of alignments */
  long unsigned int aliDx ;  /* cumulated aligned read length */
  long unsigned int aliDa ;  /* cumulated genome coverage */

  char *gzBuffer ;
  
  /*   BitSet isAligned ; */
  BigArray aligns ; /* final alignments */  
  RunSTAT runStat ;
  Array cpuStats ;
  Array errors ;
  Array wiggles ;
  BOOL isGenome ;
  int step, skips0, skips1, skips2, skips3, skips4, skipsFound, skipsNotFound ;
  long int nIntronSupportPlus ;
  long int nIntronSupportMinus ;    
  vTXT txt1, txt2 ; /* a pair of reusable txt buffer */ 
} BB ;  

typedef struct pStruct {
  AC_HANDLE h ;
  BOOL debug, gzi, gzo ;
  BOOL createIndex ;
  BOOL align ;
  BOOL wiggle ;
  BOOL introns ;
  BOOL snps ;
  BOOL ignoreIntronSeeds ;
  const char *runName ;
  const char *inFileName ;
  const char *inConfigFileName ;
  const char *outFileName ;
  const char *indexName ;
  const char *tFileName ;
  const char *tConfigFileName ;
  const char *tFileBinaryCwsName ;
  const char *tFileBinaryDnaName ;
  const char *tFileBinaryDnaRName ;
  const char *tFileBinaryIdsName ;
  const char *tFileBinaryCoordsName ;

  /* Agents:
     R read parser
     G genome parser
     L load regulator,
     C code words,
     S sort words,
     M match words
     O order hits
     A align in detail
     E export ali
  */
  CHAN *fpChan ; /* RC, Run Config -> readParser */
  CHAN *npChan ; /* int, number of BB emitted for each inFile */
  CHAN *gmChan ; /* Genome is ready, signals the Matcher */
  CHAN *plChan ; /* Parser sends a BB to the load regulator */
  CHAN *lcChan ; /* Load regulator sends a BB to the word Coder */
  CHAN *csChan ; /* Coder sends words to the Sorter */
  CHAN *smChan ; /* Sorter sends words to the Matcher */
  CHAN *moChan ; /* Matcher needs seeds to be Ordered */
  CHAN *oaChan ; /* Ordered matches  sent to the Aligner */
  CHAN *aeChan ; /* Aligner results to be Exported */
  CHAN *wwChan ; /* Wiggler */
  CHAN *doneChan ; /* return to main program */

  BB bbG ;  /* genome or genes target */
  Array runStats ;
  BigArray intronSeeds ;
  BigArray exonSeeds ;
  Array confirmedIntrons ;
  BOOL fasta, fastq, fastc, raw, solid, sra ;
  BOOL sam, exportSamSequence, exportSamQuality ;
  int bonus[256] ;
  DICT *runDict ;
  DICT *targetClassDict ;
  DICT *geneDict ;
  DICT *mrnaDict ;
  DICT *wiggleFileDict ;
  const char *method ;
  int run ;
  int nFiles ;  /* number of input sequence files */
  int agent ;  /* instance of the agent */
  int nBlocks ;
  int blocMaxBases ; /* max number of bases read as one bloc */
  int iStep ;        /* default 2, take a read seed every iStep */
  int tStep ;        /* default 3, take a target seed every tStep */
  int maxTargetRepeats ;
  int tMaxTargetRepeats ;
  int seedLength ;
  int maxIntron ;
  int BMAX ; /* max number of bases in a block, default 200M */
  int errCost ;
  int errMax ;       /* (--align case) max number of errors in seed extension */
  int minScore, minAli, minAliPerCent ;
  int errRateMax ;       /* (--align case) max number of errors in seed extension */
  int OVLN ;
  BOOL splice ;
  long int nRawReads, nRawBases ; 
} PP ;

typedef struct codeWordsStruct {
  unsigned int seed ; /* 32 bits = 16 bases, 2 bits per base */
  int nam ; /* index in readDict or chromDict << 1 | (0x1 for minus words) */
  int pos ;  /* bio coordinate of first letter of seed */
  unsigned int intron ;
} __attribute__((aligned(16))) CW ;

typedef struct hitStruct {
  unsigned int read ;  /* index in readDict */
  unsigned int chrom ; /* index in chromDict << 1 | (0x1 if minus strand) */
  unsigned int a1 ;  /* bio coordinates on chrom (base 1) */
  unsigned int x1 ;  /* bio coordinate on read */
} __attribute__((aligned(16))) HIT ;

typedef struct countChromStruct {
  float weight ;
  int seeds, chrom ;
  int seed1, seed2, seed4, seed8, seed16, seed32 ;
  int i1, i2, x1, x2 ;
  unsigned int a1 ;  /* bio coordinates on chrom (base 1) */
  unsigned int a2 ;  /* bio coordinate on read */
}  COUNTCHROM ;

typedef struct alignStruct {
  int read ;
  int targetClass ;
  int chrom ;
  int a0, x0 ;  /* coordinate of the hit before extension */
  int a1, a2 ;  /* bio position in chrom, a1 < a2 on minus strand */
  int x1, x2 ;  /* bio position in read, x1 < x2 always */
  int w1, w2 ;  /* wiggle coords, rounded, flipped if read2 */
  int chain, chainX1, chainX2, chainA1, chainA2 ;
  int id, previous, next ;
  int ali, chainAli, score, chainScore ;
  int pairScore, mateChrom, mateA1, mateA2, pairLength ;
  int nN, nErr, chainErr ;
  int nTargetRepeats ;
  int nChains ;
  int readLength ;
  int errShort, errLong ; /* bb->dict */
  int leftOverhang, rightOverhang ; /* bb->dict */
  int donor, acceptor ;
  Array errors ;
} __attribute__((aligned(32))) ALIGN ;


typedef struct geneStruct {
  int gene ; /* index in pp->geneDict */
  int chrom ;
  int a1, a2 ; 
} GENE ;
		  
typedef struct mrnaStruct {
  int mrna ; /* index in pp->mrnaDict */
  int chrom ;
  int a1, a2 ; 
} MRNA ;
		  

typedef struct intronStruct {
  int mrna ;
  int chrom ;
  int a1, a2 ; 
} INTRON ;
		  
typedef struct exonStruct {
  int mrna ;
  int chrom ;
  int a1, a2 ; 
} EXON ;
		  

typedef struct cpuStatStruct {
  char nam[32] ;
  int agent ;
  int nB ; /* blocs treated */
  long int n ;
  clock_t tA ; /* time time active */
} CpuSTAT ;
		  
#define step1 256
#define step2 512
#define step3 1024
#define step4 4096
/* was 256 1024 4096 16384 */

#define NTARGETREPEATBITS 5
#define NSHIFTEDTARGETREPEATBITS 8

#define mstep1 255
#define mstep2 510
#define mstep3 1020
#define mstep4 4080
/* was 255 1020 4080 16320 */

#include <pthread.h>
#include <time.h>

typedef struct timespec TMS ;
/* bin/sortalign -t TARGET/Targets/hs.genome.fasta.gz -i titi.fastc --align -o tatou */

static atomic_int lane ;
static int cwOrder (const void *va, const void *vb) ;

/**************************************************************/
/************** utilities *************************************/
/**************************************************************/

static void cpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, long int n)
{
  int s = arrayMax (cpuStats) ;
  CpuSTAT *bs = arrayp (cpuStats, s, CpuSTAT) ;
  bs->nB++ ;
  strncpy (bs->nam, nam, 30) ;
  bs->agent = agent ;
  bs->n += n ;
  bs->tA += t2 - t1 ;  /* active time */
  return ;
} /* cpuStatRegister */

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
/**************************************************************/
/* newMsort algorithm minimizing memcpy */

/* Taquin insertion algorithm
 * works en place
 */
static BOOL newInsertionSort (char *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb))
{
  mysize_t i, j ;
  char buf[s] ;
  BOOL clean = TRUE ;

  for (i = 1 ; i < n ; i++)
    {
      j = i - 1 ;
      if ((*cmp) (b + i*s, b + j*s) >= 0)
	continue ;
      clean = FALSE ;
      memcpy (buf, b + i*s, s) ;
      memcpy (b + i*s, b + j*s, s) ;
      while (j > 0 && (*cmp) (buf, b + (j-1)*s) < 0)
	{
	  memcpy (b + j*s, b + (j-1)*s, s) ;
	  j-- ;
	}
      memcpy (b + j*s, buf, s) ;
    }
  return clean ;
} /* insertionSort */

static void showHitsDo (HIT *hit, long int iMax) ;

/* recursivelly split the table
 * the partially sorted data oscillate between b and buf
 * they end up correctly in b because for small n
 * we switch to the insertion taquin algoright
 * on correct parity, as speed is 2 persent higher with
 * insertion n>0,  relative to n==0
 * but n=8,16,32 are equivalent speeds
 */

static BOOL newMsortDo (char *b, long int nn, int s, char *buf, BOOL hitIsTarget, int (*cmp)(const void *va, const void *vb))
{
 char *up, *vp, *wp ;
  long int n1 = nn / 2 ;
  long int n2 = nn - n1 ;
  char *b1 = b ;
  char *b2 = b + n1 * s ;
  char *b01 = buf ;
  char *b02 = buf + n1 * s ;
  int n = 0 ;
  BOOL ok = FALSE ;
  BOOL clean1, clean2, clean = TRUE ;
  /* for small n,
   * sort en place using the insertion algorithm (game of taquin)
   */
  if (hitIsTarget && nn < 8)
    {
      clean = newInsertionSort (b, nn, s, cmp) ;
      return clean ;
    }

  /* otherwise: sort the 2 halves exchanging hit and buf */
  clean1 = newMsortDo (b01, n1, s, b1, ! hitIsTarget, cmp) ;
  clean2 = newMsortDo (b02, n2, s, b2, ! hitIsTarget, cmp) ;
  
  /* then merge the 2 sorted halves */
  up = b01 ;
  vp = b02 ;
  wp = b1 ;

  if  ((*cmp) (b02 - s, b02) <= 0)
    {
      /* sortmerge is not needed, copy whole blocks */
      /* do we need to copy */
      if (! clean1)
	{ /* copy n1 (or nn idf ! clean2)  records back to b */
	  clean = FALSE ;
	  memcpy(wp, up, (clean2 ? n1 : nn)  * s) ;	  
	}
      else if (! clean2)
	{ /* just copy n2 to the second part of b */
	  clean = FALSE ;
	  memcpy(wp + n1 * s, vp, n2 * s);	  
	}
      /* if clean1 && clean2, no copying is needed */
    return clean  ;
    }
  clean = FALSE ;
  
#ifdef VECTORIZED_MEM_CPY
  /* code generated by Grok, loads and stores 16bytes (128 bits) */
  if (cmp == cwOrder)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;

	  int n = (*(unsigned int*)up <= *(unsigned int*)vp) ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
      ok = TRUE ;
    }
  else if (s == 16)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;
	  
	  int n = ((*cmp)(up, vp) <= 0) ? 1 : 0 ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
      ok = TRUE ;
    }
#endif
  if (! ok) /* either we do not have _mm_store_si128, or size s is not 16 */
    { /* classic code */
      while (n1 > 0 && n2 > 0)
	{
	  n = ((*cmp) (up, vp) <= 0) ? 1 : 0 ;

	  memcpy (wp, (n<=0) ? up : vp, s) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
	}
    }

  /* I also tried to count all greater cases and bulk copy
   * but this code was more complex and not faster
   */

  /* bulk copy the reminders */
    if (n1 > 0) memcpy(wp, up, n1 * s);
    if (n2 > 0) memcpy(wp, vp, n2 * s);

    return clean ;
} /* newMsortDo */

static void newMsort (BigArray aa, int (*cmp)(const void *va, const void *vb))
{
  long int N = bigArrayMax (aa) ;
  char *cp = N ?  (char *) bigArrp (aa, 0, HIT) : 0 ;
  int s = aa->size ;
  
#ifdef VECTORIZED_MEM_CPYzzzzzz
  if (s != 16) /* code now always works */
    messcrash ("only works for aligned 16 bytes structures") ;
#endif
  if (N <= 1) return;
  else if (N < 128)
    newInsertionSort (cp, N, s, cmp) ;
  else
    {
      char *buf = malloc (N * s) ;
      memcpy (buf, cp, N * s) ;
      newMsortDo (cp, N, s, buf, TRUE, cmp) ;
      free (buf) ;
    }
}/* newMsort */

/**************************************************************/
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

static int cwOrder (const void *va, const void *vb)
{
  const CW *up = va ;
  const CW *vp = vb ;
  int n ;
  n = (up->seed > vp->seed) - (up->seed < vp->seed) ; if (n) return n ;
  n = up->nam - vp->nam ; if (n) return n ;
  n = (up->pos > vp->pos) - (up->pos < vp->pos) ; if (n) return n ;
  return 0 ;
} /* cwOrder */

/**************************************************************/

static int intronOrder (const void *va, const void *vb)
{
  const CW *up = va ;
  const CW *vp = vb ;
  int n ;

  /* chrom order */
  n = up->nam - vp->nam ; if (n) return n ;
  /* strand order */
  n = (up->pos > up->intron) - (vp->pos > vp->intron) ; if (n) return n ;
  /* pos order */
  n = (up->pos > vp->pos) - (up->pos < vp->pos) ; if (n) return n ;
  return 0 ;
} /* intronOrder */

/**************************************************************/

static void showCws (const PP *pp, BB *bb, BigArray cws)
{
  long int ii, iMax = bigArrayMax (cws), jj ;
  char buf[17], bufR[17] ;
  const long unsigned int mask26 = (1L << 26) - 1 ;
  
  for (ii = jj = 0 ; ii < iMax && jj < 50 ; ii++)
    {
      CW *cw = bigArrp (cws, ii, CW) ;
      if (1)
	{
	  if (!bb && ((cw->intron >> 31) & 0x1)==0)
	    continue ;  /* select the intron words */
	  else
	    {
	      if (cw->pos >= 0 && (cw->pos < 80 || cw->pos > 110)) continue ;
	      if (cw->pos <  0 && (-cw->pos < 80 || -cw->pos > 110)) continue ;
	    }
	}
      jj++ ;
      for (int i = 0 ; i < 16 ; i++)
	{
	  unsigned int z = (cw->seed >> (30 - 2 * i)) & 0x3 ;
	  switch (z)
	    {
	    case 0x0: buf[i] = 'A' ; break ;
	    case 0x1: buf[i] = 'C' ; break ;
	    case 0x2: buf[i] = 'G' ; break ;
	    case 0x3: buf[i] = 'T' ; break ;
	    }
	}
      buf[16] = 0 ;
      for (int i = 0 ; i < 16 ; i++)
	{
	  unsigned int z = (cw->seed >> (30 - 2 * i)) & 0x3 ;
	  switch (z)
	    {
	    case 0x0: bufR[15-i] = 'T' ; break ;
	    case 0x1: bufR[15-i] = 'G' ; break ;
	    case 0x2: bufR[15-i] = 'C' ; break ;
	    case 0x3: bufR[15-i] = 'A' ; break ;
	    }
	}
      bufR[16] = 0 ;
 
      if (1 || (cw->nam  == 44878))
	{
	  int da1 = 999, da = 0 ;
	  if (!bb && (ii % 256) < 2)
	    continue ;
	  if ((cw->intron >> 31) & 0x1)
	    {
	      da1 =  cw->intron & 0xf ;
	      da  =  ((cw->intron >> 4) & mask26) ;
	    }
	  printf (".. r=%d\t%d\t%u\t%s\t%s\t%s\tii=%ld\tda1=%d\tda=%d\n"
		  , cw->nam,cw->pos, cw->seed
		  , bb ? dictName(bb->dict, cw->nam >> 2) : dictName(pp->bbG.dict, cw->nam >> 1)
		  , buf, bufR, ii
		  , da1, da
		  ) ;
	}
    }
  printf ("......... max %ld\n", iMax) ;
} /* showCws */

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
/* create the GLOBAL DNA
 * Copy all DNAS of the BB block in a single bigArray globalDna
 * Set the individual dna array to point into the globalDna
 * Collect their coordinates 
 */
static void globalDnaCreate (BB *bb)
{
  long int ln = 0 ;
  int n, ii, iMax = arrayMax (bb->dnas) ;
  unsigned char *cp, *cq ;
  Array dna = 0 ;

  if (iMax & 0x1) /* complete the last read pair with a zero */
    {
      array (bb->dnas, iMax, Array) = 0 ;
      iMax++ ;
    }
  
  bb->dnaCoords = bigArrayHandleCreate (2 * (iMax + 1), unsigned int, bb->h) ;

  /* compute the length of the global DNA (avoid reallocations) */
  for (ln = 0, ii = 1 ; ii < iMax ; ii++)
    {
      dna = array (bb->dnas, ii, Array) ;
      if (dna)
	{
	  bigArray (bb->dnaCoords, 2*ii, unsigned int) = ln ;    /* off set of dna ii */
	  n = arrayMax (dna) ;
	  bigArray (bb->dnaCoords, 2*ii + 1, unsigned int) = ln + n ;    /* off set of dna ii */
	  ln += n + 16 ;            /* mininal 16 zeroes protection */
	  n = n % 16 ;
	  if (n)          /* align the data */
	    ln += 16 - n ;
	}
    }
  bigArray (bb->dnaCoords, 2 * ii, unsigned int) = ln ;    /* global end */
  bigArray (bb->dnaCoords, 2 * ii + 1, unsigned int) = ln ;    /* global end */ 
  
  /* construct the global DNA  */
  bb->globalDna = bigArrayHandleCreate (ln, unsigned char, bb->h) ;
  cp = bigArrayp (bb->globalDna, ln -1, unsigned char) ; /* make room */
  for (ln = 0, ii = 1 ; ii < iMax ; ii++)
    {
      /* transfer the sequence to the globalDna array */
      dna = array (bb->dnas, ii, Array) ;
      if (dna)
	{
	  n = arrayMax (dna) ;            /* mininal 16 n protection */
	  cp = bigArrp (bb->globalDna, ln, unsigned char) ;
	  cq = arrp (dna, 0, unsigned char) ;
	  memcpy (cp, cq, n) ; 
	  messfree (dna->base) ;
	  arrayLock (dna) ;
	  dna->base = (char *) cp ; 
	  ln += n ; cp += n ;
	  /* protect each sequence with a terminal A to allow constructing w1 in codeWords */
	  memset (cp, 0, 16) ; /* *cp = A_ */ ; ln += 16 ; cp += 16 ;
	  n = n % 16 ;
	  if (n)          /* align the data */
	    {  memset (cp, 0, 16 - n) ; ln += 16 - n ; cp += 16 - n ; }
	}
    }

  return ;
} /* globalDnaCreate */
  
/**************************************************************/
/**************************************************************/
#define NAMMAX 1024

static BOOL parseOneSequence (DnaFormat format, char *namBuf, ACEIN ai, Array dna, Array qual, int *linep)
{
  int n, nn ;
  char *cp ;
  
  arrayMax (dna) = 0 ;
  if (qual) arrayMax (qual) = 0 ;
  memset (namBuf, 0, NAMMAX) ;

  if (aceInCard (ai))
    {
      (*linep)++ ;
      cp = aceInWord (ai) ;
      switch (format)
	{
	case RAW:
	  if (!cp) return FALSE ;
	  n = strlen ((char *)cp) ;
	  if (!n) return FALSE ;
 	  sprintf (namBuf, "s.%d\n", *linep) ;
	  array (dna, n+1, unsigned char) = 0 ; /* make room */
	  memcpy (arrp (dna, 0, char), cp, n) ;
	  arrayMax (dna) = n ;
	  break ; /* cp already points to the DNA sequence */
	case FASTA:
	  while (! cp || *cp == '#')
	    {
	      cp = 0 ;
	      if (aceInCard (ai))
		{
		  (*linep)++ ;
		  cp = aceInWord (ai) ;
		}
	      else
		return FALSE ;
	    }
	  if (! cp) /* no > identifier in the raninder of the fasta file */
	    return FALSE ;
	  if (cp[0] != '>' || cp[1] == 0)
	    messcrash ("\nMissing identifier at line %d of fasta sequence file %s\n", *linep, aceInFileName (ai)) ; 
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = nn = 0 ;
	  while (aceInCard (ai))
	    { /* parse the sequence, it may spread over many lines */
	      cp = aceInWord (ai) ;
	      if (!cp)
		break ;
	      if (*cp == '>')
		{
		  aceInCardBack (ai) ; /* fold back the indetifier of the next sequence in the ACEIN data flow */
		  break ;
		}
	      (*linep)++ ;
	      n = strlen ((char *)cp) ;
	      if (!n && !nn) return FALSE ;
	      if (1) { char *cq = strstr (cp, "><") ; if (cq) {*cq = 0 ; n = cq - cp ; }}
	      array (dna, nn + n + 1, unsigned char) = 0 ; /* make room */
	      memcpy (arrp (dna, nn, char), cp, n) ;
	      nn += n ;
	      arrayMax (dna) = nn ;
	      }
	  if (! nn)
	    return FALSE ;
	  break ;
	case FASTQ:
	  if (!cp) return FALSE ;
	  if (cp[0] != '@')
	    messcrash ("\nMissing @ identifier, found %s\nat line %d of fastq sequence file %s\n", cp ? cp : "NULL", *linep, aceInFileName (ai)) ;
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = 0 ;
	  if (aceInCard (ai))
	    { /* parse the sequence, it resides on a single line */
	      (*linep)++ ;
	      cp = aceInWord (ai) ;
	      n = cp ? strlen (cp) : 0 ;
	      if (n)
		{
		  array (dna, n + 1, unsigned char) = 0 ; /* make room */
		  memcpy (arrp (dna, 0, char), cp, n) ;
		  arrayMax (dna) = n ;
		}
	    }
	  if (! n)
	    messcrash ("\nMissing sequence at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  if (aceInCard (ai))
	    (*linep)++ ;
	  else
	    messcrash ("\nMissing quality identifier at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  if (aceInCard (ai))
	    {
	      (*linep)++ ;
	      if (qual) 
		{
		  cp = aceInWord (ai) ;
		  n = cp ? strlen (cp) : 0 ;
		  if (n)
		    {
		      array (qual, n + 1, unsigned char) = 0 ; /* make room */
		      memcpy (arrp (qual, 0, char), cp, n) ;
		      arrayMax (qual) = n ;
		    }
		}
	    }
	  else
	    messcrash ("\nMissing quality factors at line %d of fastq sequence file %s\n", *linep, aceInFileName (ai)) ;
	  break ;
	case FASTC:
	  while (! cp || *cp == '#')
	    {
	      cp = 0 ;
	      if (aceInCard (ai))
		{
		  (*linep)++ ;
		  cp = aceInWord (ai) ;
		}		  
	      else
		return FALSE ;
	    }
	  if (! cp) /* no > identifier in the fasta file */
	    return FALSE ;
	  if (cp[0] != '>' || cp[1] == 0)
	    messcrash ("\nMissing identifier at line %d of fastc sequence file %s\n", *linep, aceInFileName (ai)) ; 
	  strncpy (namBuf, cp + 1, NAMMAX - 2) ;
	  n = 0 ;
	  if (aceInCard (ai))
	    { /* parse the sequence, it resides on a single line */
	      (*linep)++ ;
	      cp = aceInWord (ai) ;
	      n = strlen (cp) ;
	      if (n)
		{
		  array (dna, n + 1, unsigned char) = 0 ; /* make room */
		  memcpy (arrp (dna, 0, char), cp, n) ;
		  arrayMax (dna) = n ;
		}
	    }
	  if (! n)
	    messcrash ("\nMissing sequence at line %d of fastc sequence file %s\n", *linep, aceInFileName (ai)) ;
	  break ;
	case SRA:
	  messcrash ("\nParsing SRA is not yet implemented, sorry") ;
	  break ;
	default:
	    messcrash ("\nBad format passed to sequenceParser, please edit the sdource code, sorry") ;
	}
      
      /*
	switch (DnaFormat)
	{
	case RAW:
	  break ;
	case FASTA:
	  break ;
	case FASTQ:
	  break ;
	case FASTC:
	  break ;
	case SRA; 
	  break ;
	}
      */
    }
  return arrayMax (dna) ? TRUE : FALSE ;
} /* parseOneSequence */

/**************************************************************/

static BOOL parseOnePair (DnaFormat format, char *namBuf
			  , ACEIN ai1, Array dna1, Array qual1, int *linep1
			  , ACEIN ai2, Array dna2, Array qual2, int *linep2
			  )
{
  BOOL ok1 = FALSE, ok2 = TRUE ;
  int n ;
  int line1 = *linep1 ;
  int line2 = linep2 ? *linep2 : 0 ;
  
  arrayMax (dna1) = 0 ;
  if (dna2) arrayMax (dna2) = 0 ;
  if (qual1) arrayMax (qual1) = 0 ;
  if (qual2) arrayMax (qual2) = 0 ;
    
  memset (namBuf, 0, NAMMAX) ;
  if (ai1)
    ok1 = parseOneSequence (format, namBuf, ai1, dna1, qual1, linep1) ;
  if (ok1 && ai2)
    {
      char namBuf2 [NAMMAX] ;
      int ln = strlen (namBuf) ;
      ok2 = parseOneSequence (format, namBuf2, ai2, dna2, qual2, linep2) ;
      if (strncmp (namBuf, namBuf2, ln-1))
	messcrash ("\nNon matching pair of sequence identifiers %s <> %s\n line %d of file %s\nlne %d of file %s\n"
		   , namBuf, namBuf2, *linep1, aceInFileName (ai1), *linep2, aceInFileName (ai2)) ;
    }

  if (format == FASTC && ok1)
    {    /* fastc format specifies the pair on a single line with >< separator */
      unsigned char *cp = arrp (dna1, 0, unsigned char) ;
      unsigned char *cq = cp ? (unsigned char *) strstr ((char *)cp, "><") : 0 ;
      if (cq)
	{
	  n = strlen ((char *)cq +2) ;
	  if (n)
	    {
	      array (dna2, n + 2, unsigned char) = 0 ; /* make room */
	      memcpy (arrp (dna2, 0, unsigned char), cq + 2, n) ;
	      arrayMax (dna2) = n ;
	    }
	  *cq = 0 ;
	  n = cq - cp ;
	  arrayMax (dna1) = n ;
	}      
    }
  
  n = arrayMax (dna1) ;
  if (n)
    {
      dnaEncodeArray (dna1) ;
      if (0 &&  /* - (not sequenced) i encoded as 0, i do not dare to change that in w1/dnacode.c  */
	  n != strlen (arrp (dna1, 0, char)))  /* will detect an early  zero */
	messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		   , namBuf, line1, aceInFileName (ai1)) ;
      
    }
  n = dna2 ? arrayMax (dna2) : 0 ;
  if (n)
    {
      dnaEncodeArray (dna2) ;
      if (n != strlen (arrp (dna2, 0, char)))  /* will detect an early  zero */
	{
	  if (ai2)
	    messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		       , namBuf, line2, aceInFileName (ai2)) ;
	  else  /* may happen in fastc case */
	    messcrash ("\nunrecognized DNA character in sequence %s line %d of file %s\n"
		       , namBuf, line1, aceInFileName (ai1)) ;
	}
    }

  return ok1 && ok2 ;
} /* parseOnepair */

/*************************************************************************************/

static int sraDownload (const char *sraID)
{
  AC_HANDLE h = ac_new_handle () ;
  char *fName = hprintf (h, "SRA/%s.fasta", sraID) ;
  ACEOUT ao = 0 ; 

  if (mkdir("./SRA", 0755) == -1)
    {
      if (errno != EEXIST)       /* not "already exists" */
	messcrash ("\nCannot create or cannot write in the SRA cache directory ./SRA") ;
    }
  ao = aceOutCreate (fName, 0, TRUE, h) ;
  if (!ao)
    messcrash ("\nCannot create the SRA cache file %s", fName) ;
  
  SRAObj* sra = SraObjNew(sraID);
  int num_bases = 1 << 27 ; /* 128 M */
  long unsigned int nMax = (0x1L << 31) / num_bases ; /* 2Gb */
  const char *cpp ;
  
  while (nMax-- > 0 && (cpp = SraGetReadBatch(sra, num_bases)))
    aceOutf (ao, "%s", cpp) ;

  SraObjFree(sra);
  ac_free (h) ;
  return 0 ;
}

/**************************************************************/
/* aug 2
 * alternative idea
 * f = gzopen ("gilname.gz", r)
 *  cp = gzread (f, 100Mega, buffer) ;
 *  cq = strrchr (cp, '>')
 *    if (!cq || cq == cp) continue reading intil EOF
 *    else { *(cq-1)=0; pass the cp buffer to a new agent which will decode the dna ;
 *    copy cq ... (n bytes)  to a new clean buffer and gzread in buffer+n
 *  gzclose() 
 * the general idea is that parsing big buffers is fast, while decoding them is slow
 * so in this way, even facinga single large fastq, the pipeline will no longer be hanged on the parser
 */

static void newSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b ;
  int BMAX = isGenome ? 100000 : (pp->BMAX << 20) ;
  unsigned char *buffer = halloc (BMAX, h) ;
  unsigned char *buffer2 = halloc (BMAX, h) ;
  int pos = 0 ;
  BOOL done = FALSE ;
  long int nBytes = 0 ;
  int nPuts = 0 ;

  CHAN *chan = pp->plChan ;
  gzFile file = 0 ;
  BOOL debug = FALSE ;
  
  DnaFormat format = rc->format ;
  const char *fileName1 = rc ? rc->fileName1 : tc->fileName ;
  const char *fileName2 = rc ? rc->fileName2 : 0 ;
  BOOL pairedEnd = rc ? rc->pairedEnd : FALSE ;
  char tBuf[25] ;
  clock_t t1, t2 ;

  if (isGenome || (format != FASTA && format != SRA) || fileName2 || pairedEnd)
    messcrash ("Bad internal options in newSequenceParser, please edit the code, sorry") ;

  if (format == SRA)  /* check in the cache */
    {
      char tBuf[25] ;
      char *cp = hprintf (pp->h, "SRA/%s.fasta.gz", rc->fileName1) ;
      char *cr = filName (cp, 0, "r") ;
      if (!cr)
	{
	  fprintf (stderr, "%s: SRA download start %s\n", timeBufShowNow (tBuf), rc->fileName1) ;
	  sraDownload (rc->fileName1) ;
	  fprintf (stderr, "%s: SRA download done %s\n", timeBufShowNow (tBuf), cp) ;
	  cr = filName (cp, 0, "r") ;
	}
      if (! cr)
	messcrash ("Failed to download from run %s from SRA, try sortalign --help\n", rc->fileName1) ;
      else
	fprintf (stderr, "Found cached file %s\n", cp) ;
      fileName1 = rc->fileName1 = strnew (cr, pp->h) ;
      rc->format = SRACACHE ; 
    }

  
  t1 = clock () ;
  
  file = gzopen (fileName1, "r");
  if (! file)
    messcrash ("\ncannot gzopen target file %s", fileName1) ;

  while (!done)
    {
      int err = 0 ;                    
      int bytes = gzread (file, buffer + pos, BMAX - pos) ;
      unsigned char *restrict cp ;

      bytes += pos ;
      if (bytes < BMAX)
	{
	  done = TRUE ;
	  if (! gzeof (file)) 
	    messcrash ("Error %s in gzread %s", gzerror (file, & err), fileName1) ;                
        }
      else
	{
	  /* search for beginning of last probably partial sequence */
	  cp = buffer + bytes - 1 ; /* last byte read */
	  pos = 0 ;
	  while (cp > buffer && *cp != '>')
	    { pos++ ; cp-- ;}
	  if (*cp != '>' || (cp > buffer && cp[-1] != '\n'))
	    messcrash ("gzread found a read > BMAX=%d", BMAX) ;
	  pos++ ;

	  memcpy (buffer2, cp, pos) ; /* copy the remnant */
	  bytes -= pos ;
	  cp[0] = 0 ;
	}
      nBytes += bytes ;
      if (!nBytes)
	messcrash ("No sequence found in file %s\n", fileName1) ;

      if (format == SRACACHE)
	{   /* check for identifiers signalling a paired end read */
	  unsigned char *cq = buffer ;
	  int nDots = 0 ;
	  while (cq && *cq != '\n')
	    nDots += (*cq++ == '.' ? 1 : 0) ;
	  if (nDots == 3)
	    format = SRACACHE2 ;
	  else
	    format = SRACACHE1 ;
	}
      
      /* create a data block */
      bb = &b ;
      memset (bb, 0, sizeof (BB)) ;
      bb->h = ac_new_handle () ;
	  
      bb->readerAgent = pp->agent ;
      bb->run = rc ? rc->run : 0 ;
      bb->start = timeNow () ;
      bb->lane = atomic_fetch_add (rc ? &(rc->lane) : &lane, 1) + 1 ;
      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->rc.fileName1 = fileName1 ;
      /* copy the buffer */
      bb->gzBuffer = halloc (bytes + 1, bb->h) ;
      memcpy (bb->gzBuffer, buffer, bytes) ;
      bb->gzBuffer[bytes] = 0 ;
      /* position the remnant */
      if (! done) memcpy (buffer, buffer2, pos) ;
      bb->nSeqs = 100 ;  /* a guess */
      /* export the databalock to the channel */
      nPuts++ ;
      channelPut (chan, bb, BB) ;
    }
  channelPut (pp->npChan, &nPuts, int) ; /* global counting of BB blocks accross all sequenceParser agents */
  
  gzclose (file) ;
  ac_free (h) ;
  
  t2 = clock () ;
  cpuStatRegister ("2.NewSequenceParser", pp->agent, bb->cpuStats, t1, t2, nBytes) ;

  if (debug)
     printf ("--- %s: Stop newSequenceParser %d blocks %ld bytes file %s\n", timeBufShowNow (tBuf), lane, nBytes, fileName1) ;
  
  return ;
}

/**************************************************************/

static void oldSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b ;
  int BMAX = isGenome ? 100000 : (pp->BMAX << 20 ) ;
  int NMAX = BMAX / 200 ;
  int nPuts = 0 ;
  int nn = 0, nSeqs = 0 ;
  CHAN *chan = 0 ;
  ACEIN ai1 = 0 ;
  ACEIN ai2 = 0 ;
  int line1 = 0, line2 = 0, line10 = 0 ;
  Array dna1, dna2, dnas ;
  Array qual1 = 0, qual2 = 0 ;
  char namBufG [NAMMAX+12] ;
  char *namBufX, *namBuf = namBufG + 2 ;
  DnaFormat format = rc ? rc->format : tc->format ;
  const char *fileName1 = rc ? rc->fileName1 : tc->fileName ;
  const char *fileName2 = rc ? rc->fileName2 : 0 ;
  BOOL pairedEnd = rc ? rc->pairedEnd : FALSE ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  
  unsigned char natgc[256] ;
  memset (natgc, 0, sizeof(natgc)) ;
  natgc[A_] = 1 ;
  natgc[T_] = 2 ;
  natgc[G_] = 3 ;
  natgc[C_] = 4 ;
  
  t1 = clock () ;
  
  ai1 = aceInCreate (fileName1, 0, h) ;
  if (!ai1)
    messcrash ("\ncannot read target file %s", fileName1) ;
  aceInSpecial (ai1, "\n") ;
  if (fileName2)
    {
      ai2 = aceInCreate (fileName2, 0, h) ;
      if (!ai2)
	messcrash ("\ncannot read target file %s", fileName2) ;
      aceInSpecial (ai2, "\n") ;
    }
  
  if (isGenome)
    {
      chan = 0 ;
      namBufG[0] = tc->targetClass ;
      namBufG[1] = '.' ;
      namBufX = namBufG ;
    }
  else
    {
      bb = &b ;
      bb->readerAgent = pp->agent ;
      bb->start = timeNow () ;
      bb->lane = atomic_fetch_add (rc ? &(rc->lane) : &lane, 1) + 1 ;
      memset (bb, 0, sizeof (BB)) ;
      chan = pp->plChan ;
      namBufX = namBuf ;
    }
  
  if (! bb->h)
    {  /* isGenome: keep expanding the same bbG if already initialised */
      DnaFormat format = rc ? rc->format : tc->format ;
      
      bb->h = ac_new_handle () ;
      bb->txt1 = vtxtHandleCreate (bb->h) ;
      bb->txt2 = vtxtHandleCreate (bb->h) ;
      bb->errors = arrayHandleCreate (256, int, bb->h) ;
      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->run = rc ? rc->run : 0 ;
      bb->length = 0 ;
      bb->dnas = dnas = arrayHandleCreate (64, BigArray, bb->h) ;
      if (pp->exportSamQuality && format == FASTQ)
	bb->quals = arrayHandleCreate (64, Array, bb->h) ;
      bb->dict = dictHandleCreate (NMAX, bb->h) ;
      bb->errDict = dictHandleCreate (NMAX, bb->h) ;
      bb->rc.fileName1 = fileName1 ;
    }
  if (pp->debug || isGenome)
    printf ("+++ %s: Start sequence parser %s\n", timeBufShowNow (tBuf), fileName1) ;
  
  line1 = line2 = 0 ; line10 = 1 ;
  dna1 = arrayHandleCreate (isGenome ? 1 << 20 : 256, unsigned char, bb->h) ;
  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
  if (pp->exportSamQuality && format == FASTQ)
    {
      qual1 = arrayHandleCreate (isGenome ? 1 << 20 : 256, unsigned char, bb->h) ;
      qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
    }
  while (parseOnePair (format, namBuf, ai1, dna1, qual1, &line1, ai2, dna2, qual2, &line2)) 
    {
      int mult = 1 ;
      char *cr = namBuf + strlen (namBuf) ;
      if (format == FASTC && !isGenome)
	{
	  char *cp = strchr (namBuf, '#') ;
	  if (cp)
	    {
	      int k = atoi (cp+1) ;
	      if (k > 1) mult = k ;
	    }
	}
      for (int iMult = 0 ; iMult < mult ; iMult++)
	{
	  memset (cr, 0, 10) ;
	  if (mult > 1)
	    {
	      sprintf (cr, ".%d", iMult+1) ; 
	    }
	  if (pairedEnd || arrayMax(dna2)) 
	    {
	      int k = strlen (namBuf) ;
	      int nn1, n1 = arrayMax (dna1) ;
	      int nn2, n2 = arrayMax (dna2) ;
	      
	      bb->rc.pairedEnd = pairedEnd = TRUE ;
	      nSeqs += 2 ;
	      bb->nSeqs += 2 ;
	      bb->length += n1 + n2 ;
	      bb->runStat.nPairs++ ;
	      bb->runStat.nReads += 2 ;
	      bb->runStat.nBase1 += n1 ;
	      bb->runStat.nBase2 += n2 ;
	      if (namBuf[k-1] == '>') k-- ;
	      /* namBuf[k] = '>' ; namBuf[k+1] = 0 ; */
	      namBuf[k] = 0 ;
	      dictAdd (bb->dict, namBuf, &nn1) ;
	      nn1 <<= 1 ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (bb->quals && qual1)
		array (bb->quals, nn1, Array) = qual1 ;
	      
	      if (iMult < mult - 1)
		{
		  dna1 = arrayHandleCopy (dna1, bb->h) ;
		  if (bb->quals)
		    qual1 = arrayHandleCopy (qual1, bb->h) ;
		}
	      
	      else
		{
		  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		  if (bb->quals)
		    qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		}
		  
	      /*
		namBuf[k] = '<' ;  namBuf[k+1] = 0 ;
		dictAdd (bb->dict, namBuf, &nn2) ;
	      */
	      nn2 = nn1 | 0x1 ;
	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++)
		    bb->runStat.NATGC[natgc[(int)*cp]]++ ;
		}
	      if (arrayMax (dna2))
		{
		  int i, iMax = arrayMax (dna2) ;
		  unsigned char *cp = arrp (dna2, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++)
		    bb->runStat.NATGC[natgc[(int)*cp]]++ ;
		}
	      if (arrayMax (dna2))
		array (bb->dnas, nn2, Array) = dna2 ;
		  if (bb->quals && qual2)
		    array (bb->quals, nn2, Array) = qual2 ;
		  if (iMult < mult - 1)
		    {
		      dna2 = arrayHandleCopy (dna2, bb->h) ;
		      if (bb->quals)
			qual2 = arrayHandleCopy (qual2, bb->h) ;
		    }
		  else
		    {
		      dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
		      if (bb->quals)
			qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
		    }
	    }
	  else
	    {
	      int nn1, n1 = arrayMax (dna1) ;
	      int new = FALSE ;
	      
	      nSeqs++ ;
	      bb->nSeqs++ ;
	      bb->length += n1 ;
	      bb->runStat.nReads++ ;
	      bb->runStat.nBase1 += n1 ;
	      new = dictAdd (bb->dict, namBufX, &nn1) ;
	      nn1 <<= (isGenome ? 0 : 1) ;
	      if (isGenome && ! new)
		messcrash ("\nDuplicate target sequence identifier %s at line %d of file %s\nThe doublon may occur in this file or in a previous target file\n"
			   , namBufX
			   , line10
			   , fileName1
			   ) ;
	      line10 = line1 + 1 ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (arrayMax (dna1))
		{
		  int i, iMax = arrayMax (dna1) ;
		  unsigned char *cp = arrp (dna1, 0, unsigned char) ;
		  for (i = 0 ; i < iMax ; i++)
		    bb->runStat.NATGC[natgc[(int)*cp]]++ ;
		}
	      if (bb->quals && qual1)
		array (bb->quals, nn1, Array) = qual1 ;
	      if (iMult < mult - 1)
		{
		  dna1 = arrayHandleCopy (dna1, bb->h) ;
		  if (qual1)
		    qual1 = arrayHandleCopy (qual1, bb->h) ;
		}
	      else
		{
		  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		  if (pp->sam)
		    qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
		}
	    }
	}
      if (chan && (bb->nSeqs > NMAX || bb->length > BMAX))
	{
	  t2 = clock () ;
	  	  
	  cpuStatRegister ("2.ReadParser", pp->agent, bb->cpuStats, t1, t2, arrayMax (dnas) - 1) ;
	  
	  globalDnaCreate (bb) ;
	  nPuts++ ;
	  channelPut (chan, bb, BB) ;
	  
	  memset (bb, 0, sizeof (BB)) ;
	  t1 = clock () ;
	  
	  nn = 0 ;
	  bb->readerAgent = pp->agent ;
	  bb->start = timeNow () ;
	  bb->lane = atomic_fetch_add (rc ? &(rc->lane) : &lane, 1) + 1 ;
	  bb->h = ac_new_handle () ;
	  bb->txt1 = vtxtHandleCreate (bb->h) ;
	  bb->txt2 = vtxtHandleCreate (bb->h) ;
	  bb->errors = arrayHandleCreate (256, int, bb->h) ;
	  bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
	  bb->run = rc ? rc->run : 0 ;
	  bb->isGenome = isGenome ;
	  bb->nSeqs = 0 ;
	  bb->length = 0 ;
	  bb->dnas = dnas = arrayHandleCreate (nn, BigArray, bb->h) ;
	  if (pp->exportSamQuality && format == FASTQ)
	    bb->quals = arrayHandleCreate (64, Array, bb->h) ;
	  bb->dict = dictHandleCreate (NMAX, bb->h) ;
	  bb->errDict = dictHandleCreate (NMAX, bb->h) ;
	  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
	  if (pp->exportSamQuality && format == FASTQ)
	    {
	      qual1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	      qual2 = arrayHandleCreate (256, unsigned char, bb->h) ;
	    }
	}
    }

  if (! bb->nSeqs)
    messcrash ("\nEmpty sequence file %s\n", fileName1) ;
  

  if (chan)
    {
      globalDnaCreate (bb) ;
      t2 = clock () ;
      
      cpuStatRegister (isGenome ? "1.GParser" : "2.ReadParser", pp->agent, bb->cpuStats, t1, t2, arrayMax (dnas) - 1) ;
      nPuts++  ;
      channelPut (chan, bb, BB) ;
      channelPut (pp->npChan, &nPuts, int) ; /* global counting of BB blocks accross all sequenceParser agents */
    }
  
  /* create the REVERSE COMPLEMENT of the GENOME */
  if (isGenome == 2)
    {                       /* we need the minus strand of the genome */
      int iMax = bb->nSeqs ;
      globalDnaCreate (bb) ;
      bb->globalDnaR = bigArrayHandleCopy (bb->globalDna, bb->h) ;
      
      bb->dnasR = arrayHandleCreate (iMax, Array, bb->h) ;
      unsigned char *cp0 = bigArrayp (bb->globalDnaR, 0, unsigned char) ;
      for (int ii = 1 ; ii <= iMax ; ii++)
	{
	  Array dnaR = arrayHandleCreate (8, unsigned char, bb->h) ;
	  unsigned int x1 = bigArr (bb->dnaCoords, 2*ii, unsigned int) ;      /* offset of this DNA */
	  unsigned int x2 = bigArr (bb->dnaCoords, 2*ii + 1, unsigned int) ;
	  messfree (dnaR->base) ;
	  arrayLock (dnaR) ;
	  dnaR->base = (char *) cp0 + x1 ;
	  dnaR->max = dnaR->dim = x2 - x1 ;
	  reverseComplement (dnaR)  ;             /* complement in place */
	  array (bb->dnasR, ii, Array) = dnaR ;
	}
    }
  
  if (1 || pp->debug) printf ("--- %s: Stop sequence parser %d sequences, %d lines from file %s\n", timeBufShowNow (tBuf), nSeqs, line1, fileName1) ;
  
  ac_free (h) ;
  return ;
} /* oldSequenceParser */

/**************************************************************/

static void sequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  DnaFormat format = rc ? rc->format : tc->format ;
  if (1 && 
      ! isGenome && ! rc->pairedEnd &&  (format == FASTA  || format == SRA))
    return newSequenceParser (pp, rc, tc, bb, isGenome) ;
  else
    return oldSequenceParser (pp, rc, tc, bb, isGenome) ;
} /* sequenceParser */

/**************************************************************/

static void readParser (const void *vp)
{
  const PP *pp = vp ;
  RC rc ;
  
  while (channelGet (pp->fpChan, &rc, RC))
    sequenceParser (pp, &rc, 0, 0, 0) ;
  return ;
} /* Readparser */

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
  channelCloseAt (pp->plChan, nn) ;
  return ;
} /* npCounter */

/**************************************************************/

static BigArray GenomeAddSkips (const PP *pp, BigArray cws, BB *bb, int kk)
{
  long int iMax ; 
  long int jMax ; 
  long int i, j, k ;
  AC_HANDLE h = bb->h ;
  int maxRepeats = pp->maxTargetRepeats ;
    
  BigArray aa ;
  CW *up, *vp, *wp, *upMax ;
  unsigned int wordMax = 0xffffffff ;

  if (0)
    {
      int k = 0 ;
      vp = bigArrp (cws, 0, CW) ;
      for (j = 0 ; k < 6 && j < bigArrayMax(cws) ; j++, vp++)
	{
	  if (vp->nam == 44878)
	    {
	      k++ ;
	      fprintf (stderr, "==++ %d %d %u\n", vp->nam, vp->pos, vp->seed) ;
	    }
	}
    }
  /* remove highly repeated words and register number of repeats */
  if (1)
    {
      long int ks[21], cumul = 0 ;
      up = bigArrp (cws, 0, CW) ;
      vp = bigArrp (cws, 0, CW) ;

      iMax = bigArrayMax (cws) ;
      upMax = up + iMax ;
      memset (ks, 0, sizeof(ks)) ;
      for (i = 0, j = 0 ; i < iMax ; up++, i++)
	{
	  /* int tc = *dictName(pp->bbG.dict,up->nam >> 1) ; */
	  int m, n = 0, nI = 0, nX = 0 ;
	  wp = up ;

	  if (maxRepeats && pp->exonSeeds)
	    {
	      while (wp < upMax && wp->seed == up->seed)
		{
		  nI += ((wp->intron >> 31) & 0x1) ? 1 : 0 ;
		  nX += (((wp->intron >> 30) & 0x3) == 0x1) ? 1 : 0 ;
		  wp++ ;
		}
	    }
	  else
	    {
	      while (wp < upMax && wp->seed == up->seed)
		wp++ ;
	    }
	  n = wp - up ;
	  
	  int nnX = n < (1<<16) ? n : (1<<16) - 1 ;
	  nnX += nX ;
	  if (0 && kk == 1 && up->seed == 185667857)
	    fprintf (stderr, "seed 185667857 kk=%d n=%d, nX=%d, nI=%d\n", kk, n, nX, nI) ;

	  if (!maxRepeats || n < maxRepeats)
	    {
	      for (wp = up, m = 0 ; m < n ; wp++, m++)
		{
		  if (((wp->intron >> 31) & 0x1) == 0x0)
		    wp->intron = n ;
		  if (vp < wp)
		    *vp = *wp ;
		  j++ ; vp++ ;
		}
	    }
	  else if (nX < maxRepeats)
	    {
	      int nnX = n < (1<<16) ? n : (1<<16) - 1 ;
	      nnX += nX ;
	      for (wp = up, m = 0 ; m < n ; wp++, m++)
		{
		  nI = ((wp->intron >> 31) & 0x1) ? 1 : 0 ;
		  nX = (((wp->intron >> 30) & 0x3) == 0x1) ? 1 : 0 ;
		  if (nI + nX)
		    {
		      if (((wp->intron >> 31) & 0x1) == 0x0)
			wp->intron = n ;
		      if (vp < wp)
			*vp = *wp ;
		      j++ ; vp++ ;
		    }
		}
	    }
	  up += n - 1 ; i += n - 1 ;
	  
	  if (n > 20) n = 20 ;
	  ks[n]++ ;
	}
      bigArrayMax (cws) = j ;

      if (0)
	{
	  AC_HANDLE h = ac_new_handle () ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, ".repeated_16_mers_in_target", FALSE, h) ;
	  aceOutf (ao, "#N\tWord\tInstances\tCumul");
	  for (i = 1 ; i <= 20 ; i++)
	    {
	      cumul += i * ks[i] ;
	      aceOutf (ao, "\n%ld\t%ld\t%ld\t%ld", i, ks[i], i * ks[i], cumul) ;
	    }
	  aceOut (ao, "\n") ;
	  ac_free (h) ;
	}
    }
  iMax = bigArrayMax (cws) ;
  if (! iMax) iMax = 1 ; /* insure non void */
  long int jMax0 = iMax + iMax/mstep1 + 1 ;
  aa = bigArrayHandleCreate (jMax0, CW, h) ;
  /* add skipping info */
  up = bigArrp (cws, 0, CW) ;
  vp = bigArrayp (aa, jMax0 - 1, CW) ; 
  vp = bigArrp (aa, 0, CW) ; jMax = 0 ;
  for (long int ii = 0 ; ii < iMax ; ii += mstep1)
    {
      vp->intron = ii + mstep4 < iMax ? (up + mstep4)->seed : wordMax ;
      vp->nam = ii + mstep3 < iMax ? (up + mstep3)->seed : wordMax ;
      vp->pos = ii + mstep2 < iMax ? (up + mstep2)->seed : wordMax ;
      vp->seed = ii + mstep1 < iMax ? (up + mstep1)->seed : wordMax ;

      k = iMax - ii ;
      if (k > mstep1)
	k = mstep1 ;
      vp++ ; 
      memcpy (vp, up, k * sizeof (CW)) ;
      vp += k ; up += k ;
      jMax += k + 1 ; /* k + 1_for_the_jumper */
      if (jMax > jMax0)
	messcrash ("add skipps error ") ;
    }
  bigArrayMax (aa) = jMax ;

  if (0)
    {
      int k = 0 ;
      vp = bigArrp (aa, 0, CW) ;
      for (j = 0 ; k < 6 &&j < jMax ; j++, vp++)
	{
	  if (vp->nam == 44878)
	    {
	      k++ ;
	      fprintf (stderr, "==== %d %d %u\n", vp->nam, vp->pos, vp->seed) ;
	    }
	}
    }
  return aa ;
} /* GenomeAddSkips */

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
/* parse a fasta buffer into an array of DNA */ 
static void parseReadsDo (const PP *pp, BB *bb, int step, BOOL isTarget)
{
  if (bb->gzBuffer) /* fasta buffer */
    {
      AC_HANDLE h = ac_new_handle () ;
      Array qual1 = 0, dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
      int line1 = 0 ; 
      char namBuf [NAMMAX + 12] ;
      ACEIN ai = aceInCreateFromText (bb->gzBuffer, 0, h) ;

      unsigned char natgc[256] ;
      memset (natgc, 0, sizeof(natgc)) ;
      natgc[A_] = 1 ;
      natgc[T_] = 2 ;
      natgc[G_] = 3 ;
      natgc[C_] = 4 ;
      
      memset (namBuf, 0, NAMMAX) ;
      bb->errors = arrayHandleCreate (256, int, bb->h) ;
      bb->txt1 = vtxtHandleCreate (bb->h) ;
      bb->txt2 = vtxtHandleCreate (bb->h) ;
      bb->length = 0 ;
      bb->dnas = arrayHandleCreate (bb->nSeqs, BigArray, bb->h) ;
      bb->dict = dictHandleCreate (bb->nSeqs, bb->h) ;
      bb->nSeqs = 0 ;
      bb->errDict = dictHandleCreate (100000, bb->h) ;
      
      while (parseOnePair (FASTA, namBuf, ai, dna1, qual1, &line1, 0, 0, 0, 0)) 
	{
	  int nn1, n1 = arrayMax (dna1) ;
	  int isRead2 = 0 ;
	  char *cp ;
	  
	  bb->nSeqs++ ;
	  bb->length += n1 ;
	  bb->runStat.nReads++ ;
	  bb->runStat.nBase1 += n1 ;

	  switch (bb->rc.format)
	    {
	    case SRACACHE2:
	      cp = namBuf + strlen (namBuf) - 3 ;
	      if (! strcmp (cp, ".2"))
		{
		  isRead2 = 1 ;
		  *cp = 0 ;
		}
	      if (! strcmp (cp, ".1"))
		{
		  isRead2 = 0 ;
		  *cp = 0 ;
		}
	      break ;
	    default:
	      break ;
	    }
	  
	  dictAdd (bb->dict, namBuf, &nn1) ;
	  nn1 = (nn1 << 1) | isRead2 ;
	  array (bb->dnas, nn1, Array) = dna1 ;
	  
	  if (arrayMax (dna1))
	    {
	      int i, iMax = arrayMax (dna1) ;
	      unsigned char *cp = arrp (dna1, 0, unsigned char) ;
	      for (i = 0 ; i < iMax ; i++)
		bb->runStat.NATGC[natgc[(int)*cp]]++ ;
	    }

	  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	}
      ac_free (h) ;
      globalDnaCreate (bb) ;
    }
  
  return ;
} /* parseReadsDo */

/**************************************************************/
/* if we create N index, suppose N=16, we can mask log_4(NN), i.e. 2, letters
 * effectively indexing 18-mers in 32 bits words
 * using N index also accelerates sorting and searching
 * and consumes less memory for the sorting
 */
/*
actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
 actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
  actctatcacccaggctggagtgcagtggtgccatctcggctcactgcaacctccacctcccaggttcaatcgattctcctgcctcagcctcccgagtagctgggattataggcacccgccaccatgcccggctaatttttatattt
*/
static void codeWordsDo (const PP *pp, BB *bb, int step, BOOL isTarget) 
{
  BigArray cwsN[NN] ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bb->dnas ;
  int k, ia, iaMax = arrayMax (dnas) ;
  long int dMax = bigArrayMax (bb->globalDna) / (NN * step) ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int shiftB = 64 - 2 * wLen ;
  const long unsigned int maskN = NN - 1 ;
  const long unsigned int mask32 = 0xffffffff ; /* 4 bytes integer */
  const long unsigned int maskSeedLn = (1L << 2*wLen) - 1 ;
  long int icwx = 0, icwxMax = pp->exonSeeds ? bigArrayMax (pp->exonSeeds) : 0 ;
  CW *cwX = pp->exonSeeds ? bigArrp (pp->exonSeeds, 0, CW) : 0 ;
  const int seedLength = pp->seedLength ;
  
  memset (cwsN, 0, sizeof (cwsN)) ;
  
  if (step < 1)
    messcrash ("codeWordsDo received step = %d < 1", step) ;
  if (dMax >> 31)
    messcrash ("codeWordsDo received ln=%ls, step=%d, ln/step=%ld > 2G", bb->length, step, dMax) ;

  if (NN > 1)  /* Avoid reallocation of marginal size fluctuations */
    dMax *= 1.2 ;
  for (k = 0 ; k < NN ; k++)
    cwsN[k] = bigArrayHandleCreate (dMax, CW, bb->h) ;
  
  bb->step = step ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;
      int ii, jj, p = 0, cStep, robin ;
      int iMax = arrayMax (dna) ;
      long unsigned int w, wr ;
      const int nShift = 62 ;
      if (1) step = iMax < 60 ? 1 : bb->step ;
      
      w = wr = 0 ;
      /* each q vector is used round-robin style, so no copying is needed */
      BOOL qminus[step] ;
      long unsigned int qw[step], qwr[step], qx[step], qxr[step], qz[step], qhz[step] ;
      memset (qminus, 0, sizeof (qminus)) ;
      memset (qw, 0, sizeof (qw)) ;
      memset (qwr, 0, sizeof (qwr)) ;
      memset (qx, 0, sizeof (qx)) ;
      memset (qxr, 0, sizeof (qxr)) ;
      memset (qz, 0, sizeof (qz)) ;
      memset (qhz, 0, sizeof (qhz)) ;

      cp = arrp (dna, 0, unsigned char) ;

      for (ii = 0, jj = -wLen + 1, cStep = -wLen + 1, p = 0, robin = 0 ; ii < iMax ; ii++, jj++, cp++, cStep++, robin = (robin + 1) % step) 
	{
	  /* construct a 2*wLen bits long integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base */
	  w <<= 2 ; wr >>= 2 ;
	  p++ ;
	  switch ((int)*cp)
	    {  /* alphabetic order and XOR complement */
	    case 0: p = 0 ; break ;
	    case A_:             wr |= 0x3L << nShift ; break ;
	    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
	    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
	    case T_: w |= 0x3L ;                   break ;
	      
	    default: p = 0 ;
	    }

	  qw[robin] = w ;
	  qwr[robin] = wr ;
	  qx[robin] = qw[robin] & maskSeedLn ;
	  qxr[robin] =  (qwr[robin] >> shiftB) & maskSeedLn ;
	  qminus[robin]  = (qx[robin] > qxr[robin]) ;
	  qz[robin] = qminus[robin] ? qxr[robin] : qx[robin] ;
	  qhz[robin] = (qz[robin] ^ (qz[robin] << 14)) ^ (qz[robin] << 23) ; 
	  if (
	      p < wLen ||
	      (((  w >> 2) << 2) == (w << 2) ) || /* avoid homopolymers */
	      (((  w >> 2) << 2) == (w << 2) )    /* avoid homodimers */
	      )
	    qz[robin] = 0 ;

	  
	  if (cStep == step) cStep = 0 ;
	  if (!cStep)
	    {
	      long int word = 0 ;
	      BOOL minus, qok = FALSE ;
	      int dx = 0 ;
	      
	      for (int i = step - 1 ; !qok && i >= 0 ; i--)
		if (qz[robin])
		  qok = TRUE ;
	      if (! qok)
		{ cStep = -1 ; continue ; } /* we hope to export the next word */
	      long int besthz = qhz[robin] ;
	      int bestk = robin ; dx = 0 ;
	      int myMax = (step > jj ? jj : step) ;
	      for (int i = 1 ; i < myMax ; i++)
		{
		  int k = (robin - i + step ) % step ;
		  if (qz[k] && qhz[k] < besthz)
		    { besthz = qhz[k]; bestk = k ; dx = i ; }
		}
	      word = qz[bestk] ;
	      minus = qminus[bestk] ;

	      k = word & maskN ;   /* right bits are used to select the relevant table */
	      word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
	      word &= mask32 ;     /* masked down to 32 bits integer */
	      cw = bigArrayp (cwsN[k], bigArrayMax(cwsN[k]), CW) ;
	      /* store the sign in the name, semantically implying that we treat the 2 strands as independent chromosomes */
	      cw->nam = (ia << 1) | (minus ? 0x1 : 0) ;
	      cw->seed = word ;
	      cw->pos = jj - dx + 1 ;  /* bio coordinates of the first base of the seed */
	      cw->intron = 0 ;
	      if (icwxMax)
		{
		  for ( ; icwx < icwxMax && cwX->nam < ia ; icwx++, cwX++) ;
		  for ( ; icwx < icwxMax && cwX->nam == ia && cwX->intron < cw->pos + seedLength - 1 ; icwx++, cwX++) ;
		  ;
		  if (cwX->nam == ia && cwX->intron >= cw->pos + seedLength - 1 && cwX->pos <= cw->pos)
		    cw->intron = 0x1 << 30 ;		  
		}
	      cStep = dx ;
	      if (0)
		fprintf (stderr, "codeWordsDo p=%d k=%d pos=%d seed=%d\n", p, k, cw->pos,cw->seed) ;
	    }
	}
    }
  bb->cwsN = halloc (NN * sizeof(BigArray), bb->h) ;

  for (k = 0 ; k < NN ; k++)
    bb->cwsN[k] = cwsN[k] ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== CodeWordDo finally allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  return ;
}  /*  codeWordsDo  */

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
      parseReadsDo (pp, &bb, pp->iStep, FALSE) ;
      codeWordsDo (pp, &bb, pp->iStep, FALSE) ;
      t2 = clock () ;

      for (int i = 0 ; i < NN ; i++)
	nn += bigArrayMax (bb.cwsN[i]) ;
      nnn += nn ;
      cpuStatRegister ("3.CodeWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
      if (pp->debug) printf ("--- %s: Stop code words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;

      channelPut (pp->csChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close csChan at %d,  coded %ld words\n", n, nnn) ;
      channelCloseAt (pp->csChan, n) ;
    }

  return ;
} /* codeWords */

/**************************************************************/
/**************************************************************/

static long int genomeParseBinary (const PP *pp, BB *bbG)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  BigArray seqIds = 0;
  DICT *dict = 0 ;
  long int ii, iMax, nn = 0 ;
  
  clock_t t1, t2 ;

  t1 = clock () ;
  /* initialise the bbG (genome) block */
  bbG->h = ac_new_handle () ;
  bbG->cpuStats = arrayHandleCreate (128, CpuSTAT, bbG->h) ;
  bbG->dict = dict = dictHandleCreate (128, bbG->h) ;
  
  /* memory map the target DNA, seed Index, and Identifiers */
  BOOL READONLY = TRUE ;
  /* TRUE: memory mapping,
   * FALSE: read the data from disk into memory,
   *        100s slower in human
   */
  bbG->cwsN = halloc (NN * sizeof (BigArray), bbG->h) ;
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s.%d", pp->tFileBinaryCwsName, k) ;
      bbG->cwsN[k] = bigArrayMapRead (fNam, CW, READONLY, bbG->h) ; /* memory map the seed index */
      nn += bigArrayMax (bbG->cwsN[k]) ;
    }
  
  fNam = pp->tFileBinaryDnaName ;
  bbG->globalDna = bigArrayMapRead (fNam, unsigned char, READONLY, bbG->h) ; /* memory map the DNA */

  fNam = pp->tFileBinaryDnaRName ;
  bbG->globalDnaR = bigArrayMapRead (fNam, unsigned char, READONLY, bbG->h) ; /* memory map the reversed complemented DNA */

  fNam = pp->tFileBinaryCoordsName ;
  bbG->dnaCoords = bigArrayMapRead (fNam, unsigned int, READONLY, bbG->h) ; /* memory map the shared coordinates of the individual chromosomes in the globalDna/globalDnaR arrays */

  fNam = pp->tFileBinaryIdsName ;
   seqIds = bigArrayMapRead (fNam, char, FALSE, h) ; /* memory map the words */

  /* seqids is a char array, we need to transfer it to a dict */
  iMax = bigArrayMax (bbG->dnaCoords)/2 - 2 ;
  char *cp = bigArrp (seqIds, 0, char) ;
  int tStep = bigArr (seqIds, 255, char) ;
  /* check for common divisors */
  if (0 && /* no longer mandatory, since we select the best phasing */
      tStep > 1 && pp->iStep > 1)
    {
      int k ;
      for (int i = 2 ; i <= tStep ; i++)
	{
	  k = tStep/i ;
	  if (k * i == tStep) /* i divides tStep */
	    {
	      k = pp->iStep/i ;
	      if (k * i == pp->iStep) /* i divides iStep */
		messcrash ("\nThe target is indexed with step=%d, the requested read step is step=%d, these number are not relative primes, there will be systematic false negatives,\n please set the argument --istep %d (default 2) to a different value", tStep, pp->iStep, pp->iStep) ;
	    }
	}
    }

  bbG->step = tStep ;
  /* check the version */
  char *signature = hprintf (h, "Sort Align version 1 step %d", tStep) ;
  if (strcmp (cp, signature))
    messcrash ("\nCould not read the correct signature in the index files\n\texpected %s\n\treceived %s\nPlease destroy the index files %s.*.sortali and rerun sortalign --createIndex"
	       , signature
	       , cp
	       , fNam
	       ) ;
  /* create the dictionary of the chromosome identifiers */
  for (ii = 1 ; ii <= iMax ; ii++)
    {
      int n ;
      
      cp += 256 ;
      dictAdd (dict, cp, &n) ;
      if (n != ii)
	messcrash ("\nIndexing error, sorry, please rerun sortalign --createIndex\n") ;
    }

  /* create ancilary target dna arrays,
   * their memory is shared with the globalDna array
   * the coordinate of the seeds refer to the globalDna array
   */
  bbG->length = 0 ;
  bbG->nSeqs = iMax - 1 ;
  bbG->dnas = arrayHandleCreate (iMax + 1, Array, bbG->h) ;
  bbG->dnasR = arrayHandleCreate (iMax + 1, Array, bbG->h) ;
  /* entry zero is fake, because we index via a dictionary */
  array (bbG->dnas, 0, Array) = 0 ;
  array (bbG->dnasR, 0, Array) = 0 ;
  for (ii = 1 ; ii <= iMax ; ii++)
    {
      /* coords[iMax+1] is valid and initialised to bbG->length */
      unsigned int x1 = bigArr (bbG->dnaCoords, 2*ii, unsigned int) ;
      unsigned int x2 = bigArr (bbG->dnaCoords, 2*ii + 1, unsigned int) ;
      Array dna = arrayHandleCreate (8, unsigned char, bbG->h) ;
      Array dnaR = arrayHandleCreate (8, unsigned char, bbG->h) ;

      array (bbG->dnas, ii, Array) = dna ;
      array (bbG->dnasR, ii, Array) = dnaR ;
      /* manipulate the ancilary dna arrays */
      arrayLock (dna) ; /* protect dna->base. It must not be freed */
      messfree (dna->base) ;
      dna->base = bigArrp(bbG->globalDna, x1, char) ; 
      dna->max = dna->dim = x2 - x1 ;
      bbG->length += dna->max ;
      arrayLock (dnaR) ; /* protect dnaR->base */
      messfree (dnaR->base) ;
      dnaR->base = bigArrp(bbG->globalDnaR, x1, char) ; 
      dnaR->max = dnaR->dim = x2 - x1 ;
    }
  if (pp->wiggle)
    {
      bbG->wiggles = arrayHandleCreate (4*iMax, Array, bbG->h) ;
      for (ii = 1 ; ii <= iMax ; ii++)
	{
	  Array dna = array (bbG->dnas, ii, Array) ;
	  int ln = arrayMax (dna)/WIGGLE_STEP + 1 ;
	  array (bbG->wiggles, 4*ii, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 1, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 2, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 4*ii + 3, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	}
    }
  /*  Get thread CPU time at end */
  t2 = clock () ;
  cpuStatRegister ("1.memMapTargets" , pp->agent, bbG->cpuStats, t1, t2, nn) ; 
 
  ac_free (h) ;
  return nn ; 
} /* genomeParseBinary */

/**************************************************************/

static void genomeParser (const void *vp)
{
  const PP *pp = vp ;
  BB bbG = pp->bbG ;
  char tBuf[25] ;
  
  clock_t t2 =0,       t1 = clock () ;
  printf ("+++ %s: Start genome parser\n", timeBufShowNow (tBuf)) ;
  memset (&bbG, 0, sizeof (BB)) ;
  long int nn = genomeParseBinary (pp, &bbG) ;
  t2 = clock () ;
  cpuStatRegister ("1.GParserDone" , pp->agent, bbG.cpuStats, t1, t2, nn) ;
  channelPut (pp->gmChan, &bbG, BB) ;
  channelClose (pp->gmChan) ;
  printf ("--- %s: Stop binary genome parser\n", timeBufShowNow (tBuf)) ;

  return ;
} /* genomeParser */

/**************************************************************/
/**************************************************************/

static void storeTargetIndex (const PP *pp, BB *bbG, int tStep) 
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  long int nn = 0 ;
  /* export the code words */
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s/cws.sortali.%d", pp->indexName, k) ;
      bigArrayMapWrite (bbG->cwsN[k], fNam) ;
      nn += bigArrayMax (bbG->cwsN[k]) ;
    }
  fprintf (stderr, "genomeCreateBinary exported %ld seed records\n", nn) ;

  if (0) showCws (pp, 0, bbG->cwsN[0]) ;
    
  /* export the global DNA */
  fNam = hprintf (h, "%s/dna.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->globalDna, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld target bases\n", bigArrayMax (bbG->globalDna)) ;

  /* export the complement of the global DNA */
  fNam = hprintf (h, "%s/dnaR.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->globalDnaR, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld complemented bases\n", bigArrayMax (bbG->globalDnaR)) ;

  /* memory map the coordinates */
  fNam = hprintf (h, "%s/coords.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->dnaCoords, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld coordinates\n", bigArrayMax (bbG->dnaCoords)) ;

  /* memory map the sequence identifiers (chromosome names) */
  /* transfer the dict to the seqids char array */
  int iMax = dictMax (bbG->dict) ;
  BigArray seqIds = bigArrayHandleCreate (256 * (iMax + 1), char, h) ;
  char *cp = bigArrayp (seqIds, 256 * (iMax + 1) - 1, char) ; /* make room */
  char buf[256] ;

  /* global title */
  cp = bigArrayp (seqIds, 0, char) ;
  memset (buf, 0, 256) ;
  char *signature = hprintf (h, "Sort Align version 1 step %d", tStep) ;
  strcpy (buf, signature) ;
  memcpy (cp, buf, 255) ;
  cp[255] = tStep ;
  /* list of names */
  for (int i = 1 ; i <= iMax ; i++)
    {
      cp += 256 ;
      memset (buf, 0, 256) ;
      strncpy (buf, dictName (bbG->dict, i), 255) ;
      memcpy (cp, buf, 256) ;
    }
  fNam = hprintf (h, "%s/ids.sortali", pp->indexName) ;
  bigArrayMapWrite (seqIds, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %d identifiers\n", dictMax (bbG->dict)) ;

  ac_free (h) ;
} /* storeTargetIndex */

/**************************************************************/
/**************************************************************/

static long int gffParser (PP *pp, BB *bbG, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nnE = 0, nnI = 0 ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom= 0 ;
  int a1, a2 ;
  int line = 0 ;
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

      if (isIntrons) fprintf (stderr, "Parsing gff file %s\n", cp) ;
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

      if (! isIntrons)
	{
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* method, drop it */

	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;  /* type */
	  type = 0 ;
	  if (! strcasecmp (cp, "exon")) type = 1 ;
	  else if (! strcasecmp (cp, "intron")) type = 2 ;
	  if (! type) continue ;
	}
      else
	type = 2 ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;

      if (! isIntrons)
	{
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
	}

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
      bigArraySort (pp->exonSeeds, cwOrder) ;

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

static long int intronParser (PP *pp, BB *bbG, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nn = 0 ;
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
      up->nam = chrom ;
      up->pos = a1 ;
      up->intron = a2 ;
    }
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, aceInFileName (ai)) ;
  ac_free (h) ;

  return nn ;
} /* intronParser */

/**************************************************************/

static int intronCodeWords (PP *pp, BB *bbG)
{
  BOOL debug = FALSE ;
  CW *restrict cw  ;  
  const unsigned char *restrict cp ;
  Array dnas = bbG->dnas ;
  int k, p ;
  long int nn = 0 ;
  const int wLen = pp->seedLength ;
  const int nHidden = wLen > 16 ? wLen - 16 : 0 ;
  const int nHidden2 = nHidden << 1 ;
  const int nShift = 62 ;
  const long unsigned int maskN = NN - 1 ;
  const long unsigned int mask32 = 0xffffffff ; /* 4 bytes integer */
  const long unsigned int maskSeedLn = (1L << 2*wLen) - 1 ;
  const unsigned int mask26 = (1L << 26) - 1 ;
  BigArray aa = pp->intronSeeds ; 
  long int ii, iMax = bigArrayMax (aa) ;
  CW *restrict up = 0 ;
  int chrom = 0, a1, a2, da, v1, v2, dv, w1, w2, dw ;
  Array dna = 0 ;
  BOOL isIntronDown ;
  
  bigArraySort (aa, intronOrder) ;
  
  up = iMax ? bigArrp (aa, 0, CW) : 0 ;
  for (ii = 0 ; ii < iMax ; ii++, up++)
    {
      long unsigned int w, wr ;
      
      if (up->nam != chrom)
	{
	  chrom = up->nam ;
	  dna = array (dnas, chrom, Array) ;
	}

      a1 = up->pos ;           /* 1-based Gt_ag position */
      a2 = up->intron ;        /* 1-based gt_aG position */
      if (a1 > a2)
	{ int a0 = a1 ; a1 = a2 ; a2 = a0 ; isIntronDown = FALSE ; }
      else
	isIntronDown = TRUE ;
      da = a2 - a1 + 1 ;
      if (da >= (0x1 << 26))
	continue ; /* we use only 25 bits to code the intron length */

      /* check if left exon is shorter than 16 */
      v1 = v2 = 0 ; dv = 16 ;
      if (ii > 0)
	{
	  CW *restrict vp = up - 1 ;
	  if (vp->nam == up->nam)
	    {
	      if (isIntronDown && vp->pos < vp->intron && vp->intron > a1 - 17)
		{
		  v1 = vp->pos ; v2 = vp->intron ;
		  dv = a1 - v2 - 1 ; /* length of left exon */
		}
	      else if (! isIntronDown && vp->pos > vp->intron && vp->pos > a1 - 17)
		{
		  v2 = vp->pos ; v1 = vp->intron ;
		  dv = a1 - v2 - 1 ; /* length of left exon */
		}
	    }
	}
      
      /* check if right exon is shorter than 16 */
      w1 = w2 = 0 ; dw = 16 ;
      if (ii < iMax - 1)
	{
	  CW *restrict wp = up + 1 ;
	  if (wp->nam == up->nam)
	    {
	      if (isIntronDown && wp->pos < wp->intron && wp->pos < a2 + 17)
		{
		  w1 = wp->pos ; w2 = wp->intron ;
		  dw = w1 - a2 - 1 ; /* length of right exon */
		}
	      else if (! isIntronDown && wp->pos > wp->intron && wp->intron < a2 + 17)
		{
		  w2 = wp->pos ; w1 = wp->intron ;
		  dw = w1 - a2 - 1 ; /* length of right exon */
		}
	    }
	}
      
      /* Construct a 32 bases, 64 bits long unsigned integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base 
       * first 15 bases belong to the left exon, upstream of Gt_ag
       * It does not matter if the annotated left or rigth exon is shorter than 15 bases
       */
      int dw0 = dw, dv0 = dv ;
      for (int pass = 0 ; pass < 4 ; pass++)
	{
	  /* In case of short exon  use either the dv and dw, or one of them or none */
	  switch (pass)
	    {
	    case 0: /* consider the up intron in isolation, forget the v and w introns */
	      dv = dw = 16 ;
	      break ;
	    case 1: /* consider the vp intron if it exists */
	      dv = dv0 ; dw = 16 ;
	      if (dv == 16) continue ;
	      break ;
	    case 2: /* consider the wp intron if it exists */
	      dv = 16 ; dw = dw0 ;
	      if (dw == 16) continue ;
	      break ;
	    case 3: /* consider both the vp and the wp intron if they both exit */
	      dv = dv0 ; dw = dw0 ;
	      if (dv == 16 || dw == 16) continue ;
	      break ;
	    } 
	  p = 0 ;
	  w = wr = 0 ;
	  
	  if (dv < 16) /* import some bases from 2 exons above */
	    {
	      cp = arrp (dna, v1 - (16 - dv) - 1, unsigned char) ;
	      for (int i = 0 ; i < 16 - dv ; i++, cp++)
		{
		  w <<= 2 ; wr >>= 2 ;
		  p++ ;
		  switch ((int)*cp)
		    {  /* alphabetic order and XOR complement */
		    case 0: p = 0 ; break ;
		    case A_:             wr |= 0x3L << nShift ; break ;
		    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		    case T_: w |= 0x3L ;                   break ;
		      
		    default: p = 0 ;
		    }
		}
	    }
	  
	  /* import dv bases from the left exon */
	  cp = arrp (dna, a1 - dv - 1, unsigned char) ;
	  for (int i = 0 ; i < dv ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      switch ((int)*cp)
		{  /* alphabetic order and XOR complement */
		case 0: p = 0 ; break ;
		case A_:             wr |= 0x3L << nShift ; break ;
		case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		case T_: w |= 0x3L ;                   break ;
		  
		default: p = 0 ;
		}
	    }
	  
	  /* import dw bases from the right exon */
	  cp = arrp (dna, a2, unsigned char) ;
	  for (int i = 0 ; i < dw ; i++, cp++)
	    {
	      w <<= 2 ; wr >>= 2 ;
	      p++ ;
	      switch ((int)*cp)
		{  /* alphabetic order and XOR complement */
		case 0: p = 0 ; break ;
		case A_:             wr |= 0x3L << nShift ; break ;
		case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		case T_: w |= 0x3L ;                   break ;
		  
		default: p = 0 ;
		}
	    }
	  /* import 16 - dw  bases from the following exon */
	  if (dw < 16)
	    {
	      cp = arrp (dna, w2, unsigned char) ; /* a2 is G in bio coords, so cp is the first base of the second exon */
	      for (int i = 0 ; i < 16 - dw  ; i++, cp++)
		{
		  w <<= 2 ; wr >>= 2 ;
		  p++ ;
		  switch ((int)*cp)
		    {  /* alphabetic order and XOR complement */
		    case 0: p = 0 ; break ;
		    case A_:             wr |= 0x3L << nShift ; break ;
		    case C_: w |= 0x1L ; wr |= 0x2L << nShift ; break ;
		    case G_: w |= 0x2L ; wr |= 0x1L << nShift ; break ;
		    case T_: w |= 0x3L ;                   break ;
		      
		    default: p = 0 ;
		    }
		}
	    }
	  
	  if (p == 32)   /* 4,5,6,7,8,9,10,11 = 8 possibilites, codes on 3 bits 21-24 */
	    for (int i1 = 4 - (wLen - 16)  ; i1 < 13 ; i1++) /* search every word with at least 4 bases on each side */
	      {
		if (i1 < 0) continue ;
		long unsigned int x =   (w >> (2 * (32 - wLen - i1))) & maskSeedLn ;
		long unsigned int xr =  (wr >> 2 * (i1)) & maskSeedLn ;
		int dx1 = 16 - i1 ;     /* number of letters in first exon */
		BOOL  minus  = (x > xr) ; /* word is antistrand to the exon-exon junction */
		long unsigned int word = ( minus ? xr : x ) ;
		k = word & maskN ;   /* right bits are used to select the relevant table */
		word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
		word &= mask32 ;     /* masked down to 32 bits integer */
		cw = bigArrayp (bbG->cwsN[k], bigArrayMax(bbG->cwsN[k]), CW) ;
		cw->nam = (chrom << 1) | (minus ? 0x1 : 0) ;
		cw->seed = word ;
		cw->pos = a1 ; /* in introns, a1 cannot be zero */
		cw->intron =                  /* bit 30 alone, isExon */
		  (0x1 << 31)                 /* bit 31: isIntron */
		  | (isIntronDown ? (0x1 << 30) : 0) /* bit 30: intron strand */
		  | ((da & mask26)  << 4)            /* bits 5-30 : intron length */
		  | (dx1 & 0xf)          /* bits 0-4 : number of seed letters in the left exon, dx1 < 16 */
		  
		  ;
		nn++ ;
	      }
	}
    }
  fprintf (stderr, "+++++++ Coded %ld intron seeds\n", nn) ;
  if (debug && NN == 1)
    {
      BigArray uu = bigArrayCreate (nn, CW) ;
      CW *up, *vp = bigArrayp (bbG->cwsN[0], bigArrayMax (bbG->cwsN[0]) - nn, CW) ;
      for (int i = 0 ; i < nn ; i++, vp++)
	{
	  up = bigArrayp (uu, i, CW) ;
	  *up = *vp ;
	}
      showCws (pp, 0, uu) ;
    }
  return nn ;
} /* intronCodeWords */

/**************************************************************/
/* parse, code, sort the genome and create the index on disk
 * the human index takes around 18 GigaBytes
 */
static long int createTargetIndex (PP *pp, BB *bbG, Array tArray)
{
  AC_HANDLE h = ac_new_handle () ;
  BigArray cwsN[NN] ;
  int nMax = arrayMax (tArray) ;
  TC *tc = 0 ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  long int nn = 0 ;
  int nTc = 0 ;
  
  memset (bbG, 0, sizeof (BB)) ;
  bbG->isGenome = TRUE ;
  t1 = clock () ;
  printf ("+++ %s: Parse the target files\n", timeBufShowNow (tBuf)) ;

  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  
  /* parse all targets into a single bbG->seqs array, with targetClass prefix in the sequence name */
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;
      if (tc->targetClass == 'I')
	continue ; /* we need to parse thge genome before the introns */
      nTc++ ;
    }

  for (int nn = 0, ntc = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;
      RC rc ;
      int step ;
      
      if (tc->targetClass == 'I')
	continue ; /* we need to parse the genome before the introns */
      ntc++ ;
      
      memset (&rc, 0, sizeof (RC)) ;
      rc.fileName1 = tc->fileName ;
      rc.format = tc->format ;
      rc.run = nn + 1 ;
      sequenceParser (pp, 0, tc, bbG, ntc == nTc ? 2 :1 ) ; /* 2 for last non-intron target */
      step = (bbG->length < 1<<20) ? 1 : 3 ;
      if (pp->tStep)
	step = pp->tStep ;
      if (step > pp->tStep)
	pp->tStep = step ;
      if (0)
	{
	  int mem = 0, mx = 0 ; /* megaBytes */
	  messAllocStatus (&mem) ;
	  messAllocMaxStatus (&mx) ;
	  fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
	}
    }
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;

      if (tc->targetClass == '*')
	intronParser (pp, bbG, tc) ;
    }
  
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;

      if (tc->targetClass == 'I')
	gffParser (pp, bbG, tc) ;
    }
  
  t2 = clock () ;
  cpuStatRegister ("1.Parse targets" , pp->agent, bbG->cpuStats, t1, t2, arrayMax (bbG->dnas) - 1) ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "===== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }

  t1 = clock () ;
  printf ("%s : extract the target seeds\n" , timeBufShowNow (tBuf)) ;
  codeWordsDo (pp, bbG, pp->tStep, TRUE) ;
  if (pp->intronSeeds)
    intronCodeWords (pp, bbG) ;
  for (int k = 0 ; k < NN ; k++)
    {
      long int n1 = bigArrayMax (bbG->cwsN[k]) ;
      cwsN[k] = bbG->cwsN[k] ;
      bbG->cwsN[k] = 0 ;
      nn += n1 ; 
    }
  for (int k = 0 ; k < NN ; k++)
    {
      long int n1 = bigArrayMax (cwsN[k]) ;
      if (1)
	{
	  int mem = 0, mx = 0 ; /* megaBytes */
	  messAllocStatus (&mem) ;
	  messAllocMaxStatus (&mx) ;
	  fprintf (stderr, "=== k=%d , %ld/%ld words %.1f %%,  Allocated %d Mb, max %d Mb\n", k, n1, nn, 100.0*n1/nn,  mem, mx) ;
	}
  }
  t2 = clock () ;
  cpuStatRegister ("2.Extract target seeds" , pp->agent, bbG->cpuStats, t1, t2, bbG->nSeqs) ;
  t1 = clock () ;

  printf ("%s : sort the target seeds\n" , timeBufShowNow (tBuf)) ;
  for (int k = 0 ; k < NN ; k++)
    newMsort (cwsN[k], cwOrder) ;

  t2 = clock () ;
  cpuStatRegister ("3.Sort seeds" , pp->agent, bbG->cpuStats, t1, t2, nn) ;
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }
  t1 = clock () ;

  printf ("%s : write the index to disk\n" , timeBufShowNow (tBuf)) ;
  bbG->cwsN = halloc (NN * sizeof(BigArray), bbG->h) ;
  for (int kk = 0 ; kk < NN ; kk++)
    {
      bbG->cwsN[kk] = GenomeAddSkips (pp, cwsN[kk], bbG, kk) ;
      bigArrayDestroy (cwsN[kk]) ;
    }

  storeTargetIndex (pp, bbG, pp->tStep) ;
  
  t2 = clock () ;
  cpuStatRegister ("4.Write the index to disk" , pp->agent, bbG->cpuStats, t1, t2, nn) ;
  if (1)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }

  ac_free (h) ;
  return nn ;
} /* createTargetIndex */

/**************************************************************/
/**************************************************************/

static void sortWords (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  long int nnn = 0 ;
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->csChan, &bb, BB))
    {
      if (1)
	{
	  long int nn = 0 ;
	  if (pp->debug) printf ("+++ %s: Start sort words\n", timeBufShowNow (tBuf)) ;
	  
	  t1 = clock () ;
	  for (int k = 0 ; k < NN ; k++)
	    if (bb.cwsN[k])
	      {
		newMsort (bb.cwsN[k], cwOrder) ;
		nn += bigArrayMax (bb.cwsN[k]) ;
	      }
	  t2 = clock () ;
	
	  cpuStatRegister ("4.SortWords", pp->agent, bb.cpuStats, t1, t2, nn) ;
	  if (pp->debug) printf ("--- %s: Stop sort words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;
	}
      channelPut (pp->smChan, &bb, BB) ;

    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close smChan at %d,  sorted %ld words\n", n, nnn) ;
      channelCloseAt (pp->smChan, n) ;
    }

  return ;
} /* sortWords */

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
	  cpuStatRegister ("5.MatchHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
	}
      channelPut (pp->moChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close moChan at %d,  found %ld hits\n", n, nnn) ;
      channelCloseAt (pp->moChan, n) ;
    }

  return ;
} /* matchHits */

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

static int wiggleOrder (const void *va, const void *vb)
{
  const ALIGN *up = va ;
  const ALIGN *vp = vb ;
  int n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->w1 - vp->w1  ; if (n) return n ;

  return 0 ;
} /* wiggleOrder */

/**************************************************************/
/* a0 = a1 - x1 is the putative position of base 1 of the read 
 * It also works for the negative strand (a1 < 0, x1 > 0).
 */
static int hitReadOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n ;

  n = ((up->read > vp->read) - (up->read < vp->read)) ; if (n) return n ;
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n = ((up->a1 > vp->a1) - (up->a1 < vp->a1)) ;  if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; 
  return n ;
} /* hitReadOrder */

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
static int hitPairOrder (const void *va, const void *vb)
{
  const HIT *up = va ;
  const HIT *vp = vb ;
  int n, n1, n2 ;
  
  n = ((up->read >> 1) > (vp->read >> 1)) -  ((up->read >> 1) < (vp->read >> 1)) ; if (n) return n ; 
  n = ((up->chrom > vp->chrom) - (up->chrom < vp->chrom)) ; if (n) return n ; 
  n1 = up->a1 + (up->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n2 = vp->a1 + (vp->x1 >> NSHIFTEDTARGETREPEATBITS) ;
  n = n1 - n2 ; if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; 

  return n ;
} /* hitPairOrder */

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
	      newMsort (bb.hits, hitPairOrder) ;
	      t2 = clock () ;

	      long int nn = bigArrayMax (bb.hits) ;
	      nnn += nn ;
	      cpuStatRegister ("6.SortHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
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

      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close oaChan at %d,  sorted %ld hits\n", n, nnn) ;
      channelCloseAt (pp->oaChan, n) ;
    }
  return ;
} /* sortHits */

/**************************************************************/
#ifdef __SSE2__
#define EP128
#include <emmintrin.h> // SSE2
/*  aaaa aaaa aaaa aaaa  / aaaa aaaa aaaa aaga   -> 14 (number of exact matches)
static int first_non_equal_byte(unsigned char *cp, unsigned char *cq) {
    __m128i v1 = _mm_loadu_si128((__m128i*)cp); // aaaaaaaaaaaaaaaa
    __m128i v2 = _mm_loadu_si128((__m128i*)cq); // aaaaaaaaaaaaaaga
    __m128i cmp = _mm_cmpeq_epi8(v1, v2);       // 0xff for equal, 0x00 for non-equal
    int mask = _mm_movemask_epi8(cmp);          // 0xffff for equal, bit 0 for non-equal
    return mask == 0xffff ? 16 : __builtin_ctz(~mask); // First 0 bit
}
*/
#endif
/**************************************************************/
#include <stdint.h>

void extendExact(Array dna, int *x1p, int *x2p, Array dnaG, int *a1p, int *a2p)
{
  int dnaMax = arrayMax (dna) ;
  int dnaGMax = arrayMax (dnaG) ;
  unsigned char *cp, *cq ;
  BOOL ok = TRUE ;
  int a1 = *a1p, a2 = *a2p, x1 = *x1p, x2 = *x2p ;
  
  cp = arrp (dna , x2, unsigned char) ; /* first base beyond exact match */
  cq = arrp (dnaG, a2, unsigned char) ;

  /* Forward extension 16 bases steps */
#ifdef VECTORIZED_MEM_CPY
  while (x2 < dnaMax - 16 && a2 < dnaGMax - 16)
    {
      __m128i v1 = _mm_loadu_si128((__m128i*)cp);
      __m128i v2 = _mm_loadu_si128((__m128i*)cq);
      __m128i cmp = _mm_cmpeq_epi8(v1, v2);
      int mask = _mm_movemask_epi8(cmp);
      if (mask != 0xffff)
	{
	  int pos = __builtin_ctz(~mask);
	  x2 += pos; a2 += pos;
	  ok = FALSE ;
	  break; // Mismatch
        }
      else
	{
	  cp += 16; cq += 16; x2 += 16; a2 += 16;
	}
    }
#endif
  // Byte-wise to end
  if (ok)
    while (x2 < dnaMax && a2 < dnaGMax && *cp == *cq)
      { cp++; cq++; x2++ ; a2++ ;}

  // Backward extension
  cp = arrp (dna , x1 - 1, unsigned char) ; /* first matching base */
  cq = arrp (dnaG, a1 - 1, unsigned char) ;
  ok = TRUE ;
#ifdef VECTORIZED_MEM_CPY
  while (x1 > 16 && a1 > 16)
    {
      __m128i v1 = _mm_loadu_si128((__m128i*)(cp - 16));
      __m128i v2 = _mm_loadu_si128((__m128i*)(cq - 16));
      __m128i cmp = _mm_cmpeq_epi8(v1, v2);
      int mask = _mm_movemask_epi8(cmp);
      if (mask != 0xffff)
	break ;
      else
        { cp -= 16; cq -= 16; x1 -= 16; a1 -= 16; }
    }
#endif
  // Byte-wise to end
  while (x1 > 1 && a1 > 1 && *cp == *cq)
    { cp--; cq--; x1-- ; a1-- ; }
  *a1p = a1 ; *a2p = a2 ; *x1p = x1 ; *x2p = x2 ;    
} /* extendExact */

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
  
  if (0 && errMax == 0) /* faster code, actually slower oniRefSeq and gives very poor number of bases aligned */
    {
      if (isDown)
	extendExact(dna, x1p, x2p, dnaG, a1p, a2p) ;
      else
	{
	  *a1p = chromLength - *a1p + 1 ; *a2p = chromLength - *a2p + 1 ;
	  extendExact(dna, x1p, x2p, dnaGR, a1p, a2p) ;
	  *a1p = chromLength - *a1p + 1 ; *a2p = chromLength - *a2p + 1 ;
	}
      dx = *x2p - *x1p + 1 ;
      goto done ;
    }
  

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
  
 done:

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
/* sliding intron, clip errors in overlap r636 r452 */
static void alignOptimizeIntron (BB *bb, ALIGN *vp, ALIGN *wp, Array dnaG)
{
  A_ERR *epX, *epY ;
  int x2 = vp->x2 ;
  int y1 = wp->x1 ;
  int nEx = vp->nErr ;
  int nEy = wp->nErr ;
  int i, j, nE ;
  int bestN, bestI, bestJ ;
  int dy = vp->x2 - wp->x1 + 1 ;
  int da = vp->a2 - wp->a1 ;

  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
      memset (wp, 0, sizeof (ALIGN)) ;
      wp->chain = -1 ;
      return ;
    }
  if (vp->x1 > wp->x1)
    {
      int dx = vp->x1 - wp->x1 ;
      wp->x1 = vp->x1 ;
      if (wp->a1 < wp->a2) wp->a1 += dx ;
      else wp->a1 -= dx ;
    }
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && nEx + nEy > 0)
    {  
      for (nE = 0 ; nE < nEx ; nE++)
	{  /* count the vp errors that cannot be clipped */
	  epX = arrp (vp->errors, nE, A_ERR) ;
	  int zX = epX->iShort ; /* first bad base */
	  if (zX + 1 >= y1) /* bio coords */
	    break ;
	}
      /* if we clip both ali at y1, we keep nE errors
       * let us try to find a better position
       */
      i = nE ; j = 0 ; /* next error in vp and wp */
      nE = nE + nEy ;
      bestN = nE ; bestI = i ; bestJ = j ; 
      epX = epY = 0 ;
      while (bestN && (i < nEx || j < nEy))
	{
	  int zX = x2 ;
	  if (i < nEx)
	    {
	      epX = arrp (vp->errors, i, A_ERR) ;
	      zX = epX->iShort ; /* first bad base starting from the left */
	      zX = zX  < x2 - 1 ? zX : x2 - 1 ;
	    }
	  int zY = x2 ;
	  if (j < nEy)
	    {
	      epY = arrp (wp->errors, j, A_ERR) ;
	      zY = epY->iShort  ; /* last bad base starting from the right */
	      switch (epY->type)
 		{   /* we need to adjust the coordinates */
		case INSERTION:
		  break ;
		case INSERTION_DOUBLE:
		  zY+= 1 ;
		  break ;
		case INSERTION_TRIPLE:
		  zY += 2 ;
		  break ;
		case TROU:
		case TROU_DOUBLE:
		case TROU_TRIPLE:
		  zY -- ;
		  break ;
		default:
		  break ;
		}
	      zY = zY  < x2  ? zY : x2 ;
	    }
	  if (epX && epY && epX->iLong == epY->iLong && epX->iShort == epY->iShort)
	    {
	      /* merge the 2 alignments */
	      vp->a2 = wp->a2 ;
	      vp->x2 = wp->x2 ;
	      vp->nErr += wp->nErr ;
	      vp->errors = vp->errors ? vp->errors : wp->errors ;
	      if (vp->errors && vp->nErr)
		epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
	      memset (wp, 0, sizeof (ALIGN)) ;
	      wp->chain = -1 ;
	      return ;
	    }

	  if (zX < zY)
	    { nE++ ; i++ ; }
	  else if (zX >= zY && zY < x2)
	    {
	      nE-- ;
	      if (nE < bestN)
		{ bestN = nE ; bestI = i ; bestJ = j + 1 ; }
	      j++ ;
	    }
	  else
	    break ;
	}
      if (bestJ > 0) /* if bestJ == 0, wp->x1/a1 is well positioned */
	{ /* we must clip these errors */
	  epY = arrp (wp->errors, bestJ - 1, A_ERR) ;
	  int zA = epY->iLong + ((wp->chrom & 0x1) ? -1 : 1) ;
	  int zY = epY->iShort + 1 ; /* we need to find the first exact */
	  switch (epY->type)
	    {   
	    case INSERTION:
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_DOUBLE:
	      zY += 1 ;
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case INSERTION_TRIPLE:
	      zY += 2 ;
	      zA += ((wp->chrom & 0x1) ? -1 : 1) ;
	      break ;
	    case TROU:
	      zY -- ;
	      break ;
	    case TROU_DOUBLE:
	      zY -- ;
	      zA += ((wp->chrom & 0x1) ? -1: 1) ;
	      break ;
	    case TROU_TRIPLE:
	      zY -- ;
	      zA += ((wp->chrom & 0x1) ? -2 : 2) ;
	      break ;
	    default:
	      break ;
	    }
	  wp->x1 = zY + 1 ;  /* bio coords */
	  wp->a1 = (wp->chrom & 0x1 ?  arrayMax(dnaG) - zA : zA + 1) ;
	  wp->nErr -= bestJ ;
	  for (j = 0, epY = arrp (wp->errors, 0, A_ERR) ; j < wp->nErr ; epY++, j++)
	    *epY = *(epY + bestJ) ;
	  arrayMax (wp->errors) = wp->nErr ;
	}
      if (bestI < nEx) /* we can clip the vp errors */
	{
	  epX = arrp (vp->errors, bestI, A_ERR) ;
	  int zA = epX->iLong  + ((wp->chrom & 0x1) ? +1 : -1) ;
	  int zX = epX->iShort - 1 ; /* we need to find the last exact */
	  vp->nErr -= nEx - bestI ;
	  vp->x2 = zX + 1 ;  /* bio coords */
	  vp->a2 = (vp->chrom & 0x1 ?  arrayMax(dnaG) - zA : zA + 1) ;
	  arrayMax (vp->errors) = vp->nErr ;
	}
    }

  
  BOOL isIntronDown = TRUE ;
  BOOL foundDonor = FALSE ;
  BOOL foundAcceptor = FALSE ;
  int donor = 0, acceptor = 0 ;
  BOOL isReadDown = TRUE ;
  dy = vp->x2 - wp->x1 + 1 ;
  
  /* trim if possible the vp->x2 end of the first exon on a known 'donor' */
  if (vp->a1 < vp->a2)
    {
      isReadDown = TRUE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;
      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; isIntronDown = FALSE ; }
      else if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ;  isIntronDown = FALSE ; }
      if (donor > acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 >= donor && vp->a2 <= donor + dy - 1) 
	{  /* move back to the canonical donor site */
	  dy = vp->a2 - donor + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;
	}
      if (donor && vp->a2 == donor - 1)
	foundDonor = TRUE ;
    }
  else
    {
      isReadDown = FALSE ;
      donor = vp->donor ;
      acceptor = wp->acceptor ;

      if (donor < 0)
	{ donor = - donor ; acceptor = -acceptor ; isIntronDown = FALSE ; }
      if (donor == 0 && acceptor < 0)
	{ acceptor = - acceptor ;  isIntronDown = FALSE ; }
      if (donor && donor < acceptor)
	donor = acceptor = 0 ;
      
      if (dy > 0 && donor && vp->a2 <= donor && vp->a2 >= donor - dy + 1)
	
	{  /* move back to the canonical donor site */
	  dy = donor - vp->a2 + 1 ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	}      
      if (donor && vp->a2 == donor + 1)
	foundDonor = TRUE ;
    }

  /* alternativelly  trim if possible the wp->x2 start of the second  exon on a known 'acceptor' */
  dy = vp->x2 - wp->x1 + 1 ;
  if (isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 < acceptor + 1 && wp->a1 + dy >= acceptor + 1)
	{
	  dy = acceptor + 1 - wp->a1 ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	}
      if (wp->a1 == acceptor + 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 -= dy ;
	    }
	}
    }
  else if (!isReadDown && acceptor)
    {
      if (dy > 0 && wp->a1 > acceptor - 1 && wp->a1 - dy <= acceptor - 1)
	{
	  dy = wp->a1 - acceptor + 1 ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;
	}
      if (wp->a1 == acceptor - 1)
	{
	  foundAcceptor = TRUE ;
	  dy = vp->x2 - wp->x1 + 1 ;
	  if (dy > 0) /* trim vp->x2 */
	    {
	      vp->x2 -= dy ;
	      vp->a2 += dy ;
	    }
	}
    }

  BOOL gt_ag = FALSE ;
  BOOL gc_ag = FALSE ;
  BOOL ct_ac = FALSE ;
  BOOL ct_gc = FALSE ;

  /* no annotated junction, look for gt_ag */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0 && isReadDown && ! foundDonor && !foundAcceptor  && wp->a1 > 3 && vp->a2 > dy)  
    {
      unsigned char *cp = arrp (dnaG, vp->a2 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, wp->a1 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;

      if (bb->nIntronSupportPlus >= bb->nIntronSupportMinus)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { gt_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { ct_ac = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { gc_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { ct_ac = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { gt_ag = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { ct_gc = TRUE ; bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok1 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok1 ; }
	}
    ok1:
      if (bestI >= 0)
	{
	  bestI = (bestI ? bestI : (dy+1) / 2) ; /* avod zero to keep vp oriented */
	  dy -= bestI ;
	  vp->x2 -= dy ;
	  vp->a2 -= dy ;

	  dy = bestI ;
	  wp->x1 += dy ;
	  wp->a1 += dy ;
	}	    
    }
  else if (dy > 0 && !isReadDown && ! foundDonor && !foundAcceptor && wp->a1 > dy && vp->a2 > 3) 
    {
      /* move backwards on the genome */
      unsigned char *cp = arrp (dnaG, wp->a1 - dy, unsigned char) ; /* the base just after vp->a2 - dy */
      unsigned char *cq = arrp (dnaG, vp->a2 - 3, unsigned char) ; /* the base 2 bases before wp->a1 - dy */
      int bestI = -1 ;
      if (bb->nIntronSupportPlus < bb->nIntronSupportMinus)
	{     /* favor gt_ag over ct_ac */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	}
      else
	{ /* favor ct_ac over gt_ag */
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == T_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_ && cq[i] == G_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == G_ && cp[i+1] == C_ && cq[i] == A_ && cq[i+1] == G_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cp[i] == C_ && cp[i+1] == T_)
	      { bestI = i ; goto ok2 ; }
	  for (int i = 0 ; i <= dy ; i++)
	    if (cq[i] == A_ && cq[i+1] == C_)
	      { bestI = i ; goto ok2 ; }
	}
    ok2:
      if (bestI >= 0)
	{
	  dy -= bestI ;
	  wp->x1 += dy ;
	  wp->a1 -= dy ;

	  dy = bestI ;
	  vp->x2 -= dy ;
	  vp->a2 += dy ;
	}	    
    }
  
  /* In any case trim the second exon a bouts francs */
  dy = vp->x2 - wp->x1 + 1 ;
  if (dy > 0)
    { 
      wp->x1 += dy ;
      wp->a1 += (wp->a1 < wp->a2 ? dy : -dy) ;
    }

  dy = vp->x2 - wp->x1 + 1 ;
  da = vp->a2 - wp->a1 ;

  if (da < 5 && da > -5 && dy < 5 && dy > -5 && vp->chrom == wp->chrom)
    {
      /* merge the 2 alignments */
      vp->a2 = wp->a2 ;
      vp->x2 = wp->x2 ;
      vp->nErr += wp->nErr ;
      vp->errors = vp->errors ? vp->errors : wp->errors ;
      if (vp->errors && vp->nErr)
	epX = arrayp (vp->errors, vp->nErr - 1, A_ERR) ;
      memset (wp, 0, sizeof (ALIGN)) ;
      wp->chain = -1 ;
      return ;
    }

  /* register exact intron support */
  if (isReadDown && donor == vp->a2 + 1 && acceptor == wp->a1 - 1)
    {  /* this intron is confirmed */
      INTRON *zp = arrayp (bb->confirmedIntrons, arrayMax (bb->confirmedIntrons), INTRON) ;
      zp->chrom = vp->chrom ;
      zp->a1 = isIntronDown ? donor : acceptor ;
      zp->a2 = isIntronDown ? acceptor : donor ;
      bb->nIntronSupportPlus += (isIntronDown ? 1 : 0) ;
      bb->nIntronSupportMinus += (isIntronDown ? 0 : 1) ;
    }
  else if (! isReadDown && donor == vp->a2 - 1 && acceptor == wp->a1 + 1)
    {  /* this intron is confirmed */
      INTRON *zp = arrayp (bb->confirmedIntrons, arrayMax (bb->confirmedIntrons), INTRON) ;
      zp->chrom = vp->chrom ;
      zp->a1 = isIntronDown ? acceptor : donor ;
      zp->a1 = isIntronDown ? donor : acceptor ;
      bb->nIntronSupportPlus += (isIntronDown ? 0 : 1) ;
      bb->nIntronSupportMinus += (isIntronDown ? 1 : 0) ;
    }
  gt_ag = gc_ag ^ gt_ag ^ ct_ac ^ ct_gc ; /* for compiler happiness */

} /* alignOptimizeIntron */

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
	      alignOptimizeIntron (bb, up, wp, dnaG) ;
	    }
	  else
	    alignOptimizeIntron (bb, vp, wp, dnaG) ;
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
  int step = WIGGLE_STEP ;  /* examples s=10, 5, 1 */
  int demiStep = step/2 ;
  char allTc[256] ;
  Array dna = 0, dna1 = 0, dna2 = 0 ;
  int chromA = 0 ;
  Array dnaG = 0, dnaGR = 0 ;
  int read1 = read & (~0x1) ;
  int read2 = read | 0x1 ;

  dna1 = arr (bb->dnas, read1, Array) ;
  dna2 = arr (bb->dnas, read2, Array) ;
  
  if (2*demiStep == step) demiStep-- ; /* examples d=4, 2, 0 */
  /* create chains */

  /* register the wiggle boundaries */
  iMax = arrayMax (aa) ;
  if (iMax)
    for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
      if (read == up->read)
	{
	  int a1 = up->a1 ;
	  int a2 = up->a2 ;
	  if (up->read & 0x1) { int a0 = a1 ; a1 = a2 ; a2 = a0 ;}
	  up->w1 = (a1 + demiStep)/step ;
	  up->w2 = (a2 + demiStep)/step ;
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
		  newMsort (hits2, hitReadOrder) ;
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

#ifdef JUNK

dna2dna -i Fasta/iRefSeq38/iRefSeq38.fasta.gz -get 'NM_001385438.1|Gene|NBPF15|GeneId|284565' -o AMELY/NM_001385438.1 -rightClipAt 462
run -x ../Aligners/013_SortAlignG3R1/IDX.GRCh38.18 -i NM_001385438.1.fasta --numactl --nB 1 --step 3
zcat Fasta/iRefSeq38/TARGET.human.GRCh38.2025/hs.RefSeq.GRCh38.mrnaRemap.gz | grep NM_001385438.1

(gdb) p showCountChroms (countChroms)
.. 0:  weight=1356.00	seeds=1356	1:120/2:96/4:344/8:540/16:206/32:50	index 13401 15214	chrom=31	pos=82349478 82348881	da=4294966699	x= 1 4435
.. 1:  weight=1203.00	seeds=1203	1:62/2:192/4:404/8:352/16:146/32:47	index 3171 6239	chrom=30	pos=82429982 82529588	da=99606	x= 1 4433
.. 2:  weight=1187.00	seeds=1187	1:78/2:99/4:337/8:479/16:173/32:21	index 15215 18354	chrom=31	pos=84437849 84516001	da=78152	x= 1 4435
.. 3:  weight=1182.00	seeds=1182	1:120/2:234/4:344/8:288/16:140/32:56	index 6240 9456	chrom=30	pos=84235736 84166271	da=4294897831	x= 1 4271
.. 4:  weight=983.00	seeds=983	1:163/2:225/4:240/8:200/16:91/32:64	index 9841 11435	chrom=30	pos=85244038 85249952	da=5914	x= 532 4434
.. 5:  weight=633.00	seeds=633	1:4/2:62/4:185/8:208/16:119/32:55	index 19797 20510	chrom=49	pos=25502067 25694944	da=192877	x= 7 4435
.. 6:  weight=609.00	seeds=609	1:0/2:61/4:189/8:201/16:101/32:57	index 18740 19406	chrom=48	pos=24168143 24177212	da=9069	x= 56 4433
.. 7:  weight=256.00	seeds=256	1:0/2:10/4:27/8:103/16:115/32:1	index 2559 3148	chrom=30	pos=75258593 75790118	da=531525	x= 94 4402
.. 8:  weight=215.00	seeds=215	1:4/2:0/4:15/8:76/16:101/32:19	index 12072 12513	chrom=31	pos=28358224 30409394	da=2051170	x= 44 4404
.. 9:  weight=210.00	seeds=210	1:0/2:0/4:24/8:90/16:90/32:6	index 13000 13209	chrom=31	pos=74082705 74076496	da=4294961087	x= 94 4404
.. 10:  weight=209.00	seeds=209	1:0/2:0/4:22/8:108/16:75/32:4	index 2338 2552	chrom=30	pos=72605180 72669569	da=64389	x= 96 4402
#endif
/**************************************************************/

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
	  
	  cpuStatRegister ("7.Align_r", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  if (pp->debug) printf ("--- %s: Stop align %lu ali, %lu mismatches\n", timeBufShowNow (tBuf), bb.nAli, bb.nerr) ;
	}
      channelPut (pp->aeChan, &bb, BB) ;
    }
  
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close aeChan at %d,  found %ld ali\n", n, nnn) ;
      channelCloseAt (pp->aeChan, n) ;
    }
  return ;
} /* align */

/**************************************************************/
/**************************************************************/

static void wiggleDo (const PP *pp, BB *bb)
{
  ALIGN *ap ;
  long int ii, iMax = bigArrayMax (bb->aligns) ;
  int chrom = 0 ;
  Array wig = 0 ;
  Array wiggles = pp->bbG.wiggles ;
    
  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int w1 = ap->w1, w2 = ap->w2 ;
      int mult = ap->nTargetRepeats ;
      int weight = 720/mult ;

      if (! weight)
	continue ;
      if (ap->chrom != chrom)
	{
	  chrom = ap->chrom ;
	  wig = array (wiggles, chrom + 1, Array) ;
	  if (strchr (dictName (pp->bbG.dict, chrom >> 1), '|'))
	    wig = 0 ;
	}
      if (!wig)
	continue ;
      if (w1 < w2)
	for (int i = w1 ; i <= w2 ; i++)
	  array (wig, i, int) += weight ;
    }
  return ;
} /* wiggle do */

/**************************************************************/

static void wiggleExport (const PP *pp)
{
  Array wiggles = pp->bbG.wiggles ;
  int ii, iMax = arrayMax (wiggles) ;
  int step = WIGGLE_STEP ;
  DICT *dict = pp->bbG.dict ;

  for (ii = 1 ; ii < iMax ; ii ++)
    {
      Array wig = array (wiggles, ii, Array) ;

      if (wig && arrayMax (wig))
	{
	  AC_HANDLE h = ac_new_handle () ;
	  const char *chrom = dictName (dict, ii >> 1) ;
	  char *fNam = hprintf (h, ".any.%s.u.%s.BF", chrom, (ii & 0x1 ? "r" : "f")) ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
	  aceOutDate (ao, "##", "wiggle") ;
	  aceOutf (ao, "track type=wiggle_0\n") ;
	  aceOutf (ao, "fixedStep chrom=%s start=%d step=%d\n", chrom, step, step) ;
	  
	  for (int jj = 1 ; jj < arrayMax (wig) ; jj++)
	    aceOutf (ao, "%d\n", array (wig, jj, int)) ;
	  ac_free (h) ;
	}
    }
  return ;
} /* wiggleExport */

/**************************************************************/

static void wiggle (const void *vp)
{
  const PP *pp = vp ;
  char tBuf[25] ;
  int nn ;
  BB bb ;
  
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->wwChan, &nn, int))
    {
      if (pp->wiggle && bb.aligns)
	{
	  if (pp->debug) printf ("--- %s: Start wiggle %lu seeds\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	  t1 = clock () ;

	  int k = 0 ;
	  /* 	  int k = wiggleDo (pp, &nn) ; */

	  t2 = clock () ;
	  cpuStatRegister ("8.Wiggle", pp->agent, bb.cpuStats, t1, t2, k) ;
	  if (pp->debug) printf ("--- %s: Stop wiggle\n", timeBufShowNow (tBuf)) ;
	}
      channelPut (pp->doneChan, &bb, BB) ;
    }

  if (pp->wiggle)
    {
      printf ("--- %s: Start exporting %d wiggle files\n"
	      , timeBufShowNow (tBuf)
	      , arrayMax (pp->bbG.wiggles)
	      ) ;
      wiggleExport (pp) ;
      printf ("--- %s: Stop exporting wiggles\n", timeBufShowNow (tBuf)) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      channelCloseAt (pp->doneChan, n) ;
    }
  return ;
} /* wiggle */

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
/*********************************************************************/

/* export  SAM
 * hits is  in pExportSamOrder
 * so organized by  read/target/chain in the orientation of the target 
 * as befits a good Cuban CIGAR
 */

static int exportOneSamExon (BB *bb, BOOL isDown, vTXT cigar, ALIGN *ap, int *nMp, int *nSubp, int *nInsp, int *nDelp)
{
  A_ERR *ep ;
  Array errors = ap->errors ;
  int lSam = 0, dx, ii, iMax = errors ? arrayMax (errors) : 0 ;
  int x1 = ap->x1 ;
  const int L = 16 ;
  char segs[2*iMax + 4][L] ;
  int iSeg = 0 ;
  KEYSET ks = keySetCreate () ;
  
  if (iMax)
    for (ii = 0 ; ii < iMax ; ii++)
      {
	int j = isDown ? ii : iMax - 1 - ii ;
	ep = arrp (errors, j, A_ERR)  ;
	
	int xShort = ep->iShort + 1 ;
	int xLong = ep->iLong + 1 ;
	if (xShort <= 0 || xLong <= 0)
	  continue ;
	
	dx = ep->iShort - x1 + 1 ; /* number of exact matches, not including the problem */
	if (dx <= 0)
	  continue ;
	switch (ep->type)
	  {
	  case AMBIGUE:
	  case TYPE80:
	    break ;
	  case ERREUR:
	    *nSubp += 1 ;
	    break ;
	  case INSERTION:
	    *nInsp  += 1 ;
	    x1 = ep->iShort + 2 ;
	    if (x1 >= ap->x2) dx-- ;
	    keySet (ks, iSeg) = dx + 1 ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "1I") ;
	    *nMp += dx ;
	    break ;
	  case INSERTION_DOUBLE:     
	    *nInsp += 2 ;
	    x1 = ep->iShort + 3 ;
	    if (x1 >= ap->x2) dx-- ;
	    keySet (ks, iSeg) = dx + 2 ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "2I") ;
	    *nMp += dx ;
	    break ;
	  case INSERTION_TRIPLE:     
	    *nInsp += 3 ;
	    x1 = ep->iShort + 4 ;
	    if (x1 >= ap->x2) dx-- ;
	    keySet (ks, iSeg) = dx + 3 ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "3I") ;
	    *nMp += dx ;
	    break ;
	  case  TROU:
	    *nDelp += 1 ;
	    x1 = ep->iShort + 1 ;
	    keySet (ks, iSeg) = dx ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "1D") ;
	    *nMp += dx ;
	    break ;
	  case TROU_DOUBLE:
	    *nDelp += 2 ;
	    x1 = ep->iShort + 1 ;
	    keySet (ks, iSeg) = dx ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "2D") ;
	    *nMp += dx ;
	    break ;
	  case TROU_TRIPLE:
	    *nDelp += 3 ;
	    x1 = ep->iShort + 1 ;
	    keySet (ks, iSeg) = dx ;
	    snprintf(segs[iSeg++], L, "%dM", dx) ;
	    snprintf(segs[iSeg++], L, "3D") ;
	    *nMp += dx ;
	    break ;
	  }
      }
  
  dx = ap->x2 - x1 + 1 ;
  if (dx > 0)
    {
      keySet (ks, iSeg) = dx ;
      snprintf(segs[iSeg++], L, "%dM", dx) ;
      *nMp += dx ;
    }

  dx = ap->x2 - ap->x1 + 1 ;
  for (int i = 0 ; i < iSeg ; i++)
    lSam += keySet (ks, i) ;

  for (int i = 0 ; i < iSeg ; i++)    
    vtxtPrintf (cigar, segs[isDown ? i : iSeg - i - 1]) ;

  if (1 && lSam != dx)
    fprintf (stderr,  "Bad length in sam cigar lSam=%d dx=%d ap->x1=%d ap->x2=%d\t%s\t%s\n", lSam, dx, ap->x1, ap->x2, vtxtPtr (cigar), dictName (bb->dict, ap->read >> 1)) ;

  keySetDestroy (ks) ;
  return lSam - dx ;
} /*  exportOneSamExon */

/*********************************************************************/

static char *exportOneSamCigar (BB *bb, vTXT cigar, ALIGN *ap0, int iMax, Array dna, int *nMp, int *nSubp, int *nInsp, int *nDelp, int *nGapp)
{
  int ii ;
  ALIGN *ap = ap0 ;
  int x1, x2 ;
  int dnaMax = dna ? arrayMax (dna) : 0 ;
  BOOL isDown = ap->a1 < ap->a2 ? TRUE : FALSE ;
  
  vtxtClear (cigar) ;
  if (isDown)
    {
      /* gap en tete */
      x1 = ap->x1 ;
      if (x1 > 1)
	vtxtPrintf (cigar, "%dS", x1 - 1) ;

      for (ii = 0 ; ii < iMax - 1 ; ii++, ap++)
	{
	  exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
	  int da = ap[1].a1 - ap[0].a2 - 1 ;
	  if (da > 0) vtxtPrintf (cigar, "%dN", da) ;
	  *nGapp += da ;
	  int dx = ap[1].x1 - ap[0].x2 - 1 ;
	  if (dx > 0) vtxtPrintf (cigar, "%dS", dx) ;
	}
      exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
      
      /* gap en queue */
      x2 = ap->x2 ;
      if (x2 < dnaMax)
	vtxtPrintf (cigar, "%dS", dnaMax - x2) ; 
    }
  else
    {
      ap += iMax - 1 ;

      /* gap en tete */
      x2 = ap->x2 ;
      if (x2 < dnaMax)
	vtxtPrintf (cigar, "%dS", dnaMax - x2) ; 
      
      for (ii = 0 ; ii < iMax - 1 ; ii++, ap--)
	{
	  exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
	  int da = ap[-1].a2 - ap[0].a1 - 1 ;
	  vtxtPrintf (cigar, "%dN", da) ;
	  *nGapp += da ;
	  int dx = ap[0].x1 - ap[-1].x2 - 1 ;
	  if (dx > 0) vtxtPrintf (cigar, "%dS", dx) ;
	}
      exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
      
      /* gap en queue */
      x1 = ap->x1 ;
      if (x1 > 1)
	vtxtPrintf (cigar, "%dS", x1 - 1) ;
    }
  return vtxtPtr (cigar) ;
} /* exportOneSamCigar */

/*********************************************************************/

static int exportOneSam (ACEOUT ao, const PP *pp, BB *bb, vTXT cigar, Array cigarettes, Array errors, ALIGN *ap0, long int i0, long int aMax)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *run = dictName (pp->runDict, bb->run) ;
  DICT *dictG = pp->bbG.dict ;
  DICT *dict = bb->dict ;
  ALIGN *ap = ap0 ;
  long int ii ;
  int kMax = 0 ;
  int flag = 0 ;
  int read = ap->read ;
  int chain = ap->chain ;
  int chrom = ap->chrom >> 1 ;
  int mate = 0 ;
  BOOL isDown = (ap->a1 <= ap->a2 ? TRUE : FALSE)  ;
  int chainA1 = ap->a1 ;
  BOOL pairedEnd = bb->rc.pairedEnd ;
  BOOL mateIsDown = TRUE ;
  BOOL isPrimary = ap->chain == 1 ? TRUE : FALSE ;
  BOOL isSecondary = FALSE ;
  BOOL isSupplementary = FALSE ;
  BOOL isOrphan = FALSE ; 
  Array dna = arr (bb->dnas, read, Array) ;
  int nM = 0, nSub = 0, nIns = 0, nDel = 0, nGap = 0 ;

  /* a chain [i0, iMax[, is reported as a single CIGAR */
  chainA1 = ap->chainA1 ;
  for (ii = i0, ap = ap0 ; ii < aMax && ap->chain == chain && ap->read == read ; ii++, ap++)
    kMax++ ;

  ap = ap0 ;
  if (ap->chain == 1)
    isPrimary = TRUE ;
  else
    {
      isPrimary = FALSE ;
      if (0) /* ATTENTION: meaning not clear */
	isSupplementary = TRUE ;
      else
	isSecondary = TRUE ;
    }

  /* if PAIRS locate mate */
  if (isPrimary && pairedEnd && ! mate)  
    isOrphan = TRUE ;
  
  
  /* there are 11 mandatory columns */
  aceOutf (ao, "%s/%s%s", run, dictName (dict, read >> 1), (pairedEnd ? ((read & 0x1 )? "<" : ">") : "")) ; 
  flag = 0 ;
        /* SAM flag
       * 0x4  0x8 unaligned 
       * 0x10 0x20  target strand
       * 0x80   second read of a pair
       * 0x100  secondary mappings
       * 1   	000000000001 	0x1   template having multiple templates in sequencing (read is paired)
       * 2 	000000000010 	0x2   each segment properly aligned according to the aligner (read mapped in proper pair)
       * 4 	000000000100 	0x4   this segment is unmapped (read1 unmapped)
       * 8 	000000001000 	0x8   next segment in the template unmapped (read2 unmapped)
       * 16 	000000010000 	0x10  this read is minus strand: SEQ being reverse complemented (read1 reverse complemented)
       * 32 	000000100000 	0x20  minus strand read2 : SEQ of the next segment in the template being reverse complemented (read2 reverse complemented)
       * 64 	000001000000 	0x40  the first segment in the template (is read1)
       * 128 	000010000000 	0x80  the last segment in the template (is read2)
       * 256 	000100000000 	0x100 not primary alignment
       * 512 	001000000000 	0x200 alignment fails quality checks
       * 1024 	010000000000 	0x400 PCR or optical duplicate
       * 2048 	100000000000 	0x800 supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region)
       */
if (1) 
    {   /* SAM FLAGS */
      /* this read belongs to a pair */
      if (pairedEnd) flag |= 0x1 ;
      /* all segments are well and compatibly aligned */
      if (pairedEnd && mate) flag |= 0x2 ;
      /* unmapped read */
      if (! ap->chrom) flag |= 0x4 ;  
      /* orphan pair-mate is unmapped */
      if (pairedEnd && isOrphan) flag |= 0x8 ;
      /* the seq is given in the orientation of the genome
       * the flags say if the shown read is the complement of 
       * read found on the machine fastq file 16 == 0x10,  32 == 0x20
       */
      if (! isDown) flag |= 0x10 ;
      if (mate && ! mateIsDown) flag |= 0x20 ;
      /* i am read 1 of a pair or first or central part of a composite */
      if (pairedEnd && ((read & 0x1) == 0)) flag |= 0x40 ;
      /* i am read 2 of a pair, or last or central part of a composite 128 == 0x80*/
      /* notice that central read parts are 0x4 | 0x8 == 0xb0 */
      if (pairedEnd && ((read & 0x1) == 1)) flag |= 0x80 ;
      /* multi-mapping : all have this flag except first one 256 == 0x100 */
      if (isSecondary) flag |= 0x100 ;
      /* read not passing quality filters 512 = 0x200 */
      if (ap->chainScore < 10) flag |= 0x200 ;
      /* i am a PCR duplicate: never used by magic 1024 == 0x400 */
      if (0) flag |= 0x400 ;
      /*  supplementary alignment  2048 == 0x800
       *  divers portion du meme read sont incompatibles
       *  all have this flag except the first one
       */
      if (isSupplementary) flag |= 0x800 ;
    }


  aceOutf (ao, "\t%d", flag) ;
  aceOutf (ao, "\t%s", dictName (dictG, chrom)) ;
  aceOutf (ao, "\t%d", chainA1) ;

  ap = ap0 ;
  if (1)
    {
      /* mapping quality, copied from magicBlast */
      int mapq = 255 ; /* 255 means MAPQ value unavailable */
      /* for single alignements, report 60 (like HISAT2) */
      
      if (ap->nTargetRepeats == 1)
	mapq = 60;
      else
	{
	  /* > 1 alignment (like TopHat2 and STAR) */
	  mapq = (int)((-10.0 * log10(1.0 - 1.0 / (double) ap->nTargetRepeats)) + 0.5);
	}
      aceOutf (ao, "\t%d", mapq) ;
    }
    
  /* CIGAR */
  if (exportOneSamCigar (bb, cigar, ap0, kMax, dna, &nM, &nSub, &nIns, &nDel, &nGap))
    aceOutf (ao, "\t%s", vtxtPtr (cigar)) ;
  else
    aceOutf (ao, "\t*") ;

  if (ap->mateChrom && ap->mateA1)
    {
      aceOutf (ao, "\t%s\t%d\t%d", dictName (dictG, ap->mateChrom >> 1), ap->mateA1, ap->pairLength) ;
    }
  else
    aceOutf (ao, "\t*\t0\t0") ;    

  /* export the sequence of the read */
  if (pp->exportSamSequence)
    {
      /* what should we do about clipping the polyA and adaptors ? */
      if (isDown)
	{
	  int i = arrayMax (dna) ;
	  Array dnaR = dnaCopy (dna) ;
	  char *cp = arrp (dnaR, 0, char) ;
	  dnaDecodeArray (dnaR) ;
	  if (0) while (i--)
	    { *cp = ace_upper (*cp) ; cp++; }
	  aceOutf (ao, "\t%s", arrp (dnaR, 0, char)) ;
	  ac_free (dnaR) ;
	}
      else
	{
	  int i = arrayMax (dna) ;
	  Array dnaR = dnaCopy (dna) ;
	  char *cp = arrp (dnaR, 0, char) ;
	  reverseComplement (dnaR) ;
	  dnaDecodeArray (dnaR) ;
	  if (0) while (i--)
	    { *cp = ace_upper (*cp) ; cp++; }
	  aceOutf (ao, "\t%s", arrp (dnaR, 0, char)) ;
	  ac_free (dnaR) ;
	}
    }
  else
    aceOut (ao, "\t*") ;

  if (pp->exportSamQuality && bb->quals)
    {
      Array qual = read < arrayMax (bb->quals) ? arr (bb->quals, read, Array) : 0 ;
      if (qual && arrayExists (qual) && arrayMax(qual))
	{
	  if (isDown)
	    {
	      aceOutf (ao, "\t%s", arrp (qual, 0, char)) ;
	    }
	  else
	    {
	      int i = arrayMax (qual) - 1 ;
	      char *cp = arrp (qual, i, char) ;
	      aceOut (ao, "\t") ;
	      while (i--)
		aceOutf (ao, "%c", *cp--) ;
	    }
	}
      else
	aceOut (ao, "\t*") ;
    }
  else
    aceOut (ao, "\t*") ;


#ifdef JUNK
  if (! isPrimary)
    aceOutf (ao, "\t*1") ;
  else
    {
      if (isDown)
	aceOutf (ao, "\t*2") ;  /* QUAL ascii de phread scale quaity + 33 ou * */
      else
	aceOutf (ao, "\t*3") ; /* we must reverse the order of the QUAL */
    }
#endif
  /* end of the 11 mandatory columns */

    
  aceOutf (ao, "\tNH:i:%d", ap->nChains) ; /* px->uu */ ;
  aceOutf (ao, "\tAS:i:%d", ap->pairScore ? ap->pairScore : ap->chainScore) ;
  aceOutf (ao, "\tNM:i:%d", nSub + nIns + nDel) ;
  aceOutf (ao, "\tmt:i:%d", 1) ; /* mult */
  /* aceOutf (ao, "\tMD:Z:0") ; */
  
  aceOut (ao, "\n") ;
  ac_free (h) ;
  return kMax ;
} /* exportOneSam */

/**************************************************************/
/* return the number of consumed ap (must be >=1) */
static int exportSamDo (ACEOUT ao, const PP *pp, BB *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ia, aMax = bigArrayMax (bb->aligns) ;
  vTXT cigar = vtxtHandleCreate (h) ;
  Array cigarettes = arrayHandleCreate (1024, SAMCIGAR, h) ;
  Array errors = arrayHandleCreate (1024, A_ERR, h) ;
  int read, chain ;

  for (ia = 0, read = chain = 0 ; ia < aMax ; ia++, ap++)
    {
      int k ;
      if (ap->read != read || ap->chain != chain)
	k = exportOneSam (ao, pp, bb, cigar, cigarettes, errors, ap, ia, aMax) ;
      read = ap->read ;
      chain = ap->chain ;
      ap += k - 1 ;
      ia += k - 1 ;
    }
  ac_free (h) ;
  return 1 ;
} /* exportSamDo */

/**************************************************************/

static void export (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  clock_t  t1, t2, dt = 50 * CLOCKS_PER_SEC ;
  
  char tBuf[25] ;
  int n, nn = 0 ;
  float nG = 0 ;

  if (pp->agent != 0) /* sequential agent */
    return ;
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
		{
		  char *VERSION = "0.1.1" ;
		  DICT *dictG = pp->bbG.dict ;
		  ao = aos[bb.run] = aceOutCreate (pp->outFileName, hprintf (h, ".%s.sam", dictName (pp->runDict, bb.run)), pp->gzo, h) ;
		  aceOutf (ao, "@HD VN:1.5\tSO:queryname\n") ;
		  aceOutf (ao, "@PG ID:1\tPN:Magic\tVN:%s\n", VERSION) ;
		  
		  for (int chrom = 1 ; chrom <= dictMax (dictG) ; chrom++)
		    {
		      Array dna = arr (pp->bbG.dnas, chrom, Array) ;
		      int ln = dna ? arrayMax (dna) : 0 ;
		      aceOutf (ao, "@SQ\tSN:%s\tLN:%d\n", dictName (dictG, chrom), ln) ;
		    }
		  /* aceOutf (ao, "\tCL:%s", commandBuf) ; */
		}
	      exportSamDo (ao, pp, &bb) ;
	    }
	  t2 = clock () ;
	  cpuStatRegister ("8.Export_ali", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
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
  printf ("--- %s: Export closes doneChan at %d\n", timeBufShowNow (tBuf), n) ;
  channelCloseAt (pp->doneChan, n) ;

  ac_free (h) ;
  return ;
} /* export */

/**************************************************************/

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
      parseReadsDo (pp, &bb, pp->iStep, FALSE) ;
      codeWordsDo (pp, &bb, pp->iStep, FALSE) ;

      if (pp->debug) printf ("+++ %s: Start wholeWork %ld bases againt %ld target bases\n", timeBufShowNow (tBuf), bb.length, bbG.length) ;

      /* sort words */
      for (int k = 0 ; k < NN ; k++)
	if (bb.cwsN[k])
	  newMsort (bb.cwsN[k], cwOrder) ;
      
      /* match hits */
      if (bb.length)
	nn = matchHitsDo (pp, &bbG, &bb) ;
      nnn += nn ;

      /* sorthits */
      if (pp->align && bb.hits)
	{
	  sortHitsFuse (pp, &bb) ;
	  newMsort (bb.hits, hitPairOrder) ;
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
      cpuStatRegister ("5.WholeWork", pp->agent, bb.cpuStats, t1, t2, nnn) ;
      channelPut (pp->aeChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      channelCloseAt (pp->aeChan, n) ;
    }

  return ;
} /* wholeWork */

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
/* check the existence of the input files
 * identify their absolute file names
 * associate each file to a run name and to its optional parameters
 */
static Array parseInConfig (PP *pp, Array runStats)
{
  Array rcs = arrayHandleCreate (64, RC, pp->h) ;
  RC *rc = 0 ;
  int nn = 0 ;
  
  if (pp->inFileName)
    {   /* Split the individual file names, they are coma separated 
	 * Make the names absolute, and load them in the file-parser fp channel
	 *  all names must be loaded at once, otherwise we may have a deadlock
	 *  because the agents needing genome parser completion are not yet launched 
	 */
      char *buf = strnew (pp->inFileName, pp->h) ;
      char *cp, *cq, *cr, *cs ;
      int run = 0 ;
      char *filName2 = 0 ;
      if (! pp->runName)
	pp->runName = "runX" ;
      dictAdd (pp->runDict, pp->runName, &run) ;
      cp = buf ; *(cp + strlen (cp)) = ',' ;
      while ((cq = strchr (cp, ',')))
	{
	  cs = strchr (cp, ':') ;
	  if (cs)
	    { /* switch runName */
	      *cs = 0 ;
	      dictAdd (pp->runDict, cp, &run) ;
	      cp = cs + 1 ;
	    }
	  *cq = 0 ;
	  if (! *cp)
	    continue ;
	  cs = strchr (cp, '+') ;
	  if (cs)
	    {
	      *cs = 0 ;
	      cr = filName (cs+1, 0, "r") ;
	      if (! cr)
		messcrash ("\nCannot open input file %s\n", cs+1) ;
	      filName2 = strnew (cr, pp->h) ;
	    }
	  

	  rc = arrayp (rcs, nn++, RC) ; 
	  rc->fileName1 = strnew (cp, pp->h) ;
	  rc->fileName2 = filName2 ;
	  rc->run = run ;

	  rc->format = FASTA ; /* default */
	  rc->format = FASTA ; /* default */
	  if (strstr (rc->fileName1, ".fasta")) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fna")) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fa")) rc->format = FASTA ;
	  else if (strstr (rc->fileName1, ".fastq")) rc->format = FASTQ ;
	  else if (strstr (rc->fileName1, ".fastc")) rc->format = FASTC ;
	  else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	  
	  /* user can override the defaults */
	  if (pp->raw) rc->format = RAW ;
	  if (pp->fasta) rc->format = FASTA ;
	  if (pp->fastq) rc->format = FASTQ ;
	  if (pp->fastc) rc->format = FASTC ;
	  if (pp->sra) rc->format = SRA ;

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		messcrash ("\nCannot open input file %s\n", cp) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    {
	      rc->pairedEnd = TRUE ;
	      rs->nFiles++ ;
	    }
	  if (rc->format == FASTC)
	      rc->pairedEnd = TRUE ;
	  
	  cp = cq + 1 ;
	}
    }
  else if (pp->inConfigFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *fDict = dictHandleCreate (64, h) ;
      ACEIN ai = aceInCreate (pp->inConfigFileName, 0, h) ;
      int nRuns = 0, run = 0 ;
      int line = 0 ;

      if (! ai)
	messcrash ("\nCannot find file -I %s specified on the command line", pp->inConfigFileName) ;
      while (aceInCard (ai))
	{
	  char *cq, *cr, *cp = aceInWord (ai) ;
	  
	  line++ ;
	  if (! cp || ! *cp || *cp == '#')
	    continue ;
	  /* file names */
	  rc = arrayp (rcs, nn++, RC) ; 
	  nRuns++ ;
	      /* file pairs */
	  cq = strchr (cp, '+') ;
	  if (cq) 
	    {
	      *cq++ = 0 ;
	      if (*cq)
		{
		  cr = filName (cq, 0, "r") ;
		  if (! cr)
		    messcrash ("\nCannot open input file %s\n", cq) ;
		  rc->pairedEnd = TRUE ;
		  rc->fileName2 = strnew (cr, pp->h) ;
		}
	    }
	  rc->fileName1 = strnew (cp, pp->h) ;
	  if (! dictAdd (fDict, cp, 0))
	    messcrash ("\nDuplicate file name %s\n at line %d of file -T %s\n try sortalign --help\n"
		       , cr
		       , line
		       , pp->inConfigFileName
		       ) ;

	  /* run name */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (cp && *cp && *cp != '#')
	    dictAdd (pp->runDict, cp, &run) ;
	  else
	    dictAdd (pp->runDict, hprintf (h, "r.%d", nRuns), &run) ;
	  rc->run = run ;
	  rc->RNA = TRUE ; /* default */

	  /* options */
	  cp = aceInWord (ai) ;
	  while (cp)
	    {
	      if (*cp == '#')
		break ;
	      cq = strchr (cp, ',') ;
	      if (cq)
		*cq++ = 0 ;
	      if (! strcasecmp (cp, "fasta")) rc->format = FASTA ;
	      if (! strcasecmp (cp, "fastq")) rc->format = FASTQ ;
	      if (! strcasecmp (cp, "fastc")) rc->format = FASTC ;
	      if (! strcasecmp (cp, "raw")) rc->format = RAW ;
	      if (! strcasecmp (cp, "SRA")) rc->format = SRA ;

	      if (! strcasecmp (cp, "rna")) rc->RNA = TRUE ;
	      if (! strcasecmp (cp, "dna")) rc->RNA = FALSE ;
	      
	      cp = cq ;
	    }

	  if (! rc->format)
	    {
	      rc->format = FASTA ; /* default */
	      if (strstr (rc->fileName1, ".fasta")) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fna")) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fa")) rc->format = FASTA ;
	      else if (strstr (rc->fileName1, ".fastq")) rc->format = FASTQ ;
	      else if (strstr (rc->fileName1, ".fastc")) rc->format = FASTC ;
	      else if (! strncmp (rc->fileName1, "SRR", 3)) rc->format = SRA ;
	    }

	  if (rc->format == SRA)  /* check in the cache */
	    {
	    }
	  else
	    {
	      cr = filName (rc->fileName1, 0, "r") ;
	      if (! cr)
		messcrash ("\nCannot open sequence file %s listed in configuration file -I %s\n\tDid you mean -i %s\n\tPlease try sortalign --help\n", cp, pp->inConfigFileName, pp->inConfigFileName) ;
	      rc->fileName1 = strnew (cr, pp->h) ;
	    }
	  
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    rs->nFiles++ ;
	  if (rc->fileName2 || rc->format == FASTC)
	      rc->pairedEnd = TRUE ;
	} 
      
      ac_free (h) ;
    }
  pp->nFiles = arrayMax (rcs) ;
  printf ("Found %d sequence file%s\n", arrayMax (rcs), arrayMax (rcs) > 1 ? "s" : "") ;
  return rcs ;
} /* parseInConfig */    

/*************************************************************************************/

static int parseWiggleFileList (PP *pp)
{
  ACEIN ai ;
  AC_HANDLE wh = ac_new_handle () ;
  char *cp, cutter ;
  
  pp->wiggleFileDict = dictHandleCreate (1024, pp->h) ;
  if (pp->inFileName)
    ai = aceInCreateFromText (pp->inFileName, 0, wh) ;
  else
    ai = aceInCreate (pp->inConfigFileName, 0, wh) ;
  
  while (aceInCard (ai))
    {
      while ((cp = aceInWordCut (ai, ",", &cutter)))
	{
	  char *cr = filName (cp, 0, "r") ;
	  if (! cr)
	    messcrash ("\nCannot open wiggle file %s\n", cp) ;
	  dictAdd (pp->wiggleFileDict, cr, 0) ;
	}
    }
  ac_free (wh) ;

  return dictMax (pp->wiggleFileDict) ;
} /* parseWiggleFileList */

/*************************************************************************************/
/* check the existence of the target files
 * identify their absolute file names
 * associate each file to a target class and to its optional parameters
 */
static Array parseTargetConfig (PP *pp, Array runStats)
{
  Array tcs = arrayHandleCreate (64, TC, pp->h) ;
  TC *tc = 0 ;
  int nn = 0 ;
  
  if (pp->tFileName)
    {
      char *cr ;
      const char *cp = pp->tFileName ;
      cr = filName (cp, 0, "r") ;
      if (! cr)
	messcrash ("\nCannot open the target file -t %s\n", cp) ;
      tc = arrayp (tcs, nn++, TC) ;
      tc->fileName = strnew (cr, pp->h) ;
      tc->targetClass = 'G' ;
      
      tc->format = FASTA ; /* default */
      if (strstr (tc->fileName, ".fasta")) tc->format = FASTA ;
      if (strstr (tc->fileName, ".fna")) tc->format = FASTA ;
      if (strstr (tc->fileName, ".fa")) tc->format = FASTA ;

      /* user can override the defaults */
      if (pp->raw) tc->format = RAW ;
      if (pp->fasta) tc->format = FASTA ;
    }
  else if (pp->tConfigFileName)
    {
      AC_HANDLE h = ac_new_handle () ;
      DICT *fDict = dictHandleCreate (64, h) ;
      ACEIN ai = aceInCreate (pp->tConfigFileName, 0, h) ;
      int nn = 0, line = 0 ;

      while (aceInCard (ai))
	{
	  char cc, *cq, *cr, *cp = aceInWord (ai) ;
	  
	  line++ ;
	  if (! cp || ! *cp || *cp == '#')
	    continue ;
	  /* target class */
	  if (!cp || !*cp || *cp == '#')
	    continue ;
	  cc = *cp ;
	  if (!(cc >= 'A' && cc <= 'Z')) cc = 0 ;
	  if (cp[1]) cc = 0 ;
	  if (! cc)
	    messcrash ("\n\nThe target class must be specified as a single character (A-Z), not %s,  at line %d of -T target config file %s\n try sortalign --help\n"
		       , cp
		       , line
		       , pp->tConfigFileName
		       ) ;
	  tc = arrayp (tcs, nn++, TC) ;
	  tc->targetClass = cc ;
	  if (cc == 'M') tc->bonus = 1 ;
	  if (cc == 'R') tc->bonus = 8 ;
	  if (cc == 'C') tc->bonus = 8 ;	  
	  /* file names */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  if (!cp || !*cp || *cp == '#')
	    messcrash ("\nNo file name at line %d of file -T %s\n try sortalign --help\n"
		       , line
		       , pp->tConfigFileName
		       ) ;
	  cr = filName (cp, 0, "r") ;
	  if (! cr)
	    messcrash ("\nCannot open the target file -t %s\n", cp) ;
	  if (! dictAdd (fDict, cr, 0))
	    messcrash ("\nDuplicate target file name %s\n at line %d of file -T %s\n try sortalign --help\n"
		       , cr
		       , line
		       , pp->tConfigFileName
		       ) ;
	  tc->fileName = strnew (cr, pp->h) ;
	  tc->format = FASTA ; /* default */
	  if (cc == 'I')
	    {
	      if (strstr (cp, ".introns") || strstr (cp, ".gff") || strstr (cp, ".gtf"))
		tc->format = INTRONS ;
	      else
		messcrash ("\n\nThe Introns must be specified via a .gtf or a .gff file.\n try sortalign --help\n"
			   , cp
			   , line
			   , pp->tConfigFileName
		       ) ;
	    }
	  /* options */
	  aceInStep (ai, '\t') ;
	  cp = aceInWord (ai) ;
	  
	  while (cp)
	    {
	      if (*cp == '#')
		break ;
	      cq = strchr (cp, ',') ;
	      if (cq)
		*cq++ = 0 ;
	      if (! strcasecmp (cp, "fasta")) tc->format = FASTA ;
	      if (! strcasecmp (cp, "raw")) tc->format = RAW ;
	      if (!strncasecmp (cp, "bonus=", 6))
		{
		  int bonus = atoi(cp+6) ;
		  if (bonus > -50 && bonus < 50)
		    tc->bonus = bonus ;
		  else
		    messcrash ("\nBonus %d too large at line %d of file -T %s\n try sortalign --help\n"
			       , bonus
			       , line
			       , pp->tConfigFileName
			       ) ;
		}
	      cp = cq ;
	    }
	} 
      ac_free (h) ;
    }
  printf ("Found %d target files\n", arrayMax (tcs)) ;
  return tcs ;
} /* parseTargetConfig */    
  

// Check if executable exists in PATH
static BOOL isExecutableInPath (const char *name)
{
  const char *path = getenv("PATH") ;
  if (!path)
    return FALSE ;
  char *path_copy = strdup(path) ;
  if (path_copy)
    {
      char *dir = strtok(path_copy, ":") ;
      while (dir)
	{
	  char full_path[1024] ;
	  snprintf (full_path, sizeof(full_path), "%s/%s", dir, name) ;
	  if (access(full_path, X_OK) == 0)
	    {
	      free (path_copy) ;
	      return TRUE ;
	    }
	  dir = strtok (NULL, ":") ;
	}
      free (path_copy) ;
    }
  return FALSE ;
}

/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message, int argc, const char **argv)
{
  int i ;

  if (message)
    {
      fprintf (stderr, "########## Error : %s\n", message) ;
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
	       "// --gzi\n"
	       "//    Forces decompression of the input files, useful when piping into sortalign\n"
	       "//    All files named *.gz are automatically decompressed\n"
	       "//  --sraDownload SRR123,SRR456,SRR999\n"
	       "//    Download a set of SRR entries into the cache directory ./SRA/SRR123.fasta.gz ... \n"
	       "//    This command is optional, since direct SRA access via (-I config) will also cache the sequence files\n"
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
} /* usage */

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

  if (argc < 2)
    usage (0, 0, argv) ;
  if (getCmdLineBool (&argc, argv, "-help")  ||
      getCmdLineBool (&argc, argv, "--help") ||
      getCmdLineBool (&argc, argv, "-h")
      )
    usage (0, 0, argv) ;

  if (getCmdLineBool (&argc, argv, "--version"))
    {
      fprintf (stderr, "sortalign version 0.0.12, august 2025") ;
      exit (0) ;
    }     

  {{
      const char *sraID = 0 ;
      if (getCmdLineText (h, &argc, argv, "--sraDownload", &sraID))
	{
	  char *buf = strnew (sraID, h) ;
	  char *cp = buf ;

	  while (cp)
	    {
	      char *cq = strchr (cp, ',') ;
	      if (cq) *cq = 0 ;
	      sraDownload (cp) ;
	      cp = cq ? cq + 1 : 0 ;
	    }
	  exit (0) ;
	}
    }}

  /*
  */
  /***************** pin the threads to the processors ***/

  if (! getCmdLineBool (&argc, argv, "--numactl")  &&
      !  (getenv("INVOCATION_NOTIFICATIONS") && strstr(getenv("INVOCATION_NOTIFICATIONS"), "numactl")) &&
      isExecutableInPath ("numactl")
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

      
      if (0)
	{
	  fprintf (stderr, "nodes=%d\n", nodes) ;
	  exit (0) ;
	}



      vtxtPrintf (txt, "/usr/bin/numactl  --cpunodebind=%d --membind=%d ", nodes, nodes) ;
      for (int i = 0 ; i < argc ; i++)
	vtxtPrintf (txt, " %s " , argv[i]) ;
      vtxtPrintf (txt, " --numactl ") ;

      fprintf (stderr, "%s\n", vtxtPtr (txt)) ;
      return system (vtxtPtr (txt)) ;
    }
  
  p.debug  = getCmdLineText (h, &argc, argv, "--debug", 0) ;
  p.debug |= getCmdLineText (h, &argc, argv, "--verbose", 0) ;
  
  /**************************  debugging tools, ignore *********************************/

  if ( getCmdLineInt (&argc, argv, "-n", &n))
    {
      BigArray b, a = bigArrayHandleCreate (n, HIT, 0) ;

      
      for (int i = 0 ; i < n ; i++)
	bigArray (a, i, HIT).read = ( !(n &0x1) ? (n - i)%7 :  (unsigned int)randint()) ;  

      b = bigArrayHandleCopy (a, 0) ;
      newMsort (b, hitReadOrder) ;
      
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

  /* another debugging tool, ignore */
  
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

  /***************** select the requested analyses, not used if --createIndex *************/
  p.align = getCmdLineBool (&argc, argv, "--align") ;          /* left for back compatibility, over-ridden by --do_not_align */
  p.align = ! getCmdLineBool (&argc, argv, "--do_not_align") ; /* default is to align */
  
  p.wiggle = getCmdLineBool (&argc, argv, "--wiggle") ;
  p.wiggle = FALSE ; /* code not ready */
  p.ignoreIntronSeeds = getCmdLineBool (&argc, argv, "--ignoreIntronSeeds") ;


  /* future options ***/
  p.snps = getCmdLineBool (&argc, argv, "--snp") ;
  p.introns = getCmdLineBool (&argc, argv, "--intron") ;
  
  /*****************  file names and their formats  ************************/
  
  if (! (p.createIndex = getCmdLineText (h, &argc, argv, "--createIndex", &(p.indexName))) &&
      ! getCmdLineText (h, &argc, argv, "-x", &(p.indexName))
      )
    getCmdLineText (h, &argc, argv, "--index", &(p.indexName)) ;
  getCmdLineText (h, &argc, argv, "-i", &(p.inFileName)) ;
  getCmdLineText (h, &argc, argv, "-I", &(p.inConfigFileName)) ;
  getCmdLineText (h, &argc, argv, "-t", &(p.tFileName)) ;
  getCmdLineText (h, &argc, argv, "-T", &(p.tConfigFileName)) ;
  getCmdLineText (h, &argc, argv, "-o", &(p.outFileName)) ;
  p.runName = 0  ; /* default */

  /* formats are only used in -t case, in -T case provide the format in the 3rd columnh of the tConfig file */
  p.fasta = getCmdLineBool (&argc, argv, "--fasta");
  p.fastq = getCmdLineBool (&argc, argv, "--fastq");
  p.fastc = getCmdLineBool (&argc, argv, "--fastc");
  p.raw = getCmdLineBool (&argc, argv, "--raw");
  p.sra = getCmdLineBool (&argc, argv, "--sra");

  p.gzi = getCmdLineBool (&argc, argv, "-gzi") ||
    getCmdLineBool (&argc, argv, "--gzi");
  p.gzo = getCmdLineBool (&argc, argv, "-gzo") ||
    getCmdLineBool (&argc, argv, "--gzo");
  p.sam = getCmdLineBool (&argc, argv, "--sam") ;
  if (p.sam)
    {
      p.exportSamSequence = TRUE ;
      p.exportSamQuality = TRUE ;
    }
  
  /***************** run name, only use in conjunction with -t, not use with -T ***********/

  p.run = 0 ;
  getCmdLineText (h, &argc, argv, "-r", &(p.runName)) ;
  getCmdLineText (h, &argc, argv, "--run", &(p.runName)) ;

  /***************** method name, only used in some output files, convenient when optimizing parameters */
  
  getCmdLineText (h, &argc, argv, "--method", &(p.method)) ;
  
  /*****************  seed length, only used when createIndex  ************************/

  NN = 16 ; /* default */
  getCmdLineInt (&argc, argv, "--NN", &(NN));
  p.seedLength = 18 ; /* default */
  getCmdLineInt (&argc, argv, "--seedLength", &(p.seedLength));
  
  p.maxTargetRepeats = 31 ;  /* was 81  31 12 */
  getCmdLineInt (&argc, argv, "--maxTargetRepeats", &p.maxTargetRepeats) ;

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

  if (sizeof (long unsigned int) != 8)
    messcrash ("The source code assumes that long unsigned ints use 64 bits not %d, sorry", 8 * sizeof (long unsigned int)) ;

  /*****************  Optional technical parameters ************************/
  nAgents = 40 ;
  if (! getCmdLineInt (&argc, argv, "--nAgents", &(nAgents)))
    getCmdLineInt (&argc, argv, "--nA", &(nAgents)) ;

  p.nBlocks = 40 ;  /* max number of BB blocks processed in parallel */
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
  
  /****************** Check the existence of all file names ****************************************/ 

  if (p.outFileName)
    {
      BOOL ok = FALSE ;
      int n = randint (), n2 = 0 ;
      AC_HANDLE h = ac_new_handle () ;
      char *cp = strrchr (p.outFileName, '/') ;
      if (cp)
	{
	  *cp = 0 ;
	  char *cmd = hprintf (h, "mkdir -p %s", p.outFileName) ;
	  system (cmd) ;
	  *cp = '/' ;
	}
      ACEOUT ao = aceOutCreate (p.outFileName, ".test", 0, h) ;
      if (ao)
	{
	  char * nam = strnew (aceOutFileName (ao), h) ;
	  aceOutf (ao, "%d\n", n ) ;
	  ac_free (ao) ;
	  
	  ACEIN ai = aceInCreate (nam, 0, h) ;
	  if (ai)
	    {
	      if (aceInCard (ai))
		aceInInt (ai, &n2) ;
	      if (n == n2)
		ok = TRUE ;
	      ac_free (ai) ;
	    }
	  char *cmd = hprintf (h, "rm %s", nam) ;
	  system (cmd) ;
	}
      ac_free (h) ;
      if (! ok)
	messcrash ("\nUnable to open a test file  in %s*, please check the -o parameter\n", p.outFileName) ;
    }
  
  /* check that the index directory is accessible */
  if (! p.wiggle)
    {
      if (! p.indexName)
	usage ("Missing parameter --index", argc, argv) ;
      else
	{
	  AC_HANDLE h1 = ac_new_handle () ;
	  const char *cp = p.indexName ;
	  char *cq ;
	  if (p.createIndex)
	    {
	      cq = hprintf (h1, "mkdir %s ; \\rm %s/* ; echo test > %s/sortalign.testFile", cp, cp, cp) ;
	      system (cq) ;
	    }
	  cq = hprintf (h1, "%s/sortalign.testFile", p.indexName) ;
	  ACEIN ai = aceInCreate (cq, FALSE, h1) ;
	  if (! ai)
	    usage (hprintf (p.h, "Cannot create and access the index directory: %s", p.indexName), argc, argv) ;
	  ac_free (h1) ;
	}
    }
  
  /*****************  Check tat all parameters have been parsed *******************/
  if (argc > 1)
    usage (0, argc, argv) ;

  showAli (0) ; /* for compiler happiness */
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
      /* check that input files were provided */
      if (! p.tFileName && ! p.tConfigFileName)
	usage ("--createIndex requires providing a target parameter -t or -T", argc, argv) ;
      if ( p.tFileName &&  p.tConfigFileName)
	usage ("conflicting parameters -t and -T, both define the targets", argc, argv) ;

      Array tArray = parseTargetConfig (&p, runStats) ;
      p.tMaxTargetRepeats = p.maxTargetRepeats ;
      createTargetIndex (&p, &p.bbG, tArray) ;
      if (p.tConfigFileName)
	system (hprintf(h, "\\cp %s %s\n", p.tConfigFileName, p.indexName)) ;
      ACEOUT ao = aceOutCreate (filName (p.indexName, "/seedLength", "w") , 0, 0, p.h) ;
      aceOutf (ao, "%d\t%d\t%d\n# SeedLength\ttStep\tmaxTargetRepeats\n", p.seedLength, p.tStep, p.maxTargetRepeats) ;
      ac_free (ao) ;
      
      goto done ;
    }

  /*******************  otherwise verify the existence of the indexes ********************/
  
  /* check that input files were provided */
  if (! p.wiggle)
    {
      if (! p.indexName)
	usage ("missing parameter -x or --createIndex", argc, argv) ;
    }
  if (! p.inFileName && ! p.inConfigFileName)
    usage ("missing parameter -i or -I inputFileName(s)", argc, argv) ;
  if (p.inFileName && p.inConfigFileName)
    usage ("conflicting parameter -i and -I, both define the input files", argc, argv) ;

  if (! p.wiggle)
    {
      BOOL ok = TRUE ;
      int k = 0 ;
      char *cp ;
      NN = 0 ;
      while (1)
	{
	  char *fNam = hprintf (h, "/cws.sortali.%d", k) ;
	  k++ ;
	  cp = filName (p.indexName, fNam, "rb") ;
	  if (cp)
	    {
	      NN++;
	      if (NN == 1)
		{
		  char *cq = cp + strlen(cp) - 2 ;
		  *cq = 0 ;
		  p.tFileBinaryCwsName = strnew (cp, p.h) ;
		}
	    }
	  else
	    break ;
	}
      if (NN != 1 && NN != 2 && NN != 4 && NN != 8 && NN != 16 && NN != 32 && NN != 64)
	messcrash ("The number of indexes NN=%d must be  a power of 2, say 1 2 4 8 ...", NN) ;

      ACEIN ai = aceInCreate (filName (p.indexName, "/seedLength", "r") , 0, p.h) ;
      if (ai)
	{
	  int x = 0, y = 0, z = 0 ;
	  if (aceInCard (ai))
	    {
	      if (aceInInt (ai, &x))
		p.seedLength = x ;
	      aceInStep (ai, '\t') ;
	      if (aceInInt (ai, &y))
		p.tStep = y ;
	      aceInStep (ai, '\t') ;
	      if (aceInInt (ai, &z))
		p.tMaxTargetRepeats = z ;
	    }
	  ac_free (ai) ;
	}
      else
	ok = FALSE ;
      cp = filName (p.indexName, "/dna.sortali", "rb") ;
      if (cp)
	p.tFileBinaryDnaName = strnew (cp, p.h) ;
      else
	ok = FALSE ;
      cp = filName (p.indexName, "/dnaR.sortali", "rb") ;
      if (cp)
	p.tFileBinaryDnaRName = strnew (cp, p.h) ;
      else
	ok = FALSE ;
      cp = filName (p.indexName, "/ids.sortali", "rb") ;
      if (cp)
	p.tFileBinaryIdsName = strnew (cp, p.h) ;
      else
	ok = FALSE ;
      cp = filName (p.indexName, "/coords.sortali", "rb") ;
      if (cp)
	p.tFileBinaryCoordsName = strnew (cp, p.h) ;
      else
	ok = FALSE ;
      if (!ok)
	usage ("Some of the target binary index files are missing, please run sortalign --createIndex <indexName> ", argc, argv) ;
  
      if (p.iStep > p.tStep)      /* impose phasing */
	p.iStep -= (p.iStep % p.tStep) ;
      
      /* check the existence of the input sequence files */
      inArray = parseInConfig (&p, runStats) ;
    }
  
  /***************************/
  /* launch the multiprocessing */

  wego_max_threads (maxThreads) ;

  if (! p.wiggle)
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
      wego_go (genomeParser, &p, PP) ;
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
      wego_go (npCounter, &p, PP) ;
      channelCloseAt (p.npChan, p.nFiles) ; /* close the counter of BB blocks */
      
      /* Read preprocessing agents, they do not require the genome */
      for (int i = 0 ; i < p.nFiles && i < nAgents && i < 10 ; i++)
	{
	  fprintf (stderr, "Launch readParser %d\n", i) ;
	  p.agent = i ;
	  wego_go (readParser, &p, PP) ;
	}
      for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
	{
	  p.agent = i ;

#ifndef YANN
	  wego_go (codeWords, &p, PP) ;
	  wego_go (sortWords, &p, PP) ;
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
#else

	  wego_go (matchHits, &p, PP) ;
	  wego_go (sortHits, &p, PP) ;
	  if (!i || p.align) /* at least 1 agent */
	    {
	      p.agent = 3*i ;
	      wego_go (align, &p, PP) ;
	      p.agent = 3*i + 1 ;
	      wego_go (align, &p, PP) ;
	      p.agent = 3*i + 2 ;
	      wego_go (align, &p, PP) ;
	      p.agent = i ;
	    }
#endif
	  if (!i) /* only 1 export agent */
	    wego_go (export, &p, PP) ;
	}
    }

  if (p.wiggle)
    {
      int nn = parseWiggleFileList (&p) ;

      if (nn)
	{
	  p.wwChan = channelCreate (nn+1, int, p.h) ;
	  channelDebug (p.wwChan, debug, "wwChan") ;

	  for (int ii = 0 ; ii < 4 ; ii++)
	    {
	      p.agent = ii ;
	      wego_go (wiggle, &p, PP) ;
	    }
	  
	  for (int ii = 1 ; ii <= nn ; ii++)
	    channelPut (p.wwChan, &ii, int) ;
	  channelClose (p.wwChan) ;
	}
      else
	nn = 1 ;
      channelCloseAt (p.doneChan, nn - 1) ;
    }

  while (channelGet (p.doneChan, &bb, BB))
    {
      long int n = (bb.hits ? bigArrayMax (bb.hits) : 0) ;

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
	}
      
      /* recycle the unaligned reads */
#ifdef JUNK
      if (bb->step > 1)
	{
	  BB b, *bbR ;
	  bbR = &b ;
	  memset (bbR, 0, sizeof (BB)) ;
	  bbR->h = ac_new_handle () ;
	  bbR->readerAgent = bb->agent ;
	  bbR->lane = 1000 + bb->lane ;
	  bbR->txt1 = vtxtHandleCreate (bbR->h) ;
	  bbR->txt2 = vtxtHandleCreate (bbR->h) ;
	  bbR->errors = arrayHandleCreate (256, int, bbR->h) ;
	  bbR->cpuStats = arrayHandleCreate (128, CpuSTAT, bbR->h) ;
	  bbR->run = bb->run ;
	  bbR->length = 0 ;
	  bbR->dnas = arrayHandleCreate (64, BigArray, bbR->h) ;
	  bbR->dict = dictHandleCreate (NMAX, bbR->h) ;
	  bbR->errDict = dictHandleCreate (NMAX, bbR->h) ;

	  for (int ii = 1 ; ii < bb->nSeqs ; ii++)
	    if (! bitt(bb->isAligned, ii))
	      {
		int k = 0 ;
		bbR->nSeqs++ ;
		dictAdd (bbR->dict, dictName (bb->dict, ii >>1), &k) ;
		k = k << 1 | (ii & 0x1) ;
		array (bbR->dnas, k, Array) = arrayHandleCopy (array (bb->dnas, ii, Array), bbR->h) ;
	      }
	  if (! bbR->nSeqs)
	    ac_free (bbR->h) ;
	  else
	    {


	    }

	}
#endif
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
      if (1)
	{
	  int ns = 0 ;
	  char tBuf[25], tBuf2[25] ;
	  bb.stop = timeNow () ;
	  timeDiffSecs (bb.start, bb.stop, &ns) ;
	  printf ("%s: run %d / lane %d done start %s elapsed %d s, nSeqs %ld nBases %.1g\n",  timeBufShowNow (tBuf), bb.run, bb.lane, timeShow (bb.start, tBuf2, 25), ns, bb.nSeqs, (double)bb.length) ; 
	}
      ac_free (bb.h) ;
    }
  cpuStatExport (&p, cpuStats) ;
  runStatExport (&p, runStats) ;

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
      if (1) /* valgrind, bbG->globalDnaR->base is double freed */
	ac_free (p.bbG.h) ;
    }
  /* wego_log is the thread-safe way to pass messages to stderr */
  
  if (0)   ac_free (h) ; /* blocks on channel cond destroy */
  return 0 ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
 
