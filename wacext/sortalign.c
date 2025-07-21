/*
 * DNA aligner

 * Created April 18, 2025
 * In collaboration with Greg Boratyn, NCBI

 * A new DNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality

 * The algorithm extract words from the traget and from the reads creating 2 long tables
 * Both tables are sorted
 *  sorting the target is a preprocessing, so the time does not count
 *  sorting the reads should be fast, we will try pytorch.sort as a GPU i
 * The two sorted tables are 'merged' which is local in memory
 *  exporting a sequential table, including local bloom bit maps indicating the genome sections
 * That table is sorted per genome section and per read
 *  alignments are extended using the clipalign algorithm
 * Finally auxiliary tables, like introns are exported
 */

/* #define ARRAY_CHECK  */

/*
  cd worm_2024_RSMagic1

\rm -rf WG WGR
sortalign -t _a1.GR.fasta --createIndex WGR
sortalign -t _a1.G.fasta --createIndex WG

run -i _a1.fasta -x WGR --minAli 30 -o _a3.R.hits --align --max_threads 1

\rm -rf WX.1
bin/sortalign -T tConfig --createIndex WX.1 --NN 1

time sortalign -i _unc32.1.fastc -x WX.1 --minAli 30 -o _unc32.sort.hits --align
time sortalign -i _unc32.fastc -x WX.1 --minAli 30 -o _unc32.sort.hits --align

 */


#include "ac.h"
#include "channel.h"
#define VECTORIZED_MEM_CPY
#ifdef VECTORIZED_MEM_CPY
#include <emmintrin.h> // SSE2
#endif

#define WIGGLE_STEP 10
#define METHOD_AL_JUNK

typedef struct nodeStruct { double x ; CHAN *cx, *cy, *cu, *cv, *done ; int k ; } NODE ;

typedef enum {FASTA=1, FASTQ, FASTC, RAW, SRA, INTRON} DnaFormat ;
typedef struct targetClassStruct {
  char targetClass ; /* single char a-z, A-Z */
  int step ;
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
} RC ;
		  
typedef struct runStatStruct {
  int run ;   /* index in p.runDict */
  int nFiles ;
  long int nPairs, nReads ;
  long int nBase1, nBase2 ;
  long int nMultiAligned[11] ;
  long int nAlignedPerTargetClass[256] ;
  long int nBaseAligned1, nBaseAligned2 ;
  long int compatiblePairs ;
  long int nClippedPolyA ;
  long int nClippedVectors ;
  long int nSupportedIntrons ;
  long int nIntronSupports ;
  long int nN, nErr ;
  int GF[256], GR[256] ; /* number of reads aligned per target_class on Forward/Reverse strand */
  Array errors ;  /* substitutions, insertions, deletions counts */
  /* coverage of long transcripts ? */
} RunSTAT ;
		  
static int NN = 1 ;
typedef struct bStruct {
  AC_HANDLE h ;
  int run, lane, readerAgent ;
  RC rc ;
  DICT *dict ;
  BigArray dnaCoords ;   /* offSets of the dna in the globalDna array */
  Array dnas ;           /* Array of const char Arrays */
  Array dnasR ;          /* Their reverse complement, only computed for the genome */
  BigArray globalDna ;   /* concatenation of all sequence DNAs separated by blocks of nnn */
  BigArray globalDnaR ;  /* concatenation of all reverse sequences, in the same order */
  BigArray hits ;        /* BigArray of read<->genome hits */
  BigArray *cwsN ;         /* BigArray of codeWords */
  BigArray *cwsLN ;        /* BigArray of lines in codeWords */
  BigArray *cwsAN ;        /* BigArray of target offsets */
  long unsigned int nSeqs ;  /* number of sequences in bloc */
  long unsigned int length ; /* cumulated number of bases */
  long unsigned int nerr ;   /* cumulated number of errors */
  long unsigned int nAli ;   /* cumulated number of alignments */
  long unsigned int aliDx ;  /* cumulated aligned read length */
  long unsigned int aliDa ;  /* cumulated genome coverage */
  
  BigArray aligns ; /* final alignments */  
  RunSTAT runStat ;
  Array cpuStats ;
  Array errors ;
  Array wiggles ;
  BOOL isGenome ;
  BOOL isRead2 ;
  int step, skips0, skips1, skips2, skips3, skips4, skipsFound, skipsNotFound ;
  vTXT txt1, txt2 ; /* a pair of reusable txt buffer */ 
} BB ;  

typedef struct pStruct {
  AC_HANDLE h ;
  BOOL debug, gzi, gzo ;
  BOOL createIndex ;
  BOOL align ;
  BOOL wiggle ;
  const char *runName ;
  const char *inFileName ;
  const char *inConfigFileName ;
  const char *outFileName ;
  const char *indexName ;
  const char *tFileName ;
  const char *tConfigFileName ;
  const char *tFileBinaryCwsName ;
  const char *tFileBinaryCwsAName ;
  const char *tFileBinaryCwsLName ;
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
  CHAN *moChan ; /* Matcher needs matches to be Ordered */
  CHAN *oaChan ; /* Ordered words (seeds) sent to the Aligner */
  CHAN *awChan ; /* Aligner to Wiggler */
  CHAN *weChan ; /* Aligner results to be Exported */
  CHAN *doneChan ; /* return to main program */

  BB bbG ;  /* genome or genes target */
  Array runStats ;
  BigArray intronSeeds ;
  BOOL fasta, fastq, fastc, raw, solid ;
  BOOL sam ;
  BOOL hasBonus ;
  int bonus[256] ;
  DICT *runDict ;
  DICT *targetClassDict ;
  int run ;
  int nFiles ;  /* number of input sequence files */
  int agent ;  /* instance of the agent */
  int nBlocks ;
  int blocMaxBases ; /* max number of bases read as one bloc */
  int iStep ;        /* default 2, take a read seed every iStep */
  int tStep ;        /* default 3, take a target seed every tStep */
  int maxTargetRepeats ;
  int seedLength ;
  int errMax ;       /* (--align case) max number of errors in seed extension */
  int minAli ;
  int errRateMax ;       /* (--align case) max number of errors in seed extension */
  int OVLN ;
} PP ;

typedef struct codeWordsStruct {
  unsigned int word ; /* 32 bits = 16 bases, 2 bits per base */
  int nam ; /* index in readDict or chromDict */
  int pos ;  /* position, negative on minus strand */
  unsigned int targetClass ;
} __attribute__((aligned(16))) CW ;

#ifdef METHOD_AL
typedef struct codeWordsLineStruct {
  unsigned int line ; /* offsett in target CW table */
} __attribute__((aligned(4))) CWL ;

typedef struct targetAddressStruct {
  int nam ; /* targetClass | index in chromDict */
  int pos ;  /* position, negative on minus strand */
} __attribute__((aligned(8))) CWA ;
#endif

typedef struct hitStruct {
  unsigned int read ;  /* index in readDict */
  unsigned int chrom ; /* index in chromDict */
  unsigned int a1 ;  /* 1 << 31 + (chrom pos, negative on minus strand) */
  unsigned int x1 ;  /* position in read, negative  minus strand */
} __attribute__((aligned(16))) HIT ;

typedef struct alignStruct {
  int read ;
  int targetClass ;
  int chrom ;
  int a1, a2 ;  /* bio position in chrom, a1 < a2 on minus strand */
  int x1, x2 ;  /* bio position in read, x1 < x2 always */
  int w1, w2 ;  /* wiggle coords, rounded, flipped if bb->isRead2 */
  int chain ;
  int score, chainScore ;
  int nN, nErr ;
  int readLength ;
  int errShort, errLong ; /* bb->dict */
  int leftOverhang, rightOverhang ; /* bb->dict */
  Array errors ;
} __attribute__((aligned(32))) ALIGN ;

#
typedef struct cpuStatStruct {
  char nam[32] ;
  int agent ;
  int nB ; /* blocs treated */
  long int n ;
  clock_t tA ; /* time time active */
} CpuSTAT ;
		  
#define step1 256
#define step2 1024
#define step3 4096
#define step4 16384

#define mstep1 255
#define mstep2 1020
#define mstep3 4080
#define mstep4 16320

#include <pthread.h>
#include <time.h>

typedef struct timespec TMS ;
/* bin/sortalign -t TARGET/Targets/hs.genome.fasta.gz -i titi.fastc --align -o tatou */


static int cwOrder (const void *va, const void *vb) ;

/**************************************************************/
/************** utilities *************************************/

static void cpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, int n)
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
static void newInsertionSort (char *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb))
{
  mysize_t i, j ;
  char buf[s] ;
  for (i = 1 ; i < n ; i++)
    {
      j = i - 1 ;
      if ((*cmp) (b + i*s, b + j*s) >= 0)
	continue ;
      memcpy (buf, b + i*s, s) ;
      memcpy (b + i*s, b + j*s, s) ;
      while (j > 0 && (*cmp) (buf, b + (j-1)*s) < 0)
	{
	  memcpy (b + j*s, b + (j-1)*s, s) ;
	  j-- ;
	}
      memcpy (b + j*s, buf, s) ;
    }
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
static void newMsortDo (char *b, long int nn, int s, char *buf, BOOL hitIsTarget, int (*cmp)(const void *va, const void *vb))
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
  
  /* for small n,
   * sort en place using the insertion algorithm (game of taquin)
   */
  if (hitIsTarget && nn < 8)
    {
      newInsertionSort (b, nn, s, cmp) ;
      return ;
    }

  /* otherwise: sort the 2 halves exchanging hit and buf */
  newMsortDo (b01, n1, s, b1, ! hitIsTarget, cmp) ;
  newMsortDo (b02, n2, s, b2, ! hitIsTarget, cmp) ;
  
  /* then merge the 2 sorted halves */
  up = b01 ;
  vp = b02 ;
  wp = b1 ;

#ifdef VECTORIZED_MEM_CPY
  /* code generated by Grok, loads and stores 16bytes (128 bits) */
  if (cmp == cwOrder)
    {
      while (n1 > 0 && n2 > 0)
	{
	  __m128i u = _mm_load_si128((__m128i*)up) ;
	  __m128i v = _mm_load_si128((__m128i*)vp) ;

#ifdef JUNK	  
	  int n = (*(unsigned int*)up > *(unsigned int*)vp) - (*(unsigned int*)up < *(unsigned int*)vp) ;
	  
	  _mm_store_si128((__m128i*)wp, n <= 0 ? u : v) ;
	  wp += s ;
	  up += (n <= 0) ? s : 0 ;
	  vp += (n > 0) ? s : 0 ;
	  n1 -= (n <= 0) ;
	  n2 -= (n > 0) ;
#else
	  int n = (*(unsigned int*)up <= *(unsigned int*)vp) ;
	  
	  _mm_store_si128((__m128i*)wp, n  ? u : v) ;
	  wp += s ;
	  up += n * s ;
	  vp += (1 - n) * s ;
	  n1 -= n ;
	  n2 -= 1 - n ;
#endif
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
} /* newMsortDo */

static void newMsort (BigArray aa, int (*cmp)(const void *va, const void *vb))
{
  long int N = bigArrayMax (aa) ;
  char *cp = (char *) bigArrp (aa, 0, HIT) ;
  int s = aa->size ;
#ifdef VECTORIZED_MEM_CPY
  if (s != 16)
    messcrash ("only works for aligned 16 bytes structures") ;
#endif
  if (N <= 1) return;

  char *buf = malloc (N * s) ;
  newMsortDo (cp, N, s, buf, TRUE, cmp) ;
  free (buf) ;
} /* newMsort */

/**************************************************************/
/**************************************************************/
/* cumulate int global runStats the content of bb->runStats */
static void runStatsCumulate (int run, Array aa, RunSTAT *vp)
{
  RunSTAT *up = arrayp (aa, run, RunSTAT) ;
    
  up->run = run ;
  up->nPairs += vp->nPairs ;
  up->nReads += vp->nReads ;
  up->nBase1 += vp->nBase1 ;
  up->nBase2 += vp->nBase2 ;
    
  for (int i = 0 ; i < 11 ; i++)
    up->nMultiAligned[i] += vp->nMultiAligned[i] ;
  for (int i = 0 ; i < 256 ; i++)
    up->nAlignedPerTargetClass[i] += vp->nAlignedPerTargetClass[i] ;
  up->nBaseAligned1 += vp->nBaseAligned1 ;
  up->nBaseAligned2 += vp->nBaseAligned2 ;
  up->compatiblePairs += vp->compatiblePairs ;
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
  n = (up->word > vp->word) - (up->word < vp->word) ; if (n) return n ;
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
  n = (up->pos > up->targetClass) - (vp->pos > vp->targetClass) ; if (n) return n ;
  /* pos order */
  n = (up->pos > vp->pos) - (up->pos < vp->pos) ; if (n) return n ;
  return 0 ;
} /* intronOrder */

/**************************************************************/

static void showCws (BigArray cws)
{
  long int ii, iMax = bigArrayMax (cws) ;

  for (ii = 0 ; ii < iMax && ii < 30 ; ii++)
    {
      CW *cw = bigArrp (cws, ii, CW) ;
      printf (".. r=%d\t%d\t%u\n", cw->nam,cw->pos, cw->word) ;
    }
  printf (".........\n") ;
} /* showCws */

/**************************************************************/

static void showHits (BigArray hits)
{
  long int ii, iMax = hits ? bigArrayMax (hits) : 0 ;

  for (ii = 0 ; ii < iMax && ii < 30 ; ii++)
    {
      HIT *hit = bigArrp (hits, ii, HIT) ;
      printf (".. r=%u\t%u\t%u\t%u\n"
	      , hit->read, hit->chrom, hit->a1, hit->x1
	) ;
    }
  printf (".........\n") ;
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
  
  bb->dnaCoords = bigArrayHandleCreate (2 * (iMax + 1), unsigned int, bb->h) ;

  /* compute the length of the global DNA (avoid reallocations) */
  for (ln = 0, ii = 1 ; ii < iMax ; ii++)
    {
      dna = array (bb->dnas, ii, Array) ;
      bigArray (bb->dnaCoords, 2*ii, unsigned int) = ln ;    /* off set of dna ii */
      n = arrayMax (dna) ;
      bigArray (bb->dnaCoords, 2*ii + 1, unsigned int) = ln + n ;    /* off set of dna ii */
      ln += n + 16 ;            /* mininal 16 zeroes protection */
      n = n % 16 ;
      if (n)          /* align the data */
	ln += 16 - n ;
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
      n = arrayMax (dna) ;            /* mininal 16 n protection */
      cp = bigArrp (bb->globalDna, ln, unsigned char) ;
      cq = arrp (dna, 0, unsigned char) ;
      memcpy (cp, cq, n) ; 
      messfree (dna->base) ;
      arrayLock (dna) ;
      dna->base = (char *) cp ; 
      ln += n ; cp += n ;
      /* protect each sequence with a terminal A to allow constructing w1 in codeWords */
      memset (cp, 0, 16) ; *cp = A_ ; ln += 16 ; cp += 16 ;
      n = n % 16 ;
      if (n)          /* align the data */
	{  memset (cp, 0, 16 - n) ; ln += 16 - n ; cp += 16 - n ; }
    }

  return ;
} /* globalDnaCreate */
  
/**************************************************************/
/**************************************************************/
#define NAMMAX 1024

static BOOL parseOneSequence (DnaFormat format, char *namBuf, ACEIN ai, Array dna, int *linep)
{
  int n, nn ;
  char *cp ;
  
  arrayMax (dna) = 0 ;
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
	      n = strlen (cp) ;
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
	    (*linep)++ ;
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
	  if (! cp) /* no > identifier in the raninder of the fasta file */
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
			  , ACEIN ai1, Array dna1, int *linep1
			  , ACEIN ai2, Array dna2, int *linep2
			  )
{
  BOOL ok1 = FALSE, ok2 = TRUE ;
  int n ;
  int line1 = *linep1 ;
  int line2 = *linep2 ;
  
  arrayMax (dna1) = arrayMax (dna2) = 0 ;
  memset (namBuf, 0, NAMMAX) ;
  if (ai1)
    ok1 = parseOneSequence (format, namBuf, ai1, dna1, linep1) ;
  if (ok1 && ai2)
    {
      char namBuf2 [NAMMAX] ;
      ok2 = parseOneSequence (format, namBuf2, ai2, dna2, linep2) ;
      if (strcmp (namBuf, namBuf2))
	messcrash ("\nNon matching pair of sequence identifiers %d <> %s\n line %d of file %s\nlne %d of file %s\n"
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
  n = arrayMax (dna2) ;
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

/**************************************************************/

static void sequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome)
{
  AC_HANDLE h = ac_new_handle () ;
  BB b ;
  int NMAX = isGenome ? 500 : 100000 ;  /* was 100000 for 300s on 6.5 GigaBase human 2025_05_23 */
  int BMAX = 200 * NMAX ;
  int nPuts = 0 ;
  int nn = 0, nSeqs = 0 ;
  int lane = 1 ;
  CHAN *chan = 0 ;
  ACEIN ai1 = 0 ;
  ACEIN ai2 = 0 ;
  int line1 = 0, line2 = 0, line10 = 0 ;
  Array dna1, dna2, dnas ;
  char namBufG [NAMMAX+2] ;
  char *namBufX, *namBuf = namBufG + 2 ;
  DnaFormat format = rc ? rc->format : tc->format ;
  const char *fileName1 = rc ? rc->fileName1 : tc->fileName ;
  const char *fileName2 = rc ? rc->fileName2 : 0 ;
  BOOL pairedEnd = rc ? rc->pairedEnd : FALSE ;
  char tBuf[25] ;
  clock_t t1, t2 ;

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
      bb->lane = lane++ ;
      memset (bb, 0, sizeof (BB)) ;
      chan = pp->plChan ;
      namBufX = namBuf ;
    }

  if (! bb->h)
    {  /* isGenome: keep expanding the same bbG if already initialised */
      bb->h = ac_new_handle () ;
      bb->txt1 = vtxtHandleCreate (bb->h) ;
      bb->txt2 = vtxtHandleCreate (bb->h) ;
      bb->errors = arrayHandleCreate (256, int, bb->h) ;
      bb->cpuStats = arrayHandleCreate (128, CpuSTAT, bb->h) ;
      bb->run = rc ? rc->run : 0 ;
      bb->length = 0 ;
      bb->dnas = dnas = arrayHandleCreate (64, BigArray, bb->h) ;
      bb->dict = dictHandleCreate (NMAX, bb->h) ;
    }
  if (pp->debug || isGenome)
    printf ("+++ %s: Start sequence parser %s\n", timeBufShowNow (tBuf), fileName1) ;

  line1 = line2 = 0 ; line10 = 1 ;
  dna1 = arrayHandleCreate (isGenome ? 1 << 20 : 256, unsigned char, bb->h) ;
  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
  while (parseOnePair (format, namBuf, ai1, dna1, &line1, ai2, dna2, &line2)) 
    {
      int mult = 1 ;
      char *cr = 0 ;
      if (format == FASTC && !isGenome)
	{
	  char *cp = strchr (namBuf, '#') ;
	  int k = atoi (cp+1) ;
	  if (k > 1) mult = k ;
	  cr = namBuf + strlen (namBuf) ;
	}
      for (int iMult = 0 ; iMult < mult ; iMult++)
	{
	  if (mult > 1)
	    {
	      sprintf (cr, ".%d", iMult+1) ; 
	    }
	  if (pairedEnd || arrayMax(dna2)) 
	    {
	      int k = strlen (namBuf) ;
	      int nn1, n1 = arrayMax (dna1) ;
	      int nn2, n2 = arrayMax (dna2) ;
	      
	      pairedEnd = TRUE ;
	      nSeqs += 2 ;
	      bb->nSeqs += 2 ;
	      bb->length += n1 + n2 ;
	      bb->runStat.nPairs++ ;
	      bb->runStat.nReads += 2 ;
	      bb->runStat.nBase1 += n1 ;
	      bb->runStat.nBase2 += n2 ;
	      namBuf[k] = '>' ;
	      dictAdd (bb->dict, namBuf, &nn1) ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (iMult < iMult - 1)
		dna1 = arrayHandleCopy (dna1, bb->h) ;
	      else
		dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	      namBuf[k] = '<' ;
	      dictAdd (bb->dict, namBuf, &nn2) ;
	      if (arrayMax (dna2))
		{
		  array (bb->dnas, nn2, Array) = dna2 ;
		  if (iMult < mult - 1)
		    dna2 = arrayHandleCopy (dna2, bb->h) ;
		  else
		    dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
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
	      if (isGenome && ! new)
		messcrash ("\nDuplicate target sequence identifier %s at line %d of file %s\nThe doublon may occur in this file or in a previous target file\n"
			   , namBufX
			   , line10
			   , fileName1
			   ) ;
	      line10 = line1 + 1 ;
	      array (bb->dnas, nn1, Array) = dna1 ;
	      if (iMult < mult - 1)
		dna1 = arrayHandleCopy (dna1, bb->h) ;
	      else
		dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
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
	  bb->lane = lane++ ;
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
	  bb->dict = dictHandleCreate (NMAX, bb->h) ;
	  dna1 = arrayHandleCreate (256, unsigned char, bb->h) ;
	  dna2 = arrayHandleCreate (256, unsigned char, bb->h) ;
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
#ifdef METHOD_AL
static void GenomeCreateAN_LN (const PP *pp, BigArray cws, BB *bb, int kk)
{
  long int iMax ; 
  long int i ;
  int maxRepeats = pp->maxTargetRepeats ;
  unsigned int mask24 = (0x1 << 25) - 1 ;  
  BigArray cwsLN ;
  BigArray cwsAN ;
  CWA *restrict cwa ;
  CWL *restrict cwl ;
  CW *up, *wp, *upMax ;
  long int line = 1 ;
  int sNN = 0 ;
  while (1 << sNN < NN)
    sNN++ ;
  cwsLN = bb->cwsLN[kk] = bigArrayHandleCreate ((0x1L << (32- sNN)), CWL, bb->h) ;
  cwsAN = bb->cwsAN[kk] = bigArrayHandleCreate (bigArrayMax (cws), CWA, bb->h) ;
    /* store all words, except if highly repeated */
  if (1)
    {
      long int ks[21], cumul = 0 ;
      up = bigArrp (cws, 0, CW) ;
      iMax = bigArrayMax (cws) ;
      upMax = up + iMax ;
      memset (ks, 0, sizeof(ks)) ;
      for (i = 0 ; i < iMax ; )
	{
	  int m, n = 1 ;
	  wp = up + 1 ;
	  while (wp < upMax && wp->word == up->word)
	    wp++ ;
	  n = wp - up ;
	  if (n < maxRepeats)
	    {
	      long int wordS = up->word >> sNN ;
	      cwl = bigArrayp (cwsLN, wordS, CWL) ;
	      cwl->line = line ;
	      (cwl + 1)->line = line + n ;
	      cwa = bigArrayp (cwsAN, line, CWA) ;
	      for (wp = up, m = 0 ; m < n ; wp++, m++)
		{
		  cwa = bigArrayp (cwsAN, line + m, CWA) ;
		  cwa->nam = (wp->nam & mask24) | (wp->targetClass << 24) ;
		  cwa->pos = wp->pos ;
		}
	      line += n ;
	    }
	  up += n ; i += n ;

	  if (n > 20) n = 20 ;
	  ks[n]++ ;
	}


      if (kk == 0)
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

  return ;
} /* GenomeCreateAN_LN */


#else
static BigArray GenomeAddSkips (const PP *pp, BigArray cws, BB *bb)
{
  long int iMax ; 
  long int jMax ; 
  long int i, j, k ;
  AC_HANDLE h = bb->h ;
  int maxRepeats = pp->maxTargetRepeats ;
  
  BigArray aa ;
  /*
  const CW *restrict up ;
  const CW *restrict upMax ; 
  CW *restrict vp ;
  */
  CW *up, *vp, *wp, *upMax ;
  unsigned int wordMax = 0xffffffff ;

  /* remove highly repeated words */
  if (1)
    {
      long int ks[21], cumul = 0 ;
      up = bigArrp (cws, 0, CW) ;
      vp = bigArrp (cws, 0, CW) ;
      iMax = bigArrayMax (cws) ;
      upMax = up + iMax ;
      memset (ks, 0, sizeof(ks)) ;
      for (i = 0, j = 0 ; i < iMax ; )
	{
	  int m, n = 1 ;
	  wp = up + 1 ;
	  while (wp < upMax && wp->word == up->word)
	    wp++ ;
	  n = wp - up ;
	  if (n < maxRepeats)
	    {
	      for (m = 0 ; m < n ; m++)
		{
		  if (j < i)
		    *vp = *up ;
		  i++ ; j++ ; up++ ; vp++ ;
		}
	    }
	  else
	    {
	      up += n ; i += n ;
	    }
	  if (n > 20) n = 20 ;
	  ks[n]++ ;
	}
      bigArrayMax (cws) = j ;

      if (1)
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
      vp->targetClass = ii + mstep4 < iMax ? (up + mstep4)->word : wordMax ;
      vp->nam = ii + mstep3 < iMax ? (up + mstep3)->word : wordMax ;
      vp->pos = ii + mstep2 < iMax ? (up + mstep2)->word : wordMax ;
      vp->word = ii + mstep1 < iMax ? (up + mstep1)->word : wordMax ;

      k = iMax - ii ;
      if (k > mstep1)
	k = mstep1 ;
      vp++ ; 
      memcpy (vp, up, k * sizeof (CW)) ;
      vp += k ; up += k ;
      jMax += k ;
      if (jMax > jMax0)
	messcrash ("add skipps error ") ;
    }
  bigArrayMax (aa) = jMax ;
  return aa ;
} /* GenomeAddSkips */

#endif

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
/* if we create N index, suppose N=16, we can mask log_4(NN), i.e. 2, letters
 * effectively indexing 18-mers in 32 bits words
 * using N index also accelerates sorting and searching
 * and consumes less memory for the sorting
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

  memset (cwsN, 0, sizeof (cwsN)) ;
  
  if (step < 1)
    messcrash ("codeWordsDo received step = %d < 1", step) ;
  if (dMax >> 31)
    messcrash ("codeWordsDo received ln=%ls, step=%d, ln/step=%ld > 2G", bb->length, step, dMax) ;

  if (NN > 1)  /* Avoid reallocation of marginal size fluctuations */
    dMax *= 1.2 ;
  for (k = 0 ; k < NN ; k++)
    cwsN[k] = bigArrayHandleCreate (dMax, CW, bb->h) ;
  
  if (0)
    {
      int mem = 0, mx = 0 ; /* megaBytes */
      messAllocStatus (&mem) ;
      messAllocMaxStatus (&mx) ;
      fprintf (stderr, "=== CodeWordDo Allocated %d Mb, max %d Mb\n", mem, mx) ;
    }

  bb->step = step ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;
      int ii, jj, p = 0, cStep ;
      int iMax = arrayMax (dna) ;
      long unsigned int w, wr, w1, w2, wr1 ;
      const int nShift = 62 ;
      unsigned int targetClass = isTarget ? *(dictName (bb->dict, ia)) : 0 ;
      step = iMax < 30 ? 1 : bb->step ;
      w = w1 = w2 = wr = wr1 = 0 ;
      cp = arrp (dna, 0, unsigned char) ;
      for (ii = 0, jj = -wLen, cStep = -wLen, p = 0 ; ii < iMax + 1 ; ii++, jj++, cp++, cStep++) 
	{
	  /* construct a 2*wLen bits long integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base */
	  w2 = w1 ; w1 = w ; wr1 = wr ;
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
	  if (cStep == step) cStep = 0 ;
	  if (p >= wLen + 1 && !cStep && w1 && wr1 && w != w1 && w != w2 && w1 != w2)
	    {
	      long unsigned int x1 =  w1 & maskSeedLn ;
	      long unsigned int xr1 =  (wr1 >> shiftB) & maskSeedLn ;
	      BOOL minus  = (x1 > xr1) ;
	      long unsigned int word = ( minus ? xr1 : x1 ) ;
	      k = word & maskN ;   /* right bits are used to select the relevant table */
	      word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
	      word &= mask32 ;     /* masked down to 32 bits integer */
	      cw = bigArrayp (cwsN[k], bigArrayMax(cwsN[k]), CW) ;
	      cw->nam = ia ;
	      cw->word = word ;
	      cw->pos = 1 + jj ; /* add 1 to avoid 0 = -0 */
	      if (minus) cw->pos = - cw->pos ;
	      /* x = minus ? seedLngth - 1 - j : j - (seedLength -1) ;  signed bio coordinates of base 1 of the w1 seed */
	      cw->targetClass = targetClass ;
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
  clock_t  t1, t2 ;
  

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->lcChan, &bb, BB))
    {
      t1 = clock () ;
      if (pp->debug) printf ("+++ %s: Start code words\n", timeBufShowNow (tBuf)) ;

      codeWordsDo (pp, &bb, pp->iStep, FALSE) ;

      t2 = clock () ;
      cpuStatRegister ("3.CodeWords", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.cwsN[0])) ;
      if (pp->debug) printf ("--- %s: Stop code words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;

      channelPut (pp->csChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
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
#ifdef METHOD_AL
  bbG->cwsAN = halloc (NN * sizeof (BigArray), bbG->h) ;
  bbG->cwsLN = halloc (NN * sizeof (BigArray), bbG->h) ;
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s.%d", pp->tFileBinaryCwsLName, k) ;
      bbG->cwsLN[k] = bigArrayMapRead (fNam, CWL, READONLY, bbG->h) ;
      fNam = hprintf (h, "%s.%d", pp->tFileBinaryCwsAName, k) ;
      bbG->cwsAN[k] = bigArrayMapRead (fNam, CWA, READONLY, bbG->h) ; /* memory map the seed index */
      nn += bigArrayMax (bbG->cwsAN[k]) ;
    }
#else
  bbG->cwsN = halloc (NN * sizeof (BigArray), bbG->h) ;
  for (int k = 0 ; k < NN ; k++)
    {
      fNam = hprintf (h, "%s.%d", pp->tFileBinaryCwsName, k) ;
      bbG->cwsN[k] = bigArrayMapRead (fNam, CW, READONLY, bbG->h) ; /* memory map the seed index */
      nn += bigArrayMax (bbG->cwsN[k]) ;
    }
#endif
  
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
  if (tStep > 1 && pp->iStep > 1)
    {
      int k ;
      for (int i = 2 ; i <= tStep ; i++)
	{
	  k = tStep/i ;
	  if (k * i == tStep) /* i divides tStep */
	    {
	      k = pp->iStep/i ;
	      if (k * i == pp->iStep) /* i divides iStep */
		messcrash ("\nThe target is indexed with step=%d, the requested read step=%d, these number are not relative primes, there will be systematic false negatives,\n please set the argument --iStep %d (default 2) to a lower value", tStep, pp->iStep, pp->iStep) ;
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
      bbG->wiggles = arrayHandleCreate (2*iMax, Array, bbG->h) ;
      for (ii = 1 ; ii <= iMax ; ii++)
	{
	  Array dna = array (bbG->dnas, ii, Array) ;
	  int ln = arrayMax (dna)/WIGGLE_STEP + 1 ;
	  array (bbG->wiggles, 2*ii, Array) = arrayHandleCreate (ln, int, bbG->h) ;
	  array (bbG->wiggles, 2*ii + 1, Array) = arrayHandleCreate (ln, int, bbG->h) ;
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
#ifdef METHOD_AL
      fNam = hprintf (h, "%s/cws.sortali.L.%d", pp->indexName, k) ;
      bigArrayMapWrite (bbG->cwsLN[k], fNam) ;
      fNam = hprintf (h, "%s/cws.sortali.A.%d", pp->indexName, k) ;
      bigArrayMapWrite (bbG->cwsAN[k], fNam) ;
      nn += bigArrayMax (bbG->cwsAN[k]) ;
#else  
      fNam = hprintf (h, "%s/cws.sortali.%d", pp->indexName, k) ;
      bigArrayMapWrite (bbG->cwsN[k], fNam) ;
      nn += bigArrayMax (bbG->cwsN[k]) ;
#endif
    }
  fprintf (stderr, "genomeCreateBinary exported %ld seed records\n", nn) ;

  if (0) showCws (bbG->cwsN[0]) ;
    
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

  tc->step = 1 ;
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
      up->targetClass = a2 ;
    }
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, aceInFileName (ai)) ;
  ac_free (h) ;

  return nn ;
} /* intronParser */

/**************************************************************/

static int intronCodeWords (PP *pp, BB *bbG)
{
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
  const unsigned int mask21 = (0x1 << 21) - 1 ;
  BigArray aa = pp->intronSeeds ; 
  long int ii, iMax = bigArrayMax (aa) ;
  CW *restrict up = 0 ;
  int chrom = 0, a1, a2, da ;
  Array dna = 0 ;
  BOOL isIntronUp ;
  int tc = 'G' - 'A' ;
  
  bigArraySort (aa, intronOrder) ;
  
  up = bigArrp (aa, 0, CW) ;
  for (ii = 0 ; ii < iMax ; ii++, up++)
    {
      long unsigned int w, wr ;
      
      if (up->nam != chrom)
	{
	  chrom = up->nam ;
	  dna = array (dnas, chrom, Array) ;
	}

      a1 = up->pos ;           /* 1-based Gt_ag position */
      a2 = up->targetClass ;   /* 1-based gt_aG position */
      da = (a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1) ;
      if (a1 > a2)
	{ int a0 = a1 ; a1 = a2 ; a2 = a0 ; isIntronUp = TRUE ; }
      else
	isIntronUp = FALSE ;
      
      /* Construct a 30 bases, 60 bits long integer representing wLen bases (A_, T_, G_, C_) using 2 bits per base 
       * first 15 bases belong to the left exon, upstream of Gt_ag
       * It does not matter if the annotated left or rigth exon is shorter than 15 bases
       */
      p = 0 ;
      w = wr = 0 ;
      cp = arrp (dna, a1 - 15, unsigned char) ;
      for (int i = 0 ; i < 15 ; i++, cp++)
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
      /* next wLen - 8 bases belong to the right exon, downstream of gt_aG */
      cp = arrp (dna, a2, unsigned char) ; /* a2 is G in bio coords, so cp is the first base of the second exon */
      for (int i = 0 ; i < 15  ; i++, cp++)
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
      
      if (p == 30)
	for (int i1 = 4 - (wLen - 15)  ; i1 < 11 ; i1++) /* search every word with at least 4 bases on each side */
	  {
	    if (i1 < 0) continue ;
	    long unsigned int x =   (w >> 2 * (30 - wLen - i1)) & maskSeedLn ;
	    long unsigned int xr =  (wr >> 2 * (i1)) & maskSeedLn ;
	    int dx1 = 15 - i1 ;     /* number of letters in first exon */
	    BOOL  minus  = (x > xr) ;
	    long unsigned int word = ( minus ? xr : x ) ;
	    k = word & maskN ;   /* right bits are used to select the relevant table */
	    word >>= nHidden2 ;  /* remove excessive bits, by construction of NN, they are stored in k */
	    word &= mask32 ;     /* masked down to 32 bits integer */
	    cw = bigArrayp (bbG->cwsN[k], bigArrayMax(bbG->cwsN[k]), CW) ;
	    cw->nam = chrom ;
	    cw->word = word ;
	    cw->pos = a1 ; /* in introns, a1 cannot be zero */
	    if (minus ^ isIntronUp) cw->pos = - cw->pos ; /* word is antistrand to the exon-exon junction */
	    cw->targetClass =
	      (0x1 << 31)                 /* bit 31: isIntron */
	      | (a1 > a2 ? 0x1 << 30 : 0) /* bit 30: intron strand */
	      | (tc & 0x1f) << 25         /* bits 25-29 : target class - 'A', must be < 'Z' */
	      | (dx1 & 0xf) << 21         /* bits 22-24 : number of seed letters in the left exon, dx1 < 16 */
	      | (da & mask21)             /* bits 0-21 : intron length */
	      ;
	  nn++ ;
	}
    }
  fprintf (stderr, "+++++++ Coded %ld intron seeds\n", nn) ;
  return nn ;
} /* intronParser */

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

      if (tc->targetClass == 'I')
	continue ; /* we need to parse thge genome before the introns */
      ntc++ ;
      
      memset (&rc, 0, sizeof (RC)) ;
      rc.fileName1 = tc->fileName ;
      rc.format = tc->format ;
      rc.run = nn + 1 ;
      sequenceParser (pp, 0, tc, bbG, ntc == nTc ? 2 :1 ) ; /* 2 for last non-intron target */
      tc->step = (bbG->length < 1<<20) ? 1 : 3 ;
      if (tc->step > 1 && pp->tStep)
	tc->step = pp->tStep ;
      if (! pp->tStep && tc->step > 1)
	pp->tStep = tc->step ;
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

      if (tc->targetClass == 'I')
	intronParser (pp, bbG, tc) ;
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
  printf ("%s : extract the target seeds\n" , timeShowNow (tBuf)) ;
  codeWordsDo (pp, bbG, tc->step, TRUE) ;
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

  printf ("%s : sort the target seeds\n" , timeShowNow (tBuf)) ;
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

  printf ("%s : write the index to disk\n" , timeShowNow (tBuf)) ;
#ifdef METHOD_AL
  bbG->cwsAN = halloc (NN * sizeof(BigArray), bbG->h) ;
  bbG->cwsLN = halloc (NN * sizeof(BigArray), bbG->h) ;
  for (int kk = 0 ; kk < NN ; kk++)
    {
      /* create bbG->cwsLN[k] and bbG->cwsAN[k] */
      GenomeCreateAN_LN (pp, cwsN[kk], bbG, kk) ;
    }
#else
  bbG->cwsN = halloc (NN * sizeof(BigArray), bbG->h) ;
  for (int kk = 0 ; kk < NN ; kk++)
    {
      bbG->cwsN[kk] = GenomeAddSkips (pp, cwsN[kk], bbG) ;
      bigArrayDestroy (cwsN[kk]) ;
    }
#endif
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
  
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->csChan, &bb, BB))
    {
      t1 = clock () ;
      if (pp->debug) printf ("+++ %s: Start sort words\n", timeBufShowNow (tBuf)) ;
       for (int k = 0 ; k < NN ; k++)
	 if (bb.cwsN[k])
	   newMsort (bb.cwsN[k], cwOrder) ;
      t2 = clock () ;
      cpuStatRegister ("4.SortWords", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.cwsN[0])) ;
      channelPut (pp->smChan, &bb, BB) ;
      if (pp->debug) printf ("--- %s: Stop sort words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cwsN[0])) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      channelCloseAt (pp->smChan, n) ;
    }

  return ;
} /* sortWords */

/**************************************************************/
/**************************************************************/
#ifdef METHOD_AL
static long int  matchHitsLineDo (const PP *pp, BB *bbG, BB *bb)
{
  BigArray hitsArray = bigArrayHandleCreate(64, BigArray, bb->h);
  long int nn = 0, k = 0 ;
  int nHA = 0, hMax = 100000 ;
  BigArray hits = bigArrayHandleCreate(hMax, HIT, bb->h);
  bigArray (hitsArray, nHA++, BigArray) = hits ;
  int seedL1 = pp->seedLength - 1;
  int sNN = 0 ;
  while (1 << sNN < NN)
    sNN++ ;
  
  for (int kk = 0; kk < NN ; kk++)
    {
      long int iMax = bigArrayMax (bbG->cwsAN[kk]);
      long int j = 0, jMax = bigArrayMax (bb->cwsN[kk]);
      BigArray cwsLN = bbG->cwsLN[kk] ;
      BigArray cwsAN = bbG->cwsAN[kk] ;
      long int lnMax = bigArrayMax (cwsLN) ;
      const CWA *restrict cwa ;
      const CW  *restrict rw ;
      HIT *restrict hit;
        
      if (0)
	{
	  if (1) showCws (bb->cwsN[kk]) ;
	}

      for  (j = 0, rw = bigArrp(bb->cwsN[kk], 0, CW) ; j < jMax ; j++, rw++)  /* for all read-seeds */
	{
	  unsigned int word = rw->word ;
	  unsigned int wordS = word >> sNN ;
	  if (wordS >= lnMax)
	    continue ;
	  unsigned int line = bigArrayp (cwsLN, wordS, CWL)->line ;
	  unsigned int line2 = (wordS < lnMax -1) ?  bigArrayp (cwsLN, wordS + 1, CWL)->line : iMax ;
	  if (!line || !line2)
	    continue ;
	  for (cwa = bigArrp(cwsAN, line, CWA) ; line < line2 ; line++, cwa++)  
	    {
	      /* all lines corresponding to actual hits to wordS */
	      long int a1, x1 ;

	      nn++ ;
	      hit = bigArrayp (hits, k++, HIT) ;
	      hit->read = rw->nam ;
	      hit->chrom = cwa->nam ;
	      /* i.e. cwl->nam | (cw1->targetClass << 24) */
	      a1 = cwa->pos ;
	      x1 = rw->pos ;
	      
	      if (x1 > 0)
		{
		  x1-- ;
		  if (a1 > 0)
		    a1-- ;
		  if (a1 < 0)
		    { a1 = a1 + 1 - seedL1 ; }
		}
	      else /* x1 < 0 but strand info must reside on the a1 coordinate */
		{
		  
		  if (a1 < 0)
		    { a1 = - a1 - 1 ; x1 = -x1 - 1 ; }
		  else
		    { a1 = - a1 - seedL1 + 1 ; x1 = -x1 - 1 ; }
		}
	      
	      hit->a1 = (unsigned int) ((a1 - x1) ^ 0x80000000) ;
	      /* to simplify the sorting */
	      hit->x1 = x1 ;
	    }
	      
	  if (k == hMax)
	    {
	      bigArrayMax (hits) = k  ;
	      hits = bigArrayHandleCreate(hMax, HIT, bb->h);
	      bigArray (hitsArray, nHA++, BigArray) = hits ;
	      k = 0 ;
	    }
	}
      bigArrayDestroy (bb->cwsN[kk]) ; /* the read words are no longer needed */
      bigArrayMax (hits) = k ;
    }

  bb->hits = hitsArray ;
  
  return nn ;
}  /* matchHitsDo */

#else

static long int  matchHitsDo (const PP *pp, BB *bbG, BB *bb)
{
  BigArray hitsArray = bigArrayHandleCreate(64, BigArray, bb->h);
  long int nn = 0, k = 0 ;
  int nHA = 0, hMax = 100000 ;
  BigArray hits = bigArrayHandleCreate(hMax, HIT, bb->h);
  bigArray (hitsArray, nHA++, BigArray) = hits ;
  int seedL1 = pp->seedLength - 1;
  
  for (int kk = 0; kk < NN ; kk++)
    {
      long int i = 0, iMax = bigArrayMax (bbG->cwsN[kk]);
      long int j = 0, jMax = bigArrayMax (bb->cwsN[kk]);
      unsigned int mask = step1 - 1 ;
      const CW *restrict cw = bigArrp(bbG->cwsN[kk], 0, CW) ;
      const CW *restrict rw = bigArrp(bb->cwsN[kk], 0, CW) ;
      const CW *restrict cw1;
      const CW *restrict cwMax = cw + iMax ;
      HIT *restrict hit;
        
      if (0)
	{
	  if (1) showCws (bb->cwsN[kk]) ;
	  if (1) showCws (bbG->cwsN[kk]) ;
	}
      while  (i < iMax && j < jMax)
	{
	  if ((i & mask) == 0)
	    {
	      if (rw->word < (unsigned int) cw->word)
		{
		  cw++ ;
		  i++ ;
		  bb->skips0++;
		  continue;
		}
	      else if (rw->word < (unsigned int) cw->pos)
		{
		  cw += step1 ;
		  i += step1 ;
		  bb->skips1++;
		  continue;
		}
	      else if (rw->word < (unsigned int) cw->nam)
		{
		  cw += step2 ;
		  i += step2 ;
		  bb->skips2++;
		  continue;
		}
	      else if (rw->word < (unsigned int) cw->targetClass)
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
	  if (cw->word < rw->word)
	    { i++; cw++; bb->skipsNotFound++ ; }
	  else if (cw->word > rw->word)
	    { j++ ; rw++ ; }
	  else
	    {
	      long int a1, x1, i1 = i ;
	      bb->skipsFound++ ;
	      for (cw1 = cw ; cw1 < cwMax ; i1++, cw1++)
		{
		  if ((i1 & mask) == 0)
		    continue ;
		  if (cw1->word != rw->word)
		    break ;
		  if (1) /* success */
		    {
		      nn++ ;
		      hit = bigArrayp (hits, k++, HIT) ;
		      hit->read = rw->nam ;
		      hit->chrom = cw1->nam  | (cw1->targetClass << 24) ;
		      a1 = cw1->pos ;
		      x1 = rw->pos ;
		      
		      if (x1 > 0)
			{
			  x1-- ;
			  if (a1 > 0)
			    a1-- ;
			  if (a1 < 0)
			    { a1 = a1 + 1 - seedL1 ; }
			}
		      else /* x1 < 0 but strand info must reside on the a1 coordinate */
			{
			  
			  if (a1 < 0)
			    { a1 = - a1 - 1 ; x1 = -x1 - 1 ; }
			  else
			    { a1 = - a1 - seedL1 + 1 ; x1 = -x1 - 1 ; }
			}
		      
		      hit->a1 = (unsigned int) ((a1 - x1) ^ 0x80000000) ;
		      /* to simplify the sorting */
		      hit->x1 = x1 ;
		    }
		  
		  if (k == hMax)
		    {
		      bigArrayMax (hits) = k  ;
		      hits = bigArrayHandleCreate(hMax, HIT, bb->h);
		      bigArray (hitsArray, nHA++, BigArray) = hits ;
		      k = 0 ;
		    }
		}
	      j++ ; rw++ ;
	    }
	}
      bigArrayDestroy (bb->cwsN[kk]) ; /* the read words are no longer needed */
      bigArrayMax (hits) = k ;
    }
  bb->hits = hitsArray ;

  return nn ;
}  /* matchHitsDo */
#endif  /* METHOD_AL */
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
      long int nn = 0 ;
      t1 = clock () ;

      if (bb.length)
	{
	  if (pp->debug) printf ("+++ %s: Start match %ld bases againt %ld sarget bases\n", timeBufShowNow (tBuf), bb.length, bbG.length) ;
#ifdef METHOD_AL
	  nn = matchHitsLineDo (pp, &bbG, &bb) ;
#else
	  nn = matchHitsDo (pp, &bbG, &bb) ;
#endif
	  if (pp->debug) printf ("--- %s: Stop match hits constructed %ld arrays\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	}
      nnn += nn ;
      t2 = clock () ;
      cpuStatRegister ("5.MatchHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
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

static int alignOrder (const void *va, const void *vb)
{
  const ALIGN *up = va ;
  const ALIGN *vp = vb ;
  int n ;
  n = up->read - vp->read ; if (n) return n ;
  n = up->score - vp->score ; if (n) return n ;
  n = up->chain - vp->chain ; if (n) return n ;
  n = up->chrom - vp->chrom ; if (n) return n ;
  n = up->x1 - vp->x1  ; if (n) return n ;
  n = up->x2 - vp->x2  ; if (n) return n ;
  n = up->a1 - vp->a1  ; if (n) return n ;
  n = up->a2 - vp->a2  ; if (n) return n ;

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
  n = ((up->a1 > vp->a1) - (up->a1 < vp->a1)) ; if (n) return n ;
  n = ((up->x1 > vp->x1) - (up->x1 < vp->x1)) ; if (n) return n ;
    ;
  return n ;
} /* hitReadOrder */

/**************************************************************/
/* order the hits */
static void orderHits (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  memset (&bb, 0, sizeof (BB)) ;
  char tBuf[25] ;
  
  clock_t  t1, t2 ;

  while (channelGet (pp->moChan, &bb, BB))
    {
      BigArray aa = bb.hits ;
      long int iMax = aa ? bigArrayMax (aa) : 0 ;
      if (pp->debug) printf ("+++ %s: Start sort hits merging %ld arrays\n", timeBufShowNow (tBuf), iMax) ;
      t1 = clock () ;

      if (iMax == 0)
	bb.hits = 0 ;
      else if (iMax == 1)
	bb.hits = bigArray (aa, 0, BigArray) ;
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
	  bb.hits = bigArrayHandleCreate (n, HIT, bb.h) ;
	  up = bigArrayp (bb.hits, n - 1, HIT) ; /* make room */
	  up = bigArrayp (bb.hits, 0, HIT) ;
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
      t2 = clock () ;
      if (pp->debug) printf ("....... %s: merge done, arraySort start %ld\n", timeBufShowNow (tBuf), bb.hits ? bigArrayMax (bb.hits) : 0) ;

      if (pp->align && bb.hits)
	{
	  newMsort (bb.hits, hitReadOrder) ;
	  
	  t2 = clock () ;
	  cpuStatRegister ("6.OrderHits", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.hits)) ;
	  
	  if (pp->debug) printf ("--- %s: Stop sort hits %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	}
      channelPut (pp->oaChan, &bb, BB) ;
    }

  int n = channelCount (pp->plChan) ;
  memset (&bb, 0, sizeof (BB)) ;
  bb.isGenome = TRUE ;
  channelPut (pp->doneChan, &bb, BB) ; /* destroy bbG.cws, all Matches are already computed */ 
  channelCloseAt (pp->oaChan, n) ;
  
  return ;
} /* orderHits */

/**************************************************************/
/**************************************************************/

static BOOL alignExtendHit (Array dna, Array dnaG, Array dnaGR, Array err
			    , BOOL isDown, int chromLength
			    , int *a1p, int *a2p, int *x1p, int *x2p
			    , int errCost
			    , int errMax, int errRateMax, int minAli)
{
  int nN = 0, dx = 0, nerr = 0 ;
  int x1 = *x1p, x2 = *x2p ;
  int a1 = *a1p, a2 = *a2p ;
#define MAXJUMP 8
  arrayMax (err) = 0 ;

  if (! isDown)
    {
      Array dna = dnaG ; dnaG = dnaGR ; dnaGR = dna ;
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
    }
  
  aceDnaDoubleTrackErrors (dna, &x1, &x2, TRUE /* bio coordinates, extend = TRUE */
			   , dnaG, dnaGR, &a1, &a2
			   , &nN, err, MAXJUMP, errMax, TRUE, 0) ;
  if (x1 > *x1p || x2 < *x2p)
    return FALSE ;
  dx = x2 - x1 + 1 ;
  if (dx < arrayMax (dna) && dx < minAli)
    return FALSE ;
  
  /* reclip */
  /* clip left errors */
  nerr = arrayMax (err) ;
  if (nerr)
    {
      int i, j, y1 = x1 - 1 ;  /* natural coordinates */
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
	    }
	  else  /* register all remaining errors */
	    {	      
	      wonder = FALSE ;
	      if (j < i) *vp = *up ;
	      j++ ; vp++ ;
	    }
	}
      nerr =  arrayMax (err) = j ; /* some 5' errors are dropped */
    }
  
  /* clip right errors */
  if (nerr)
    {
      int i, y2 = x2 - 1 ;  /* natural coordinates */
      A_ERR *up = arrp (err, nerr - 1, A_ERR) ;

      for (i = nerr - 1 ; i >= 0 ; i--, up--)
	{
	  if (errCost > y2 - up->iShort) /* clip this error */
	    {
	      nerr-- ; /* this 3' error is dropped */
	      x2 = up->iShort ; /* -1 to drop the error, +1 for bio coordinates */
	      a2 = up->iLong  ;
	      y2 = x2 - 1 ; /* natural coordinates */
	    }
	  else
	    break ;
	}
      arrayMax (err) = nerr ;
    }

  dx = x2 - x1 + 1 ;
  if (dx < minAli)
    return FALSE ;
  if (100 * nerr > errRateMax * dx)
    return FALSE ;
  if (errMax >= 0 &&  nerr > errMax)
    return FALSE ;

  if (! isDown)
    {
      a1 = chromLength - a1 + 1 ; a2 = chromLength - a2 + 1 ;
    }
  
  *x1p = x1 ; *x2p = x2 ;
  *a1p = a1 ; *a2p = a2 ;
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
  
  vtxtClear (txt1) ;
  vtxtClear (txt2) ;
  for (ii = 0, sep = "" ; ii < nerr ; ii++, ep++, sep = ",")
    {
      xShort = isUp ? arrayMax (dna) - ep->iShort : ep->iShort + 1 ;
      xLong = ep->iLong + 1 ;
      switch (ep->type)
	{
	  case TYPE80:
	    {
	      char cc1a, cc2a ;
	      int xLongR = arrayMax (dna) - xLong + 1 ;
	      
	      cc1a = arr (dnaG, (isUp ? xLongR - 0 : xLong - 2), unsigned char) ;
	      cc2a = arr (dnaG, (isUp ? xLongR - 1 : xLong - 1), unsigned char) ;
	      
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
	    char cc1o, cc2o, cc1oc, cc2oc, cc1a, cc2a, cc1ac, cc2ac ;
	    int xLongR ;
	    
	    xLongR = arrayMax (dnaG) - xLong + 1 ;
	    
	    cc1a = arr (dna, xShort - 1, char) ;
	    cc1o = arr(dnaG, isUp ? xLongR - 1 : xLong - 1,unsigned char) ;
	    
	    cc2o = isUp ? complementBase[(int)cc1o] : cc1o ; 
	    cc2a = isUp ? complementBase[(int)cc1a] : cc1a ; 
	    cc1oc = dnaDecodeChar[(int)cc1o] ;
	    cc2oc = dnaDecodeChar[(int)cc2o] ;
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    
	    vtxtPrintf(txt2, "%s%d:%c>%c"
		       , sep
		       , xLong 
		       , cc2oc, cc2ac
		       ) ;
	    
	    vtxtPrintf(txt1, "%s%d:%c>%c"
		       , sep
		       , xShort 
		       , cc1oc, cc1ac
		       ) ;
	    
	  }  
	  break ;
	case TROU: 
	  {
	    char *ss = "-", cc1a, cc2a, cc1ac, cc2ac ;
	    int xLongR = arrayMax (dnaG) - xLong + 1 ;
	    
	    cc1a = arr(dnaG, (isUp ? xLongR - 1 + 0 : xLong - 1 + 0), unsigned char) ;
	    cc2a = isUp ? complementBase[(int)cc1a] : cc1a ; 
	    
	    cc1ac = dnaDecodeChar[(int)cc1a] ;
	    cc2ac = dnaDecodeChar[(int)cc2a] ;
	    
	    /*
	      if (chenillette (arrp(dnaG, 0, unsigned char), isUp ? xLongR - 1 - (pp->solid ? 1 : 0) : xLong - 1 + (pp->solid ? 1 : 0), 1, &dxL, &dxR))
	      ss = "*-" ;
	    */
	    
	    vtxtPrintf (txt2,"%s%d:%s%c"
			, sep
			, xLong 
			, ss 
			, cc2ac
			) ;
	    
	    vtxtPrintf (txt1, "%s%d:%s%c"
			, sep
			, xShort
			, ss
			, cc1ac
			) ;
	  }
	  break ;

	case TROU_DOUBLE:
	  {
	    char *ss = "--", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;
	    int xLongR = arrayMax (dnaG) - xLong + 1 ;
	    /*
	      if (chenillette (arrp(pp->dna0, 0, unsigned char), isUp ? xLongR - 2 : xLong - 1, 2, &dxL, &dxR))
	      {
	      ss = "*--" ; 
	      xShort += (pp->solid ? 0 : 0) ;
	      }
	    */
	    cc1a = arr(dnaG, (isUp ? xLongR - 2 + 0 : xLong -1 + 0), unsigned char) ;
	    cc1b = arr(dnaG, (isUp ? xLongR - 2 + 1 : xLong -1 + 1), unsigned char) ;
	    
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
			, xShort + 1
			, ss
			, cc1ac, cc1bc
			) ;
	  }
	  break ;
	case TROU_TRIPLE:
	  {
	    char *ss = "---", cc1a, cc1b, cc1c, cc2a, cc2b, cc2c, cc1ac, cc1bc, cc1cc, cc2ac, cc2bc, cc2cc ;
	    /*
	      if (ii == 0)
	      {
	      lastShort-- ;
	      lastShort-- ;
	      lastShort-- ;
	      }
	    */
	    /*
	      if (chenillette (arrp(pp->dna0, 0, unsigned char), isUp ? xLongR - 3 : xLong - 1, 3, &dxL, &dxR))
	      {
	      ss = "*---" ;
	      xShort += (pp->solid ? 0 : 0 ) ;
	      }
	    */
	    int xLongR = arrayMax (dnaG) - xLong + 1 ;
	    cc1a = arr(dnaG, (isUp ? xLongR - 3 + 0 : xLong -1 + 0), unsigned char) ;
	    cc1b = arr(dnaG, (isUp ? xLongR - 3 + 1 : xLong -1 + 1), unsigned char) ;
	    cc1c = arr(dnaG, (isUp ? xLongR - 3 + 2 : xLong -1 + 2), unsigned char) ;
	    
	    
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
			, xShort + (isUp ? 1 : 0)
			, ss
			, cc1ac, cc1bc, cc1cc
			) ;
	  }
	  break ;
	case INSERTION: 
	  {
	    char *ss = "+", cc1a, cc2a, cc1ac, cc2ac ;
	    /*
	      if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 1, &dxL, &dxR))
	      ss = "*+" ;
	    */
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
			, xShort + (isUp ? -1 : 0)
			, ss
			, cc1ac
			) ;
	  }
	  break ;

	case INSERTION_DOUBLE:
	  {
	    char *ss = "++", cc1a, cc1b, cc2a, cc2b, cc1ac, cc1bc, cc2ac, cc2bc ;
	    /*	    
		    if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 2, &dxL, &dxR))
	      ss = "*++" ;
	    */
	    cc1a = arr (dna, xShort - 1 + (isUp ? -1 : 0), char) ;
	    cc1b = arr (dna, xShort - 1 + (isUp ? -0 : 1), char) ;

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

	    /*
	    if ((ii == 0 || ep->iShort > (ep-1)->iShort + 1) && chenillette (probeDna, xShort - 1, 3, &dxL, &dxR))
	      ss = "*+++" ;
	    */
	    cc1a = arr (dna, xShort - 1 + (isUp ? -2 : 0), char) ;
	    cc1b = arr (dna, xShort - 1 + (isUp ? -1 : 1), char) ;
	    cc1c = arr (dna, xShort - 1 + (isUp ? -0 : 2), char) ;
	    
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
			, xShort + (isUp ? - 1 : 0)
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
  dictAdd (bb->dict, vtxtPtr (txt1), &up->errShort) ;
  dictAdd (bb->dict, vtxtPtr (txt2), &up->errLong) ;
  
  return ;
} /* alignFormatErrors */

/**************************************************************/
/* Establish chain scores, select best */
static void  alignSelectBestChain (const PP *pp, BB *bb, BigArray aaa, Array aa, int errCost, Array dna, int chromA, Array dnaG, Array dnaGR)
{
  ALIGN *up, *vp ;

  int ii, jj ;
  int iMax = arrayMax (aa) ;
  int bestScore = 0, bestDx = 0 ;
  int chain = 0 ;
  int maxIntron = 100000 ;
  BOOL isRead2 = bb->isRead2 ;
  int step = WIGGLE_STEP ;  /* examples s=10, 5, 1 */
  int demiStep = step/2 ;
  if (2*demiStep == step) demiStep-- ; /* examples d=4, 2, 0 */
  /* create chains */
  arraySort (aa, alignOrder) ;
  for (ii = 0, up = arrp (aa, 0, ALIGN) ; ii < iMax ; ii++, up++)
    {
      if (! up->chain)
	up->chain = ++chain ;
      up->score = up->x2 - up->x1 + 1 - up->nErr * errCost ;
      if (ii < iMax - 1)
	{
	  vp = up + 1 ;
	  int da = up->a2 - up->a1 ;
	  int db = vp->a2 - vp->a1 ;
	  int dy = vp->x1 - up->x2 - 1 ;
	  int dc = vp->a1 - up->a2 ;
	  if (da < 0) dc = - dc ;
	  dc-- ;
	  if (vp->read == up->read && vp->chrom == up->chrom &&
	      ((da > 0 && db > 0) || (da < 0 && db < 0)) &&
	      dy < 3 && dc > -3 && dc < maxIntron
	      )
	    {
	      vp->chain = chain ;
	      if (dy < 0)
		up->score += dy  ;
	    }
	}
    }

  
  /* create chain scores */
  for (ii = 0, up = arrp (aa, 0, ALIGN) ; ii < iMax ; ii++, up++)
    {
      int score = up->score, chain = up->chain, dx = up->x2 - up->x1 + 1 ;
      for (jj = ii + 1, vp = up + 1 ; jj < iMax && vp->chain == chain ; jj++, vp++)
	{ score += vp->score ; dx += up->x2 - up->x1 + 1 ; }
      if (pp->hasBonus)
	score += pp->bonus[up->targetClass] ;
      for (jj = ii, vp = up ; jj < iMax && vp->chain == chain ; jj++, vp++)
	vp->score = score ;
      ii = jj - 1 ; up = vp - 1 ;
      if (bestScore < score)
	{ bestDx = dx ; bestScore = score ; }
    }
  /* keep happy few */
  int tc = 0, nChains = 0 ;
  chain = 0 ;
  for (ii = jj = 0, up = vp = arrp (aa, 0, ALIGN) ; ii < iMax ; ii++, up++)
    if (up->score >= bestScore)
      {
	if (ii > jj) *vp = *up ;
	vp++ ; jj++ ;
	if (chain != up->chain)
	  { chain = up->chain ; nChains++ ; 
	    if (tc != up->targetClass)
	      {
		tc = up->targetClass ;
		if (tc)  /* do not double count targetClass 0 */
		  bb->runStat.nAlignedPerTargetClass[tc]++ ;
	      }
	  }
      }  
    else
      arrayDestroy (up->errors) ;
  iMax = arrayMax (aa) = jj ;

  if (iMax)
    {
      bb->runStat.nAlignedPerTargetClass[0]++ ;
      bb->runStat.nMultiAligned[0]++ ;
      if (bb->runStat.nPairs && ((vp-1)->read & 0x1))
	bb->runStat.nBaseAligned2 += bestDx ;
      else
      	bb->runStat.nBaseAligned1 += bestDx ;
      if (nChains > 10) nChains = 10 ;
      bb->runStat.nMultiAligned[nChains]++ ;
    }
  /* sort by score */
  arraySort (aa, alignOrder) ;
  
  /* copy the best chains, increase the block stats */
  long int kMax = bigArrayMax (aaa) ;

  /* stranding */
  for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax && ii < 1 ; ii++, up++)
    {
      vp = bigArrayp (aaa, kMax++, ALIGN) ;
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      int tc = up->targetClass ;
      BOOL s = a1 < a2 ;
      if (a1 != a2)
	{
	  if (bb->runStat.nBase2 && (up->read & 0x1))
	    s = !s ;
	  if (s)
	    bb->runStat.GF[tc]++ ;
	  else
	    bb->runStat.GR[tc]++ ;
	}
    }

  /* register the best alignmenmts */
  for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
    {
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      if (isRead2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ;}
      up->w1 = (a1 + demiStep)/step ;
      up->w2 = (a2 + demiStep)/step ;
      bb->runStat.nErr += up->nErr ;
      bb->nAli++ ;
      bb->aliDx += up->x2 - up->x1 + 1 ;
      bb->aliDa += (up->a1 < up->a2 ? up->a2 - up->a1 + 1 : up->a1 - up->a2 + 1) ;
      if (ii < iMax - 1 && up->chain == (up+1)->chain)
	{
	  /* remove overlap */
	  int dx = (up+1)->x1 - up->x2 - 1 ;
	  if (dx < 0)
	    { bb->aliDx += dx ; bb->aliDa += dx ;}
	}
      if (ii == 0 || up->chain != (up-1)->chain)
	alignFormatLeftOverhang (pp, bb, up, dna, dnaG, dnaGR) ;
      if (ii == iMax - 1 || up->chain != (up+1)->chain)
	alignFormatRightOverhang (pp, bb, up, dna, dnaG, dnaGR) ;
      if (up->errors)
	{
	  unsigned int flip = 0 ;
	  if (bb->runStat.nPairs && (up->read & 0x1))
	    flip = 0x0f ; /* will flip last 4 bits */
	  if (1) mergeErrors (bb->errors, up->errors, flip) ;
	  if (up->chrom == chromA)
	    alignFormatErrors (pp, bb, up, dna, dnaG, dnaGR) ;	  
	}
	  
      arrayDestroy (up->errors) ;
      vp = bigArrayp (aaa, kMax++, ALIGN) ;
      *vp = *up ;
    }
  return ;
} /* alignSelectBestChain */

/**************************************************************/

static void alignDo (const PP *pp, BB *bb)
{
  BOOL debug = FALSE  ;
  AC_HANDLE h = ac_new_handle () ;
  HIT * restrict hit ;
  HIT * restrict h1 ;
  ALIGN *ap ;
  long int ii, jj, iMax = bigArrayMax (bb->hits), kMax = 0 ;
  Array err = arrayHandleCreate (256, A_ERR, h) ;
  BigArray aaa = bigArrayHandleCreate (iMax, ALIGN, h) ;
  Array aa = arrayHandleCreate (128, ALIGN, h) ;
  int a1, a2, x1, x2 ;
  int b1, b2, y1, y2, readOld = 0, chromOld = 0, readA = 0, chromA = 0, read1 = 0 ;
  BOOL isDownOld = TRUE ;
  Array dna = 0, dnaG = 0, dnaGR = 0 ;
  int errMax = pp->errMax ; /* 999999 ; */
  int errRateMax = pp->errRateMax ;
  int minAli = pp->minAli ;
  int chromLength = 0 ;
  int r1 = 0, nh1 = 0 ;
  int errCost = 8 ;
  unsigned int uu = 0, u1 = 0 ;
  unsigned int mask24 = (1 << 24) - 1 ;
  
  for (ii = 0, hit = bigArrp (bb->hits, 0, HIT) ; ii < iMax ; ii++, hit++)
    {
      int read = hit->read ;
      int tc = hit->chrom >> 24 ;
      int chrom = hit->chrom & mask24 ;
      BOOL isDown = TRUE ;

      if (! read || ! chrom)
	continue ;
      uu = hit->a1 ;  /* meaning 1<<31 + a1 - x1 */
      if (read != read1)
	{
	  read1 = read ;
	  if (arrayMax (aa))
	    alignSelectBestChain (pp, bb, aaa, aa, errCost, dna, chromA, dnaG, dnaGR) ;
	  arrayMax (aa) = kMax = 0 ;
	}
      /* count all matching hits of that read */
      if (read != r1 || uu != u1)
	{
	  r1 = read ;
	  nh1 = 1 ;
	  u1 = uu ;
	  for (jj = ii + 1, h1 = hit + 1 ; jj < iMax && h1->read == read && h1->a1 == u1 ; jj++, h1++)
	    nh1++ ;
	}	      
      if (nh1 < 2) /* we may loose exons below 22, do loose all exons below 20 */
	continue ;
      x1 = hit->x1 + 1 ; x2 = x1 + 15 ;    /* bio coordinates */   
      a1 = ((int) (hit->a1 ^ 0x80000000)) + hit->x1 ;  /* signed bio coordinates */
      if (a1 < 0)     /* x1 matches a1 */
	{ a1 = -a1 + 1 ; a2 = a1 - 15 ; isDown = FALSE ; }
      else
	{ a1++ ; a2 = a1 + 15 ; isDown = TRUE ; }

      if (read == readOld && chrom == chromOld && isDown == isDownOld &&
	  x1 >= y1 && x2 <= y2 &&
	  (
	   (isDown && a1 >= b1 && a2 <= b2) ||
	   (! isDown && a1 <= b1 && a2 >= b2)
	   )
	  )
	{
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tDoublet\n", ii, read, x1, x2, chrom, a1, a2) ;
	  hit->read = 0 ; /* remove doublet */
	}
      else 
	{
	  if (read != readA)
	    { readA = read ; dna = arr (bb->dnas, read, Array) ; }
	  if (chrom != chromA)
	    {
	      chromA = chrom ;
	      dnaG = arr (pp->bbG.dnas, chrom, Array) ;
	      dnaGR = arr (pp->bbG.dnasR, chrom, Array) ;
	      chromLength = arrayMax (dnaG) ;
	    }
	  if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tbefore align\n", ii, read, x1, x2, chrom, a1, a2) ;
	  if (alignExtendHit (dna, dnaG, dnaGR, err, isDown, chromLength, &a1, &a2, &x1, &x2, errCost, errMax, errRateMax, minAli))
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tAccepted\n", ii, read, x1, x2, chrom, a1, a2) ;
	      ap = arrayp (aa, kMax++, ALIGN) ;
	      memset (ap, 0, sizeof (ALIGN)) ;
	      ap->read = read ;
	      ap->targetClass = tc ;
	      ap->chrom = chrom ;
	      ap->a1 = a1 ;
	      ap->a2 = a2 ;
	      ap->x1 = x1 ;
	      ap->x2 = x2 ;
	      ap->readLength = arrayMax (dna) ;
	      ap->nErr = arrayMax (err) ;
	      if (ap->nErr)
		ap->errors = arrayCopy (err) ;
	    }
	  else
	    {
	      if (debug) fprintf (stderr, "Hit %ld\tr=%d\t%d\t%d\tc=%d\t%d\t%d\tRejected\n", ii, read, x1, x2, chrom, a1, a2) ;
	      hit->read = 0 ; /* remove false positive */
	    }
	  readOld = read ;
	  chromOld = chrom ;
	  isDownOld = isDown ;
	  b1 = a1 ; b2 = a2 ; y1 = x1 ; y2 = x2 ;
	}
    }
  if (arrayMax (aa))
    alignSelectBestChain (pp, bb, aaa, aa, errCost, dna, chromA, dnaG, dnaGR) ;
  bb->aligns = bigArrayHandleCopy (aaa, bb->h) ; /* resize */

  ac_free (h) ;
  return ;
} /* alignDo */

/**************************************************************/

static void align (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  
  clock_t  t1, t2 ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->oaChan, &bb, BB))
    {
      if (pp->align && bb.hits)
	{
	  if (pp->debug) printf ("--- %s: Start align %lu seeds\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	  t1 = clock () ;

	  alignDo (pp, &bb) ;

	  t2 = clock () ;
	  cpuStatRegister ("7.Align_r", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  if (pp->debug) printf ("--- %s: Stop align %lu ali, %lu mismatches\n", timeBufShowNow (tBuf), bb.nAli, bb.nerr) ;
	}
      channelPut (pp->awChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      channelCloseAt (pp->awChan, n) ;
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
  Array wigF = 0, wigR = 0 ;
  Array wiggles = pp->bbG.wiggles ;
    
  for (ii = 0, ap = bigArrp (bb->aligns, 0, ALIGN) ; ii < iMax ; ap++, ii++)
    {
      int w1 = ap->w1, w2 = ap->w2 ;
      if (ap->chrom != chrom)
	{
	  chrom = ap->chrom ;
	  wigF = array (wiggles, 2*chrom, Array) ;
	  wigR = array (wiggles, 2*chrom + 1, Array) ;
	  if (strchr (dictName (pp->bbG.dict, chrom), '|'))
	    wigF = 0 ;
	}
      if (!wigF)
	continue ;
      if (w1 < w2)
	for (int i = w1 ; i <= w2 ; i++)
	  array (wigF, i, int)++ ;
      if (w1 > w2)
	for (int i = w2 ; i <= w1 ; i++)
	  array (wigR, i, int)++ ;
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
  
  for (ii = 1 ; ii < iMax ; ii += 2)
    {
      Array wigF = array (wiggles, 2*ii, Array) ;
      Array wigR = array (wiggles, 2*ii + 1, Array) ;

      if (wigF && arrayMax (wigF))
	{
	  AC_HANDLE h = ac_new_handle () ;
	  const char *chrom = dictName (dict, ii) ;
	  char *fNam = hprintf (h, ".any.%s.u.f.BF", chrom) ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
	  aceOutDate (ao, "##", "wiggle") ;
	  aceOutf (ao, "track type=wiggle_0\n") ;
	  aceOutf (ao, "fixedStep chrom=%s start=10 step=10\n", chrom) ;

	  for (int jj = 1 ; jj < arrayMax (wigF) ; jj++)
	    aceOutf (ao, "%d\n", array (wigF, jj, int)) ;
	  ac_free (h) ;
	}
      if (wigR && arrayMax (wigR))
	{
	  AC_HANDLE h = ac_new_handle () ;
	  const char *chrom = dictName (dict, ii) ;
	  char *fNam = hprintf (h, ".any.%s.u.r.BF", chrom) ;
	  ACEOUT ao = aceOutCreate (pp->outFileName, fNam, pp->gzo, h) ;
	  aceOutDate (ao, "##", "wiggle") ;
	  aceOutf (ao, "track type=wiggle_0\n") ;
	  aceOutf (ao, "fixedStep chrom=%s start=10 step=10\n", chrom) ;
	  
	  for (int jj = 0 ; jj < arrayMax (wigR) ; jj++)
	    aceOutf (ao, "%d\t%d\n", step * jj, array (wigR, jj, int)) ;
	  ac_free (h) ;
	}
    }
  return ;
} /* wiggleExport */

/**************************************************************/

static void wiggle (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  char tBuf[25] ;
  
  clock_t  t1, t2 ;

  if (pp->agent != 0) /* sequential agent */
    return ;  /* accumulates in  pp->wiggles */
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->awChan, &bb, BB))
    {
      if (pp->wiggle && bb.aligns)
	{
	  if (pp->debug) printf ("--- %s: Start wiggle %lu seeds\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	  t1 = clock () ;

	  bigArraySort (bb.aligns, wiggleOrder) ;
	  wiggleDo (pp, &bb) ;

	  t2 = clock () ;
	  cpuStatRegister ("8.Wiggle", pp->agent, bb.cpuStats, t1, t2, bb.aligns ? bigArrayMax (bb.aligns) : 0) ;
	  if (pp->debug) printf ("--- %s: Stop wiggle %lu ali, %lu mismatches\n", timeBufShowNow (tBuf), bb.nAli, bb.nerr) ;
	}
      channelPut (pp->weChan, &bb, BB) ;
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
      channelCloseAt (pp->weChan, n) ;
    }
  return ;
} /* wiggle */

/********************************************************************/

static void sortAlignTableCaption (const PP *pp, ACEOUT ao)
{
  aceOutf (ao, "## %s Magic aligner\n", timeShowNow()) ;
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
  DICT *dictG = pp->bbG.dict ;
  const char *run = dictName (pp->runDict, bb->run) ;
  AC_HANDLE h = ac_new_handle () ;
  char *runNam = hprintf (h, ".%s.f2.%d.%d.hits", run, bb->readerAgent, bb->lane) ;
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
      int dx = x2 - x1 + 1 ;
      int nerr = ap->nErr ;
      int nN = ap->nN ;

      if (read)
	{
	  aceOutf (ao, "\n%s/%s", run, dictName (dict, read)) ; 
	  aceOutf (ao, "\t%d", ap->score) ;
	  aceOutf (ao, "\t%d", 1) ; /* ap->multiplicity */
	  aceOutf (ao, "\t%d", ap->readLength) ;
	  aceOutf (ao, "\t%d\t%d\t%d", dx, x1, x2) ;
	  
	  aceOutf (ao, "\t%c\t-", ap->targetClass) ;
	  aceOutf (ao, "\t%d", 1) ; /* ap->targetMultiplicity */
	  aceOutf (ao, "\t%s", dictName (dictG, chrom)) ; 
	  aceOutf (ao, "\t%d\t%d\t%d", a1, a2) ;
	  aceOutf (ao, "\t%d\t%d", nN, nerr) ;
	  aceOutf (ao, "\t%s", ap->errShort ? dictName (dict, ap->errShort) : "-") ;
	  aceOutf (ao, "\t%s", ap->errLong ? dictName (dict, ap->errLong) : "-") ;
	  aceOut  (ao, "\t-\t-") ;  /* prefix suffix in genome */
	  aceOutf (ao, "\t%s", ap->leftOverhang ? dictName (dict, ap->leftOverhang) : "-") ;
	  aceOutf (ao, "\t%s", ap->rightOverhang ? dictName (dict, ap->rightOverhang) : "-") ;
	  aceOutf (ao, "\tchain %d 1", ap->chain) ;
	  aceOutf (ao, "\t%d", ap->chainScore) ;
	  aceOut  (ao, "\n") ;
	}
    }
  aceOut (ao, "\n") ;

  ac_free (h) ;
  return ;
} /* exportDo */

/**************************************************************/
/*********************************************************************/
#ifdef JUNK

/* export  SAM
 * hits is  in pExportSamOrder
 * so organized by  read/target/chain in the orientation of the target 
 * as befits a good Cuban CIGAR
 */
static char *clipAlignExportOneSamCigar (CLIPALIGN *pp, PEXPORT *px, int di, vTXT cigar, Array cigarettes, Array dnaShort, Array err, int *nErrp, AC_HANDLE h)
{
  int ii, nErr = 0 ;
  PEXPORT *py ;
  int a1, a2, b2 = 0, x1, x2, y2 ;
  const char *ccp, *errors ;
  int isDown = 1 ;
  int ln = 0 ;

  vtxtClear (cigar) ;

  for (py = px, ii = 0 ; ii < di ; ii++, py++)
    {
      BOOL isGtAg = FALSE ;
      BOOL isGcAg = FALSE ;
      int oldA, ddN = 0, ddA = 0 ;
      int da1 = 0, da2 = 0, intron = 0, intronLn = 0 ;
      a1 = py->a1 ;
      a2 = py->a2 ;
      x1 = py->x1 ;
      x2 = py->x2 ;

      if (a1 < a2) 
	{
	  intron = py->intron ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      isGtAg = !strcmp (pi->foot, "gt_ag") ;
	      isGcAg = !strcmp (pi->foot, "gc_ag") ;
	      if (pi->a2 > pi->a1)
		{
		  da1 = pi->a2 - a1 ;
		  x1 += da1 ;
		  intronLn = pi->a2 - pi->a1 - 1 ;
		}
	      else
		{
		  da1 = pi->a1 - a1 ;
		  x1 += da1 ;
		  intronLn = pi->a1 - pi->a2 - 1 ;
		}
	    }
	  intron = (ii < di - 1) ? py[1].intron : 0 ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      if (pi->a2 > pi->a1)
		{
		  da2 = pi->a1 - a2 ;
		  x2 += da2 ;
		}
	      else
		{
		  da2 = pi->a2 - a2 ;
		  x2 += da2 ;
		}
	    }
	}
      else
	{
	  intron = (ii > 0) ? py[-1].intron : 0 ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      isGtAg = !strcmp (pi->foot, "gt_ag") ;
	      isGcAg = !strcmp (pi->foot, "gc_ag") ;
	      if (pi->a2 > pi->a1)
		{
		  da2 = pi->a2 - a2 ;
		  x2 -= da2 ;
		  intronLn = pi->a2 - pi->a1 - 1 ;
		}
	      else
		{
		  da2 = pi->a1 - a2 ;
		  x2 -= da2 ;
		  intronLn = pi->a1 - pi->a2 - 1 ;
		}
	    }
	  intron = py->intron ;
	  if (intron)
	    {
	      PEXPORT *pi = bigArrp (pp->exportIntrons, intron - 1, PEXPORT) ;
	      if (pi->a2 > pi->a1)
		{
		  da1 = pi->a1 - a1 ;
		  x1 -= da1 ;
		}
	      else
		{
		  da1 = pi->a2 - a1 ;
		  x1 -= da1 ;
		}
	    }
	}
      a1 += da1 ;
      a2 += da2 ;

      errors = py->errLeftVal ? dictName (pp->exportDict, py->errLeftVal) : 0 ;   /* error list in target coordinates and orientation */
      
      if (a1 > a2)
	{
	  int a0 ;
	  a0 = a1 ; a1 = a2 ; a2 = a0 ; a0 = x1 ; x1 = x2 ; x2 = a0 ;
	  isDown = -1 ;
	  if (ii > 0) /* look for an overlap on the read */
	    {
	      ddN = y2 - 1 - x1 ;
	    }
	}
      else if (ii > 0) /* look for an overlap on the read */
	{
	  ddN = x1 - y2 - 1 ;
	}
      b2 += intronLn ;
      ddA = a1 - b2 - 1 ; 


      if (ii == 0) /* gap en tete de sequence */
	{
	  ddN = 0 ; ddA = 0 ;
	  if (isDown == 1 && x1 > 1)
	    { vtxtPrintf (cigar, "%dS", x1 - 1) ;  ln += x1 - 1 ; }
	  else if (isDown == -1 && x1 < arrayMax (dnaShort))
	    { vtxtPrintf (cigar, "%dS", arrayMax (dnaShort) - x1) ; ;  ln += arrayMax (dnaShort) - x1 ; }
	}
      if (intronLn)
	{
	  vtxtPrintf (cigar, "%d%c", intronLn, isGtAg ? 'N' :(isGcAg ? 'n' : 'D')) ;
	}

      if (ddN < 0) /* overlap in x , on detruit cet overlap */ 
	{
	  if (1) { a1 -= ddN ; x1 -= ddN ; ddA -= ddN ; }
	  ddN = 0 ;
	}
      if (ddN > 0 && ddA >= ddN)
	{  /* on alonge l'exon precedent avec des erreurs */ 
	  vtxtPrintf (cigar, "%dX", ddN) ; ln += ddN ;
	  ddA -= ddN ; 
	  b2 += ddN ;
	  y2 += ddN * isDown ;
	  ddN = 0 ;
	}
      if (ddN > 0 && ddA < 0)
	{
	  ddN -= ddA ;
	  a1 -= ddA ;
	  x1 -= ddA * isDown ;
	  ddA = 0 ;
	}
      if (ddN > 0 && ddA > 0 && ddN >= ddA)
	{  /* on alonge l'exon precedent avec des erreurs */ 
	  vtxtPrintf (cigar, "%dX", ddA) ; ln += ddA ;
	  b2 += ddA ;
	  y2 += ddA * isDown ;
	  ddN -= ddA ;
	  ddA = 0 ;
	}
      if (ddN > 0 && ddA == 0)
	{  /* insertion */
	  vtxtPrintf (cigar, "%dI", ddN) ; ln += ddN ;
	  y2 += ddN * isDown ;
	  ddN = 0 ;
	}

      if (ddA < 0) /* duplication */
	{
	  vtxtPrintf (cigar, "%dI", -ddA) ; ln += -ddA ;
	  a1 -= ddA ; x1 -= ddA * isDown ; ddA = 0 ;
	}
      else if (0 && ddA > 30 && pp->strategy == STRATEGY_RNA_SEQ) /* intron */
	{
	  vtxtPrintf (cigar, "%dN", ddA) ; ddA = 0 ;
	}
      else if (ddA > 0)  /* deletion */
	{
	  vtxtPrintf (cigar, "%dD", ddA) ; ddA = 0 ;
	}
 

      oldA = a1 - 1 ;

      if (errors)
	{
	  int a, n, nSub = 0 ;
	  ccp = errors - 1 ;
	  
	  while (ccp[1] && ccp[2])
	    {
	      a = 0 ; while (*++ccp != ':') { a = 10 * a + (*ccp - '0') ; } /* coord of the error */ 
	      if (a < oldA) /* clip this error */
		{
		  while (ccp[1] && *ccp != ',')
		    ccp++ ;
		  continue ;		  
		}
	      
	      n = 0 ; ccp++ ; /* jump the : */
	      if (*ccp == '*') ccp++ ;
	      if (a - oldA  + (*ccp == '+' ? 1 : 0) > 1)
		{
		  int k ;
		  if (nSub)
		    { vtxtPrintf (cigar, "%dX", nSub) ;  ln += nSub ; }
		  nSub = 0 ;
		  k =  a - oldA - 1 + (*ccp == '+' ? 0 : 0) ;
		  if (k)
		    {
		      vtxtPrintf (cigar, "%d=",k) ;
		      ln +=  k ;
		    }
		  oldA  = a ;
		}
	      if (n == 0)
		{
		  while (*ccp == '+') { n++ ; ccp++; }
		  if (n) 
		    { 
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ;  ln += nSub ; }
		      nSub = 0 ; 	
		      oldA = a ;
		      vtxtPrintf (cigar, "%dI", n) ; oldA += 0*n - 1 ; ccp += n ;
		      ln += n ;
		    }
		}
	      if (n == 0)
		{
		  while (*ccp == '-') { n++ ; ccp++; }
		  if (n) 
		    {
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ; ln += nSub ; }
		      nSub = 0 ;
		      oldA = a ;
		      vtxtPrintf (cigar, "%dD", n) ; oldA += n - 1 ; ccp += n ; 
		    }
		}
	      if (n == 0) /* substitution */
		{
		  if (a == oldA  + 1)
		    ;
		  else
		    {
		      if (nSub)
			{ vtxtPrintf (cigar, "%dX", nSub) ; ln +=  nSub ; } 
		      nSub = 0 ;
		      if (a - oldA) 
			{ vtxtPrintf (cigar, "%d=", a - oldA - 1) ; ln += a - oldA - 1 ; }
		    }
		  oldA = a ; 
		  nSub++ ;
		  ccp += 3 ; n = 1 ;
		}
	      nErr += n ; if (! *ccp) break ; 
	    }
	  if (nSub)
	    { vtxtPrintf (cigar, "%dX", nSub) ; ln +=  nSub ; }
	  nSub = 0 ;
	}
      if (a2 > oldA)
	{ vtxtPrintf (cigar, "%d=", a2 - oldA) ; ln +=  a2 - oldA ; }

      y2 = x2 ;
      b2 = a2 ;
    }  

  /* gap en queue */
  if (isDown == 1 && x2 <  arrayMax (dnaShort))
    { vtxtPrintf (cigar, "%dS", arrayMax (dnaShort) - x2) ;  ln += arrayMax (dnaShort) - x2 ; }
  else if (isDown == -1 && x2 > 1)
    { 
      int k = arrayMax (dnaShort) - (ln + x2 - 1) ;
      x2 += k ;
      if (x2 > 1)
	vtxtPrintf (cigar, "%dS", x2 - 1) ;  
      ln += x2 - 1 ; 
    }

  samCheckCigar(0, vtxtPtr(cigar), cigarettes, px->a1, arrayMax (dnaShort)) ;

  *nErrp = nErr ;
  return vtxtPtr (cigar) ;
} /* clipAlignExportOneSamCigar */
#endif

/*********************************************************************/
#ifdef JUNK
/* count the hit lines summarized in this sam line */ 
static int clipAlignExportOneSam (ACEOUT ao, CLIPALIGN *pp
				  , vTXT cigar, Array cigarettes
				  , Array dnaShort, Array err
				  , BigArray hits
				  , long int ii0, long int iMax)
{
  AC_HANDLE h = 0 ;
  int flag, di = 1 ;					
  long int j ;
  PEXPORT *pz = 0, *py, *px = bigArrayp (hits, ii0, PEXPORT) ;
  int nerr = 0 ;
  int a1 = px->a1, a2 = px->a2 ;
  int b1 = 0, b2 = 0 ;
  int target = px->target ;
  int chain = px->chain ;
  BOOL isDown = TRUE ;
  BOOL mateIsDown = TRUE ;
  BOOL isPrimary = TRUE ;
  BOOL isSecondary = FALSE ;
  BOOL isSupplementary = FALSE ;
  BOOL isOrphan = FALSE ; 
  char *seq ;
  MM *mm = px->mm ;
  int fragment = mm ? mm->fragment : 0 ;
  MM *mate = 0 ;

  if (!px->target || ! px->score || ! px->pairScore)
    return 1 ;
  h = ac_new_handle () ;

  if ((a1 < a2 && px->x1 < px->x2) || (a1 > a2 && px->x1 > px->x2))
    isDown = TRUE ;
  else
    isDown = FALSE ;

  if (a1 > a2)
    { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }

  /* the chain starting at px including di lines will be reported as a single CIGAR */
  for (di = 0, j = ii0, py = px ; j < iMax && py->mm == mm && py->target == target && py->chain == chain ; j++, py++)
    {
      if (py->a1 < a1) a1 = py->a1 ;
      if (py->a2 < a1) a1 = py->a2 ;
      di++ ;
    }

  /* check f I am a primary alignment, i.e. first alignment of the current read */
  for (j = ii0, py = px ; j >= 0 && py->mm == mm && py->target == target && py->chain == chain ; j--, py--)
    ;
  if (j >= 0 && py->mm == mm) /* i changed chain or target */
    {
      isPrimary = FALSE ;
      if (px->uu == 1)
	isSupplementary = TRUE ;
      else
	isSecondary = TRUE ;
    }

  /* proceed to the next read of the same fragment */
  if (! strcmp (stackText (pp->probeStack, mm->probeName), "seq.144399>"))
    invokeDebugger () ;
  if (isPrimary && mm->pair)
    {
      /* try upwards */
      mate = 0 ; pz = 0 ;
      for (j = ii0, py = px ; j >= 0 && py->mm->fragment == fragment ; j--, py--)
	;
      j++ ; py++ ;
      if (j >= 0 && py->mm->fragment == fragment && py->mm != mm)
	{ pz = py ; mate = pz->mm ; } /* first hit of other read */
      
      if (mate == 0)
	{ /* try downwards */
	  for (j = ii0, py = px ; j < iMax && py->mm == mm ; j++, py++)
	    ;
	  if (j < iMax && py->mm && py->mm->fragment == fragment && py->mm != mm)
	    { pz = py ; mate = pz->mm ; } /* first hit of other read */
	}
      if (mate) /* select mateTarget and matePosition and fragment length */
	{
	  b1 = b2 = pz->a1 ;
	  if (pz->x1 < pz->x2 && pz->a1 < pz->a2)
	    mateIsDown = TRUE ;
	  else
	    mateIsDown = FALSE ;
	  for (py = pz ; j < iMax && py->mm == pz->mm && py->target == pz->target && py->chain == pz->chain ; j++, py++)
	    {
	      if (py->a1 < b1) b1 = py->a1 ;
	      if (py->a2 < b1) b1 = py->a2 ;
	      if (py->a1 > b2) b2 = py->a1 ;
	      if (py->a2 > b2) b2 = py->a2 ;
	    }
	}
    }
 
   if (isPrimary && mm->pair && ! mate)  
     isOrphan = TRUE ;

  /* there are 11 mandatory columns */
  aceOutf (ao, "%s", stackText (pp->probeStack, mm->probeName)) ;
  flag = 0 ;
  if (1) 
    {   /* SAM FLAGS */
      /* this read belongs to a pair */
      if (pp->hasPairs && mm->pair) flag |= 0x1 ;
      /* all segments are well and compatibly aligned */
      if (pp->hasPairs && mm->pair && mate) flag |= 0x2 ;
      /* unmapped read */
      if (! px->target) flag |= 0x4 ;  
      /* orphan pair-mate is unmapped */
      if (pp->hasPairs && mm->pair && isOrphan) flag |= 0x8 ;
      /* the seq is given in the orientation of the genome
       * the flags say if the shown read is the complement of 
       * read found on the machine fastq file 16 == 0x10,  32 == 0x20
       */
      if (! isDown) flag |= 0x10 ;
      if (pp->hasPairs && mm->pair && mate && ! mateIsDown) flag |= 0x20 ;
      /* i am read 1 of a pair or first or central part of a composite */
      if (pp->hasPairs && mm->pair > 1) flag |= 0x40 ;
      /* i am read 2 of a pair, or last or central part of a composite 128 == 0x80*/
      /* notice that central read parts are 0x4 | 0x8 == 0xb0 */ 
      if (pp->hasPairs && mm->pair < 0) flag |= 0x80 ;
      /* multi-mapping : all have this flag except first one 256 == 0x100 */
      if (isSecondary) flag |= 0x100 ;
      /* read not passing quality filters 512 = 0x200 */
      if (px->score < 10) flag |= 0x200 ;
      /* i am a PCR duplicate: never used by magic 1024 == 0x400 */
      if (0) flag |= 0x400 ;
      /*  supplementary alignment  2048 == 0x800
       *  divers portion du meme read sont incompatibles
       *  all have this flag except the first one
       */
      if (isSupplementary) flag |= 0x800 ;
    }

  aceOutf (ao, "\t%d", flag) ;
  aceOutf (ao, "\t%s", dictName (pp->targetDict, px->target)) ;
  aceOutf (ao, "\t%d", a1) ;
  { /* do we trust the ali */
    int z = 0 ;
    if (
	(! mm->pair || px->pairScore > px->score) &&
	10*px->ali >= 9*px->ln && 100 * px->nErr < px->ali
	)
      z = 33 ;
    else if (px->nErr < 4)
      z = 20 ;
    else
      z = 10 ;
    z = px->pairScore ;
    aceOutf (ao, "\t%d", z) ;
  }
  
  seq =  stackText (pp->probeStack, mm->probeSequence) ;
  array (dnaShort, strlen (seq) + mm->probeLeftClip + mm->probeRightClip, char) = 0 ; /* zero terminate */
  arrayMax (dnaShort) = strlen (seq) + mm->probeLeftClip + mm->probeRightClip  ;
  if (mm->probeLeftClip)
    memset (arrp (dnaShort, 0, char), T_,  mm->probeLeftClip) ;
  memcpy (arrp (dnaShort, mm->probeLeftClip, char), seq, strlen (seq)) ;
  if (mm->probeRightClip)
    memset (arrp (dnaShort, mm->probeLeftClip + strlen (seq), char), A_,  mm->probeRightClip) ;

  if (! isDown)
    reverseComplement (dnaShort) ;
  
/* CIGAR */
  aceOutf (ao, "\t%s", clipAlignExportOneSamCigar (pp, px, di, cigar, cigarettes, dnaShort, err, &nerr, h)) ;
  if (1)
    aceOutf (ao, "\t%s", isPrimary && pz && pz->target ? (px->target == pz->target ? "=" :dictName (pp->targetDict, pz->target)) : "0") ; /* RNEXT name of mate or next chain of same read */
  else
    aceOutf (ao, "\t%s", isPrimary && mate ? (mm == mate ? "=" : stackText (pp->probeStack, mate->probeName)) : "0") ; /* RNEXT name of mate or next chain of same read */
  aceOutf (ao, "\t%d", b1) ; /* PNEXT pmate position or next string if on same target */
  aceOutf (ao, "\t%d", b2 ? (a1 < b2 ? b2 - a1 + 1 : b1 - a2 - 1) : 0) ; /* TLEN observed template length */
  
  if (0 && ! isPrimary)
    aceOutf (ao, "\t*") ;
  else
    {
      dnaDecodeArray (dnaShort) ;
      aceOutf (ao, "\t%s", arrp (dnaShort, 0, char)) ;
    }
  
  if (! isPrimary)
    aceOutf (ao, "\t*") ;
  else
    {
      if (isDown)
	aceOutf (ao, "\t*") ;  /* QUAL ascii de phread scale quaity + 33 ou * */
      else
	aceOutf (ao, "\t*") ; /* we must reverse the order of the QUALs */
    }
  /* end of the 11 mandatory columns */
  
  aceOutf (ao, "\tNH:i:%d", px->uu) ;
  aceOutf (ao, "\tAS:i:%d", px->score) ;
  aceOutf (ao, "\tNM:i:%d", nerr) ;
  aceOutf (ao, "\tmt:i:%d", mm->mult) ;
  /* aceOutf (ao, "\tMD:Z:0") ; */
  
  aceOut (ao, "\n") ;
  ac_free (h) ;
  return di ;
}  /* clipAlignExportOneSam */
#endif
/*********************************************************************/

static void exportOneSam (ACEOUT ao, const PP *pp, BB *bb, vTXT cigar, Array cigarettes, ALIGN *ap)
{
  const char *run = dictName (pp->runDict, bb->run) ;
  DICT *dictG = pp->bbG.dict ;
  DICT *dict = bb->dict ;

  int read = ap->read ;
  int chrom = ap->chrom ;
  int x1 = ap->x1 ;
  int x2 = ap->x2 ;
  int a1 = ap->a1 ;
  int a2 = ap->a2 ;
  int dx = x2 - x1 + 1 ;
  int nerr = ap->nErr ;
  int nN = ap->nN ;
  
  if (read)
    {
      aceOutf (ao, "\n%s/%s", run, dictName (dict, read)) ; 
      aceOutf (ao, "\t%d", ap->score) ;
      aceOutf (ao, "\t%d", 1) ; /* ap->multiplicity */
      aceOutf (ao, "\t%d", ap->readLength) ;
      aceOutf (ao, "\t%d\t%d\t%d", dx, x1, x2) ;
      
      aceOutf (ao, "\t%c\t-", ap->targetClass) ;
      aceOutf (ao, "\t%d", 1) ; /* ap->targetMultiplicity */
      aceOutf (ao, "\t%s", dictName (dictG, chrom)) ; 
      aceOutf (ao, "\t%d\t%d\t%d", a1, a2) ;
      aceOutf (ao, "\t%d\t%d", nN, nerr) ;
      aceOutf (ao, "\t%s", ap->errShort ? dictName (dict, ap->errShort) : "-") ;
      aceOutf (ao, "\t%s", ap->errLong ? dictName (dict, ap->errLong) : "-") ;
      aceOut  (ao, "\t-\t-") ;  /* prefix suffix in genome */
      aceOutf (ao, "\t%s", ap->leftOverhang ? dictName (dict, ap->leftOverhang) : "-") ;
      aceOutf (ao, "\t%s", ap->rightOverhang ? dictName (dict, ap->rightOverhang) : "-") ;
      aceOutf (ao, "\tchain %d 1", ap->chain) ;
      aceOutf (ao, "\t%d", ap->chainScore) ;
      aceOut  (ao, "\n") ;
    }
} /* exportOneSam */

/**************************************************************/
   
static void exportSamDo (ACEOUT ao, const PP *pp, BB *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ii, aMax = bigArrayMax (bb->aligns) ;
  vTXT cigar = vtxtHandleCreate (h) ;
  Array cigarettes = arrayHandleCreate (1024, SAMCIGAR, h) ;
  char *VERSION = "1.1" ;
  
  aceOutf (ao, "@HD VN:1.5\tSO:queryname\n") ;
  aceOutf (ao, "@PG ID:1\tPN:Magic\tVN:%s", VERSION) ;
  /* aceOutf (ao, "\tCL:%s", commandBuf) ; */
  aceOut (ao, "\n") ;

  for (ii = 0 ; ii < aMax ; ii++, ap++)
    exportOneSam (ao, pp, bb, cigar, cigarettes, ap) ;
   
  ac_free (h) ;
return ;
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
  ACEOUT ao ;

  if (!pp->sam)
    {
    }
  else
    {
      ao = aceOutCreate (pp->outFileName, ".sam", pp->gzo, h) ;
    }
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->weChan, &bb, BB))
    {
      t1 = clock () ;
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  bigArraySort (bb.aligns, alignOrder) ;
	  if (! pp->sam)
	    exportDo (pp, &bb) ;
	  else
	    exportSamDo (ao, pp, &bb) ;

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

/*************************************************************************************/

static void reportRunStats (PP *pp, Array runStats)
  {
    Array aa = runStats ;
    int ii, iMax = arrayMax (aa) ;
    RunSTAT *s0 = arrayp (aa, 0, RunSTAT) ;
	
    printf ("\n####### Run Statistics") ;

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
    
    printf ("\nLn1") ;
    for (ii = 0 ; ii < iMax ; ii++) 
      {	
	RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	if (! s->nReads)
	  continue ;
	float nR1 = s->nPairs ? s->nPairs : s->nReads ;
	printf ("\t%.1f", s->nBase1/nR1) ;
      }
    
    printf ("\nLn1") ;
    for (ii = 0 ; ii < iMax ; ii++) 
      {	
	RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	if (! s->nReads)
	  continue ;
	float nR2 = s->nPairs ? s->nPairs : 1 ;
	printf ("\t%.1f", s->nBase2/nR2) ;
      }
    
    printf ("\nReads aligned") ;
    for (ii = 0 ; ii < iMax ; ii++) 
      {	
	RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	if (! s->nReads)
	  continue ;
	printf ("\t%ld", s->nMultiAligned[0]) ;
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
    printf ("\n") ;
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
	pp->runName = "f.1" ;
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
	  
	  cr = filName (cp, 0, "r") ;
	  if (! cr)
	    messcrash ("\nCannot open input file %s\n", cp) ;
	  
	  rc = arrayp (rcs, nn++, RC) ; 
	  rc->fileName1 = strnew (cr, pp->h) ;
	  rc->fileName2 = filName2 ;
	  rc->run = run ;

	  rc->format = RAW ; /* default */
	  if (strstr (rc->fileName1, ".fasta")) rc->format = FASTA ;
	  if (strstr (rc->fileName1, ".fna")) rc->format = FASTA ;
	  if (strstr (rc->fileName1, ".fa")) rc->format = FASTA ;
	  if (strstr (rc->fileName1, ".fastq")) rc->format = FASTQ ;
	  if (strstr (rc->fileName1, ".fastc")) rc->format = FASTC ;

	  /* user can override the defaults */
	  if (pp->raw) rc->format = RAW ;
	  if (pp->fasta) rc->format = FASTA ;
	  if (pp->fastq) rc->format = FASTQ ;
	  if (pp->fastc) rc->format = FASTC ;

	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    {
	      rc->pairedEnd = TRUE ;
	      rs->nFiles++ ;
	    }
	  
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
	  cq = strchr (cp, ',') ;
	  if (!cq) cq = strchr (cp, '+') ;
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
	  cr = filName (cp, 0, "r") ;
	  if (! cr)
	    messcrash ("\nCannot open input file %s\n", cp) ;
	  rc->fileName1 = strnew (cr, pp->h) ;
	  if (! dictAdd (fDict, cr, 0))
	    messcrash ("\nDuplicate target file name %s\n at line %d of file -T %s\n try sortalign --help\n"
		       , cr
		       , line
		       , pp->inConfigFileName
		       ) ;

	  /* run name */
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

	      if (! strcasecmp (cp, "rna")) rc->RNA = TRUE ;
	      if (! strcasecmp (cp, "dna")) rc->RNA = FALSE ;
	      
	      cp = cq ;
	    }

	  if (! rc->format)
	    {
	      rc->format = RAW ; /* default */
	      if (strstr (rc->fileName1, ".fasta")) rc->format = FASTA ;
	      if (strstr (rc->fileName1, ".fa")) rc->format = FASTA ;
	      if (strstr (rc->fileName1, ".fastq")) rc->format = FASTQ ;
	      if (strstr (rc->fileName1, ".fastc")) rc->format = FASTC ;
	    }
	 
	  RunSTAT *rs = arrayp (runStats, run, RunSTAT) ;
	  rs->nFiles++ ;
	  if (rc->fileName2)
	    rs->nFiles++ ;
	} 
      
      ac_free (h) ;
    }
  printf ("Found %d sequence files\n", arrayMax (rcs)) ;
  return rcs ;
} /* parseInConfig */    

/*************************************************************************************/
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
	  if (!(cc >= 'A' && cc <= 'Z') && !(cc >= 'A' && cc <= 'Z')) cc = 0 ;
	  if (cp[1]) cc = 0 ;
	  if (! cc)
	    messcrash ("\n\nThe target class must be specified asa single character (a-z, A-Z)  at line %d of -T target config file %s\n try sortalign --help\n"
		       , *cp
		       , line
		       , pp->tConfigFileName
		       ) ;
	  tc = arrayp (tcs, nn++, TC) ;
	  tc->targetClass = cc ;
	  
	  /* file names */
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
	      if (strstr (cp, ".introns"))
		tc->format = INTRON ;
	      else
		messcrash ("\n\nThe Introns must be specified via a .gtf file named *.intronss (not ,%s), with 3 columns tab delimited: chromosome,a1,a2. If (as happens most frequently) the intron is GT_AG, a1 and a2 are the (1-based) coordinates of these 2 G. Thus a1<a2 on the top strand, a1>a2 on the bottom strand. Line %d of -T target config file %s\n try sortalign --help\n"
			   , cp
			   , line
			   , pp->tConfigFileName
		       ) ;
	    }
	  /* options */
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
	       "//      try: -h --help\n"
	       "// EXAMPLES:\n"
	       "//      sortalign --createIndex XYZ -t target.fasta (needed once)\n"
	       "//      sortalign --index XYZ -i f.fastq.gz --wiggle -o results/xxx \n"
	       "// OBJECTIVE:\n"
	       "//    Sort align maps deep-sequencing files and can export at once alignments, coverage plots, introns, SNPs.\n"
	       "//      On the first pass, sortalign analyses the (-t or -T) target(s) and creates a directory of binary files.\n"
	       "//      On subsequent calls, sortalign memory maps the indexes and analyses the (-i or -I) sequence files.\n"
	       "//      All output files use the -o parameter (we recommend to give the name of a directory) as prefix.\n"
	       "//    The amount of data to be analysed can be large (say 100 Gbases), as it is processed in batches.\n"
	       "//      Depending on the  hardware, using 10 to 20 CPUs, it may process about 1 Gigabases of human RNA-seq per minute\n"
	       "//      On human or mouse, the program requires around (20 + nB) gigabytes of RAM (see --nB --nA options).\n"
	       "// TARGET CREATION:\n"
	       "// --createIndex <directory_name>\n"
	       "// --seedLength <int> : default 16, min 10, max 19\n"
	       "// --maxTargetRepeats <int> : default[12], do not index highly repeated target words\n"
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


	       "//         I: Introns, in three columns the coordinates on this genome of the known introns:\n"
	       "//               Chromosome_name  i1 i2 : i1 and i2 give the position of the 2 G of the Gt_aG motif.\n"
	       "//               On the top strand i1 < i2, on the bottom strand i1 > i2, use 1 for the first base.\n"
	       "//         T: Transcripts and gene coordinates, used to compute gene expression.\n"
	       "//         B: Frequent bacteria contaminants fasta file.\n"
	       "//         V: Frequent virus contaminants fasta file.\n"
	       "//      Column 2: f1[,f2] : one or several coma separated file names.\n"
	       "//      Column 3: Coma separated options, no blank/spaces/tabs please\n"
	       "//           fast[default],raw,bonus=<int>  (the bonus is added to the alignment score)\n" 
	       "// SEQUENCE FILES  TO BE ANALYZED:\n"
	       "// -x, --index <directory_name> : a directory of index files created  previously using --createIndex\n" 
	       "// -r, --run <runName>   [default -]  : global run name\n"
	       "//    the run name may also be attached to the file names as shown now.\n"
	       "// -i <sequence_fileName[s]> \n"
	       "//     example:   -i f0,r1:f1,f2,f3,r2:f4,f5\n"
	       "//     Align files f0 to f5, optionally the run name changes to r1, r2 when : is found\n"


	       "//     The file format is implied by the file name, or may be provided explicitly using\n"
	       "//       --raw   | --fasta  | --fastq | --fastc\n"
	       "//     Use '-i - ' to pipe sequence files into the pipeline, --gzi to for unzipping\n"
	       "//        example: zcat fx.*.fastq.gz | sortalign -x XYZ -i - --gzi --fastq -o results/fx\n"
	       "// -I <config_fileName>\n"
	       "//     A more precise definition of a set of sequencing files to be analysed\n"
	       "//     Eeach line contains 1 to 3 tab separated columns\n"
	       "//          FileName[s]  RunName Descriptors\n"
	       "//     Example:\n"
	       "//          f1.fasta run1 RNA,Nanopore\n"
	       "//          f.R1.gz,f.R2.gz  run2  DNA,fastq,Illumina\n"
	       "//        A table with a single column of sequence files is acceptable, the format will be deduced from the file names\n"
	       "//     1: file name, mandatory\n"
	       "//        For paired end sequencing, provide, as in the second example, two coma separated file names\n"
	       "//        A file named (.fa, .fasta, .fastq, .fastc) implies that format, the file may be gzipped (.gz)\n"
	       "//     2: Run name, optional,the name must not contain spaces or special characters, as it will be used to name subdirectories\n"
	       "//     3: Descriptors, optional, the options are coma separated, possibilities are\n"
	       "//        DNA/RNA : [default RNA], if DNA sortalign will not clip polyA or jump introns\n"
	       "//        Machine : [default Illumina], one of Illumina, PacBio, Nanopore,..\n"
	       "//        Adaptors=atagg,cctg   run specific adaptor sequences\n"
	       "//        Format : [default raw], one of raw, fasta, fastq, fastc, SRA, only needed if not implied by the file name\n"
	       "// --gzi\n"
	       "//    Forces decompression of the input files, useful when piping into sortalign\n"
	       "//    All files named *.gz are automatically decompressed\n"
	       "// OUTPUT:\n"
	       "// -o <outFileNamePrefix> [default stdout]\n"
	       "//   All output files will share this prefix, large outputs are split\n"
	       "// --oMax  <int> : max number of lines per output file, [default 100M]\n"
	       "// --gzo : the output files will be gziped\n"      
	       "// REQUEST\n"
	       "// --align : [default] Extend the seed alignments hopefully to full read length\n"
	       "// --do_not_align : do not align, just test the validity of the index directory\n"
	       "//   --minAli <int> : minimal length of aech aligned fragmant (exon)\n"
	       "//   --errMax <int> : [default NA] maximal number of mismatches in any (partial) alignment\n"
	       "//   --errRateMax <int> : [default 10] maximal percentage of mismatches in any (partial) alignment\n"
	       "// --wiggle  : Report target coverage wiggles in UCSC BF (fixed) format\n"
	       "// --intron  : (not yet ready) Report intron support\n"
	       "// --snp : (not yet ready) Report candidate SNP counts (substitutions and short indels)\n"
	       "// STEPPING\n"
	       "// --step <int>, take a seed every <int> base\n"
	       "//   while createIndex, the default is 1 for targets < 1Mb, 3 for larger targets, 5 would accelerate the code\n"
	       "//   while analysing the reads, the default is 2. --step 1 would increase sensitivity\n"
	       "//   The programs checks that the read and target steppinngs are realtive prime\n"
	       "// OPTIONAL TECHNICAL PARAMETERS\n"
	       "//   The program sortalign is parallelized and synchronized by GO-language like channels\n"
	       "//   All code layers execute at the same time, the work-load among layers is self balancing\n"
	       "//   There is no limit to the size of the input, as data continuously flow in and out\n"
	       "// --nA or --nAgents <int> : [default 10] number of agent in each code layer\n"
	       "// --nB or --nBlocks <int> : [default 20] number of simultaneous data blocks circulating in the pipeline\n"
	       "// --NN <int> : [default 16] split the seed files in NN parts, allowed values (1,2,4,8,16,32,64)\n"
	       "//              seedLength > 16 imply NN >= 4^(sedLength - 16)\n"
	       "// --max_threads <int>  : [default 128] maximal number of simultaneous UNIX threads\n"
	       "//    Possible values are \n"
	       "//        --nAgents 1  --nBlocks 1 --max_threads 8 : for a small test\n"
	       "//        --nAgents 10 --nBlocks 10 --max_threads 128 : default, uses less than 32 Gbytes of RAM\n"
	       "//        --nAgents 20 --nBlocks 30 --max_threads 512 : on a large machine\n"
	       "//    These parameters do not limit the ammount of data to be processed, they affect the speed\n"
	       "//    If no message is emited within the first minute, the computer is overloaded,\n"
	       "//    please try to increase --max_threads to 256 or 512. alternativelly lower --nAgents and --nBlocks\n"
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
  int maxThreads = 128 ;
  long unsigned int nHits = 0, nSeqs = 0, nBaseAligned = 0 ;
  long unsigned int nerr = 0 ;   /* cumulated number of errors */
  long unsigned int aliDx = 0 ; /* cumulated aligned read length */
  long unsigned int aliDa = 0 ;  /* cumulated genome coverage */
  long unsigned int nAli = 0 ;
  BOOL debug = FALSE ;
  int n = 0 ;
  long int skips0 = 0, skips1 = 0, skips2 = 0, skips3 = 0, skips4 = 0, skipsFound = 0, skipsNotFound = 0 ;
  char tBuf[25] ;
  char tBuf0[25] ;
  AC_HANDLE h ;
  int nAgents = 10 ;
  int channelDepth = 1 ;
  mytime_t t0, t1 ;
  
  freeinit () ; 
  messErrorInit (argv[0]) ;

  memset (&p, 0, sizeof (PP)) ;
  h = ac_new_handle () ;
  p.h = h ;

  if (argc < 2)
    usage (0, 0, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (0, 0, argv) ;


  p.debug = getCmdLineOption (&argc, argv, "--debug", 0) ;
  p.debug |= getCmdLineOption (&argc, argv, "--verbose", 0) ;


  /**************************  debugging tools, ignore *********************************/

  if ( getCmdLineInt (&argc, argv, "-n", &n))
    {
      BigArray b, a = bigArrayHandleCreate (n, HIT, 0) ;

      
      for (int i = 0 ; i < n ; i++)
	bigArray (a, i, HIT).read = ( !(n &0x1) ? (n - i)%7 :  (unsigned int)randint()) ; 
      b = bigArrayHandleCopy (a, 0) ;
      newMsort (b, hitReadOrder) ;

      for (int i = 0 ; i < 12 && i < n ; i++)
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
  p.align = getCmdLineBool (&argc, argv, "--align") ;          /* over-ridden by --do_not_align */
  p.align = ! getCmdLineBool (&argc, argv, "--do_not_align") ; /* default is to align */
  
  p.wiggle = getCmdLineBool (&argc, argv, "--wiggle") ;


  /* future options ***
     p.snps = getCmdLineBool (&argc, argv, "--snp") ;
     p.introns = getCmdLineBool (&argc, argv, "--intron") ;
  */
  
  /*****************  file names and their formats  ************************/
  
  if (! (p.createIndex = getCmdLineOption (&argc, argv, "--createIndex", &(p.indexName))) &&
      ! getCmdLineOption (&argc, argv, "-x", &(p.indexName))
      )
    getCmdLineOption (&argc, argv, "--index", &(p.indexName)) ;
  getCmdLineOption (&argc, argv, "-i", &(p.inFileName)) ;
  getCmdLineOption (&argc, argv, "-I", &(p.inConfigFileName)) ;
  getCmdLineOption (&argc, argv, "-t", &(p.tFileName)) ;
  getCmdLineOption (&argc, argv, "-T", &(p.tConfigFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(p.outFileName)) ;
  p.runName = 0  ; /* default */

  /* formats are only used in -t case, in -T case provide the format in the 3rd columnh of the tConfig file */
  p.fasta = getCmdLineBool (&argc, argv, "--fasta");
  p.fastq = getCmdLineBool (&argc, argv, "--fastq");
  p.fastc = getCmdLineBool (&argc, argv, "--fastc");
  p.raw = getCmdLineBool (&argc, argv, "--raw");

  p.gzi = getCmdLineBool (&argc, argv, "-gzi") ||
    getCmdLineBool (&argc, argv, "--gzi");
  p.gzo = getCmdLineBool (&argc, argv, "-gzo") ||
    getCmdLineBool (&argc, argv, "--gzo");
  p.sam = getCmdLineBool (&argc, argv, "--sam") ;
  
  /***************** run name, only use in conjunction with -t, not use with -T ***********/

  p.run = 0 ;
  getCmdLineOption (&argc, argv, "-r", &(p.runName)) ;
  getCmdLineOption (&argc, argv, "--run", &(p.runName)) ;

  /*****************  seed length, only used when createIndex  ************************/

  NN = 16 ; /* default */
  getCmdLineInt (&argc, argv, "--NN", &(NN));
  p.seedLength = 16 ; /* default */
  getCmdLineInt (&argc, argv, "--seedLength", &(p.seedLength));
#ifdef METHOD_AL
  if (p.seedLength > 15)
    p.seedLength = 15 ; /* default */
#endif
  
  p.maxTargetRepeats = 12 ;
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
  nAgents = 10 ;
  if (! getCmdLineInt (&argc, argv, "--nAgents", &(nAgents)))
    getCmdLineInt (&argc, argv, "--nA", &(nAgents)) ;

  p.nBlocks = 20 ;  /* max number of BB blocks processed in parallel */
  if (!getCmdLineInt (&argc, argv, "--nBlocks", &(p.nBlocks)))
    getCmdLineInt (&argc, argv, "--nB", &(p.nBlocks));
  if (p.nBlocks == 1)
    { nAgents = 1 ; }

  maxThreads = 128 ;  /* UNIX  max om lmem12 machine */
  getCmdLineInt (&argc, argv, "--max_threads", &maxThreads) ;
  if (maxThreads < 24)
    maxThreads = 24 ;

  if (p.createIndex)
    { nAgents = 1 ; maxThreads = 1 ; p.nBlocks = 1 ; }

  p.iStep = 2 ;   /* read default */
  p.tStep = 0 ;   /* default 1 or 3 set in createIndex or read in existing index */
  if (p.createIndex)
    getCmdLineInt (&argc, argv, "--step", &(p.tStep)) ;
  else
    getCmdLineInt (&argc, argv, "--step", &(p.iStep)) ;
  
  /*****************  Aligner filters  ************************/
  p.minAli = 30 ;
  p.errMax = 10 ;
  p.errRateMax = 10 ;
  p.OVLN = 30 ;
  getCmdLineInt (&argc, argv, "--errMax", &(p.errMax)) ;
  getCmdLineInt (&argc, argv, "--minAli", &(p.minAli)) ;
  getCmdLineInt (&argc, argv, "--errRatMax", &(p.errRateMax)) ;

  /****************** Check the existence of all file names ****************************************/ 

  /* check that the index directory is accessible */
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
  
  /*****************  Check tat all parameters have been parsed *******************/
  if (argc > 1)
    usage (0, argc, argv) ;

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
      createTargetIndex (&p, &p.bbG, tArray) ;
      if (p.tConfigFileName)
	system (hprintf(h, "\\cp %s %s\n", p.tConfigFileName, p.indexName)) ;
      ACEOUT ao = aceOutCreate (filName (p.indexName, "/seedLength", "w") , 0, 0, p.h) ;
      aceOutf (ao, "%d\n", p.seedLength) ;
      ac_free (ao) ;
      
      goto done ;
    }

  /*******************  otherwise verify the existence of the indexes ********************/
  
  /* check that input files were provided */
  if (! p.indexName)
    usage ("missing parameter -x or --createIndex", argc, argv) ;
  if (! p.inFileName && ! p.inConfigFileName)
     usage ("missing parameter -i or -I inputFileName(s)", argc, argv) ;
  if (p.inFileName && p.inConfigFileName)
     usage ("conflicting parameter -i and -I, both define the input files", argc, argv) ;

  if (1)
    {
      BOOL ok = TRUE ;
      int k = 0 ;
      char *cp ;
      NN = 0 ;
      while (1)
	{
#ifdef METHOD_AL
	  char *fNam = hprintf (h, "/cws.sortali.A.%d", k) ;
#else
	  char *fNam = hprintf (h, "/cws.sortali.%d", k) ;
#endif
	  k++ ;
	  cp = filName (p.indexName, fNam, "rb") ;
	  if (cp)
	    {
	      NN++;
	      if (NN == 1)
		{
		  char *cq = cp + strlen(cp) - 2 ;
		  *cq = 0 ;
#ifdef METHOD_AL
		  p.tFileBinaryCwsAName = strnew (cp, p.h) ;
		  cq = cp + strlen(cp) - 1 ;
		  *cq = 'L' ;
		  p.tFileBinaryCwsLName = strnew (cp, p.h) ;
#else
		  p.tFileBinaryCwsName = strnew (cp, p.h) ;
#endif
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
	  int x = 0 ;
	  if (aceInCard (ai) && aceInInt (ai, &x))
	    p.seedLength = x ;
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
    }

  /* check the existence of the input sequence files */
  Array inArray = parseInConfig (&p, runStats) ;

  /***************************/
  /* launch the multiprocessing */

  wego_max_threads (maxThreads) ;
  
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
  p.csChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.csChan, debug, "csChan") ;
  p.smChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.smChan, debug, "smChan") ;
  p.moChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.moChan, debug, "moChan") ;
  p.oaChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.oaChan, debug, "oaChan") ;
  p.awChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.awChan, debug, "awChan") ;
  p.weChan = channelCreate (channelDepth, BB, p.h) ;
  channelDebug (p.weChan, debug, "weChan") ;
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
  /* Read preprocessing agents, they do not require the genome */
  for (int i = 0 ; i < nAgents && i < 5 ; i++)
    {
      p.agent = i ;
      wego_go (readParser, &p, PP) ;
    }
  for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
    {
      p.agent = i ;
      wego_go (codeWords, &p, PP) ;
      wego_go (sortWords, &p, PP) ;
    }
  
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
#ifdef METHOD_AL
  if (! p.bbG.cwsLN[0])
    messcrash ("matchHits received no target lines") ;
#else
  if (! p.bbG.cwsN[0])
    messcrash ("matchHits received no target words") ;
#endif
  
  /* map the reads to the genome in parallel */
  for (int i = 0 ; i < nAgents && i < p.nBlocks ; i++)
    {
      p.agent = i ;
      wego_go (matchHits, &p, PP) ;
      wego_go (orderHits, &p, PP) ;
      if (!i || p.align) /* at least 1 agent */
	{
	  wego_go (align, &p, PP) ;
	  wego_go (align, &p, PP) ;
	}
      if (!i) /* only 1 wiggle agent */
	wego_go (wiggle, &p, PP) ;
      if (!i) /* only 1 export agent */
	wego_go (export, &p, PP) ;
    }
  
  while (channelGet (p.doneChan, &bb, BB))
    {
      long int n = (bb.hits ? bigArrayMax (bb.hits) : 0) ;

      if (bb.isGenome && (p.bbG.cwsN || p.bbG.cwsAN))
	{
	  cpuStatCumulate (cpuStats, p.bbG.cpuStats) ;
	  for (int k = 0 ; k < NN ; k++)
	    {
	      if (p.bbG.cwsN) bigArrayDestroy (p.bbG.cwsN[k]) ;
	      if (p.bbG.cwsAN) bigArrayDestroy (p.bbG.cwsAN[k]) ;
	      if (p.bbG.cwsLN) bigArrayDestroy (p.bbG.cwsLN[k]) ;
	    }
	  ac_free (p.bbG.cwsN) ;
	  ac_free (p.bbG.cwsAN) ;
	  ac_free (p.bbG.cwsLN) ;
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

      if (bb.run)
	{
	  runStatsCumulate (0, runStats, &(bb.runStat)) ;
	  runStatsCumulate (bb.run, runStats, &(bb.runStat)) ;

	  runErrorsCumulate (0, runErrors, bb.errors) ;
	  runErrorsCumulate (bb.run, runErrors, bb.errors) ;
	}
      
      /* release block  memory */
      if (bb.dnas)
	{
	  Array dna = 0 ;
	  int iMax = arrayMax (bb.dnas) ;
	  for (int i = 1 ; i < iMax ; i++)
	    {  /* the dna array points into the memory globalDna array */
	      dna = arr (bb.dnas, i, Array) ;
	      if (dna->lock)
		{
		  dna->base = 0 ;
		  arrayUnlock (dna) ;
		}
	    }
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
  if (1 || p.debug) printf ("Skips: 0=%ld, 256=%ld, 1024=%ld, 4096=%ld, 16384=%ld, found=%ld, notFound=%ld\n",
			    skips0, skips1, skips2, skips3, skips4, skipsFound, skipsNotFound);
  if (1 || p.debug) printf ("SeedLength %d, nAgents=%d nBlocks=%d NN=%d\n",
			    p.seedLength, nAgents, p.nBlocks, NN) ;
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
 
