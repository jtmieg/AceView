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


#include "ac.h"
#include "channel.h"
#define WIGGLE_STEP 10

typedef struct nodeStruct { double x ; CHAN *cx, *cy, *cu, *cv, *done ; int k ; } NODE ;

typedef enum {FASTA=1, FASTQ, FASTC, RAW, SRA} DnaFormat ;
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
  long int nErr ;
  Array errors ;  /* substitutions, insertions, deletions counts */
  /* coverage of long transcripts ? */
} RunSTAT ;
		  
typedef struct bStruct {
  AC_HANDLE h ;
  int run ;
  RC rc ;
  DICT *dict ;
  BigArray dnaCoords ;   /* offSets of the dna in the globalDna array */
  Array dnas ;           /* Array of const char Arrays */
  Array dnasR ;          /* Their reverse complement, only computed for the genome */
  BigArray globalDna ;   /* concatenation of all sequence DNAs separated by blocks of nnn */
  BigArray globalDnaR ;  /* concatenation of all reverse sequences, in the same order */
  BigArray cws ;         /* BigArray of codeWords */
  BigArray hits ;        /* BigArray of read<->genome hits */
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
  int step, skips0, skips1, skips2, skips3, skips4 ;
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
  BOOL fasta, fastq, fastc, raw ;
  BOOL hasBonus ;
  int bonus[256] ;
  DICT *runDict ;
  DICT *targetClassDict ;
  int run ;
  int nFiles ;  /* number of input sequence files */
  int agent ;  /* instance of the agent */
  int maxBB ;
  int blocMaxBases ; /* max number of bases read as one bloc */
  int step ;         /* take a word every step (read default 2, target default 3) */
  int errMax ;       /* (--align case) max number of errors in seed extension */
  int minAli ;
  int errRateMax ;       /* (--align case) max number of errors in seed extension */ 
} PP ;

typedef struct codeWordsStruct {
  unsigned int word ; /* 32 bits = 16 bases, 2 bits per base */
  int nam ; /* index in readDict or chromDict */
  int pos ;  /* position, negative on minus strand */
  unsigned int targetClass ;
} __attribute__((aligned(16))) CW ;

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
  int score ;
  int nerr ;
  int readLength ;
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
/* Dedicated radixSort on cw->word */
static void cwRadixSort (BigArray cws)
{
  long int N = bigArrayMax (cws) ;
  CW *cw = bigArrp (cws, 0, CW) ;
  CW *temp = malloc(N * sizeof(CW)) ;
  if (!temp)
    messcrash ("cwRadixSort: temp malloc failure\n") ;
  unsigned int *keys = malloc(N * sizeof(unsigned int)) ;
  if (!keys)
    messcrash ("cwEadixSort: keys malloc failure\n") ;

  /* Extract keys */
  for (size_t i = 0 ; i < N ; i++)
    keys[i] = cw[i].word ;

  /* 8-bit radix sort, 4 passes + 1 for counts */
  unsigned int *temp_keys = malloc(N * sizeof(unsigned int)) ;
  if (!temp_keys)
    messcrash ("radixSort: temp_keys malloc failure\n") ;
  
  for (int shift = 0 ; shift < 32 ; shift += 8)
    {
      size_t count[256] = {0} ;
      // Count occurrences
      for (size_t i = 0 ; i < N ; i++) 
	count[(keys[i] >> shift) & 0xFF]++ ;
      
      // Compute offsets
      for (int i = 1 ; i < 256 ; i++) 
	count[i] += count[i - 1] ;
      
      // Scatter
      for (size_t i = N ; i > 0 ; i--)
	{
	  size_t idx = --count[(keys[i - 1] >> shift) & 0xFF] ;
	  temp_keys[idx] = keys[i - 1] ;
	  temp[idx] = cw[i - 1] ;
        }
      // Swap buffers
      unsigned int *swap = keys ;
      keys = temp_keys ;
      temp_keys = swap ;
      memcpy(cw, temp, N * sizeof(CW)) ;
    }
  
  free(temp_keys) ;
  free(keys) ;
  free(temp) ;
} /* cwRadixSort */

/**************************************************************/
/* Dedicated hitSort on all 4 fields */
static void hitRadixSort (BigArray hits)
{
  long int N = bigArrayMax (hits)  ;
  HIT *hp = bigArrp (hits, 0, HIT) ;
  
  if (N <= 1) return ;

  // Allocate temporary arrays
  HIT *temp = malloc(N * sizeof(HIT)) ;
  if (!temp) 
    messcrash ("hitRadixSort: temp malloc failure\n") ;
  
  unsigned int *keys = malloc(N * 4 * sizeof(unsigned int)) ;
  if (!keys)
    messcrash ("hitRadixSort: keys malloc failure\n") ;



  unsigned int *temp_keys = malloc(N * 4 * sizeof(unsigned int));
  if (!temp_keys) 
    messcrash ("hitRadixSort: temp_keys malloc failure\n") ;


  // Extract keys
  for (int i = 0; i < N; i++) {
    keys[i * 4 + 0] = hp[i].x1;    // x1
    keys[i * 4 + 1] = hp[i].a1;    // a1
    keys[i * 4 + 2] = hp[i].chrom; // chrom
    keys[i * 4 + 3] = hp[i].read;  // read
  }
  
  HIT *src = hp;
  HIT *dst = temp;
  
  // Process fields in reverse order: x1, a1, chrom, read
  for (int field = 0; field < 4; field++) {
    for (int shift = 0; shift < 32; shift += 8) {
      size_t count[256] = {0};
      
      // Count occurrences
      for (int i = 0; i < N; i++) {
	unsigned int val = keys[i * 4 + field];
	count[(val >> shift) & 0xFF]++;
      }
      
      // Compute offsets
      for (int i = 1; i < 256; i++) {
	count[i] += count[i - 1];
      }
      
      // Scatter
      for (size_t i = N; i > 0; i--) {
	unsigned int val = keys[(i - 1) * 4 + field];
	size_t idx = --count[(val >> shift) & 0xFF];
	dst[idx] = src[i - 1];
	temp_keys[idx * 4 + 0] = keys[(i - 1) * 4 + 0];
	temp_keys[idx * 4 + 1] = keys[(i - 1) * 4 + 1];
	temp_keys[idx * 4 + 2] = keys[(i - 1) * 4 + 2];
	temp_keys[idx * 4 + 3] = keys[(i - 1) * 4 + 3];
      }
      
      // Swap buffers
      HIT *swap = src;
      src = dst;
      dst = swap;
      unsigned int *key_swap = keys;
      keys = temp_keys;
      temp_keys = key_swap;
    }
  }
  
  // Ensure final result is in hp
  if (src != hp) 
    memcpy(hp, src, N * sizeof(HIT));
  
  free(temp_keys);
  free(keys);
  free(temp);
} /* hitRadixSort */

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
  while (n1 > 0 && n2 > 0)
    {
      n = (*cmp) (up, vp) ;
      
      if (n <= 0)
	{ n1-- ; memcpy (wp, up, s) ; wp += s ; up += s ; }
      else
	{ n2-- ; memcpy (wp, vp, s) ; wp += s ; vp += s ; }
    }

  while (n1 > 0 && n2 > 0)
    {
      n = (*cmp) (up, vp) ;
      if (n <= 0)
	{ n1-- ; memcpy (wp, up, s) ; wp += s ; up += s ; }
      else
	{ n2-- ; memcpy (wp, vp, s) ; wp += s ; vp += s ; }
    }

  #ifdef JUNK
  /* more complex equivalent code
   * we do less memcpy of larger size, on average k=2
   * but this code is more complex and not faster
   * i keep just to avoid trying the same idea again later
   */

  n = (*cmp) (up, vp) ;
  while (n1 > 0 && n2 > 0)
    {
      if (n <= 0)
	{
	  int k = 1 ; n1-- ;
	  while (n1 > 0 && n2 > 0 &&  (*cmp) (up + k, vp) <= 0)
	    { k++ ; n1-- ;}
	  memcpy (wp, up, k * s) ;
	  wp += k * s ; up += k * s ;
	  n = 1 ;
	}
      else
	{
	  int k = 1 ; n2-- ;
	  while (n1 > 0 && n2 > 0 &&  (*cmp) (up, vp + k) > 0)
	    { k++ ; n2-- ;}
	  memcpy (wp, vp, k * s) ;
	  wp += k * s ; vp += k * s ;
	  n = -1 ;
	}
    }
#endif
  /* bulk copy the reminders */
  if (n1 > 0)
    {
      memcpy (wp, up, n1 * s) ;
      wp += n1 ;
    }
  if (n2 > 0)
    {
      memcpy (wp, vp, n2 * s) ;
      wp += n2 ;
    }
} /* newMsortDo */

static void newMsort (BigArray aa, int (*cmp)(const void *va, const void *vb))
{
  long int N = bigArrayMax (aa) ;
  char *cp = (char *) bigArrp (aa, 0, HIT) ;
  int s = aa->size ;
  
  if (N <= 1) return;
  /* mescrash ("this sort might be faster but is in development) ; */
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
  int NMAX = 50000 ;  /* was 100000 for 300s on 6.5 GigaBase human 2025_05_23 */
  int BMAX = 200 * NMAX ;
  int nPuts = 0 ;
  int nn = 0, nSeqs = 0 ;
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
      memset (bb, 0, sizeof (BB)) ;
      chan = pp->plChan ;
      namBufX = namBuf ;
    }

  if (! bb->h)
    {  /* isGenome: keep expanding the same bbG if already initialised */
      bb->h = ac_new_handle () ;
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
  dna1 = arrayHandleCreate (isGenome ? 1 << 30 : 256, unsigned char, bb->h) ;
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
	  bb->h = ac_new_handle () ;
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

static BigArray GenomeAddSkips (const PP *pp, BigArray cws, AC_HANDLE h)
{
  long int iMax ; 
  long int jMax ; 
  long int i, j, k ;
  
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
	  if (n > 20) n = 20 ;
	  ks[n]++ ;
	  if (n < 12)
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

/**************************************************************/
/**************************************************************/

static void loadRegulator (const void *vp)
{
  BB bb ;
  const PP *pp = vp ;
  int nn = 0, nMax = pp->maxBB ;
  
  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->plChan, &bb, BB))
    {
      nn++ ;
      while  (nn > nMax)
	{
	  int nd = channelCount (pp->doneChan) ;
	  nMax = pp->maxBB - nd ;
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

static BigArray codeWordsDo (const PP *pp, BB *bb, int step, BOOL isTarget) 
{
  BigArray cws = 0 ;
  CW *restrict cw0 ;
  CW *restrict cw ;  
  const unsigned char *restrict cp ;
  Array dnas = bb->dnas ;
  int ia, iaMax = arrayMax (dnas) ;
  long int kMax = bigArrayMax (bb->globalDna) / step ;
			     
  if (step < 1)
    messcrash ("codeWordsDo received step = %d < 1", step) ;
  if (kMax >> 31)
    messcrash ("codeWordsDo received ln=%ls, step=%d, ln/step=%ld > 2G", bb->length, step, kMax) ;
    
  cws = bigArrayHandleCreate (kMax, CW, bb->h) ;
  cw = bigArrayp (cws, kMax - 1, CW) ; /* make room */
  cw = cw0 = bigArrayp (cws, 0, CW) ; /* initialize */

  bb->step = step ;
  for (ia = 1 ; ia < iaMax ; ia++)
    {
      Array dna = arr (dnas, ia, Array) ;
      if (! dna)
	continue ;
      int ii, jj, p = 0, cStep ;
      int iMax = arrayMax (dna) ;
      unsigned int w, wr, w1, w2, wr1 ;
      unsigned int targetClass = isTarget ? *(dictName (bb->dict, ia)) : 0 ;
      step = iMax < 30 ? 1 : bb->step ;
      w = w1 = w2 = wr = wr1 = 0 ;
      cp = arrp (dna, 0, unsigned char) ;
      for (ii = 0, jj = -16, cStep = -16, p = 0 ; ii < iMax + 1 ; ii++, jj++, cp++, cStep++) 
	{
	  /* construct a 32 bits integer representing 16 bases (A_, T_, G_, C_) using 2 bits per base */
	  w2 = w1 ; w1 = w ; wr1 = wr ;
	  w <<= 2 ; wr >>= 2 ;
	  p++ ;
	  switch ((int)*cp)
	    {  /* alphabetic order and XOR complement */
	    case 0: p = 0 ; break ;
	    case A_:            wr |= 0x3 << 30 ; break ;
	    case C_: w |= 0x1 ; wr |= 0x2 << 30 ; break ;
	    case G_: w |= 0x2 ; wr |= 0x1 << 30 ; break ;
	    case T_: w |= 0x3 ;                   break ;
	      
	    default: p = 0 ;
	    }
	  if (cStep == step) cStep = 0 ;
	  if (p >= 17 && !cStep && w1 && wr1 && w != w1 && w != w2 && w1 != w2)
	    {
	      char minus  = (w1 > wr1 ? 0x1 : 0x0) ;
	      cw->nam = ia ;
	      cw->word = ( minus ? wr1 : w1 ) ;
	      cw->pos = 1 + jj ; /* add 1 to avoid 0 = -0 */
	      if (minus) cw->pos = - cw->pos ;
	      /* x = minus ? 15 - j : j - 15 ;  signed bio coordinates of base 1 of the w1 seed */
	      cw->targetClass = targetClass ;
	      cw++ ;		    
	    }
	}
    }
  if (cw - cw0 >= kMax)
    messcrash ("exceeded cws allocated array length, pllease edit the source code") ;
  bigArrayMax (cws) = cw - cw0 ;
  if (0) showCws (cws) ;
  
  return cws ;
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

      bb.cws = codeWordsDo (pp, &bb, pp->step, FALSE) ;

      t2 = clock () ;
      cpuStatRegister ("3.CodeWords", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.cws)) ;
      if (pp->debug) printf ("--- %s: Stop code words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cws)) ;

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

static BOOL genomeParseBinary (const PP *pp, BB *bbG)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  BigArray seqIds = 0;
  DICT *dict = 0 ;
  long int ii, iMax ;
  
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
  fNam = pp->tFileBinaryCwsName ;
  bbG->cws = bigArrayMapRead (fNam, CW, READONLY, bbG->h) ; /* memory map the seed index */

  if (0) showCws (bbG->cws) ;
  
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
  int step = bigArr (seqIds, 255, char) ;
  if (
      (step > 1 && pp->step > 1) &&
      (
       (step >= pp->step && (step % pp->step) != 1)  ||
       (step < pp->step && (pp->step % step) != 1)
       )
      )
    messcrash ("\nThe target is indexed with step=%d, the requested read step=%d, these number are not relative primes, there will be systematic false negatives,\n please set the argu7ment --step %d (default 2) to a lower value", step, pp->step, pp->step) ;

  bbG->step = step ;
  /* check the version */
  char *signature = hprintf (h, "Sort Align version 1 step %d", step) ;
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
  cpuStatRegister ("1.memMapTargets" , pp->agent, bbG->cpuStats, t1, t2, bigArrayMax (bbG->cws)) ;
 
  ac_free (h) ;
  return TRUE ;
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
  genomeParseBinary (pp, &bbG) ;
  t2 = clock () ;
  cpuStatRegister ("1.GParserDone" , pp->agent, bbG.cpuStats, t1, t2, bigArrayMax (bbG.cws) - 1) ;
  channelPut (pp->gmChan, &bbG, BB) ;
  channelClose (pp->gmChan) ;
  printf ("--- %s: Stop binary genome parser\n", timeBufShowNow (tBuf)) ;

  return ;
} /* genomeParser */

/**************************************************************/
/**************************************************************/

static void storeTargetIndex (const PP *pp, BB *bbG, int step) 
{
  AC_HANDLE h = ac_new_handle () ;
  const char *fNam = 0 ;
  
  /* export the code words */
  fNam = hprintf (h, "%s/cws.sortali", pp->indexName) ;
  bigArrayMapWrite (bbG->cws, fNam) ;
  fprintf (stderr, "genomeCreateBinary exported %ld seed records\n", bigArrayMax (bbG->cws)) ;

  if (0) showCws (bbG->cws) ;
    
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
  char *signature = hprintf (h, "Sort Align version 1 step %d", step) ;
  strcpy (buf, signature) ;
  memcpy (cp, buf, 255) ;
  cp[255] = step ;
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
/* parse, code, sort the genome and create the index on disk
 * the human index takes around 18 GigaBytes
 */
static void createTargetIndex (const PP *pp, BB *bbG, Array tArray)
{
  AC_HANDLE h = ac_new_handle () ;
  BigArray cws = 0 ;
  int nMax = arrayMax (tArray) ;
  TC *tc = 0 ;
  char tBuf[25] ;
  clock_t t1, t2 ;
  int step = 1 ;
  
  memset (bbG, 0, sizeof (BB)) ;
  bbG->isGenome = TRUE ;
  t1 = clock () ;
  printf ("+++ %s: Parse the target files\n", timeBufShowNow (tBuf)) ;

  /* parse all targets into a single bbG->seqs array, with targetClass prefix in the sequnce name */
  for (int nn = 0 ; nn < nMax ; nn++)
    {
      tc = arrayp (tArray, nn, TC) ;
      RC rc ;

      memset (&rc, 0, sizeof (RC)) ;
      rc.fileName1 = tc->fileName ;
      rc.format = tc->format ;
      rc.run = nn + 1 ;
      sequenceParser (pp, 0, tc, bbG, nn == nMax - 1 ? 2 : 1 ) ; /* 2 for last target */
      tc->step = (bbG->length < 1<<29) ? 1 : 3 ;
    }
  t2 = clock () ;
  cpuStatRegister ("1.Parse targets" , pp->agent, bbG->cpuStats, t1, t2, arrayMax (bbG->dnas) - 1) ;

  t1 = clock () ;
  printf ("%s : extract the target seeds\n" , timeShowNow (tBuf)) ;
  cws = codeWordsDo (pp, bbG, tc->step, TRUE) ;
  t2 = clock () ;
  cpuStatRegister ("1.Create target seeds" , pp->agent, bbG->cpuStats, t1, t2, bbG->nSeqs) ;
  t1 = clock () ;

  printf ("%s : sort the target seeds\n" , timeShowNow (tBuf)) ;

  if (0)
    {
      AC_HANDLE h = ac_new_handle () ;
      BigArray cws1 = bigArrayHandleCopy (cws, h) ;
      BigArray cws2 = bigArrayHandleCopy (cws, h) ;
      BigArray cws3 = bigArrayHandleCopy (cws, h) ;
      
      t1 = clock () ;
      bigArraySort (cws1, cwOrder) ;
      t2 = clock () ;
      printf ("\ncws bigArraySort %ld \n" , t2 - t1) ;
      
      t1 = clock () ;
      hitRadixSort (cws2) ;
      t2 = clock () ;
      printf ("cws radixSort %ld \n" , t2 - t1) ;

      t1 = clock () ;
      newMsort (cws3, cwOrder) ;
      t2 = clock () ;
      printf ("\ncws newMsort %ld \n" , t2 - t1) ;
      exit (0) ;
    }
  newMsort (cws, cwOrder) ;

  t2 = clock () ;
  cpuStatRegister ("2.Sort seeds" , pp->agent, bbG->cpuStats, t1, t2, bigArrayMax (cws) - 1) ;
  t1 = clock () ;

  printf ("%s : write the index to disk\n" , timeShowNow (tBuf)) ;
  bbG->cws = GenomeAddSkips (pp, cws, bbG->h) ;
  storeTargetIndex (pp, bbG, step) ;
  
  t2 = clock () ;
  cpuStatRegister ("3.Write the index to disk" , pp->agent, bbG->cpuStats, t1, t2, bigArrayMax (bbG->cws) - 1) ;

  ac_free (h) ;
  return ;
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
       if (0)
	 {
	   AC_HANDLE h = ac_new_handle () ;
	   BigArray cws1 = bigArrayHandleCopy (bb.cws, h) ;
	   BigArray cws2 = bigArrayHandleCopy (bb.cws, h) ;
	   BigArray cws3 = bigArrayHandleCopy (bb.cws, h) ;
	   
	   t1 = clock () ;
	   bigArraySort (cws1, cwOrder) ;
	   t2 = clock () ;
	   printf ("\ncws bigArraySort %ld \n" , t2 - t1) ;
	   
	   t1 = clock () ;
	   hitRadixSort (cws2) ;
	   t2 = clock () ;
	   printf ("cws radixSort %ld \n" , t2 - t1) ;
	   
	   t1 = clock () ;
	   newMsort (cws3, cwOrder) ;
	   t2 = clock () ;
	   printf ("\ncws newMsort %ld \n" , t2 - t1) ;
	   exit (0) ;
	 }
      if (bb.cws)
	newMsort (bb.cws, cwOrder) ;
      t2 = clock () ;
      cpuStatRegister ("4.SortWords", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.cws)) ;
      channelPut (pp->smChan, &bb, BB) ;
      if (pp->debug) printf ("--- %s: Stop sort words %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.cws)) ;
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

static long int  matchHitsDo (BB *bbG, BB *bb)
{
  long int i = 0, iMax = bigArrayMax (bbG->cws);
  long int j = 0, jMax = bigArrayMax (bb->cws);
  long int k = 0, nn = 0 ;
  int nHA = 0, hMax = 100000 ;
  unsigned int mask = step1 - 1 ;
  const CW *restrict cw = bigArrp(bbG->cws, 0, CW) ;
  const CW *restrict rw = bigArrp(bb->cws, 0, CW) ;
  const CW *restrict cw1;
  const CW *restrict cwMax = cw + iMax ;
  HIT *restrict hit;
  BigArray hitsArray = bigArrayHandleCreate(64, BigArray, bb->h);
  BigArray hits = bigArrayHandleCreate(hMax, HIT, bb->h);

  if (0)
    {
      showCws (bb->cws) ;
      showCws (bbG->cws) ;
    }
  bigArray (hitsArray, nHA++, BigArray) = hits ;
  while  (i < iMax && j < jMax)
    {
      if ((i & mask) == 0)
	{
#ifdef JUNK	  
	  if (rw->word > (unsigned int) cw->targetClass)
	    {
	      cw += step4 ;
	      i += step4 ;
	      bb->skips4++;
	      continue;
	    }
	  else if (rw->word > (unsigned int) cw->nam)
	    {
	      cw += step3 ;
	      i += step3 ;
	      bb->skips3++;
	      continue;
	    }
	  else if (rw->word > (unsigned int) cw->pos)
	    {
	      cw += step2 ;
	      i += step2 ;
	      bb->skips2++;
	      continue;
	    }
	  else if (rw->word > (unsigned int) cw->word)
	    {
	      cw += step1 ;
	      i += step1 ;
	      bb->skips1++;
	      continue;
	    }
	  else
	    {
	      cw++ ;
	      i++ ;
	      bb->skips0++;
	      continue;
	    }
#endif
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
	{ i++; cw++; }
      else if (cw->word > rw->word)
	{ j++ ; rw++ ; }
      else
	{
	  long int a1, x1, i1 = i ;
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
			{ a1 = a1 + 1 - 15 ; }
		    }
		  else /* x1 < 0 but strand info must reside on the a1 coordinate */
		    {
		      
		      if (a1 < 0)
			{ a1 = -a1 - 1; x1 = -x1 - 1 ; }
		      else
			{ a1 = -a1 - 15 + 1 ; x1 = -x1 - 1 ; }
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
  bigArrayDestroy (bb->cws) ; /* the read words are no longer needed */
  bigArrayMax (hits) = k ;
  bb->hits = hitsArray ;

  return nn ;
}  /* matchHitsDo */

/**************************************************************/
static void matchHits (const void *vp)
{
  const PP *pp = vp ;
  BB bb ;
  BB bbG = pp->bbG;
  char tBuf[25] ;
  
  clock_t  t1, t2 ;
	    
  memset (&bb, 0, sizeof (BB)) ;
  /* grab and match a block of reads */
  while (channelGet (pp->smChan, &bb, BB))
    {
      long int nn = 0 ;
      long int iMax = bbG.cws ? bigArrayMax (bbG.cws) : 0 ;
      long int jMax = bb.cws ? bigArrayMax (bb.cws) : 0 ;
      t1 = clock () ;

      if (iMax && jMax)
	{
	  if (pp->debug) printf ("+++ %s: Start match %ld target %ld read words\n", timeBufShowNow (tBuf), iMax, jMax) ;
	  nn = matchHitsDo (&bbG, &bb) ;
	  if (pp->debug) printf ("--- %s: Stop match hits constructed %ld arrays\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
	}
      t2 = clock () ;
      cpuStatRegister ("5.MatchHits", pp->agent, bb.cpuStats, t1, t2, nn) ;
      channelPut (pp->moChan, &bb, BB) ;
    }
  if (1)
    {
      int n = channelCount (pp->plChan) ;
      if (pp->debug) printf ("..... close moChan at %d\n", n) ;
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
      if (0 && pp->align && bb.hits)
	{ /* speed tests */
	  
	  AC_HANDLE h = ac_new_handle () ;
	  BigArray b1, b2 , b3, b4 ;
	  b1 =bigArrayHandleCopy (bb.hits, h) ;
	  b2 =bigArrayHandleCopy (bb.hits, h) ;
	  b3 =bigArrayHandleCopy (bb.hits, h) ;

	  t1 = clock () ;
	  bigArraySort (b1, hitReadOrder) ;
	  t2 = clock () ;
	  printf ("\nbigArraySort %ld \n" , t2 - t1) ;

	  t1 = clock () ;
	  hitRadixSort (b2) ;
	  t2 = clock () ;
	  printf ("radixSort %ld \n" , t2 - t1) ;

	  t1 = clock () ;
	  newMsort (b3, hitReadOrder) ;
	  t2 = clock () ;
	  printf ("newMsort %ld \n" , t2 - t1) ;
 

	  for (int i = 0 ; i < 20 ; i++)
	    printf("%d\t%d\t%d\n"
		   , bigArr(b1,i,HIT).read
		   , bigArr(b2,i,HIT).read
		   , bigArr(b3,i,HIT).read
		   ) ;
	  for (int i = 0 ; i < bigArrayMax (b1) ; i++)
	    if (bigArr(b1,i,HIT).read != bigArr(b2,i,HIT).read)
	      messcrash ("Error in bigArraySort or in radixSort line %d\n", i) ;
	  for (int i = 0 ; i < bigArrayMax (b1) ; i++)
	    if (bigArr(b1,i,HIT).read != bigArr(b3,i,HIT).read)
	      messcrash ("\nError in newMsort line %d/%ld\n", i, bigArrayMax(b1)) ;
	  exit (0) ; 
	}	   
      
      newMsort (bb.hits, hitReadOrder) ;
      
      t2 = clock () ;
      cpuStatRegister ("6.OrderHits", pp->agent, bb.cpuStats, t1, t2, bigArrayMax (bb.hits)) ;
    
      if (pp->debug) printf ("--- %s: Stop sort hits %ld\n", timeBufShowNow (tBuf), bigArrayMax (bb.hits)) ;
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
/* Establish chain scores, select best */
static void  alignSelectBestChain (const PP *pp, BB *bb, BigArray aaa, Array aa, int errCost)
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
      up->score = up->x2 - up->x1 + 1 - up->nerr * errCost ;
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
      if (bb->runStat.nPairs && (up->read & 0x1))
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
  
  for (ii = 0, up = arrp (aa, ii, ALIGN) ; ii < iMax ; ii++, up++)
    {
      vp = bigArrayp (aaa, kMax++, ALIGN) ;
      int a1 = up->a1 ;
      int a2 = up->a2 ;
      if (isRead2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ;}
      up->w1 = (a1 + demiStep)/step ;
      up->w2 = (a2 + demiStep)/step ;
      *vp = *up ;
      bb->runStat.nErr += up->nerr ;
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
      if (up->errors)
	{
	  unsigned int flip = 0 ;
	  if (bb->runStat.nPairs && (up->read & 0x1))
	    flip = 0x0f ; /* will flip last 4 bits */
	  if (1) mergeErrors (bb->errors, up->errors, flip) ;
	  arrayDestroy (up->errors) ;
	}
    }
  return ;
} /* alignSelectBestChain */

/**************************************************************/

static void alignDo (const PP *pp, BB *bb)
{
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
  BOOL debug = FALSE ;
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
	    alignSelectBestChain (pp, bb, aaa, aa, errCost) ;
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
	      ap->nerr = arrayMax (err) ;
	      if (ap->nerr)
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
    alignSelectBestChain (pp, bb, aaa, aa, errCost) ;
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

/**************************************************************/

static void exportDo (ACEOUT ao, const PP *pp, BB *bb)
{
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ii, aMax = bigArrayMax (bb->aligns) ;
  DICT *dict = bb->dict ;
  DICT *dictG = pp->bbG.dict ;
  const char *run = dictName (pp->runDict, bb->run) ;
  
  for (ii = 0 ; ii < aMax ; ii++, ap++)
    {
      int read = ap->read ;
      int chrom = ap->chrom ;
      int x1 = ap->x1 ;
      int x2 = ap->x2 ;
      int a1 = ap->a1 ;
      int a2 = ap->a2 ;
      int dx = x2 - x1 + 1 ;
      int nerr = ap->nerr ;

      aceOutf (ao, "\n%s\t%s", run, dictName (dict, read)) ; 
      aceOutf (ao, "\t%d", ap->score) ;
      aceOutf (ao, "\t%d", ap->chain) ;
      aceOutf (ao, "\t%d", ap->readLength) ;
      aceOutf (ao, "\t%d\t%d\t%d", x1, x2, dx) ;
      aceOutf (ao, "\t%s", dictName (dictG, chrom)) ; 
      aceOutf (ao, "\t%d\t%d", a1, a2) ;
      aceOutf (ao, "\t%d", nerr) ;
    }
  aceOut  (ao, "\n") ;
  return ;
} /* exportDo */

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
  ACEOUT ao = aceOutCreate (pp->outFileName, ".hits", pp->gzo, h) ;
  aceOutDate (ao, "##", "sortaling test") ;
  aceOutf (ao, "#Run\tRead\tScore\tChain\tLength\tx1\tx2\tdx\tTarget\ta1\ta2\tMismatch") ;

  memset (&bb, 0, sizeof (BB)) ;
  while (channelGet (pp->weChan, &bb, BB))
    {
      t1 = clock () ;
      if (bb.aligns && bigArrayMax (bb.aligns))
	{
	  bigArraySort (bb.aligns, alignOrder) ;
	  exportDo (ao, pp, &bb) ;

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
  aceOut (ao, "\n") ;

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
	
    printf ("\n####### Run Statistics\n") ;
    printf ("# Run\tFiles\tPairs\tReads\tBases1\tBases2\tLn1\tLn2") ;
    printf ("\tReads aligned\tBase aligned\tBase2 aligned\tMismatches") ;
    for (int j = 0 ; j < 256 ; j++)
      if (s0->nAlignedPerTargetClass[j])
	printf ("\tReads Aligned in %c", j) ;
    for (int j = 1 ; j < 11 ; j++)
      printf ("\t%d aligmnents", j) ;
    for (ii = 0 ; ii < iMax ; ii++) 
      {	
	RunSTAT *s = arrayp (aa, ii, RunSTAT) ;
	if (! s->nReads)
	  continue ;
	printf ("\n%s", s->run ? dictName (pp->runDict, s->run) : "Any") ;
	printf ("\t%d", s->nFiles) ;
	printf ("\t%ld", s->nPairs) ;
	printf ("\t%ld", s->nReads) ;
	printf ("\t%ld\t%ld", s->nBase1, s->nBase2) ;
	float nR1 = s->nPairs ? s->nPairs : s->nReads ;
	float nR2 = s->nPairs ? s->nPairs : 1 ;
	printf ("\t%.1f\t%.1f", s->nBase1/nR1, s->nBase2/nR2) ;
	printf ("\t%ld", s->nMultiAligned[0]) ;
	printf ("\t%ld", s->nBaseAligned1) ;
	printf ("\t%ld", s->nBaseAligned2) ;
	printf ("\t%ld", s->nErr) ;
	for (int j = 0 ; j < 256 ; j++)
	  if (s0->nAlignedPerTargetClass[j])
	    printf ("\t%ld", s->nAlignedPerTargetClass[j]) ;
	for (int j = 1 ; j < 11 ; j++)
	  printf ("\t%ld", s->nMultiAligned[j]) ;
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
	       "//      sortalign --index XYZ -i f.fastq.gz --align --wiggle -o results/xxx \n"
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
	       "//        example: zcat fx.*.fastq.gz | sortalign -x XYZ -i - --gzi --fastq --align -o results/fx\n"
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
	       "// --align : Extend the seed alignments hopefully to full read length\n"
	       "//   --minAli <int> : minimal length of aech aligned fragmant (exon)\n"
	       "//   --errMax <int> : [default NA] maximal number of mismatches in any (partial) alignment\n"
	       "//   --errRateMax <int> : [default 10] maximal percentage of mismatches in any (partial) alignment\n"
	       "// --wiggle  : Report target coverage wiggles in UCSC BF (fixed) format\n"
	       "// --intron  : (not yet ready) Report intron support\n"
	       "// --snp : (not yet ready) Report candidate SNP counts (substitutions and short indels)\n"
	       "// OPTIONAL TECHNICAL PARAMETERS\n"
	       "//   The program sortalign is parallelized and synchronized by GO-language like channels\n"
	       "//   All code layers execute at the same time, the work-load among layers is self balancing\n"
	       "//   There is no limit to the size of the input, as data continuously flow in and out\n"
	       "// --nA or --nAgents <int> : [default 10] number of agent in each code layer\n"
	       "// --nB or --nBlocks <int> : [default 10] number of simultaneous data blocks circulating in the pipeline\n"
	       "// --max_threads <int>  : [default 128] maximal number of simultaneous UNIX threads\n"
	       "//    Possible values are \n"
	       "//        --nAgents 1 --nBlocks 1 --max_threads 8 : for a small test\n"
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
  long int skips0 = 0, skips1 = 0, skips2 = 0, skips3 = 0, skips4 = 0 ;
  char tBuf[25] ;
  char tBuf0[25] ;
  AC_HANDLE h ;
  int nAgents = 10 ;
  int channelDepth = 1 ;
  mytime_t t0, t1 ;
  
  freeinit () ; 
  messErrorInit (argv[0]) ;

  if (0)
    {
      int n = 4 ;
      getCmdLineInt (&argc, argv, "-n", &n) ;
      BigArray b, a = bigArrayHandleCreate (n, int, 0) ;
  
      for (int i = 0 ; i < n ; i++)
	bigArray (a, i, int) = randint() ;    
      b = bigArrayHandleCopy (a, 0) ;
      bigArraySort (b, intOrder) ;
      
      for (int i = 0 ; i < n && i < 12 ; i++)
	printf("%d %d\n"
	       , bigArr (a,i,int)   
	       , bigArr (b,i,int)
	       ) ;
      for (int i = 0 ; i < n-1 ; i++)
	if (bigArr (b,i, int) > bigArr (b,i+1, int))
	  messcrash ("error line %d\n", i) ;
      exit (0) ;
    }
  
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
  
  memset (&p, 0, sizeof (PP)) ;
  h = ac_new_handle () ;
  p.h = h ;
  cpuStats = arrayHandleCreate (1024, CpuSTAT, h) ;
  runStats = arrayHandleCreate (1024, RunSTAT, h) ;
  runErrors = arrayHandleCreate (1024, Array, h) ;
  
  if (argc < 2)
    usage (0, 0, argv) ;
  if (getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)||
      getCmdLineOption (&argc, argv, "-h", 0)
      )
    usage (0, 0, argv) ;


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

  p.debug = getCmdLineOption (&argc, argv, "--debug", 0) ;
  p.debug |= getCmdLineOption (&argc, argv, "--verbose", 0) ;

  nAgents = 10 ;
  if (! getCmdLineInt (&argc, argv, "--nAgents", &(nAgents)))
    getCmdLineInt (&argc, argv, "--nA", &(nAgents)) ;
  maxThreads = 128 ;  /* UNIX  max om lmem12 machine */
  p.maxBB = 10 ;  /* max number of BB blocks processed in parallel */
  if (!getCmdLineInt (&argc, argv, "--nBlocks", &(p.maxBB)))
    getCmdLineInt (&argc, argv, "--nB", &(p.maxBB));
  p.align = getCmdLineBool (&argc, argv, "--align");
  p.wiggle = getCmdLineBool (&argc, argv, "--wiggle");

  p.fasta = getCmdLineBool (&argc, argv, "--fasta");
  p.fastq = getCmdLineBool (&argc, argv, "--fastq");
  p.fastc = getCmdLineBool (&argc, argv, "--fastc");
  p.raw = getCmdLineBool (&argc, argv, "--raw");

  p.gzi = getCmdLineBool (&argc, argv, "-gzi") ||
    getCmdLineBool (&argc, argv, "--gzi");
  p.gzo = getCmdLineBool (&argc, argv, "-gzo") ||
    getCmdLineBool (&argc, argv, "--gzo");

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
  p.run = 0 ;
  getCmdLineOption (&argc, argv, "-r", &(p.runName)) ;
  getCmdLineOption (&argc, argv, "--run", &(p.runName)) ;

  if  ( ! getCmdLineInt (&argc, argv, "--max_threads", &maxThreads))
    maxThreads = 128 ;
  if (p.maxBB == 1)
    { nAgents = 1 ; }
  if (maxThreads < 24)
    maxThreads = 24 ;
  if (p.createIndex) { nAgents = 1 ; maxThreads = 2 ; p.maxBB = 1 ; }

  p.minAli = 30 ;
  p.errMax = 10 ;
  p.errRateMax = 10 ;
  getCmdLineInt (&argc, argv, "--errMax", &(p.errMax)) ;
  getCmdLineInt (&argc, argv, "--minAli", &(p.minAli)) ;
  getCmdLineInt (&argc, argv, "--errRatMax", &(p.errRateMax)) ;
  p.step = 2 ;   /* read default */
  if (p.createIndex) p.step = 0 ; /* default 3 for large targets, 1 for short, set in  createIndex() */
  getCmdLineInt (&argc, argv, "--step", &(p.step)) ;

  if (argc > 1)
    usage (0, argc, argv) ;

  t0 = timeNow () ;
  printf ("%s: Start\n", timeBufShowNow (tBuf0)) ;

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
	  cq = hprintf (h1, "mkdir %s ; echo test > %s/sortalign.testFile", cp, cp) ;
	  system (cq) ;
	}
      cq = hprintf (h1, "%s/sortalign.testFile", p.indexName) ;
      ACEIN ai = aceInCreate (cq, FALSE, h1) ;
      if (! ai)
	usage (hprintf (p.h, "Cannot create and access the index directory: %s", p.indexName), argc, argv) ;
      ac_free (h1) ;
    }
  
  p.runDict = dictHandleCreate (16, p.h) ;
  p.targetClassDict = dictHandleCreate (16, p.h) ;

  dictAdd (p.targetClassDict, "rRNA", 0) ;
  dictAdd (p.targetClassDict, "Mito", 0) ;
  dictAdd (p.targetClassDict, "Genome", 0) ;
  dictAdd (p.targetClassDict, "Bacteria", 0) ;
  dictAdd (p.targetClassDict, "Virus", 0) ;

  /* create the index */
  if (p.createIndex)
    { /* The human genome index consumes around 18 Gigabytes of RAM */
      /* check that input files were provided */
      if (! p.tFileName && ! p.tConfigFileName)
	usage ("--createIndex requires providing a target parameter -t or -T", argc, argv) ;
      if ( p.tFileName &&  p.tConfigFileName)
	usage ("conflicting parameters -t and -T, both define the targets", argc, argv) ;

      Array tArray = parseTargetConfig (&p, runStats) ;
      if (p.tConfigFileName)
	system (hprintf(h, "\\cp %s %s\n", p.tConfigFileName, p.indexName)) ;
      createTargetIndex (&p, &p.bbG, tArray) ;

      goto done ;
    }
  
  /* check the existence of the index files */
  if (p.indexName)
    {
      BOOL ok = TRUE ;
      char *cp = filName (p.indexName, "/cws.sortali", "rb") ;
      if (cp)
	p.tFileBinaryCwsName = strnew (cp, p.h) ;
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

  /* check that input files were provided */
  if (! p.inFileName && ! p.inConfigFileName)
     usage ("missing parameter -i or -I inputFileName(s)", argc, argv) ;
  if (p.inFileName && p.inConfigFileName)
     usage ("conflicting parameter -i and -I, both define the input files", argc, argv) ;
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
  for (int i = 0 ; i < nAgents && i < p.maxBB ; i++)
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
  if (! p.bbG.cws)
    messcrash ("matchHits received no target words") ;
  
  /* map the reads to the genome in parallel */
  for (int i = 0 ; i < nAgents && i < p.maxBB ; i++)
    {
      p.agent = i ;
      wego_go (matchHits, &p, PP) ;
      wego_go (orderHits, &p, PP) ;
      if (!i) /* only 1 wiggle agent */
	wego_go (wiggle, &p, PP) ;
    }
  for (int i = 0 ; i < 2*nAgents && i < 2*p.maxBB ; i++)
    {
      p.agent = i ;
      if (!i || p.align)
	wego_go (align, &p, PP) ;
    }
  /* We can only declare a single export agent */
  p.agent = 0 ;
  wego_go (export, &p, PP) ; 
  
  while (channelGet (p.doneChan, &bb, BB))
    {
      long int n = (bb.hits ? bigArrayMax (bb.hits) : 0) ;

      if (bb.isGenome && p.bbG.cws)
	{
	  cpuStatCumulate (cpuStats, p.bbG.cpuStats) ;
	  bigArrayDestroy (p.bbG.cws) ;
	  continue ;
	}
      if (p.debug) printf ("%s:Block done\n", timeBufShowNow (tBuf)) ;
      if (p.debug) printf ("Found %ld hits\n", n) ; 

      skips0 += bb.skips0 ;
      skips1 += bb.skips1 ;
      skips2 += bb.skips2 ;
      skips3 += bb.skips3 ;
      skips4 += bb.skips4 ;

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
  if (1 || p.debug) printf ("Skips: 0=%ld, 256=%ld, 1024=%ld, 4096=%ld, 16384=%ld\n",
			    skips0, skips1, skips2, skips3, skips4);
  if (arrayMax (runStats))
    reportRunStats (&p, runStats) ;
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
 
