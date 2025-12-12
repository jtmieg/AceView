/*
 * sortalign : RNA aligner

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
/* Examples:  sortalign -h */

/***********************************************************************************/
/*  Greg: uncomment this lines and edit the code at the very bottom of sa.sort.c and recompile sa.*.c */
/* #define GREGGPU  */

#ifdef GREGGPU 
  #define NAGENTS 12
  #define NBLOCKS 12
#else
  #define YANN 
  #define NAGENTS 40
  #define NBLOCKS 40
#endif

/***********************************************************************************/

#include "ac.h"
#include "channel.h"
#include <dirent.h>
#include <zlib.h>
#include <stdatomic.h>
#include "../wsra/sra_read.h"




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
  atomic_int laneDone ;
  ACEOUT aoSam ;
} RC ;
		  
typedef struct runStatStruct {
  int run ;   /* index in p.runDict */
  int nFiles ;
  long int nPairs, nReads ;
  long int nCompatiblePairs, nIncompatiblePairs, n2ChromsPairs, nOrphans, nCirclePairs ;
  long int nBase1, nBase2 ;
  long int nPairsAligned, nBaseAligned1, nBaseAligned2 ;
  long int intergenic, intronic, exonic ;
  long int nMultiAligned[11] ;
  long int nAlignedPerTargetClass[256] ;
  long int nPerfectReads ;
  long int nAlignments ;
  long int nClippedPolyA ;
  long int nClippedVectors ;
  long int nSupportedIntrons ;
  long int nSupportedPolyA ;
  float intronStranding ;
  long int nIntronSupportPlus ;
  long int nIntronSupportMinus ;
  long int wiggleCumul ; /* in million bases */
  long int wiggleLCumul ; /* in million bases */
  long int wiggleRCumul ; /* in million bases */

  int minReadLength, maxReadLength ;
  char *adaptor1, *adaptor2 ;
  long int nN, nErr ;
  long int NATGC[5] ;
  int GF[256], GR[256] ; /* number of reads aligned per target_class on Forward/Reverse strand */
  Array errors ;  /* substitutions, insertions, deletions counts */
  /* coverage of long transcripts ? */
} RunSTAT ;
		  
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
  Array wigglesL ;
  Array wigglesR ;
  Array wigglesP ;
  Array wigglesNU ;
  Array confirmedIntrons ;
  BOOL isGenome ;
  int step, skips0, skips1, skips2, skips3, skips4, skipsFound, skipsNotFound ;
  long int nPolyASuppor ;
  long int nPolyASupport ;
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
  BOOL wiggleEnds ;
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
  const char *tFileSpongeFileNameF ;
  const char *tFileSpongeFileNameR ;

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
  CHAN *doneChan ; /* return to main program */
  CHAN *wwChan ; /* Wiggler */
  CHAN *wwDoneChan ; /* Wiggler */


  BB bbG ;  /* genome or genes target */
  Array runStats ;
  Array runLanes ;
  Array runLanesDone ;
  Array geneExons ;           /* gene coordinates */
  Array geneBoxes ;           /* gene coordinates */
  Array geneExonCounts ;      /* gene expression */
  Array geneBoxCounts ;      /* gene expression */
  BigArray intronSeeds ;
  BigArray exonSeeds ;
  Array wiggles ;
  Array wiggleCumuls ;
  Array wigglesL ;
  Array wiggleLCumuls ;
  Array wigglesR ;
  Array wiggleRCumuls ;
  Array wigglesP ;
  Array wigglesNU ;
  Array exonics ;
  Array intronics ;
  Array intergenics ;
  Array confirmedIntrons ;
  BOOL fasta, fastq, fastc, raw, solid, sra, sraCaching ;
  BOOL sam, exportSamSequence, exportSamQuality ;
  BOOL strand, antiStrand ;
  int bonus[256] ;
  DICT *runDict ;
  DICT *targetClassDict ;
  DICT *geneDict ;
  DICT *mrnaDict ;
  DICT *wiggleFileDict ;
  const char *method ;
  int run ;
  int nFiles ;  /* number of input sequence files */
  int agent ;   /* instance of the agent */
  int nIndex ;  /* number of parallel target indexes */
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
  int maxSraGb ; /* max number of Gigabases in each SRA download, 0 : no max */
  long int wiggleCumul ; /* in million bases */
  long int exonic, intronic, intergenic ;
  BOOL splice ;
  long int nRawReads, nRawBases ;
  int wiggle_step ;
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
}  __attribute__((aligned(32))) COUNTCHROM ;

typedef struct alignStruct {
  int read ;
  int targetClass ;
  int chrom ;
  int a0, x0 ;  /* coordinate of the hit before extension */
  int a1, a2 ;  /* bio position in chrom, a1 < a2 on minus strand */
  int x1, x2 ;  /* bio position in read, x1 < x2 always */
  int chain, chainX1, chainX2, chainA1, chainA2 ;
  int id, previous, next ;
  int ali, chainAli, score, chainScore ;
  int pairScore, mateChrom, mateA1, mateA2, pairLength ;
  int nN, nErr, chainErr ;
  int nTargetRepeats ;
  int nChains ;
  int readLength ;
  int leftClip, rightClip ;
  int pA ;
  int errShort, errLong ; /* bb->dict */
  int leftOverhang, rightOverhang ; /* bb->dict */
  int donor, acceptor ;
  Array errors ;
} __attribute__((aligned(32))) ALIGN ;


typedef struct geneStruct {
  int chrom ;
  int a1, a2 ; 
  int gene ; /* index in pp->geneDict */
} __attribute__((aligned(16))) GENE ;
		  
typedef struct mrnaStruct {
  int chrom ;
  int a1, a2 ; 
  int mrna ; /* index in pp->mrnaDict */
} __attribute__((aligned(16))) MRNA ;
		  

typedef struct intronStruct {
  int run ;
  int mrna ;
  int chrom ;
  int n, nR, a1, a2 ;
  char feet[6] ;
} __attribute__((aligned(32))) INTRON ;
		  
typedef struct exonStruct {
  int chrom ;
  int a1, a2 ; 
  int mrna ;
} __attribute__((aligned(32))) EXON ;
		  

typedef struct polyAStruct {
  int chrom ;
  int a1, a2 ; 
  int n ;
} __attribute__((aligned(32))) POLYA ;
		  

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

/**************************************************************/

/* sa.main.c */
void saUsage (char *message, int argc, const char **argv) ; 
void saCpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, long int n) ;
  
/* sa.config.c */
void saConfigIsIndexAccessible (PP *pp) ;
void saConfigIsOutDirWritable (PP *pp) ;
int saConfigCheckTargetIndex (PP *pp) ; 
Array saConfigGetRuns (PP *pp, Array runStats) ;    

/* sa.gff.c */
long int saGffParser (PP *pp, TC *tc) ;
long int saIntronParser (PP *pp, TC *tc) ; 
int saSpongeParser (PP *pp, TC *tc) ;
int saSpongeParserDirect (PP *pp) ;

/* sa.sort.c */
void saSort (BigArray aa, int type) ;
/* sa.gpusort.c */
void saGPUSort (char *cp, long int number_of_records, int size_of_record, int (*cmp)(const void *va, const void *vb)) ;

/* sa.sam.c */
int saSamExport (ACEOUT ao, const PP *pp, BB *bb) ;       
ACEOUT saSamCreateFile (const PP *pp, BB *bb, AC_HANDLE h) ; 

/* sa.introns.c */
void saIntronsOptimize (BB *bb, ALIGN *vp, ALIGN *wp, Array dnaG)  ; 
void saIntronsExport (PP *pp, Array aaa) ;
void saIntronsCumulate (PP *pp, BB *bb) ;
int saSupportedIntrons (const PP *pp, int run) ;

/* sa.targetIndex.c */
void saTargetIndexCreate (PP *pp) ;
void saTargetIndexGenomeParser (const void *vp) ;

/* sa.seeds.c */
int saCodeIntronSeeds (PP *pp, BB *bbG) ;
void saCodeSequenceSeeds (const PP *pp, BB *bb, int step, BOOL isTarget) ;
  
/* sa.sequenceParser.c */
void saSequenceParse (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome) ;
int saSequenceParseSraDownload (const char *sraID, int maxGb) ;
void saSequenceParseGzBuffer (const PP *pp, BB *bb) ;

/* sa.wiggle */
void saWiggleCumulate (const PP *pp, BB *bb) ;
void saWiggleExport (PP *pp, int nAgents) ;

/* sa.stats */
void saRunStatExport (const PP *pp, Array runStats) ;
void saCpuStatExport (const PP *pp, Array stats) ;
void saCpuStatCumulate (Array aa, Array a) ;
void saRunStatsCumulate (int run, PP *pp, BB *bb) ;

/* sa.align */
void saAlign (const void *vp) ;
void saAlignDo (const PP *pp, BB *bb) ;
int saAlignOrder (const void *va, const void *vb) ;

/**************************************************************/
/**************************************************************/
/**************************************************************/
