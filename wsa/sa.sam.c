/*
 * sa.sam.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This code handles the sam format
 * Given an alignment in internal format (ALI structures and _ERR) it exports a sam file
 */

#include "sa.h"

/*********************************************************************/
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
  int lSam = 0, da, ii, iMax = errors ? arrayMax (errors) : 0 ;
  int a1 = isDown ? ap->a1 : ap->a2 ;
  int a2 = isDown ? ap->a2 : ap->a1 ;
  const int L = 16 ;
  char segs[2*iMax + 4][L] ;
  int iSeg = 0 ;
  
  if (iMax)
    for (ii = 0 ; ii < iMax ; ii++)
      {
	ep = arrp (errors, ii, A_ERR)  ;
	if (ep->iShort < ap->x1 - 1 || ep->iShort >= ap->x2) continue ; /* may happen when fixing a duplication not authorized in sam */
	if (isDown && ep->iLong + 1 <= ap->a1)
	  continue ;
	if (! isDown && ep->iLong + 1 > ap->a1)
	  continue ;
	int xShort = ep->iShort + 1 ;
	int xLong = ep->iLong + 1 ;
	if (xShort <= 0 || xLong <= 0)
	  continue ;
	
	da = ep->iLong - a1 + 1 ; /* number of exact matches, not including the problem */
	if (da < 0)
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
	    a1 = ep->iLong + 1 ;
	    lSam += da + 1 ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "1I") ;
	    *nMp += da ;
	    break ;
	  case INSERTION_DOUBLE:     
	    *nInsp += 2 ;
	    a1 = ep->iLong + 1 ;
	    lSam += da + 2 ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "2I") ;
	    *nMp += da ;
	    break ;
	  case INSERTION_TRIPLE:     
	    *nInsp += 3 ;
	    a1 = ep->iLong + 1 ;
	    lSam += da + 3 ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "3I") ;
	    *nMp += da ;
	    break ;
	  case  TROU:
	    *nDelp += 1 ;
	    a1 = ep->iLong + 2 ;
	    lSam += da ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "1D") ;
	    *nMp += da ;
	    break ;
	  case TROU_DOUBLE:
	    *nDelp += 2 ;
	    a1 = ep->iLong + 3 ;
	    lSam += da ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "2D") ;
	    *nMp += da ;
	    break ;
	  case TROU_TRIPLE:
	    *nDelp += 3 ;
	    a1 = ep->iLong + 4 ;
	    lSam += da ;
	    if (da > 0) snprintf(segs[iSeg++], L, "%dM", da) ;
	    snprintf(segs[iSeg++], L, "3D") ;
	    *nMp += da ;
	    break ;
	  }
      }
  
  da = a2 - a1 + 1 ;
  if (da > 0)
    {
      snprintf(segs[iSeg++], L, "%dM", da) ;
      *nMp += da ;
      lSam += da ;  
    }
  for (int i = 0 ; i < iSeg ; i++)    
    vtxtPrintf (cigar, segs[i]) ;

  
  return lSam ;
} /*  exportOneSamExon */

/*********************************************************************/
/*
  samtools stats titi1/Roche.sam | gawk -f scripts/sam_stats.awk
  bad:          24S16M1D14M3I2I46M1I96M3I10N27M2I182M
  good: 	23S31M3I2I46M1I121M8N3M3I182M
  
*/
static BOOL exportOneSamCigar (BB *bb, vTXT cigar, ALIGN *ap0, int iMax, Array dna, int *nMp, int *nSubp, int *nInsp, int *nDelp, int *nGapp)
{
  int ii ;
  ALIGN *ap = ap0 ;
  int x1, x2 ;
  int dnaMax = dna ? arrayMax (dna) : 0 ;
  BOOL isDown = ap->a1 < ap->a2 ? TRUE : FALSE ;
  int lSam = 0 ;
  
  vtxtClear (cigar) ;
  if (isDown)
    {
      /* gap en tete */
      x1 = ap->x1 ;
      if (x1 > 1)
	{
	  lSam += x1 - 1 ;
	  vtxtPrintf (cigar, "%dS", x1 - 1) ;
	}
      
      for (ii = 0 ; ii < iMax - 1 ; ii++, ap++)
	{
	  lSam += exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;

	  int dx = ap[1].x1 - ap[0].x2 - 1 ;
	  if (dx > 0)
	    {
	      lSam += dx ;
	      vtxtPrintf (cigar, "%dS", dx) ;
	    }
	  if (dx < 0)
	    {
	      ap[1].x1 += dx ;
	      ap[1].a1 += dx ;
	    }

	  int da = ap[1].a1 - ap[0].a2 - 1 ;
	  if (da > 0)
	    {
	      if (da > 20)
		vtxtPrintf (cigar, "%dN", da) ;
	      else
		vtxtPrintf (cigar, "%dD", da) ;
	      *nGapp += da ;
	    }
	}

      lSam += exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
      
      /* gap en queue */
      x2 = ap->x2 ;
      if (x2 < dnaMax)
	{
	  vtxtPrintf (cigar, "%dS", dnaMax - x2) ;
	  lSam += dnaMax - x2 ;
	}
      
    }
  else
    {
      ap += iMax - 1 ;

      /* gap en tete */
      x2 = ap->x2 ;
      if (x2 < dnaMax)
	{
	  lSam += dnaMax - x2 ;
	  vtxtPrintf (cigar, "%dS", dnaMax - x2) ;
	}
      
      for (ii = 0 ; ii < iMax - 1 ; ii++, ap--)
	{
	  lSam += exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
	  int dx = ap[0].x1 - ap[-1].x2 - 1 ;
	  if (dx > 0)
	    {
	      lSam += dx ;
	      vtxtPrintf (cigar, "%dS", dx) ;
	    }
	  if (dx < 0)
	    {
	      ap[-1].x2 += dx ;
	      ap[-1].a2 -= dx ;
	    }

	  int da = ap[-1].a2 - ap[0].a1 - 1 ;
	  if (da > 0)
	    {
	      if (da > 20)
		vtxtPrintf (cigar, "%dN", da) ;
	      else
		vtxtPrintf (cigar, "%dD", da) ;
	      *nGapp += da ;
	    }
	}
      lSam += exportOneSamExon (bb, isDown, cigar, ap, nMp, nSubp, nInsp, nDelp) ;
      
      /* gap en queue */
      x1 = ap->x1 ;
      if (x1 > 1)
	{
	  vtxtPrintf (cigar, "%dS", x1 - 1) ;
	  lSam += x1 - 1 ;
	}
    }

  return lSam ;
} /* exportOneSamCigar */

/*********************************************************************/

static int exportOneSam (ACEOUT ao, ACEOUT aoe, const PP *pp, BB *bb, vTXT record, vTXT cigar, Array cigarettes, Array errors, ALIGN *ap0, long int i0, long int aMax)
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
  int chainA2 = ap->a2 ;
  BOOL pairedEnd = bb->rc.pairedEnd ;
  BOOL mateIsDown = TRUE ;
  BOOL isPrimary = ap->chain == 1 ? TRUE : FALSE ;
  BOOL isSecondary = FALSE ;
  BOOL isSupplementary = FALSE ;
  BOOL isOrphan = FALSE ; 
  Array dna = arr (bb->dnas, read, Array) ;
  int nM = 0, nSub = 0, nIns = 0, nDel = 0, nGap = 0 ;

  vtxtClear (record)  ;
  
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
  vtxtPrintf (record, "%s/%s%s", run, dictName (dict, read >> 1), (pairedEnd ? ((read & 0x1 )? "<" : ">") : "")) ; 
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
       *  diverses  portions du read sont incompatibles
       *  all have this flag except the first one
       */
      if (isSupplementary) flag |= 0x800 ;
    }


  vtxtPrintf (record, "\t%d", flag) ;
  vtxtPrintf (record, "\t%s", dictName (dictG, chrom)) ;

  if (0 && ap->read == 187746 && chrom==10)
    invokeDebugger () ;
  /* a chain [i0, iMax[, is reported as a single CIGAR */
  chainA1 = ap->chainA1 ;
  for (ii = i0, ap = ap0 ; ii < aMax && ap->chain == chain && ap->read == read ; ii++, ap++)
    chainA2 = ap->chainA2 ;
  
  vtxtPrintf (record, "\t%d", chainA1 < chainA2 ? chainA1 : chainA2) ;    
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
      vtxtPrintf (record, "\t%d", mapq) ;
    }
    
  /* CIGAR */
  int lSam = exportOneSamCigar (bb, cigar, ap0, kMax, dna, &nM, &nSub, &nIns, &nDel, &nGap) ;
  int da = arrayMax (dna) - lSam ;
  if (da == 0) /* success */
    vtxtPrintf (record, "\t%s", vtxtPtr (cigar)) ;
  else
    {
      aceOutf (aoe, "ERROR: Cigar error in lSam != dnaMax:%d != %d %s\n>%s\n"
	       , lSam, arrayMax(dna), vtxtPtr (cigar)
	       , dictName (dict, read >> 1)
	       )  ;
      
      if (da > 0 && da < 6)
	{
	  vtxtPrintf (record, "\t%s%dS", vtxtPtr (cigar), da) ;
	  aceOutf (aoe, "PADDED\n") ;
	}
      else
	{
	  Array dnaR = dnaCopy (dna) ;
	  dnaDecodeArray (dnaR) ;
	  aceOutf (aoe, "%s\n", arrp (dnaR, 0, char)) ;
	  ac_free (dnaR) ;
	  goto done ;
	}
    }
  if (ap->mateChrom && ap->mateA1)
    {
      vtxtPrintf (record, "\t%s\t%d\t%d", dictName (dictG, ap->mateChrom >> 1), ap->mateA1, ap->pairLength) ;
    }
  else
    vtxtPrintf (record, "\t*\t0\t0") ;    

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
	  vtxtPrintf (record, "\t%s", arrp (dnaR, 0, char)) ;
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
	  vtxtPrintf (record, "\t%s", arrp (dnaR, 0, char)) ;
	  ac_free (dnaR) ;
	}
    }
  else
    vtxtPrintf (record, "\t*") ;

  if (pp->exportSamQuality && bb->quals)
    {
      Array qual = read < arrayMax (bb->quals) ? arr (bb->quals, read, Array) : 0 ;
      if (qual && arrayExists (qual) && arrayMax(qual))
	{
	  if (isDown)
	    {
	      vtxtPrintf (record, "\t%s", arrp (qual, 0, char)) ;
	    }
	  else
	    {
	      int i = arrayMax (qual) - 1 ;
	      char *cp = arrp (qual, i, char) ;
	      vtxtPrintf (record, "\t") ;
	      while (i-- >= 0)
		vtxtPrintf (record, "%c", *cp--) ;
	    }
	}
      else
	vtxtPrintf (record, "\t*") ;
    }
  else
    vtxtPrintf (record, "\t*") ;


  /* end of the 11 mandatory columns */

    
  vtxtPrintf (record, "\tNH:i:%d", ap->nChains) ; /* px->uu */ ;
  vtxtPrintf (record, "\tAS:i:%d", ap->pairScore ? ap->pairScore : ap->chainScore) ;
  vtxtPrintf (record, "\tNM:i:%d", nSub + nIns + nDel) ;
  vtxtPrintf (record, "\tmt:i:%d", 1) ; /* mult */
  /* vtxtPrintf (record, "\tMD:Z:0") ; */
  
  aceOutf (ao, "%s\n", vtxtPtr (record)) ;

 done:
  ac_free (h) ;
  return kMax ;
} /* exportOneSam */

/**************************************************************/
/* return the number of consumed ap (must be >=1) */
int saSamExport (ACEOUT ao, ACEOUT aoe, const PP *pp, BB *bb)
{
  AC_HANDLE h = ac_new_handle () ;
  ALIGN *ap = bigArrp (bb->aligns, 0, ALIGN) ;
  long int ia, aMax = bigArrayMax (bb->aligns) ;
  vTXT record = vtxtHandleCreate (h) ;
  vTXT cigar = vtxtHandleCreate (h) ;
  Array cigarettes = arrayHandleCreate (1024, SAMCIGAR, h) ;
  Array errors = arrayHandleCreate (1024, A_ERR, h) ;
  int read, chain ;

  for (ia = 0, read = chain = 0 ; ia < aMax ; ia++, ap++)
    {
      int k ;
      if (ap->read != read || ap->chain != chain)
	k = exportOneSam (ao, aoe, pp, bb, record, cigar, cigarettes, errors, ap, ia, aMax) ;
      read = ap->read ;
      chain = ap->chain ;
      ap += k - 1 ;
      ia += k - 1 ;
    }
  ac_free (h) ;
  return 1 ;
} /* exportSamDo */

/**************************************************************/

ACEOUT saSamCreateFile (const PP *pp, BB *bb, BOOL isError, AC_HANDLE h)
{
  char *VERSION = "0.1.1" ;
  DICT *dictG = pp->bbG.dict ;
  
  ACEOUT ao = aceOutCreate (pp->outFileName, hprintf (h, ".%s.sam%s", dictName (pp->runDict, bb->run), isError ? ".error" : "" ), pp->gzo, h) ;
  if (! isError)
    {
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
  return ao ;
} /* saSamFile */

/**************************************************************/
/**************************************************************/
/**************************************************************/

/*
runX/rob1	272	G.chr14.T2T.NC_060938.1	24	1	12M-12N24M

runX/rob1	272	G.chr15.T2T.NC_060939.1	24	1	12M-12N24M
runX/rob1	256	G.chr15.T2T.NC_060939.1	78327069	1	7M29M


runX/rob1	36	1	36	36	1	24	G	-	6	G.chr14.T2T.NC_060938.1	24	1
runX/rob1	36	1	36	36	25	36	G	-	6	G.chr14.T2T.NC_060938.1	12	1


runX/rob1	36	1	36	36	1	24	G	-	6	G.chr15.T2T.NC_060939.1	24	1
runX/rob1	36	1	36	36	25	36	G	-	6	G.chr15.T2T.NC_060939.1	12	1


runX/rob1	36	1	36	36	1	7	G	-	6	G.chr15.T2T.NC_060939.1	78327069	78327075	0	0	-	-	-	-	-	chain 4 1	36
runX/rob1	36	1	36	36	8	36	G	-	6	G.chr15.T2T.NC_060939.1	78327070	78327098	0	0	-	-	-	-	-	chain 4 1	36
*/
