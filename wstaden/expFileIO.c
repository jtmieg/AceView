/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*
 * File: expFileIO.c
 * Version:
 *
 * Description: Routines for reading and writing to experiment files.
 *
 * 1. Opening experiment files
 * 2. Reading information from an experiment file
 * 3. Appending to experiment files
 * 4. Closing an opened experiment file
 *
 * Created:
 * Updated:
 *
 */

/*
 * Tag format:
 *
 * 0        10 
 * |----.----|-
 * TG   TYPE S position..length
 * TG        One or more comment lines starting at character position 10
 * TG        Each line represents a line of tag.
 * TG         Extra indentation is simply added to the comment.
 *
 * Where S is the strand, either "+", "-", or "=" (both).
 * Eg:
 *
 * TG   COMM = 100..110
 * TG        This comment contains
 * TG          several lines.
 *
 * So the above is a COMMent tag on both strands from bases 100 to 110
 * inclusive containing the annotation
 * "This comment contains\n  several lines.\n"
 *
 * This is written using exp_put_str giving the multi line string:
 * "COMM = 100..110\nThis comment contains\n  several lines."
 *
 * (ie the indentation is added by the experiment file format, not by the
 *  calling routines. Similarly this indentation is stripped out again when
 *  reading back.)
 */


#include <stdio.h>
#include <string.h> /* IMPORT: strdup (hopefully!) */
#include <ctype.h>

#include "stadenarray.h"
#include "expFileIO.h"
#include "xalloc.h"

/* Fixup for broken SunOS 4.x systems */
#ifndef FOPEN_MAX
#define FOPEN_MAX 20
#endif

/*************************************************************
 * Line types for experiment file
 *************************************************************/

char eflt_feature_ids[MAXIMUM_EFLTS][MAXIMUM_EFLT_LENGTH_PLUS_1] ;
#ifdef JUNK
= {
    "CF", /*  0 cloning vector sequence file */
    "CN", /*  1 clone name */
    "CS", /*  2 cloning vector sequence present in sequence */
    "CV", /*  3 cloning vector type */
    "DR", /*  4 direction of read */
    "DT", /*  5 date of experiment */
    "EN", /*  6 experiment name */
    "EX", /*  7 experimental notes */
    "FM", /*  8 sequencing vector fragmentation method */
    "LN", /*  9 local format trace file name */
    "LT", /* 10 local format trace file type */
    "MC", /* 11 machine on which experiment ran */
    "MN", /* 12 machine generated trace file name */
    "MT", /* 13 machine generated trace file type */
    "OP", /* 14 operator */
    "PN", /* 15 primer name */
    "QR", /* 16 poor quality sequence present at right (3') end */
    "SC", /* 17 sequencing vector cloning site */
    "SF", /* 18 sequencing vector sequence file */
    "SI", /* 19 sequencing vector insertion length */
    "SL", /* 20 sequencing vector present at left (5') end */
    "SP", /* 21 sequencing vector primer site (relative to cloning site) */
    "SQ", /* 22 sequence */
    "SR", /* 23 sequencing vector present at right (3') end */
    "ST", /* 24 strands */
    "SV", /* 25 sequencing vector type */
    "TN", /* 26 template name */
    "QL", /* 27 poor quality sequence present at left (5') end */
    "PS", /* 28 processing status */
    "CC", /* 29 comments */
    "SS", /* 30 sequence to screen against */
    /* added 27-May-93 */
    "TG", /* 31 gel tag line */
    "ID", /* 32 identifier */
    /* added 24-Sep-93 */
    "AQ", /* 33 average quality measure */
    /* added 15-Oct-93 */
    "PR", /* 34 primer type */
    "LI", /* 35 subclone library (mtd) */
    "LE", /* 36 subclone library entry (well) */
    /* added 19-Apr-94 */
    "TC", /* 37 contig tag line */
    "AC", /* 38 accession number */
    /* added 11-Nov-94 */
    "BC", /* 39 base calling software */
    "ON", /* 40 original base numbers (positions) */
    "AV", /* 41 accuracy (quality) values */
    "PC", /* 42 position in contig */
    "SE", /* 43 sense, whether it is complemented */
    /* added 5-4-95 */
    "CL", /* 44 cloning vector left end*/
    "CR", /* 45 cloning vector right end*/
    "AP", /* 46 assembly position */
    "CH", /* 47 special chemistry used (eg taq) */
    "PD", /* 48 primer data - the sequence of a primer */
    "WT"  /* 49 wild type trace */
};
#endif




/*************************************************************
 * Output/update lines
 *************************************************************/

int exp_print_line_(FILE *fp, char *eflt, char *entry)
/*
 * Output an experiment file line
 */
{
    return fprintf(fp,
		   "%-5s%s\n",
		   eflt,
		   entry
		   )  <  0;
}

int exp_print_line(FILE *fp, Exp_info *e, int eflt, int i)
/*
 * Output an experiment file line
 */
{
    return exp_print_line_(fp,
			   eflt_feature_ids[eflt], 
			   arr(char *,e->entries[eflt],i)
			   );
}

/*
 * Outputs a multi-line experiment file line.
 * Continuation lines are automatically added by adding 5 characters of extra
 * indentation at the start of each continuation.
 *
 * returns -1 for failure, 0 for success.
 */
int exp_print_mline(FILE *fp, Exp_info *e, int eflt, int i) {
    char *p, *c;

    p = arr(char *, e->entries[eflt], i);

    /* first line */
    if ((c = strchr(p, '\n')))
	*c = '\0';
    if (-1 == exp_print_line_(fp, eflt_feature_ids[eflt], p))
	return -1;

    while (c) {
	*c = '\n';
	p = c+1;
	
	if ((c = strchr(p, '\n'))) {
	    *c = '\0';
	}
	
	if (-1 == fprintf(fp, "%-10s%s\n", eflt_feature_ids[eflt], p))
	    return -1;
    }

    return 0;
}


int exp_print_seq(FILE *fp, Exp_info *e, int eflt, int i)
/*
 * Output an experiment file multi line
 */
{
    int j, l;
    char *seq;
    if (fprintf(fp,"%-5s",eflt_feature_ids[eflt])<0) return 1;

    l = strlen(seq = arr(char *,e->entries[eflt],i));
    for(j=0;j<l;j++) {
	if (j%60==0) if ( fprintf(fp,"\n    ") < 0 ) return 1;
	if (j%10==0) if ( fprintf(fp," ") < 0 ) return 1;
	if ( fprintf(fp,"%c",seq[j]) < 0 ) return 1;
    }
    if ( fprintf(fp,"\n//\n") < 0 ) return 1;

    return 0;
}

int exp_get_feature_index(char *e)
{
    int i;
    
    for (i = 0; i < MAXIMUM_EFLTS; i++)
	if(strcmp(eflt_feature_ids[i],e)==0) return i;
    
    return -1;
}


/*************************************************************
 * Utility routines
 *************************************************************/

/*
 * Creates a string of 'range format' from the start and end points.
 * The string (of form start..end) is also returned.
 */
char *exp_create_range(char *str, int start, int end) {
    sprintf(str, "%d..%d", start, end);
    return str;
}

/*
 * Extracts the start and end points from a range string.
 * Returns 0 for success and -1 for failure.
 */
int exp_extract_range(char *str, int *start, int *end) {
    return sscanf(str, "%d..%d", start, end) != 2;
}

Exp_info *exp_create_info(void)
/*
 * Allocate space for new experiment file information
 */
{
    Exp_info *new;
    int i;
    
    new = (Exp_info *)xmalloc(sizeof(Exp_info));
    if (new != NULL) {
	for(i=0; i< MAXIMUM_EFLTS ; i++) {
	    new->Nentries[i] = 0;
	    new->entries[i] = ArrayCreate(sizeof(char *), 1/*one entry*/);
	}
	new->fp = NULL;
    }
    
    return new;
}


void exp_destroy_info(Exp_info *e)
/*
 * Destroy experiment file information
 */
{
    int i;
    int j;
    if (e != NULL_Exp_info) {
	for (i = 0; i < MAXIMUM_EFLTS; i++) {
	    Array a = e->entries[i];
	    for(j=0;j<e->Nentries[i];j++)
		if (arr(char *,a,j) != NULL) xfree(arr(char *,a,j));
	    ArrayDestroy(a);
	}
	if (e->fp != NULL) fclose(e->fp);
	xfree(e);
    }
}






static char *exp_read_sequence(FILE *fp)
/*
 * Read from file a sequence, discarding all white space til a // is encountered
 */
{
    char *seq;
    int seql;
    char line[EXP_FILE_LINE_LENGTH+1];
    char *l;
    
    seql = 0;
    seq = (char *)xmalloc(seql+1);
    if (NULL == seq)
	return NULL;
    seq[0] = '\0';
    
    l = fgets(line,EXP_FILE_LINE_LENGTH,fp);
    while (l!= NULL && strncmp(l,"//",2)) {
	char *a, *b;
	for(a=b=line;*a;a++)
	    if (! isspace((int)*a)) *b++=*a;
	*b = '\0';
	seql = seql + b-line;
	seq = (char *)xrealloc(seq,seql+1);
	if (NULL == seq)
	    return NULL;
	strcat(seq,line);
	l = fgets(line,EXP_FILE_LINE_LENGTH,fp);
    }
    
    return seq;
}


/*
 * Converts the opos[] array into a char array.
 * In doing so this shrinks the data size by using a .. notation.
 * No check is made that buf is large enough. It is recommended that buf is
 * allocated to 5*len which covers the worst case (for sequences less that
 * 9999 bases long).
 *
 * Note that on older systems sprintf may return the first argument rather
 * than the number of characters written.
 * For this reason we have to do the counting ourselves.
 */
char *opos2str(int2 *opos, int len, char *buf) {
    int i, st, f, dir = 0;
    char *r = buf, *rs = buf;
    
    f = opos[st = 0];
    for (i = 1; i < len; f=opos[i++]) {
	if (dir == 0)
	  {
	    if (opos[i] == f+1)
	      dir=1;
	    else if (opos[i] == f-1)
	      dir=-1;
	  }
	if (dir && opos[i] != f + dir) {
	    if (st != i-1)
		sprintf(buf, "%d..%d ", opos[st], opos[i-1]);
	    else
		sprintf(buf, "%d ", opos[st]);
	    st = i;
	    dir = 0;

	    buf += strlen(buf);
		
	} else if (dir == 0) {
	    sprintf(buf, "%d ", f);

	    st = i;
	    buf += strlen(buf);
	}

	if (buf - rs > 60) {
	    *buf++ = '\n';
	    *buf = '\0';
	    rs = buf - 6;
	}	
    }
    
    if (st != i-1)
	sprintf(buf, "%d..%d", opos[st], opos[i-1]);
    else
	sprintf(buf, "%d", opos[st]);
    
    return r;
}


/*
 * Expands from the character string .. notation to the opos[] array, up to
 * a maximum of len elements in opos[].
 *
 * Returns the length of the opos array.
 */
int str2opos(int2 *opos, int len, char *buf) {
    int i, n1, n2, st, en, m, j = 0;

    while (*buf) {
	m = sscanf(buf, "%d%n..%d%n", &st, &n1, &en, &n2);

	if (m == 1) {
	    opos[j++] = st;
	    buf += n1;
	    if (j >= len)
		break;
	} else if (m == 2) {
	    if (en >= st)
		for (i = st; i <= en && j < len; i++)
		    opos[j++] = i;
	    else
		for (i = st; i >= en && j < len; i--)
		    opos[j++] = i;
	    buf += n2;
	    if (j >= len)
		break;
	} else {
	    buf++;
	}
    }

    return j;
}


/*
 * Converts the accuracy value string (AV) to the confidence array up to
 * a maximum of len elements in conf[].
 *
 * The AV string is of format:
 * "x y z ..." where x, y and z are confidence values for the first three
 * called bases. Or:
 * "a,b,c,d e,f,g,h i,j,k,l ..." where the 4-tuples represent the four
 * confidence values for each base.
 *
 * Returns: number of confidence values read, or -1 for error.
 */
int str2conf(int1 *conf, int len, char *buf) {
    int val1, ind = 0, pos;

    while (*buf && sscanf(buf, "%d%n", &val1, &pos) == 1 && ind < len) {
	if (buf[pos] == ',') {
	    fprintf(stderr, "4-tuple system is currently unsupported\n");
	    return -1;
	}

	conf[ind++] = val1;
	buf += pos;
    }

    return ind;
}

/*
 * Converts the confidence array to the accuracy value string (AV).
 *
 * Note no memory overrun checks are performed on buf. It is recommended
 * that it is allocated to 4*len (worst case of "100 " for each base).
 *
 * Returns the buf argument.
 */
char *conf2str(int1 *conf, int len, char *buf) {
    int i;
    char *ret = buf, *rs = buf;

    for (i = 0; i < len; i++) {
	sprintf(buf, "%d ", conf[i]);
	buf += strlen(buf);

	if (buf - rs > 60) {
	    *buf++ = '\n';
	    *buf = '\0';
	    rs = buf - 6;
	}
    }

    return ret;
}

/*************************************************************
 * Main C interface routines
 *************************************************************/


/*
 * Closes an experiment file (if open), but does not free it.
 */
void exp_close(Exp_info *e) {
    if (e->fp) {
	fclose(e->fp);
	e->fp = NULL;
    }
}


Exp_info *exp_read_info(char *file)
/*
 * Read in an experiment file and return handle
 */
{
    Exp_info *e;
    FILE *fp;
    
    /*
     * open for read
     */
    if ((fp = fopen(file,"r"))==NULL) {
	return NULL_Exp_info;
    }

    e = exp_fread_info(fp);
    fclose(fp);

    if (NULL_Exp_info == e) {
	return NULL_Exp_info;
    }	

    /*
     * reopen for appending
     */
    e->fp = fopen(file,"a");
    
    return e;
    
}


/*
 * Read in an experiment file and return handle
 */
Exp_info *exp_fread_info(FILE *fp)
{
    Exp_info *e;
    char line[EXP_FILE_LINE_LENGTH+1];
    
    e = exp_create_info();
    
    /*
     * open for read
     */
    e->fp = fp;
    
    if (e != NULL_Exp_info) {
	while (fgets(line,EXP_FILE_LINE_LENGTH,e->fp) != NULL) {
	    char *c;
	    int entry;

	    /*
	     * zero terminate first argument
	     * set c to point to second argument
	     *
	     * FIXME: c should point to character 6 always. Indentation is
	     * important when considering continuation lines.
	     */
	    for (c=line;*c && !isspace((int)*c); c++) ;
	    if (*c) {
		*c++ = '\0';
		for (;*c && isspace((int)*c); c++) ;
	    }
	    
	    entry = exp_get_feature_index(line);
	    if (entry >= 0) {
		/*
		 * Tag lines may be split over multiple lines. If we have no
		 * tag type then we append to the existing tag.
		 */
		if ((int)(c-line) >= 10/* continuation lines */
		    && (entry == EFLT_TG || entry == EFLT_TC ||
			entry == EFLT_ON || entry == EFLT_AV)) {
		    char *en;
		    size_t l1, l2;

		    /*
		     * Extend our current line by the appropriate amount
		     */
		    en = exp_get_entry(e,entry);
		    l1 = strlen(en);
		    l2 = strlen(&line[10]);

		    if (NULL == (en = exp_get_entry(e, entry) =
				 (char *)xrealloc(en, l1 + l2 + 1)))
			return NULL;
		    

		    /*
		     * Append the new line (without the \n char)
		     */
		    en[l1] = '\n';
		    line[l2+9] = '\0';
		    strcpy(&en[l1+1], &line[10]);
		} else {
		    /*
		     * Increment number of entries for line type entry
		     * This will force exp_get_entry() to return pointer to
		     * next free element in array
		     */
		    (void)ArrayRef(e->entries[entry],e->Nentries[entry]++);
		    
		    if (entry == EFLT_SQ)
			exp_get_entry(e,entry) = exp_read_sequence(e->fp);
		    else {
			char *eoln = strchr(c,'\n');
			int i;

			if (eoln!=NULL) *eoln='\0';

			if (entry == EFLT_LT)
			    for (i=3; isspace((int)c[i]) && i >= 0; c[i--]='\0');

			exp_get_entry(e,entry) = (char *)strdup(c);
		    }
		}
	    }
	}
    }
    
    return e;
}

static int exp_check_eid_read(Exp_info *e,int id)
/*
 * Check these are a valid combination and that
 * an entry exists for read
 */
{
    return (
	    e == NULL ||
	    id < 0 ||
	    id >= MAXIMUM_EFLTS ||
	    e->Nentries[id] == 0 ||
	    eflt_feature_ids[id][0]=='\0'
	    );
}

static int exp_check_eid_write(Exp_info *e,int id)
/*
 * Check these are a valid combination and that
 * an entry exists for write
 */
{
    return (e == NULL ||
	    id < 0 ||
	    id >= MAXIMUM_EFLTS ||
	    e->fp == NULL ||
	    eflt_feature_ids[id][0]=='\0');
}






int exp_get_int(Exp_info *e, int id, int *val)
/*
 * Get the integer for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */
{
    if ( exp_check_eid_read(e,id) ) return 1;
    *val = atoi(exp_get_entry(e,id));
    return 0;
}


int exp_get_rng(Exp_info *e, int id, int *from, int *to)
/*
 * Get the integer pair for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */
{
    if ( exp_check_eid_read(e,id) ) return 1;
    (void)exp_extract_range(exp_get_entry(e,id), from, to);

    return 0;
}



int exp_get_str(Exp_info *e, int id, char *s, f_implicit s_l)
/*
 * Get the string for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */
{
    if ( exp_check_eid_read(e,id) ) return 1;
    strncpy(s,exp_get_entry(e,id),s_l);
    
    return 0;
}


static int exp_append_str(Exp_info *e, int id, char *s, int len)
/*
 * Append the string to experiment file for entry id
 * returns:
 *    0 - success
 *    1 - no update
 */
{
    (void)ArrayRef(e->entries[id],e->Nentries[id]++);
    exp_get_entry(e,id) = (char *)xmalloc(len+1);
    strncpy(exp_get_entry(e,id), s, len);
    exp_get_entry(e,id)[len] = '\0';
    
    if ( id == EFLT_SQ )
	return exp_print_seq(e->fp,e,id,e->Nentries[id]-1);
    else if (id == EFLT_TG || id == EFLT_TC ||
	     id == EFLT_ON || id == EFLT_AV)
	return exp_print_mline(e->fp,e,id,e->Nentries[id]-1);
    else
	return exp_print_line(e->fp,e,id,e->Nentries[id]-1);
}


int exp_put_int(Exp_info *e, int id, int *val)
/*
 * Append the integer for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */
{
    char buf[EXP_FILE_LINE_LENGTH];
    if ( exp_check_eid_write(e,id) ) return 1;
    sprintf(buf,"%d",*val);
    return exp_append_str(e,id,buf,strlen(buf));
}


int exp_put_rng(Exp_info *e, int id, int *from, int *to)
/*
 * Append the integer pair for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */
{
    char buf[EXP_FILE_LINE_LENGTH];
    if ( exp_check_eid_write(e,id) ) return 1;

    (void )exp_create_range(buf, *from, *to);

    return exp_append_str(e,id,buf,strlen(buf));
}



int exp_put_str(Exp_info *e, int id, char *s, f_implicit s_l)
/*
 * Append the string for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */
{
    if ( exp_check_eid_write(e,id) ) return 1;
    return exp_append_str(e,id,s,s_l);
}


/*
 * FORTRAN INTERFACE
 */

static void c2fstr(char *c, int max_c, char *f, int max_f)
{
#define min(A,B) ( (A) < (B) ? (A) : (B) )
    int i;
    i = strlen(c);
    i = min(i,max_f);
    i = min(i,max_c);
    strncpy(f,c,i);
    for( ; i<max_f; i++) f[i]=' ';
    
}

static int fstrlen(char *f, int max_f)
{
    for (; max_f > 0 && (isspace((int)f[max_f-1]) || f[max_f-1]=='\0'); max_f--);
    return max_f;
}

static void f2cstr(char *f, int max_f, char *c, int max_c)
{
    int i;
    
    i = min(fstrlen(f,max_f),max_c);
    strncpy(c,f,i);
    c[i]='\0';
}





/*************************************************************
 * FORTRAN INTERFACE
 *************************************************************/



static int init_done = 0;
static int NHandles = 0;
static Exp_info **Handles = NULL;

static int initialise(void)
{
    int i;
    
    if (init_done) return 0;
    init_done++;
    
    NHandles = FOPEN_MAX;
    
    if ( (Handles = (Exp_info **)xmalloc(sizeof(Exp_info *) * NHandles)) == NULL) {
	NHandles = 0;
	return 1;
    }
    
    for (i=0; i<NHandles; i++) Handles[i] = NULL;
    
    return 0;
}


static int get_free_handle(void)
/*
 * find a free entry in the Exp array
 * returns -1 if there is none
 */
{
    int i;
    
    (void) initialise();
    
    if (!NHandles) return -1; /* no slots! */
    for (i=0; i<NHandles && Handles[i]!=NULL; i++) ;
    return (i==NHandles)?-1:i;
}


static int check_handle(f_int *handle)
{
    return (handle == NULL ||
	    (int) (*handle) <= 0 ||
	    (int) (*handle) > NHandles);
}



f_int expopn_(char *fn, f_implicit fn_l)
/*
 * FORTRAN interface to exp_open_file()
 */
{
    char cfn[1025];
    int handle;
    
    if ( (handle = get_free_handle()) >= 0 ) {
	f2cstr(fn,fn_l,cfn,1024);
	Handles[handle] = exp_read_info(cfn);
    }
    
    return (f_int) (handle+1);
}



f_proc_ret expkil_(f_int *handle)
/*
 * FORTRAN interface to exp_destroy_info
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) f_proc_return();
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    exp_destroy_info(e);
    
    Handles[(int)(*handle)-1] = NULL;
    *handle = 0;
    
    f_proc_return();
}

f_int expri_(f_int *handle, f_int *id, f_int *val)
/*
 * FORTRAN interface to exp_get_int
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    return exp_get_int(e, (int)*id, (int *)val);
}


f_int exprr_(f_int *handle, f_int *id, f_int *from, f_int *to)
/*
 * FORTRAN interface to exp_get_rng
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    return exp_get_rng(e,(int)*id,(int *)from,(int *)to);
    
}

/* ARGSUSED */
f_int exprsa_(f_int *handle, f_int *id, char *s, f_int *max_len, f_implicit s_l)
/*
 * FORTRAN interface to exp_get_str workalike
 * NOTE: for use with FORTRAN CHARACTER arrays instead CHARACTER strings
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    if ( exp_check_eid_read(e,*id) ) return 1;
    c2fstr(exp_get_entry(e,*id),(int)*max_len,s,(int)*max_len);
    return 0;
}


f_int exprs_(f_int *handle, f_int *id, char *s, f_implicit s_l)
/*
 * FORTRAN interface to exp_get_str workalike
 * NOTE: for use with FORTRAN CHARACTER strings instead CHARACTER arrays
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    if ( exp_check_eid_read(e,*id) ) return 1;
    c2fstr(exp_get_entry(e,*id),s_l,s,s_l);
    return 0;
}


f_int expwi_(f_int *handle, f_int *id, f_int *val)
/*
 * FORTRAN interface to exp_put_int
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    return exp_put_int(e, (int)*id, (int *)val);
}


f_int expwr_(f_int *handle, f_int *id, f_int *from, f_int *to)
/*
 * FORTRAN interface to exp_put_rng
 */
{
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    return exp_put_rng(e, (int)*id, (int *)from, (int *)to);
}


/* ARGSUSED */
f_int expwsa_(f_int *handle, f_int *id, char *s, f_int *max_len, f_implicit s_l)
/*
 * FORTRAN interface to exp_put_str workalike
 * NOTE: for use with FORTRAN CHARACTER arrays instead CHARACTER strings
 */
{
    Exp_info *e;
    char buf[EXP_FILE_LINE_LENGTH];
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    
    if ( exp_check_eid_write(e,*id) ) return 1;
    /* don't allow multi-line entries to be written */
    if (*id == EFLT_SQ ) return 1;
    f2cstr(s,(int)*max_len,buf,sizeof(buf));
    return exp_append_str(e,*id,buf,strlen(buf));
    
}

f_int expws_(f_int *handle, f_int *id, char *s, f_implicit s_l)
/*
 * FORTRAN interface to exp_put_str workalike
 * NOTE: for use with FORTRAN CHARACTER strings instead CHARACTER arrays
 */
{
    char buf[EXP_FILE_LINE_LENGTH];
    Exp_info *e;
    if ( check_handle(handle) ) return 1;
    e = (Exp_info *) Handles[(int)(*handle)-1];
    
    
    if ( exp_check_eid_write(e,*id) ) return 1;
    /* don't allow multi-line entries to be written */
    if (*id == EFLT_SQ ) return 1;
    f2cstr(s,s_l,buf,sizeof(buf));
    return exp_append_str(e,*id,buf,s_l);
}

/*
 * FORTRAN interface to exp_create_range()
 */
void expcr_(char *str, f_int *start, f_int *end, f_implicit str_l) {
    exp_create_range(str, *start, *end);
    c2fstr(str, str_l, str, str_l);

    f_proc_return();
}

/*
 * FORTRAN interface to exp_extract_range()
 */
/* ARGSUSED */
f_int exper_(char *str, f_int *start, f_int *end, f_implicit str_l) {
    return exp_extract_range(str, start, end);
}




/*************************************************************
 * Go for it!
 *************************************************************/

static void print_line(FILE *fp, Exp_info *e, int eflt, int all)
{
    if (all) {
	int i;
	for(i=0;i<e->Nentries[eflt];i++) exp_print_line(fp,e,eflt,i);
    } else if (e->Nentries[eflt] > 0) {
	exp_print_line(fp,e,eflt,e->Nentries[eflt]-1);
    }
}


static void print_mline(FILE *fp, Exp_info *e, int eflt, int all)
{
    if (all) {
	int i;
	for(i=0;i<e->Nentries[eflt];i++) exp_print_mline(fp,e,eflt,i);
    } else if (e->Nentries[eflt] > 0) {
	exp_print_mline(fp,e,eflt,e->Nentries[eflt]-1);
    }
}



static void print_seq(FILE *fp, Exp_info *e, int eflt)
{
    if (e->Nentries[eflt] > 0)
        exp_print_seq(fp,e,eflt,e->Nentries[eflt]-1);
}




void exp_print_file(FILE *fp, Exp_info *e)
{
    print_line(fp,e,EFLT_ID, 0);
    print_line(fp,e,EFLT_AC, 0);
    print_line(fp,e,EFLT_EN, 0);

    print_line(fp,e,EFLT_CC, 1);
    print_line(fp,e,EFLT_EX, 1);
    print_line(fp,e,EFLT_PS, 1);

    print_line(fp,e,EFLT_LN, 0);
    print_line(fp,e,EFLT_LT, 0);

    print_line(fp,e,EFLT_CF, 0);
    print_line(fp,e,EFLT_CV, 0);
    print_line(fp,e,EFLT_CS, 0);
    print_line(fp,e,EFLT_CL, 0);
    print_line(fp,e,EFLT_CR, 0);

    print_line(fp,e,EFLT_SF, 0);
    print_line(fp,e,EFLT_SV, 0);
    print_line(fp,e,EFLT_SI, 0);
    print_line(fp,e,EFLT_SC, 0);
    print_line(fp,e,EFLT_SP, 0);
    print_line(fp,e,EFLT_PD, 0);
    print_line(fp,e,EFLT_FM, 0);
    print_line(fp,e,EFLT_SL, 0);
    print_line(fp,e,EFLT_SR, 0);

    print_line(fp,e,EFLT_QL, 0);
    print_line(fp,e,EFLT_QR, 0);

    print_mline(fp,e,EFLT_TG,1);
    print_mline(fp,e,EFLT_TC,1);

    print_line(fp,e,EFLT_CN, 0);
    print_line(fp,e,EFLT_TN, 0);
    print_line(fp,e,EFLT_PN, 0);
    print_line(fp,e,EFLT_PR, 0);
    print_line(fp,e,EFLT_LI, 0);
    print_line(fp,e,EFLT_LE, 0);
    print_line(fp,e,EFLT_CH, 0);

    print_mline(fp,e,EFLT_ON,0);
    print_line(fp,e,EFLT_AQ, 0);
    print_mline(fp,e,EFLT_AV,0);

    print_line(fp,e,EFLT_DR, 0);
    print_line(fp,e,EFLT_SE, 0);
    print_line(fp,e,EFLT_PC, 0);
    print_line(fp,e,EFLT_AP, 0);
    print_line(fp,e,EFLT_ST, 0);

    print_line(fp,e,EFLT_DT, 0);
    print_line(fp,e,EFLT_MC, 0);
    print_line(fp,e,EFLT_MN, 0);
    print_line(fp,e,EFLT_MT, 0);
    print_line(fp,e,EFLT_OP, 1);
    print_line(fp,e,EFLT_BC, 0);
    print_line(fp,e,EFLT_SS, 0);

    print_seq (fp,e,EFLT_SQ);
}


/*
 * Allocate an set a new experiment file entry
 */
char *exp_set_entry(Exp_info *e, int eflt, char *str) {
    char *s;
    size_t l;

    if (NULL == ArrayRef(e->entries[eflt], e->Nentries[eflt]))
	return NULL;
    else
	e->Nentries[eflt]++;

    l = strlen(str);
    if (NULL == (s = exp_get_entry(e, eflt) = (char *)xmalloc(l+1))) {
	e->Nentries[eflt]--;
	return NULL;
    }
    strcpy(s, str);

    return s;
}
