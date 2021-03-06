/*  File: gnbkclientlib.c
 *  Author: Jean Thierry-Mieg (mieg@kaa.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * I started from a sample code generated by rpcgen on Solaris
 * and a first version by Peter Kocab.
 * Does not require any ACEDB library code.
 * 
 * Exported functions:
	 openServer()
	 closeServer()
	 askServer()
	 askServerBinary()
 * HISTORY:
 * Last edited: Jan 22 15:55 1996 (mieg)
 * Created: Wed Nov 25 20:02:45 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: gnbkclientlib.c,v 1.1.1.1 2002/07/19 20:23:18 sienkiew Exp $ */

#include "regular.h"
#define __malloc_h  
#include <errno.h>
#include <rpc/rpc.h>
#include "rpcgnbk.h"

/*************************************************************
Open RPC connection to server
INPUT
 char *host    hostname running server 
 int  timeOut  maximum peroid to wait for answer

OUTPUT
 return value:
 gnbk_handle *  pointer to structure containing open connection
               and client identification information
*/
static CLIENT *clnt = 0 ;


void tStatus (void) { return ; } /* missing in the client */

BOOL openServer(char *host, u_long rpc_port)
{
  struct timeval tv;

/* open rpc connection */
  clnt = clnt_create(host, RPC_GNBK, RPC_GNBK_VERS, "tcp");
  if (!clnt) return FALSE ;
  tv.tv_sec = 600 ;
  tv.tv_usec = 0;
  clnt_control(clnt, CLSET_TIMEOUT, (char *)&tv);

  return TRUE ;
}

/*************************************************************
transfer request to server, and wait for binary answer
INPUT
 char * request  string containing request
 unsigned char ** answer  ptr to char ptr, that has to be filled with answer
 clnt     pointer to structure containing open connection
OUTPUT
 unsigned char ** answer  ptr to char ptr. Pointing to allocated memory containing 
                 answer string. This memory must be XDRfreed
 return value:
 int      error condition
  0   :  no error.
  else:  a server generated error 
*/

int askServer(char *request, char **answerPtr)
{ char *loop ;
  gnbk_data question ;
  gnbk_reponse *reponse = 0 ;
  int length ;

  /* opens connection */
  
/* generate question structure */
  question.reponse.reponse_len = 0;
  question.reponse.reponse_val = "";
  question.question = request;
  reponse = gnbk_server_1(&question, clnt);
 
  if (!reponse)  /* no data was received, return error */
    return 1 ;
  
  /* answer received. allocate memory and fill with answer */
  loop = (char *) reponse->gnbk_reponse_u.res_data.reponse.reponse_val;
  length = reponse->gnbk_reponse_u.res_data.reponse.reponse_len;
  *answerPtr = length ? strnew (loop, 0) : 0 ;
    
  xdr_free((xdrproc_t )xdr_gnbk_reponse, (char *)reponse);
  return 0 ;
}

/*****************************************/
/************** end of file **************/
