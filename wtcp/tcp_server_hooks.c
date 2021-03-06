#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>


/*
 * We need at least some of these.  It is easier to just include them
 * all than to figure out which.
 */
#include "acedb.h"
#include "sysclass.h"
#include "session.h"
#include "pick.h"
#include "bs.h"
#include "systags.h"
#include "tags.h"
#include "call.h"
#include "mytime.h"
#include "aceclient.h"
#include "parse.h"
#include "freeout.h"
#include "version.h"
#include "command.h"
#include "keyset.h"
#include "query.h"
#include "a.h"
#include "dump.h"
#include "lex.h"
#include <wtcp/tcp_connect.h>

/* tracking clients */
static int next_client_serial_number = 1;
static int current_active_clients = 0;
static int current_bad_clients = 0;

extern BOOL TCPWEBSERVER ;
extern FILE *tcp_log_file ;


/*
 * dump function for debugging
 */

void ac_dump(unsigned char *cp, int n)
{
  int x;
  int addr;
  addr = 0;
  while (n > 0)
    {
      printf("%03x %4d  ",addr,addr);
      for (x=0; x< 16; x++)
	{
	  if (x >= n)
	    break;
	  if (x % 4 == 0)
	    printf(" ");
	  printf("%02x ",cp[x]);
	}
      printf("\n          ");
      for (x=0; x< 16; x++)
	{
	  if (x >= n)
	    break;
	  if (x % 4 == 0)
	    printf(" ");
	  printf(" %c ",isprint(cp[x]) ? cp[x] : '.' );
	}
      printf("\n");

      n -= 16;
      cp += 16;
      addr += 16;
    }
}

void process_incoming(char *s, int n, void *cookie);

/*
 * include code that implements acetcp transactions
 *
 * The included functions and variables are all static.  They are
 * used once here and once in the acclient library.
 */

/***********************************************************************
 * generating pseudo-random challenge of the challenge/response pair
 *
 * not much here at the moment - challenge/response auth is not implemented yet
 */

static void init_random()
{
}

static char *generate_challenge(void)
{
  char b[100];
  static int n = 0;
  sprintf(b,"%d-challenge-%d",current_active_clients,n);
  n++;
  return strnew(b, 0);
}


/***********************************************************************
 * access control by source IP of the connection.  There is a list of
 * IP numbers and matching access type.
 *
 * Sorry for the ugliness - I have other things to do, so this is just
 * exactly as neat and clean as necessary and no more.
 */

/*
 * client access types
 *
 */
#define CLIENT_READ		1
#define CLIENT_WRITE		2
#define CLIENT_SHUTDOWN		4


/*
 * a list of ip addresses allowed access.  Currently only works
 * for IPV4 systems.
 */

struct tcp_access_entry
{
  /*
   * we know int is 4 bytes, 
   */
  unsigned char address[4];
  unsigned char mask[4];
  char type;
  struct tcp_access_entry *next;
};

static struct tcp_access_entry *access_list = 0;
static struct tcp_access_entry *access_list_end = 0;


/*
 * parse a single line of access control; add it to the list
 */

static void read_tcp_access_entry(char *s)
{
  struct tcp_access_entry *t;
  char *s1;
  int n;
  BOOL debug = FALSE ;

  while (isspace((int)*s))
    s++;

  t = halloc(sizeof(struct tcp_access_entry), 0);
  n = 0;
  while (n < 4 )
    {
      if (*s == '*')
	{
	  t->address[n] = 0;
	  t->mask[n] = 0;
	}
      else
	{
	  t->address[n] = atoi(s);
	  t->mask[n] = 0xff;
	}
      n++;
      s1 = strchr(s,'.');
      if (! s1 )
	{
	  while (! isspace((int)*s))
	    s++;
	  break;
	}
      s = s1;
      s++;
    }

  if (n != 4)
    return;

  while (isspace((int)*s))
    s++;

  t->type = 0;
  while (*s)
    {
      switch (*s)
	{
	case 'r':
	  t->type |= CLIENT_READ;
	  break;
	case 'w':
	  t->type |= CLIENT_WRITE | CLIENT_READ;
	  break;
	case 'u':
	  t->type |= CLIENT_SHUTDOWN;
	  break;
	case 'x':
	  /*
	   * x by itself means no access
	   */
	  break;
	}
      s++;
    }

  if (debug)
    printf("access is %x\n",t->type);

  t->next = 0;

  if (access_list_end)
    {
      access_list_end->next = t;
      access_list_end = t;
    }
  else
    access_list = access_list_end = t;
}

/*
 * open the access control list file.  read all the entries into memory.
 */

static void read_tcp_access_list()
{
  FILE *f;
  char b[1000];
  f = filopen("wspec/acetcp_access","wrm","r");
  if (!f)
    {
      messdump("no file acetcp_access - accepting read connections from local host");
      read_tcp_access_entry("127.0.0.1 rw");
      return;
    }
  while (fgets(b,sizeof(b),f))
    if ( (b[0] != '#') && (!isspace((int)b[0])) )
      read_tcp_access_entry(b);
  fclose(f);
  /*
   * always list access for localhost.  It seems a reasonable default, but
   * if you do not want localhost to have access, you can just list it
   * in the file as "127.0.0.1 x".
   */
  read_tcp_access_entry("127.0.0.1 rw");
}

/*
 * find what access a particular host has
 */
static int check_tcp_access_list( unsigned char * there )
{
  static struct tcp_access_entry *t;
  for (t = access_list; t; t=t->next)
    {
      int n;
      for (n=0; n<4; n++)
	{
	  if ((there[n] & t->mask[n]) != t->address[n])
	    goto not_it;
	}
      return t->type;
    not_it: ;
    }
  return -1;
}


/***********************************************************************
 * tracking clients
 *
 */


/*
 * There is a struct client for each connection
 */
struct client
{
  int	fd;
  /*
   * file descriptor of the socket to the client
   */

  int	serial_number;
  /*
   * what serial number of this connection
   */

  int	last_transaction_time;
  /*
   * time of the last client transaction
   */
  int	access_type;
  /*
   * access type granted to this client
   *	'r'	may read database
   *	'w'	may read or write database
   */

  int isAttack ;
  int attackTimeOut ;
  char *challenge;
  /*
   * challenge used for authenticating requests
   * from this client.
   */

  
  AceCommand	look;
  /*
   * the magic COMMAND_LOOK structure that a server needs.
   */

  int	saved_option;
  /*
   * the "option" value saved at the end of the last
   * transaction.  In the event of a "encore" command,
   * this option has some meaning about what the command
   * to be "encored" was.
   */

  struct incoming_transaction_rec *incoming_trans;
  struct outgoing_transaction_rec *outgoing_trans;
  /*
   * acetcp transaction state blocks
   */
};

/*
 * a pointer to the struct client for the connection is stored in the
 * array clients.  The index is the socket number.
 */
static Array clients;

/*
 * the max fd and the fd_set used to select on client connections.
 * also includes the listen socket.
 */
static int tcp_server_fd_count;
static fd_set open_sockets;


/*
 * the listen socket is ready.  accept the connection and make a
 * client record for it.
 */


void accept_client (int listen_fd)
{
  struct sockaddr_in there;
  int s;
  struct client *c;
  int good = 0 ;
  int isAttack = 0 ;
  int one ;
  int access_type = 0 ;
  socklen_t n ;
  static int pass = 0, attackers[128], gooddies[128], goodAccessType[128] ;
  extern int tcp_nAttack ;
  int clientAddr = 0 ;

  if (! tcp_log_file)
    tcp_log_file = stderr ;

  n = sizeof(there);
  s = accept(listen_fd, (struct sockaddr *) &there, &n);

  if (s < 0)
    {
      perror("accept");
      return;
    }

  /* printf("accept client %d\n",s); */
  if (! pass++)
    {
      memset (gooddies, 0, sizeof (gooddies)) ;
      memset (goodAccessType, 0, sizeof (goodAccessType)) ;
      memset (attackers, 0, sizeof (attackers)) ;
    }

  memcpy (&clientAddr,  &there.sin_addr, 4) ;
  if (1)
    {
      int i, *gp ;

      for (i = 0, gp = gooddies ; *gp && !good ; i++, gp++)
	if (clientAddr == *gooddies)
	  { good = TRUE ; access_type =  goodAccessType[i] ; }
    }
  if (! good)
    {
      int *ap ;
    
      for (ap = attackers ; *ap ; ap++)
	if (clientAddr == *ap)
	  { 
	    /* Normally i could
	     *       shutdown(s, 2);
	     * but i do not
	     * I keep the attacker hanging till I time him out: 180s
	     */
	    current_bad_clients++ ;
	    tcp_nAttack++ ;
	    isAttack = 1 ;
	  }
    }
  
  
  if (isAttack)
    access_type = -1 ;
  else
    {
      if (! good)
	access_type = check_tcp_access_list( (unsigned char *) & there.sin_addr );
      
      if (access_type < 0) /* accept and reply politely once */
	{ 
	  int i, *ap ;
	  for (i = 0, ap = attackers ; i < 125 && *ap ; i++, ap++) ;
	  attackers[i] = clientAddr ;
	}
      else if (! good)
	{ 
	  int i, *gp ;
	  
	  for (i = 0, gp = gooddies ; i < 125 && *gp ; i++, gp++) ;
	  gooddies[i] = clientAddr ;
	  goodAccessType[i] = access_type ;
	  if (i < 127)
	    fprintf(tcp_log_file, "\n%s Accepting New IP - source IP %d.%d.%d.%d\n"
		    , timeShowNow()
		    , ((unsigned char *) & there.sin_addr)[0]
		    , ((unsigned char *) & there.sin_addr)[1]
		    , ((unsigned char *) & there.sin_addr)[2]
		    , ((unsigned char *) & there.sin_addr)[3]
		    ) ;
	}
    }

  /*
   * it is important to set TCP_NODELAY on the socket.  if you
   * do not, transactions will be VERY slow.  Without TCP_NODELAY,
   * the system hopes you will write more data to the socket
   * very soon, so it does NOT transmit the data immediately -- it
   * wants to put the results of your next write in the same
   * packet.
   *
   * 6 is what getprotobyname returns for "tcp" - technically, I'm 
   * cheating but TCP is _always_ protocol number 6 because that
   * is the IP protocol number for TCP.
   */
  one = 1;
  if (setsockopt( s, 6, TCP_NODELAY , &one, sizeof(int)) < 0)
    {
      /*
       * this setsockopt can never fail, so it is ok to crash if it does.
       */
      messcrash("could not set TCP_NODELAY");
    }

  /*
   * make a client record with incoming and outgoing transaction
   * states and other basic initialization.
   */
  c = halloc(sizeof(struct client), 0);
  c->fd = s;
  if (isAttack)
    {
      c->isAttack = TRUE ;
      c->attackTimeOut = time(0) + 180 ;
    }
  c->challenge = generate_challenge() ;
  c->incoming_trans = incoming_transaction_create();
  c->outgoing_trans = outgoing_transaction_create(s);
  c->saved_option = 0;
  c->serial_number = next_client_serial_number++;
  array(clients, s, struct client *) = c;
  if (s >= tcp_server_fd_count)
    tcp_server_fd_count = s+1;
  FD_SET(s, &open_sockets);

  /*
   * check if their host has access.  We do it this late so that we 
   * have everything set up to send back a properly formatted 
   * error message.  From a computer security perspective, I would
   * just close the connection without saying anything but I think
   * that is likely to cause a lot of confusion when people don't
   * have their configuration correct.
   *
   * The type of access is stored in the client block, but that can
   * be changed later by the client if it responds to the challenge/response
   * authentication.  Otherwise the access remains at the default
   * granted to the host.
   */
  /* c->access_type = check_tcp_access_list( (unsigned char *) & there.sin_addr ); */
  c->access_type = access_type ;
  if (c->access_type < 0)
    {
      /*
       * no access - bummer for you
       */
      fprintf(tcp_log_file, "\n%s Refusing New Client %d - source IP %d.%d.%d.%d\n"
	      , timeShowNow(),  c->serial_number
	      , ((unsigned char *) & there.sin_addr)[0]
	      , ((unsigned char *) & there.sin_addr)[1]
	      , ((unsigned char *) & there.sin_addr)[2]
	      , ((unsigned char *) & there.sin_addr)[3]
	      ) ;

      /* last param used to be  inet_ntoa(*(struct in_addr *) &there));
       * but that is not decoded on my linux 64
       */

      /*
       * We increment the client count because we will decrement it
       * when we close this socket.
       */
      ++current_active_clients; /* used to be on since Mark, mieg 2019_02_06, try to put it off, in robot case we seem not to decrease this number if we close brutally on client side */

      /*
       * send back a meaningful error message.
       */
      if (! TCPWEBSERVER)
	{
	  char goodMsg[] = "your computer does not have access, this message will not be repeated";
	  if (isAttack)
	    ;
	  else
	    {
	      outgoing_transaction_send (c->outgoing_trans, 'E', goodMsg, sizeof (goodMsg)) ;
	    }
	}
    }
  else
    {
      /*
       * some kind of valid access - send the challenge.
       */ 
      fprintf(tcp_log_file, "\n%s New Client %d, %d active clients %d baddies %d attacks - source IP %d.%d.%d.%d\n" 
	      , timeShowNow(), c->serial_number, ++current_active_clients
	      , current_bad_clients, tcp_nAttack
	      , ((unsigned char *) & there.sin_addr)[0]
	      , ((unsigned char *) & there.sin_addr)[1]
	      , ((unsigned char *) & there.sin_addr)[2]
	      , ((unsigned char *) & there.sin_addr)[3]
	      ) ;

      if (! TCPWEBSERVER) outgoing_transaction_send(c->outgoing_trans, 
				'C', c->challenge, strlen(c->challenge)+1);
    }
}


/*
 * This is the entry point from the main program.  Set up to receive
 * connections, receive them, and process them.
 */

extern int listen_socket(int);

void  acetcp_wait_for_client_using (u_long port, void (do_process_incoming)(char *command, int n, void *v_client_descriptor))
{
  int listen_fd;
  int x;
  int selecterror;
  int debug = 0 ;
  static fd_set read_fds;

  /*
   * make some randomness for the challenges
   */
  init_random();

  /*
   * find which hosts have any access at all to the database
   */
  read_tcp_access_list();


  /*
   * IGNORE sigpipe - return values from the write tell us
   * when we hit a broken pipe.
   */

  signal(SIGPIPE, SIG_IGN);

  /*
   * set up to listen for incoming connections
   */

  FD_ZERO(&open_sockets);

  listen_fd = listen_socket(port);
  if (listen_fd < 0)
    messcrash("unable to listen on port %d",port);
  FD_SET(listen_fd, &open_sockets);
  tcp_server_fd_count = listen_fd + 1;

  clients = arrayCreate(tcp_server_fd_count+10, struct clients *);

  printf("waiting for client\n");

  /*
   * count select errors - we get one once in a while for reasons that
   * are not clear to me.  I think it may relate to alarm signals 
   * originating somewhere else in the database.  The server only
   * quits if it gets too many consecutive select errors.
   */
  selecterror = 0;

  for (;;)
    {
      read_fds = open_sockets;

      /*
       * wait for activity from the listening socket or any client.
       *
       * If you need a timeout for server inactivity, do it in this select.
       */
      if (select(tcp_server_fd_count, &read_fds, NULL, NULL, NULL) < 0)
	{
	  perror("select");
	  if (selecterror++ > 5)
	    messcrash("too many select errors");
	}
      else
	selecterror = 0;

      /*
       * if there is a new connection, accept it
       */
      if (FD_ISSET(listen_fd, &read_fds))
	{
          FD_CLR(listen_fd, &read_fds);
	  accept_client(listen_fd);
	}

      /*
       * process transactions from any connection that has data in it.
       */
      for (x=0; x<tcp_server_fd_count; x++)
	{
	  struct client *c;

	  c = arr(clients, x, struct client *);
	  if (c && c->isAttack)
	    {
	      if (c->attackTimeOut < time(0))
		{
		  fprintf(tcp_log_file, "\n%s Closing Attack Client %d, %d active clients\n",
			  timeShowNow(), c->serial_number, --current_active_clients);
		  current_bad_clients-- ;
		  if (1)
		    {
		      static char *badMsg = "Sorry, you do not have acces to this server" ;
			outgoing_transaction_send (c->outgoing_trans, 'E', badMsg, sizeof (badMsg)) ;
		    }
		  if (c->look)
		    aceCommandDestroy(c->look);
		  incoming_transaction_free(c->incoming_trans);
		  outgoing_transaction_free(c->outgoing_trans);
		  close(x);
		  FD_CLR(x, &open_sockets);
		  /*   bug: leaking struct client *c    */
		  messfree(c->challenge) ;
		  messfree(c);
		  arr(clients, x, struct client *) = 0;
		}
	      continue ;
	    }
	  /*
	   * nothing on this fd
	   */
	  if (! FD_ISSET(x, &read_fds))
	    continue;

	  /*
	   * find the client
	   */
	  if (!c)
	    {
	      /*
	       * how odd - there is no client, but our select says
	       * we have data waiting on that socket.
	       */
	      close(x);
	      continue;
	    }

	  c->last_transaction_time = time(0);
	  /* 
	   * process client data
	   */
	  if (incoming_transaction(c->incoming_trans, x, FALSE, do_process_incoming, c))
	    {
	      /*
	       * it is seeing too many EOF - shut it down
	       */
 	      fprintf(tcp_log_file, "\n%s Closing Client %d, EOF %d active clients\n",
		      timeShowNow(), c->serial_number, --current_active_clients);
	      if (c->look)
		aceCommandDestroy(c->look);
	      incoming_transaction_free(c->incoming_trans);
	      outgoing_transaction_free(c->outgoing_trans);
	      if (debug) printf("close client %d\n",x);
	      close(x);
	      FD_CLR(x, &open_sockets);
	      /*   bug: leaking struct client *c    */
	      messfree(c->challenge) ;
	      messfree(c);
	      arr(clients, x, struct client *) = 0;
	    }
	  else
	    {
	      /*
	       * it is a valid transaction - it was handled by the callback
	       * process_incoming()
	       */
	    }
	}
    }
}

void  acetcp_wait_for_client (u_long port)
{
  acetcp_wait_for_client_using (port, process_incoming) ;
}

void closePortMap()
{
  /*
   * This function is called by the generic server main program just
   * before it exits.  We have to provide it, but there isn't anything
   * we should do.
   */
}


void server_shutdown()
{
  if (isWriteAccess())
    aceQuit (TRUE);
  else
    aceQuit (FALSE);
  fprintf(tcp_log_file,"\n\n%s #### Server normal exit\n  A bientot \n\n",
          timeShowNow()) ;
  exit(0);
}


/*
 * This is the callback for when a message arrives from the client.
 * It processes the message and sends back a response if necessary.
 */

void process_incoming(char *command, int n, void *v_client_descriptor)
{
  static Stack s = 0 ;
  int level_in, level_out;
  int option, encore;
  int nbytes;
  static int max_response_size = 49152;
  char *response;
  struct client *client;
  int command_type;

  static int transaction_count = 0;

  level_out = -1;
  level_in = -1;
  transaction_count++;

  client = v_client_descriptor;

  /*
   * skip the length bytes in the command
   */
  if (n >= 4 && TCPWEBSERVER && ! strncmp (command, "GET ", 4))
    {
      command_type = 'G' ;
       command += 4;
       n -= 4;
    }
  else if (n >= 4 && ! TCPWEBSERVER)
    {
      struct client *c = (struct client *) v_client_descriptor ;
      command += 11;
      n -= 11;
            if (0) outgoing_transaction_send(c->outgoing_trans, 
				'C', c->challenge, strlen(c->challenge)+1);


      /*
       * if you implement authentication, check it here.
       */
      
      /*
       * what type of command is this?
       */
      command_type = *command;
      command++; n--;  /* skip type byte */
    }

  if (n < 0 || (n > 0 && command[n-1] != '\0'))
    {
      fprintf (tcp_log_file,"client command not ending with null\n");
      goto closeit;
    }

  /* ac_dump(command,n); */

  s = stackReCreate(s, 20000) ;

  /*
   * the first byte of the 'R' response is the encore value.  we put something at
   * the beginning of the output buffer and write over it later.  Do not use '\0'
   * in case something in Stack decides to eat it.
   */
  catBinary(s, "X", 1);

  level_out = freeOutSetStack(s);

  freespecial ("\n\"\t\\/@%") ;  /* forbid sub shells, otherwise client might attempt to exe on the server side */

  if (! client->look )
    client->look = aceCommandDoCreate ( 0xff , 0, s );

  switch (command_type)
    {
    case 't':
      /*
       * tace command
       */

    case 'n':
      /*
       * tace command with no expected response
       */

    case 'e':
      /*
       * encore command
       */
      {
	level_in = freesettext (command, "") ; 

	if ( ! (client->access_type & CLIENT_READ) )
	  {
	    fprintf (tcp_log_file,"client does not have read access\n") ;
	    goto closeit;
	  }

#if 0
	if (client->access_type & CLIENT_WRITE )
	  client->look->noWrite = FALSE ;
	else
	  client->look->noWrite = TRUE ;
#endif
	/*
	 * logging the query
	 */
	fprintf(tcp_log_file, "%s Client %d (%d active). %c Query: %s\n", timeShowNow() ,
		client->serial_number, current_active_clients, command_type, command );


	/*
	 * if this is an "encore" command, we use the "option" value saved
	 * from last time.  otherwise, we use 0.  (You can't have an encore
	 * command with no response, but that is not meaningful anyway.)
	 */
	if (command_type == 'e')
	  option = client->saved_option;
	else
	  option = 0;

	option = aceCommandDoExecute (client->look, level_in,
				      FALSE /* no prompt */,
				      option, max_response_size ) ;

	encore = 0;

	switch (option)
	  {
	  case 'q':
	    fprintf (tcp_log_file, "client quits\n");
	    goto closeit;
	  case 'm':
	  case 'D':
	  case 'B':
	  case 'M':  /* automatic looping */
	    encore = 2;
	    client->saved_option = option;
	    break ;  
	  case 'W':  /* who */
	    catText(s, messprintf("// %d active Clients, %d clients so far, %d transactions so far\n",
				  current_active_clients, next_client_serial_number-1, transaction_count));
	    break;
	  case 'U':  /* shUtdown */
	    if (client->access_type & CLIENT_WRITE )
	      {
		/*
		 * server shutting down, but first send a valid response to the shutdown
		 * command();
		 */
		outgoing_transaction_send(client->outgoing_trans, 'R', "\0shutdown\n", 10 );
		server_shutdown();
		/* never get here */
		abort();
	      }
	    else
	      {
		outgoing_transaction_send(client->outgoing_trans, 'R', "\0Sorry\n",  7 );
	      }
	    break;
	  default:
	    client->saved_option = 0 ;
	    break ; 
	  }

	if ( (command_type == 't') || (command_type == 'e') )
	  {
	    /*
	     * For type 't' or 'e' commands, we send back the command output
	     * as a response to the user.  For type 'n' commands, the
	     * user is not expecting a response, so we don't send it.
	     */

	    /*
	     * send back the content of Stack s
	     */
	    nbytes = stackMark(s);
	    nbytes = nbytes - 1;
	    response = stackText(s, 0);

	    *response = encore;	/* remember we put an X there at the beginning */

	    /*  ac_dump(response, nbytes < 200 ? nbytes : 200 ); */
	    if (outgoing_transaction_send(client->outgoing_trans, 'R', response, nbytes ))
	      {
		/* client vanished - clean it up */
		fprintf (tcp_log_file,"client disconnect detected on write\n");
		goto closeit;
	      }
	  }
 	  
      }
      break;
#ifdef ACEMBLY
    case 'G':
      {
	vTXT txt = vtxtCreate () ;
	if (1)
	  magicProcessWebRequest (txt, command) ;
	
        vtxtPrintf (txt, "\n%c%c", 0x0D,0x0A) ;
	response = vtxtPtr (txt) ;
	nbytes = strlen (response) ;
	
	/*  ac_dump(response, nbytes < 200 ? nbytes : 200 ); */
	if (outgoing_transaction_send(client->outgoing_trans, 'G', response, nbytes ))
	  {
	    /* client vanished - clean it up */
	    fprintf (tcp_log_file,"client vanished detected on write\n");
	    goto closeit;
	  }
	ac_free (txt) ;
	goto closeit; /* to terminate communication with the web browser */
      }
      break ;
#endif
    default:
      /*
       * if the client can't follow the protocol, just hang up.
       */
      if (0) fprintf (tcp_log_file," client cant follow the protocol, just hang up\n") ;
      goto returning;
    }

  goto returning;

 closeit:
  /*
   * closing - just shut down the socket.  the next time around the loop,
   * it will have eof and we will clean it up.
   */
  shutdown(client->fd, 2);

 returning:
  if (level_out != -1)
    freeOutClose(level_out);
  if (level_in != -1)
    freeclose(level_in);
  stackDestroy (s) ;
}

