/*******************************************************************************
*  program pfs_sample
*  $Id$
*  This programs reads some data from the portable fast sampler
*
*  usage:
*  	pfs_sample -m mode [-p (print all data)] [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*	the -p option specifies to print all the data recorded
*                     (default is to print a sample every megabyte)
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log$
   Revision 3.0  2003/02/25 22:00:59  cvs
   Adapted to use Joseph Jao's byte unpacking.

   Revision 1.3  2002/04/27 20:19:17  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.2  2002/04/27 06:44:07  margot
   Updated device filename.

   Revision 1.1  2002/04/27 06:25:21  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "edtinc.h"
#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

FILE   *fpinput;		/* pointer to input file */
FILE   *fpoutput;		/* pointer to output file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_files();
void copy_cmd_line();

int main(int argc, char *argv[])
{
  EdtDev *edt_p ;
  int mode;
  int bufsize = 1048576; /* 1048576 size of read buffer, default 1 MB */
  unsigned char *buffer;		/* buffer for packed data */
  char *rcp,*lcp;	/* buffer for unpacked data */
  int smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int levels;		/* # of levels for given quantization mode */
  int open_flags;	/* flags required for open() call */
  int printall;
  int bytesw;
  int i;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&printall);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open input and output files, stdin & stdout default */
  open_files(infile,outfile,&fpinput,&fpoutput);

  switch (mode)
    {
    case -1: smpwd = 8; break;  
    case  1: smpwd = 8; break;
    case  2: smpwd = 4; break;
    case  3: smpwd = 2; break; 
    case  5: smpwd = 4; break;
    case  6: smpwd = 2; break;
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer = NULL;
  rcp = (char *) malloc(2 * nsamples * sizeof(char));
  lcp = (char *) malloc(2 * nsamples * sizeof(char));
  if (rcp == NULL || lcp == NULL)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  if ((edt_p = edt_open("edt", 0)) == NULL)
    {
      perror("edt_open") ;
      exit(1) ;
    }
  else
    fprintf(stderr,"Device opened\n");
  
  if (edt_configure_ring_buffers(edt_p, bufsize, 32, EDT_READ, NULL) == -1)
    {
      perror("edt_configure_ring_buffers") ;
      exit(1) ;
    }
  else
    fprintf(stderr,"Buffers configured\n");
  

  /* arm trigger */
  edt_reg_write(edt_p, PCD_FUNCT, 0x01 | (mode << 1));
  
  /* flush fifo */
  edt_flush_fifo(edt_p);

  /* start the transfers in free running mode */
  if (edt_start_buffers(edt_p, 0) == -1)
    {
      perror("edt_start_buffers") ;
      exit(1) ;
    }
  else
    fprintf(stderr,"Buffers started\n");
  
  for (;;)
    {
      buffer = edt_wait_for_buffers(edt_p, 1) ;

      switch (mode)
	{
	case 1:
  	  unpack_pfs_2c2b(buffer, rcp, bufsize);
	  if (printall)
	    for (i = 0; i < 2*nsamples; i+=2) 
	      fprintf(stdout,"% 4.0d % 4.0d\n",rcp[i],rcp[i+1]);
	  else
	    fprintf(stdout,"% 4.0d % 4.0d\n",rcp[0],rcp[1]);

	  break;
	case 2: 
	  unpack_pfs_2c4b(buffer, rcp, bufsize);
	  if (printall)
	    for (i = 0; i < 2*nsamples; i+=2) 
	      fprintf(stdout,"% 4.0d % 4.0d\n",rcp[i],rcp[i+1]);
	  else
	    fprintf(stdout,"% 4.0d % 4.0d\n",rcp[0],rcp[1]);

	  break;
	case 3: 
	  unpack_pfs_2c8b(buffer, rcp, bufsize); 
	  if (printall)
	    for (i = 0; i < 2*nsamples; i+=2) 
	      fprintf(stdout,"% 4.0d % 4.0d\n",rcp[i],rcp[i+1]);
	  else
	    fprintf(stdout,"% 4.0d % 4.0d\n",rcp[0],rcp[1]);

	  break;
	case 5:
	  unpack_pfs_4c2b_rcp (buffer, rcp, bufsize); 
	  unpack_pfs_4c2b_lcp (buffer, lcp, bufsize); 

	  if (printall)
	    for (i = 0; i < 2*nsamples; i+=2) 
	      fprintf(stdout,"% 4.0d % 4.0d % 4.0d % 4.0d\n",
		      rcp[i],rcp[i+1],lcp[i],lcp[i+1]);
	  else
	    fprintf(stdout,"% 4.0d % 4.0d % 4.0d % 4.0d\n",
		    rcp[0],rcp[1],lcp[0],lcp[1]);

	  break;
	case 6:
	  unpack_pfs_4c4b_rcp (buffer, rcp, bufsize);
	  unpack_pfs_4c4b_lcp (buffer, lcp, bufsize);

	  if (printall)
	    for (i = 0; i < 2*nsamples; i+=2) 
	      fprintf(stdout,"% 4.0d % 4.0d % 4.0d % 4.0d\n",
		      rcp[i],rcp[i+1],lcp[i],lcp[i+1]);
	  else
	    fprintf(stdout,"% 4.0d % 4.0d % 4.0d % 4.0d\n",
		    rcp[0],rcp[1],lcp[0],lcp[1]);

	  break;
	default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
	}

    }
  
  return 0;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,printall)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *printall;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_sample program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:p"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_sample -m mode [-p (print all data)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *printall = 0;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
    {
      case 'o':
 	       *outfile = optarg;	/* output file name */
               arg_count += 2;		/* two command line arguments */
	       break;

      case 'm':
 	       sscanf(optarg,"%d",mode);
               arg_count += 2;		/* two command line arguments */
	       break;
	    
      case 'p':
 	       *printall = 1;
               arg_count += 1;
	       break;
	    
      case '?':			 /*if not in myoptions, getopt rets ? */
               goto errout;
               break;
    }
  }

  if (arg_count < argc)		 /* 1st non-optioned param is infile */
    *infile = argv[arg_count];

  /* must specify a valid mode */
  if (*mode == 0) goto errout;
  
  return;

  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE1);
          fprintf(stderr,"%s",USAGE2);
	  exit(1);
}

/******************************************************************************/
/*	open files    							      */
/******************************************************************************/
void	open_files(infile,outfile,fpinput,fpoutput)
char	*infile;		/* input file name */
char	*outfile;		/* output file name */
FILE    **fpinput;		/* pointer to input file */
FILE    **fpoutput;		/* pointer to output file */
{
  /* opens the input and output files, sdtin and stdout are default */

  if (infile[0] == '-') 
    *fpinput=stdin;			/* use stdin */
  else
  {
    *fpinput=fopen(infile,"r");		/* open file input */
    if (*fpinput == NULL)
    {
      perror("open_files: input file open error");
      exit(1);
    }
  }

  if (outfile[0] == '-')
  {
    *fpoutput=stdout;
  }
  else
  {
    *fpoutput=fopen(outfile,"w");
     if (*fpoutput == NULL)
    {
      perror("open_files: output file open error");
      exit(1);
    }
  }

  return;
}

/******************************************************************************/
/*	copy_cmd_line    						      */
/******************************************************************************/
void	copy_cmd_line(argc,argv,command_line)
int	argc;
char	**argv;			/* command line arguements */
char	command_line[];		/* command line parameters in single string */
{
  /* copys the command line parameters in argv to the single string
     command line
  */
  char 	*result;
  int	i;

  result = strcpy(command_line,argv[0]);
  result = strcat(command_line," ");

  for (i=1; i<argc; i++)
  {
    result = strcat(command_line,argv[i]);
    result = strcat(command_line," ");
  }

  return;
}	

