/*******************************************************************************
*  program pfs_hist
*  $Id$
*  This programs unpacks some data from the portable fast sampler
*  and prints and histogram of count values for all channels
*
*  usage:
*  	pfs_hist -m mode [-a (parse all data)] [-e (parse data at eof)] 
*               [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*	the -e option specifies to parse data at the end of the file
*	the -a option specifies to parse all the data recorded
*                     (default is to parse the first megabyte)
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log$
   Revision 1.6  2002/04/27 20:22:26  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.5  2001/07/10 00:37:54  margot
   Adjusted input buffer size according to file size.

   Revision 1.4  2001/07/10 00:24:07  margot
   Added unpacking of signed bytes.

   Revision 1.3  2001/07/10 00:18:39  margot
   Added -e option for parsing at end of file.

   Revision 1.2  2000/10/30 05:21:41  margot
   Added variable nsamples.

   Revision 1.1  2000/10/30 04:46:22  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <asm/fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();

void iq_hist(float *inbuf, int nsamples, int levels);
void iq_hist_8b(float *inbuf, int nsamples, int levels);

int main(int argc, char *argv[])
{
  struct stat filestat;	/* input file status structure */
  int mode;		/* data acquisition mode */
  int bufsize = 1048576;/* size of read buffer, default 1 MB */
  char *buffer;		/* buffer for packed data */
  float *rcp,*lcp;	/* buffer for unpacked data */
  int smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int levels;		/* # of levels for given quantization mode */
  int open_flags;	/* flags required for open() call */
  int parse_all;
  int parse_end;
  int i;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&parse_all,&parse_end);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open output file, stdout default */
  open_file(outfile,&fpoutput);

  /* open file input */
#ifdef LARGEFILE
  open_flags = O_RDONLY|O_LARGEFILE;
#else
  open_flags = O_RDONLY;
#endif
  if((fdinput = open(infile, open_flags)) < 0 )
    {
      perror("open input file");
      exit(1);
    }

  /* get file status */
  if (fstat (fdinput, &filestat) < 0)
    {
      perror("input file status");
      exit(1);
    }
  
  /* adjust buffer size if needed */
  if (filestat.st_size < bufsize)
    bufsize = filestat.st_size;

  switch (mode)
    {
    case -1: smpwd = 8; levels =   4; break;  
    case  1: smpwd = 8; levels =   4; break;
    case  2: smpwd = 4; levels =  16; break;
    case  3: smpwd = 2; levels = 256; break; 
    case  5: smpwd = 4; levels =   4; break;
    case  6: smpwd = 2; levels =  16; break;
    case  8: smpwd = 2; levels = 256; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer = (char *) malloc(bufsize);
  rcp = (float *) malloc(2 * nsamples * sizeof(float));
  lcp = (float *) malloc(2 * nsamples * sizeof(float));
  if (lcp == NULL) 
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  if (parse_end)
    lseek(fdinput, -bufsize, SEEK_END);

  /* read first buffer */
  if (bufsize != read(fdinput, buffer, bufsize))
    {
      fprintf(stderr,"Read error\n");
      exit(1);
    }  

  switch (mode)
    { 
    case 1:
      unpack_pfs_2c2b(buffer, rcp, bufsize);
      iq_hist(rcp, nsamples, levels);
      break;
    case 2: 
      unpack_pfs_2c4b(buffer, rcp, bufsize);
      iq_hist(rcp, nsamples, levels);
      break;
    case 3: 
      unpack_pfs_2c8b(buffer, rcp, bufsize);
      iq_hist_8b(rcp, nsamples, levels);
      break;
    case 5:
      unpack_pfs_4c2b(buffer, rcp, lcp, bufsize);
      fprintf(fpoutput,"RCP hist\n");
      iq_hist(rcp, nsamples, levels);
      fprintf(fpoutput,"LCP hist\n");
      iq_hist(lcp, nsamples, levels);
      break;
    case 6:
      unpack_pfs_4c4b(buffer, rcp, lcp, bufsize);
      fprintf(fpoutput,"RCP hist\n");
      iq_hist(rcp, nsamples, levels);
      fprintf(fpoutput,"LCP hist\n");
      iq_hist(lcp, nsamples, levels);
      break;
    case 8: 
      unpack_pfs_signedbytes(buffer, rcp, bufsize);
      iq_hist_8b(rcp, nsamples, levels);
      break;
    default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
    }
  
  return 0;
}

/******************************************************************************/
/*	iq_hist						         	      */
/******************************************************************************/
void iq_hist(float *inbuf, int nsamples, int levels)
{
  double i=0;
  double q=0;

  int ihist[2*levels];
  int qhist[2*levels];

  int k;

  /* initialize */  
  for (k = 0; k < 2 * levels; k += 2)
    ihist[k] = qhist[k] = 0;

  /* compute histogram */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      ihist[(int)inbuf[k]   + levels - 1] += 1; 
      qhist[(int)inbuf[k+1] + levels - 1] += 1; 
    }

  /* print results */
  for (k = 0; k < 2 * levels; k += 2)
    {
      fprintf(fpoutput,"%10d %15d",k - levels + 1,ihist[k]);
      fprintf(fpoutput,"\t"); 
      fprintf(fpoutput,"%10d %15d",k - levels + 1,qhist[k]);
      fprintf(fpoutput,"\n"); 
    }

  return;
}    

/******************************************************************************/
/*	iq_hist_8b					         	      */
/******************************************************************************/
void iq_hist_8b(float *inbuf, int nsamples, int levels)
{
  double i=0;
  double q=0;

  int ihist[levels];
  int qhist[levels];

  int k;

  /* initialize */  
  for (k = 0; k < levels; k++)
    ihist[k] = qhist[k] = 0;

  /* compute histogram */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      ihist[(int)inbuf[k]   + levels/2] += 1; 
      qhist[(int)inbuf[k+1] + levels/2] += 1; 
    }

  /* print results */
  for (k = 0; k < levels; k ++)
    {
      fprintf(fpoutput,"%10d %15d",k - levels/2,ihist[k]);
      fprintf(fpoutput,"\t"); 
      fprintf(fpoutput,"%10d %15d",k - levels/2,qhist[k]);
      fprintf(fpoutput,"\n"); 
    }

  return;
}    

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,parse_all,parse_end)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *parse_all;
int     *parse_end;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_hist program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:ae"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_hist -m mode [-e (parse data at eof)] [-a (parse all data)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *parse_all = 0;
  *parse_end = 0;

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
	    
      case 'a':
 	       *parse_all = 1;
               arg_count += 1;
	       break;
	    
      case 'e':
 	       *parse_end = 1;
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
  
  /* code still in development */
  if (*parse_all)
    {
      fprintf(stderr,"-a option not implemented yet\n"); 
      exit(1);
    }
  return;

  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE1);
          fprintf(stderr,"%s",USAGE2);
	  exit(1);
}

/******************************************************************************/
/*	open file    							      */
/******************************************************************************/
void	open_file(outfile,fpoutput)
char	*outfile;		/* output file name */
FILE    **fpoutput;		/* pointer to output file */
{
  /* opens the output file, stdout is default */
  if (outfile[0] == '-')
    *fpoutput=stdout;
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

