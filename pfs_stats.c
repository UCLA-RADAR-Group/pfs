/*******************************************************************************
*  program pfs_stats
*  $Id$
*  This programs unpacks some data from the portable fast sampler
*  and prints statistics such as mean and standard deviation for 
*  all channels
*
*  usage:
*  	pfs_stats -m mode [-a (parse all data)] [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*	the -a option specifies to parse all the data recorded
*                     (default is to parse the first megabyte)
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log$
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
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

void iq_stats(float *inbuf, int nsamples);

int main(int argc, char *argv[])
{
  int mode;
  char *buffer;
  float *buf;
  float *rcp,*lcp;
  int bufsize = 1024 * 1024;
  int parseall;
  int smpwd;		/* # of single pol complex samples in a 4 byte word */
  int open_flags;	/* flags required for open() call */
  int i;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&parseall);

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
    perror("open input file");

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
  if (mode < 5)
    buf = (float *) malloc(2 * bufsize * smpwd / 4 * sizeof(float));
  else
    {
      rcp = (float *) malloc(2 * bufsize * smpwd / 4 * sizeof(float));
      lcp = (float *) malloc(2 * bufsize * smpwd / 4 * sizeof(float));
    }
  buffer = (char *) malloc(bufsize);
  if (buffer == NULL)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  if (parseall)
    {
      fprintf(stderr,"-a option not implemented yet\n"); 
      exit(1);
    }

  if (bufsize != read(fdinput, buffer, bufsize))
    fprintf(stderr,"Read error\n");
  
  switch (mode)
    { 
    case -1:
      unpack_gsb(buffer, buf, bufsize);
      iq_stats(buf, bufsize * 2);
      break;
    case 1:
      unpack_pfs_2c2b(buffer, buf, bufsize);
      iq_stats(buf, bufsize * 2);
      break;
    case 2: 
      unpack_pfs_2c4b(buffer, buf, bufsize);
      iq_stats(buf, bufsize);
      break;
    case 3: 
      unpack_pfs_2c8b(buffer, buf, bufsize);
      iq_stats(buf, bufsize / 2);
      break;
    case 5:
      unpack_pfs_4c2b(buffer, rcp, lcp, bufsize);
      fprintf(fpoutput,"RCP stats\n");
      iq_stats(rcp, bufsize);
      fprintf(fpoutput,"LCP stats\n");
      iq_stats(lcp, bufsize);
      break;
    case 6:
      unpack_pfs_4c4b(buffer, rcp, lcp, bufsize);
      fprintf(fpoutput,"RCP stats\n");
      iq_stats(rcp, bufsize / 2);
      fprintf(fpoutput,"LCP stats\n");
      iq_stats(lcp, bufsize / 2);
      break;

    default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
    }
  
  return 0;
}

/******************************************************************************/
/*	iq_stats							      */
/******************************************************************************/
void iq_stats(float *inbuf, int nsamples)
{
  double i=0;
  double q=0;
  double ii=0;
  double qq=0;
  double iq=0;

  int k;

  /* sum Is and Qs */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      i  += inbuf[k];
      q  += inbuf[k+1];
      iq += inbuf[k] * inbuf[k+1];
      ii += inbuf[k] * inbuf[k];
      qq += inbuf[k+1] * inbuf[k+1];
    }

  /* compute mean and standard deviation */
  i = i / nsamples;
  q = q / nsamples;
  
  iq = iq / nsamples;

  ii = ii / nsamples;
  qq = qq / nsamples;

  ii = sqrt(ii - i*i);
  qq = sqrt(qq - q*q);

  fprintf(fpoutput,"% 10.4f % 10.4f ",i,ii);
  fprintf(fpoutput,"% 10.4f % 10.4f ",q,qq);
  fprintf(fpoutput,"% 10.4f ",fabs(iq - i*q)/ii/qq);
  fprintf(fpoutput,"\n"); 

  return;
}    

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,parseall)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *parseall;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_stats program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:a"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_stats -m mode [-a (parse all data)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t-1: GSB\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *parseall = 0;

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
 	       *parseall = 1;
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

