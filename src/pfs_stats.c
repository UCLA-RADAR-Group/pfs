/*******************************************************************************
*  program pfs_stats
*  $Id: pfs_stats.c,v 3.2 2009/11/16 19:08:21 jlm Exp $
*  This programs unpacks some data from the portable fast sampler
*  and prints statistics such as mean and standard deviation for 
*  all channels
*
*  usage:
*  	pfs_stats -m mode [-a (parse all data)] [-e (parse data at eof)]
*                [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*       the -e option specifies to parse data at the end of the file
*	the -a option specifies to parse all the data recorded
*                     (default is to parse the first megabyte)
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log: pfs_stats.c,v $
   Revision 3.2  2009/11/16 19:08:21  jlm
   Added ifdef flag for Mac compilation.

   Revision 3.1  2003/02/26 00:42:59  margot
   Fixed first argument to floatsum()

   Revision 3.0  2003/02/25 22:40:14  cvs
   Adapted to use Joseph Jao's byte unpacking.

   Revision 2.3  2002/06/05 16:56:01  cvs
   Added > sign for easier parsing of output.

   Revision 2.2  2002/05/26 04:22:23  cvs
   Prints warning for odd file sizes.

   Revision 2.1  2002/05/26 01:01:49  cvs
   Better scheme for handling short buffers, including those at end of file.

   Revision 2.0  2002/05/13 00:01:32  cvs
   Added -a option for statistics on entire input file.

   Revision 1.7  2002/05/12 20:38:27  cvs
   Added mode for floating point values.

   Revision 1.6  2002/05/03 18:08:11  cvs
   Added dBm scale.

   Revision 1.5  2002/04/27 20:22:10  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.4  2001/07/06 19:49:26  margot
   Added -e option for statistics at end of file.
   Added unpacking of signed bytes.

   Revision 1.3  2000/10/30 22:42:32  margot
   Added voltage scale.

   Revision 1.2  2000/10/30 05:20:46  margot
   Added variable nsamples.

   Revision 1.1  2000/10/30 04:45:22  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id: pfs_stats.c,v 3.2 2009/11/16 19:08:21 jlm Exp $";

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();

void sum(char *inbuf, int nsamples, double *i, double *q, double *ii, double *qq, double *iq);
void floatsum(float *inbuf, int nsamples, double *i, double *q, double *ii, double *qq, double *iq);

int main(int argc, char *argv[])
{
  struct stat filestat;	/* input file status structure */
  int bufsize = 1000000;/* size of read buffer, default 1 MB */
  char *buffer;		/* buffer for packed data */
  char *rcp,*lcp;	/* buffer for unpacked data */
  float *fbuffer;	/* buffer for unpacked data */
  float smpwd;		/* # of single pol complex samples in a 4 byte word */
  double ri,rq,rii,rqq,riq;/* accumulators for statistics */
  double li,lq,lii,lqq,liq;/* accumulators for statistics */
  int nsamples;		/* # of complex samples in each buffer */
  int ntotal;		/* total number of samples used in computing statistics */
  int bytesread;	/* number of bytes read from input file */
  int levels;		/* # of levels for given quantization mode */
  int open_flags;	/* flags required for open() call */
  int parse_all;
  int parse_end;
  int mode;
  int k;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&parse_all,&parse_end);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open output file, stdout default */
  open_file(outfile,&fpoutput);

  /* open file input */
#ifndef __APPLE__
  open_flags = O_RDONLY|O_LARGEFILE;
#else
  open_flags = O_RDONLY;
#endif
  if((fdinput = open(infile, open_flags)) < 0 )
    {
      perror("open input file");
      exit(1);
    }

  /* check file size */
  if (fstat (fdinput, &filestat) < 0)
    {
      perror("input file status");
      exit(1);
    }
  if (filestat.st_size % 4 != 0)
    fprintf(stderr,"Warning: file size %d is not a multiple of 4\n",
	    (int) filestat.st_size);

  switch (mode)
    {
    case -1: smpwd = 8; levels =   4; break;  
    case  1: smpwd = 8; levels =   4; break;
    case  2: smpwd = 4; levels =  16; break;
    case  3: smpwd = 2; levels = 256; break; 
    case  5: smpwd = 4; levels =   4; break;
    case  6: smpwd = 2; levels =  16; break;
    /* A/D levels do not apply for packing modes below */  
    case   8: smpwd =   2; break; 
    case  32: smpwd = 0.5; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = (int) rint(bufsize * smpwd / 4.0);
  buffer = (char *) malloc(bufsize);
  rcp = (char *) malloc(2 * nsamples * sizeof(char));
  lcp = (char *) malloc(2 * nsamples * sizeof(char));
  fbuffer = (float *) malloc(2 * nsamples * sizeof(float));
  if (!lcp || !rcp || !fbuffer)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* initialize counters */
  ntotal = 0;
  ri  = li  = 0; 
  rq  = lq  = 0; 
  rii = lii = 0;
  rqq = lqq = 0;
  riq = liq = 0;

  /* go to end of file if requested */
  if (parse_end)
    lseek(fdinput, -bufsize, SEEK_END);

  /* infinite loop */
  while (1)
    {
      /* read one buffer */
      bytesread = read(fdinput, buffer, bufsize);
      /* check for end of file */
      if (bytesread == 0) break;
      /* handle small buffers */
      if (bytesread != bufsize) 
	{
	  bufsize = bytesread;
	  nsamples = (int) rint(bufsize * smpwd / 4.0);
	}

      switch (mode)
	{ 
	case 1:
	  unpack_pfs_2c2b(buffer, rcp, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  break;
	case 2: 
	  unpack_pfs_2c4b(buffer, rcp, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  break;
	case 3: 
	  unpack_pfs_2c8b(buffer, rcp, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  break;
	case 5:
	  unpack_pfs_4c2b_rcp(buffer, rcp, bufsize);
	  unpack_pfs_4c2b_lcp(buffer, lcp, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  sum(lcp, nsamples, &li, &lq, &lii, &lqq, &liq);
	  break;
	case 6:
	  unpack_pfs_4c4b_rcp(buffer, rcp, bufsize);
	  unpack_pfs_4c4b_lcp(buffer, lcp, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  sum(lcp, nsamples, &li, &lq, &lii, &lqq, &liq);
	  break;
	case 8: 
	  memcpy (rcp, buffer, bufsize);
	  sum(rcp, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  break;
	case 32: 
	  memcpy (fbuffer, buffer, bufsize);
	  floatsum(fbuffer, nsamples, &ri, &rq, &rii, &rqq, &riq);
	  break;
	default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
	}
      
      ntotal += nsamples;
      if (!parse_all) break;
    }

  /* compute mean and standard deviation (RCP) */
  ri = ri / ntotal;
  rq = rq / ntotal;
  riq = riq / ntotal;
  rii = rii / ntotal;
  rqq = rqq / ntotal;
  rii = sqrt(rii - ri*ri);
  rqq = sqrt(rqq - rq*rq);

  /* compute mean and standard deviation (LCP) */
  li = li / ntotal;
  lq = lq / ntotal;
  liq = liq / ntotal;
  lii = lii / ntotal;
  lqq = lqq / ntotal;
  lii = sqrt(lii - li*li);
  lqq = sqrt(lqq - lq*lq);

  /* print results */
  if (mode > 8) 
    fprintf(fpoutput,"Statistics on %d samples:\n",ntotal);
  else
    fprintf(fpoutput,"In digitizer counts (x2):\n");
  fprintf(fpoutput,"     DC I      RMS I       DC Q      RMS Q       rIQ\n");

  if (mode == 5 || mode == 6)
    {
      fprintf(fpoutput,"RCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",ri,rii);
      fprintf(fpoutput,"% 10.4f % 10.4f ",rq,rqq);
      fprintf(fpoutput,"% 10.4f ",fabs(riq - ri*rq)/rii/rqq);
      fprintf(fpoutput,"\n"); 
      
      fprintf(fpoutput,"LCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",li,lii);
      fprintf(fpoutput,"% 10.4f % 10.4f ",lq,lqq);
      fprintf(fpoutput,"% 10.4f ",fabs(liq - li*lq)/lii/lqq);
      fprintf(fpoutput,"\n"); 
    }
  else
    {
      fprintf(fpoutput,">");
      fprintf(fpoutput,"% 10.4f % 10.4f ",ri,rii);
      fprintf(fpoutput,"% 10.4f % 10.4f ",rq,rqq);
      fprintf(fpoutput,"% 10.4f ",fabs(riq - ri*rq)/rii/rqq);
      fprintf(fpoutput,"\n"); 
    }

  /* exit here for bytes and floats */
  if (mode >= 8) return 0;

  /* convert to volts */
  /* AD range is 1 Vpp */
  fprintf(fpoutput,"\nIn Volts:\n");
  fprintf(fpoutput,"     DC I      RMS I       DC Q      RMS Q       rIQ\n");

  if (mode == 5 || mode == 6)
    {
      fprintf(fpoutput,"RCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",ri/levels/2.0,rii/levels/2.0);
      fprintf(fpoutput,"% 10.4f % 10.4f ",rq/levels/2.0,rqq/levels/2.0);
      fprintf(fpoutput,"\n"); 

      fprintf(fpoutput,"LCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",li/levels/2.0,lii/levels/2.0);
      fprintf(fpoutput,"% 10.4f % 10.4f ",lq/levels/2.0,lqq/levels/2.0);
      fprintf(fpoutput,"\n"); 
    }
  else
    {
      fprintf(fpoutput,"% 10.4f % 10.4f ",ri/levels/2.0,rii/levels/2.0);
      fprintf(fpoutput,"% 10.4f % 10.4f ",rq/levels/2.0,rqq/levels/2.0);
      fprintf(fpoutput,"\n"); 
    }

  /* convert to dBm */
  /* AD range is 1 Vpp */
  fprintf(fpoutput,"\nIn dBm:\n");
  fprintf(fpoutput,"     DC I      RMS I       DC Q      RMS Q       rIQ\n");

  if (mode == 5 || mode == 6)
    {
      fprintf(fpoutput,"RCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(rii/levels/2.0)+13);
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(rqq/levels/2.0)+13);
      fprintf(fpoutput,"\n"); 

      fprintf(fpoutput,"LCP stats\n");
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(lii/levels/2.0)+13);
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(lqq/levels/2.0)+13);
      fprintf(fpoutput,"\n");
    } 
  else
    {
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(rii/levels/2.0)+13);
      fprintf(fpoutput,"% 10.4f % 10.4f ",0.0,20*log10(rqq/levels/2.0)+13);
      fprintf(fpoutput,"\n"); 
    }
 
  return 0;
}

/******************************************************************************/
/*	sum								      */
/******************************************************************************/
void sum(char *inbuf, int nsamples, double *i, double *q, double *ii, double *qq, double *iq)
{
  int k;

  /* sum Is and Qs */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      *i  += inbuf[k];
      *q  += inbuf[k+1];
      *iq += inbuf[k] * inbuf[k+1];
      *ii += inbuf[k] * inbuf[k];
      *qq += inbuf[k+1] * inbuf[k+1];
    }

  return;
}    

/******************************************************************************/
/* floatsum								      */
/******************************************************************************/
void floatsum(float *inbuf, int nsamples, double *i, double *q, double *ii, double *qq, double *iq)
{
  int k;

  /* sum Is and Qs */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      *i  += inbuf[k];
      *q  += inbuf[k+1];
      *iq += inbuf[k] * inbuf[k+1];
      *ii += inbuf[k] * inbuf[k];
      *qq += inbuf[k+1] * inbuf[k+1];
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
     This is a template which has been customised for the pfs_stats program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:ae"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_stats -m mode [-e (parse data at eof)] [-a (parse all data)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t32: 32bit floats\n";
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

