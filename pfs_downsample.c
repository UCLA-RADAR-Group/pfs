/*******************************************************************************
*  program pfs_downsample
*  $Id$
*  This programs downsamples data from the portable fast sampler
*  by summing coherently a specified number of consecutive samples
*
*  usage:
*  	pfs_downsample -m mode -d downsampling factor 
*                      [-s scale fudge factor]
*                      [-f output floating point numbers]
*                      [-I dcoffi] [-Q dcoffq] 
*                      [-c channel] 
*                      [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*	the -d argument specifies the downsampling factor
*       the -c argument specifies which channel (1 or 2) to process
*
*  output:
*	the -o option identifies the output file, stdout is default
*       the sums are written as signed bytes or floats if -f is used
*******************************************************************************/

/* 
   $Log$
   Revision 1.4  2002/04/27 20:23:06  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.3  2002/04/27 06:13:14  margot
   Now exits if input file cannot be opened.

   Revision 1.2  2001/07/10 01:43:35  margot
   Added check on file size and compared to downsampling factor
   and input buffer size.

   Revision 1.1  2001/07/06 19:46:23  margot
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
int verbose = 0;
int floats  = 0;

void processargs();
void open_file();
void copy_cmd_line();
void iq_downsample(float *inbuf, int nsamples, int downsample, float maxvalue, float dcoffi, float dcoffq);

int main(int argc, char *argv[])
{
  struct stat filestat;	/* input file status structure */
  int mode;		/* data acquisition mode */
  char *buffer;		/* buffer for packed data */
  float *rcp,*lcp;	/* buffer for unpacked data */
  int smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  float maxunpack;	/* maximum unpacked value from libunpack */
  float maxvalue;	/* maximum achievable value by downsampling */
  float scale;		/* scaling factor to fit in a byte */
  float fudge;		/* scale fudge factor */
  float dcoffi,dcoffq;	/* dc offsets */
  int bufsize;		/* input buffer size */
  int downsample;	/* factor by which to downsample */
  int chan;		/* channel to process (1 or 2) for dual pol data */
  int open_flags;	/* flags required for open() call */
  int i;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&downsample,&chan,&dcoffi,&dcoffq,&fudge);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

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

  /* open output file, stdout default */
  open_file(outfile,&fpoutput);


  switch (mode)
    {
    case -1:  smpwd = 8; maxunpack =   +3; break;  
    case  1:  smpwd = 8; maxunpack =   +3; break;
    case  2:  smpwd = 4; maxunpack =  +15; break;
    case  3:  smpwd = 2; maxunpack = +255; break; 
    case  5:  smpwd = 4; maxunpack =   +3; break;
    case  6:  smpwd = 2; maxunpack =  +15; break;
    case  8:  smpwd = 2; maxunpack = +255; break;
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* compute dynamic range parameters */
  fprintf(stderr,"Downsampling file of size %d by %d\n",filestat.st_size,downsample);
  maxvalue = maxunpack * sqrt(downsample);
  scale = fudge * 0.25 * 128 / maxvalue;

  if (!floats && maxvalue > 255)
    {
      fprintf(stderr,"You may have a dynamic range problem\n");
      fprintf(stderr,"Turning verbose mode on so you can detect clipping\n");
      verbose = 1;
    }

  /* compute buffer size */
  /* we need something of order a Megabyte, multiple of the downsampling factor */
  /* and ideally consistent with the input file size */
  bufsize = (int) rint(1000000/downsample) * downsample;

  if (filestat.st_size % downsample != 0)
    fprintf(stderr,"This may leave the last buffer truncated\n");
  else
    /* if the filesize is a multiple of the downsampling factor, */
    /* we try to choose a buffer size that will fit nicely as well */
    while (filestat.st_size % bufsize != 0)
      bufsize -= downsample;
  fprintf(stderr,"Using %d buffers of size of %d\n", 
	  filestat.st_size / bufsize, bufsize);
  

  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer = (char *) malloc(bufsize);
  rcp = (float *) malloc(2 * bufsize * smpwd / 4 * sizeof(float));
  lcp = (float *) malloc(2 * bufsize * smpwd / 4 * sizeof(float));
  if (!lcp) fprintf(stderr,"Malloc error\n");

  while (1)
    {
      /* read one buffer */
      if (bufsize != read(fdinput, buffer, bufsize))
	{
	  perror("read");
	  fprintf(stderr,"Read error or EOF\n");
	  exit(1);
	}

      /* unpack and downsample */
      switch (mode)
	{
	case 1:
	  unpack_pfs_2c2b(buffer, rcp, bufsize); 
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq);
	  break;
	case 2: 
	  unpack_pfs_2c4b(buffer, rcp, bufsize);
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq);
	  break;
	case 3: 
	  unpack_pfs_2c8b(buffer, rcp, bufsize);
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq);
	  break;
	case 5:
	  unpack_pfs_4c2b(buffer, rcp, lcp, bufsize);
	  if (chan == 2) memcpy(rcp, lcp, 2 * nsamples * sizeof(float));
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq); 
	  break;
	case 6:
	  unpack_pfs_4c4b(buffer, rcp, lcp, bufsize);
	  if (chan == 2) memcpy(rcp, lcp, 2 * nsamples * sizeof(float));
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq); 
	  break;
     	case 8: 
	  unpack_pfs_signedbytes(buffer, rcp, bufsize);
	  iq_downsample(rcp, nsamples, downsample, scale, dcoffi, dcoffq);
	  break;
	default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
	}
    }

  return 0;
}

/******************************************************************************/
/*	iq_downsample							      */
/******************************************************************************/
void iq_downsample(float *inbuf, int nsamples, int downsample, float scale, float dcoffi, float dcoffq)
{
  double i = 0;
  double q = 0;
  signed char *x;
  float *y;
  double is,qs;
  int j,k,l=0;
  int nbytes = 2 * nsamples / downsample;
  int nclipped = 0;

  if (floats) 
    y = (float *) malloc(4 * nbytes);
  else
    x = (signed char *) malloc(nbytes);

  for (j = 1, k = 0; j <= nsamples; j += 1, k += 2)
    {
      /* sum Is and Qs */
      i  += inbuf[k];
      q  += inbuf[k+1];

      /* finished coherent sum */
      if (j % downsample == 0)
	{
	  if (floats)
	    {
	      /* scaling is unnecessary, but it makes */
	      /* comparisons with other modes simpler */
	      is = ((i-dcoffi*downsample) * scale);
	      qs = ((q-dcoffq*downsample) * scale);

	      y[l++] = (float) is;
	      y[l++] = (float) qs;
	    }
	  else
	    {
	      /* signed bytes have limited dynamic range */
	      /* scale */
	      is = ((i-dcoffi*downsample) * scale);
	      qs = ((q-dcoffq*downsample) * scale);

	      /* compute clipped */
	      if (is >  127) {is =  127; nclipped++;}
	      if (is < -128) {is = -128; nclipped++;}
	      if (qs >  127) {qs =  127; nclipped++;}
	      if (qs < -128) {qs = -128; nclipped++;}
	      
	      x[l++] = (signed char) is;
	      x[l++] = (signed char) qs;
	    }
	  /* zero accumulator */
	  i = 0;
	  q = 0;
	}
    }
  if (l != nbytes) fprintf(stderr,"oops\n");
  
  /* print diagnostics */
  if (verbose && !floats) 
    fprintf(stderr,"this buffer: output samples %d nclipped %d\n",
	    nbytes,nclipped); 

  /* write it out */
  if (floats)
    {
      if (1 != fwrite(y, 4 * nbytes, 1, fpoutput))
	fprintf(stderr,"Write error\n");
      free(y);
    }
  else
    {
      if (1 != fwrite(x, nbytes, 1, fpoutput))
	fprintf(stderr,"Write error\n");
      free(x);
    }
    
  return;
}    

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,downsample,chan,dcoffi,dcoffq,fudge)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *downsample;
int     *chan;
float   *dcoffi;
float   *dcoffq;
float   *fudge;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_downsample program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:d:c:s:I:Q:f"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_downsample -m mode -d downsampling factor [-s scale fudge factor] [-f output floating point numbers] [-I dcoffi] [-Q dcoffq] [-c channel (1 or 2)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *downsample = 0;
  *chan = 1;
  *dcoffi = 0;
  *dcoffq = 0;
  *fudge = 1;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
    {
    case 'I':
      sscanf(optarg,"%f",dcoffi);
      arg_count += 2;		/* two command line arguments */
      break;
      
    case 'Q':
      sscanf(optarg,"%f",dcoffq);
      arg_count += 2;		/* two command line arguments */
      break;
      
    case 'o':
      *outfile = optarg;	/* output file name */
      arg_count += 2;		/* two command line arguments */
      break;
      
    case 'm':
      sscanf(optarg,"%d",mode);
      arg_count += 2;		/* two command line arguments */
      break;
      
    case 'd':
      sscanf(optarg,"%d",downsample);
      arg_count += 2;
      break;
      
    case 'c':
      sscanf(optarg,"%d",chan);
      arg_count += 2;
      break;
      
    case 's':
      sscanf(optarg,"%f",fudge);
      arg_count += 2;
      break;
      
    case 'f':
      floats = 1;
      arg_count += 1;
      break;

    case '?':			 /*if not in myoptions, getopt rets ? */
      goto errout;
      break;
    }
  }

  if (arg_count < argc)		 /* 1st non-optioned param is infile */
    *infile = argv[arg_count];

  /* must specify a valid mode and downsampling factor */
  if (*mode == 0 || *downsample < 1) goto errout;
  
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

