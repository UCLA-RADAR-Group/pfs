/*******************************************************************************
*  program pfs_dehop
*  $Id: pfs_dehop.c,v 1.7 2009/11/16 19:11:45 jlm Exp $
*  This programs dehops spectra obtained with the pfs_fft program
*  It expects a time series of one or more ffts per dwell time
*  Input data are assumed to be four byte floating point numbers
*
*  usage:
*  	pfs_dehop
*              [-f sampling frequency (KHz)] 
*              [-d dwell time (s)] 
*              [-r frequency resolution (Hz)] 
*	       [-h f0,df,n (KHz)]
*              [-b (binary output)] 
*              [-i (invert frequency axis)] 
*              [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*       the -f argument specifies the data taking sampling frequency in KHz
*	the -d argument specifies the dwell time in seconds
*	the -r argument specifies the fft frequency resolution in Hz
*       the -h argument specifies the hop parameters:
*		starting frequency, frequency increment, number of hops
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log: pfs_dehop.c,v $
   Revision 1.7  2009/11/16 19:11:45  jlm
   Added ifdef flag for Mac compilation.

   Revision 1.6  2002/04/27 20:25:31  margot
   Changed location of include file for O_LARGEFILE

   Revision 1.5  2002/04/27 06:10:57  margot
   Added -i option.

   Revision 1.4  2001/07/07 15:35:24  margot
   Added -b option for binary output.
   Fixed rounding of integer quantities.

   Revision 1.3  2001/07/06 19:34:40  margot
   Added baseline subtraction.

   Revision 1.2  2001/07/06 19:02:34  margot
   Simplified adding.

   Revision 1.1  2001/07/06 18:54:07  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef MAC
#include <fcntl.h>
#else
#include <asm/fcntl.h>
#endif
#include "unpack.h"
#include <fftw.h>

/* revision control variable */
static char const rcsid[] = 
"$Id: pfs_dehop.c,v 1.7 2009/11/16 19:11:45 jlm Exp $";

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();
int  no_comma_in_string();	
void zerofill(float *data, int len);
void cmplxmul(float *in1, float *in2, float *out, int len);

int main(int argc, char *argv[])
{
  float *fftbuf;
  float *total;
  float *baseline;
  int inbufsize;
  int outbufsize;

  double freq;		/* frequency */
  double fsamp;		/* sampling frequency, KHz */
  double freqres;	/* frequency resolution, Hz */
  double f0,df;		/* hop start frequency and increment, KHz */
  int dwell;		/* dwell time, seconds */
  int hops;		/* number of hops */
  int fftsperhop;	/* number of transforms at each hop */
  int fftlen;		/* transform length, complex samples */
  int init;		/* location of first hop in fft array */
  int shift;		/* shift by this many locations between hops */
  int open_flags;	/* flags required for open() call */
  int binary;		/* binary output */
  int inverted;		/* frequency axis inverted */
  int i,j,k,l,m;	/* indices */

  /* get the command line arguments */
  processargs(argc,argv,&infile,&outfile,&fsamp,&freqres,&dwell,&f0,&df,&hops,&binary,&inverted);

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


  /* compute transform parameters */
  fftlen     = (int) rint(fsamp * 1e3 / freqres);
  shift      = (int) rint(df    * 1e3 / freqres);
  fftsperhop = (int) rint(dwell * freqres);
  init       = (int) rint(fftlen / 2.0 + f0 * 1e3 / freqres);
  inbufsize = fftlen * sizeof(float); 
  outbufsize = shift * sizeof(float);

  fprintf(stderr,"\n%s\n\n",command_line);
  fprintf(stderr,"FFT length                     : %d\n",fftlen);
  fprintf(stderr,"Frequency resolution           : %e Hz\n",freqres);
  fprintf(stderr,"Processed bandwidth            : %e Hz\n\n",freqres*fftlen);

  fprintf(stderr,"Data required for one transform: %d bytes\n",inbufsize);
  fprintf(stderr,"Number of ffts per hop         : %d\n",fftsperhop);
  fprintf(stderr,"Data required for one hop seq  : %d bytes\n",
	  inbufsize*fftsperhop*hops);
  fprintf(stderr,"Initial location and shift     : %d,%d\n",init,shift);
  fprintf(stderr,"Dehopped output buffer size    : %d bytes\n",outbufsize);
  fprintf(stderr,"\n");
    
  /* allocate storage */
  total    = (float *) malloc(outbufsize);
  baseline = (float *) malloc(outbufsize);
  fftbuf   = (float *) malloc(inbufsize);
  if (!fftbuf)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* sum transforms */
  zerofill(total, fftlen);
  zerofill(baseline, shift);
  /* read data buffers until EOF */
  while (1)
    {
      for (i = 0; i < hops; i++)
	for (j = 0; j < fftsperhop; j++)
	  {
	    if (inbufsize != read(fdinput, fftbuf, inbufsize))
	      goto write;

	    /* now we split the data array in hop-sized chunks */
	    /* if data  (k==i), sum in array total */
	    /* if noise (k!=i), sum in array baseline */
	    for (k = 0; k < hops; k++)
	      {
		/* compute correct shift index location */
		if (inverted)
		  l = init - shift / 2 + (hops-k-1) * shift;
		else
		  l = init - shift / 2 + k * shift; 
		
		/* need to check bounds on l or have modulo here */

		/* sum data */
		if (k == i)
		  for (m = 0; m < shift; m++)
		    total[m]    += fftbuf[l+m];
		/* sum noise */
		else
		  for (m = 0; m < shift; m++)
		    baseline[m] += fftbuf[l+m];
	      }
	  }
    }

  /* scale baseline array and subtract */
  for (m = 0; m < shift; m++)
    total[m] = total[m] - baseline[m] / (hops - 1);

 write:
  /* write output */
  if (binary)
    {
      if (shift != fwrite(total,sizeof(float),shift,fpoutput))
	fprintf(stderr,"Write error\n");
    }
  else
    for (i = -shift/2; i < shift/2; i++)
      {
	freq = i*freqres;
	fprintf(fpoutput,"%.3f %.2f\n",freq,total[i+shift/2]);
      }

  return 0;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,fsamp,freqres,dwell,f0,df,hops,binary,inverted)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
double  *fsamp;
double  *freqres;
int     *dwell;
double  *f0;
double  *df;
int     *hops;
int     *binary;
int     *inverted;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_dehop program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "f:d:r:h:o:bi"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_dehop [-f sampling frequency (KHz)] [-d dwell time (s)] [-r frequency resolution (Hz)] [-h f0,df,n (KHz)] [-b (binary output)] [-i (invert frequency axis)] [-o outfile] [infile]";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *fsamp = 0;
  *freqres = 1;
  *dwell = 0;
  *f0 = 0;
  *df = 0;
  *hops = 1;
  *binary = 0;
  *inverted = 0;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
      {
      case 'o':
	*outfile = optarg;	/* output file name */
	arg_count += 2;		/* two command line arguments */
	break;
	
      case 'f':
	sscanf(optarg,"%lf",fsamp);
	arg_count += 2;
	break;
	
      case 'r':
	sscanf(optarg,"%lf",freqres);
	arg_count += 2;
	break;

      case 'd':
	sscanf(optarg,"%d",dwell);
	arg_count += 2;
	break;

      case 'b':
	*binary = 1;
	arg_count += 1;
	break;	

      case 'i':
	*inverted = 1;
	arg_count += 1;
	break;	

      case 'h':
	if ( no_comma_in_string(optarg) )
	  {
	    fprintf(stderr,"\nERROR: require comma between -h args\n");
	    goto errout;
	  }
	else
	  {
	    if (sscanf(optarg,"%lf,%lf,%d",f0,df,hops) != 3)
	      goto errout;
	    arg_count += 2;          /* two command line arguments */
	  }
	break;
	
      case '?':			 /*if not in myoptions, getopt rets ? */
	goto errout;
	break;
      }
  }
  
  if (arg_count < argc)		 /* 1st non-optioned param is infile */
    *infile = argv[arg_count];

  return;

  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE1);
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

/******************************************************************************/
/*	no_comma_in_string						      */
/******************************************************************************/
int no_comma_in_string(params)
char    params[];               /* parameter string */
{
  /* searches for a comma in the parameter string which indicates
     multiple arguements; returns TRUE if there is no comma
  */
  int no_comma;         /* flag false if comma found */
  int index;            /* index to string */

  no_comma = 1;

  for (index = 0; no_comma && params[index]!='\0' ; index++)
    no_comma = (params[index] == ',') ? 0 : no_comma;

  return(no_comma);
}

/******************************************************************************/
/*	zerofill							      */
/******************************************************************************/
void zerofill(float *data, int len)
{
  /* zero fills array of floats of size len */

  int i;

  for (i=0; i<len; i++)
    data[i] = 0.0;

  return;
}
