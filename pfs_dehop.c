/*******************************************************************************
*  program pfs_dehop
*  $Id$
*  This programs dehops spectra obtained with the pfs_fft program
*
*  usage:
*  	pfs_dehop
*              [-f sampling frequency (KHz)] 
*              [-d dwell time (s)] 
*              [-r frequency resolution (Hz)] 
*	       [-h f0,df,n (KHz)]
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
   $Log$
   Revision 1.1  2001/07/06 18:54:07  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "unpack.h"
#include <fftw.h>

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
int  no_comma_in_string();	
void zerofill(float *data, int len);
void cmplxmul(float *in1, float *in2, float *out, int len);

int main(int argc, char *argv[])
{
  float *fftbuf;
  float *total;
  int bufsize;

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
  int binary = 1;	/* binary output */
  int inverted = 1;	/* frequency axis inverted */
  int i,j,k,l;

  /* get the command line arguments */
  processargs(argc,argv,&infile,&outfile,&fsamp,&freqres,&dwell,&f0,&df,&hops);

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
  fftlen = fsamp * 1e3 / freqres;
  bufsize = fftlen * sizeof(float); 
  fftsperhop = dwell * freqres;
  init  = fftlen / 2 + f0 * 1e3 / freqres;
  shift = df * 1e3 / freqres;

  fprintf(stderr,"\n%s\n\n",command_line);
  fprintf(stderr,"FFT length                     : %d\n",fftlen);
  fprintf(stderr,"Frequency resolution           : %e Hz\n",freqres);
  fprintf(stderr,"Processed bandwidth            : %e Hz\n\n",freqres*fftlen);

  fprintf(stderr,"Data required for one transform: %d bytes\n",bufsize);
  fprintf(stderr,"Number of ffts per hop         : %d\n",fftsperhop);
  fprintf(stderr,"Data required for one hop seq  : %d bytes\n",
	  bufsize*fftsperhop*hops);
  fprintf(stderr,"Initial location and shift     : %d,%d\n",init,shift);
  fprintf(stderr,"\n");
    
  /* allocate storage */
  fftbuf = (float *) malloc(bufsize);
  total  = (float *) malloc(bufsize);
  if (!total)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* sum transforms */
  zerofill(total, fftlen);
  /* read one data buffer */
  while (1)
    {
      for (i = 0; i < hops; i++)
	for (j = 0; j < fftsperhop; j++)
	  {
	    if (bufsize != read(fdinput, fftbuf, bufsize))
	      goto write;
	    /* shift to correct location */
	    if (inverted)
	      l = init - shift / 2 + (hops-i-1) * shift;
	    else
	      l = init - shift / 2 + i * shift; 
	    /* sum */
	    for (k = 0; k < shift; k++)
	      total[k]   += fftbuf[l+k];
	  }
    }

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
	fprintf(fpoutput,"%.3f %.0f\n",freq,total[i+shift/2]);
      }

  return 0;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,fsamp,freqres,dwell,f0,df,hops)
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

  char *myoptions = "f:d:r:h:o:"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_dehop [-f sampling frequency (KHz)] [-d dwell time (s)] [-r frequency resolution (Hz)] [-h f0,df,n (KHz)] [-o outfile] [infile]";
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
