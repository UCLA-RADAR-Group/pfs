/*******************************************************************************
*  program pfs_r2c
*  $Id: $
*  This programs converts a time series of real-valued samples to 
*  IQ samples by setting the imaginary part to zero and optionally
*  applying a frequency downconversion.  If the sampling frequency 
*  is four times larger than the sampling frequency (fs = 4 fc),
*  we apply an efficient downconversion.
*
*  usage:
*  	pfs_r2c
*               -f sampling frequency (MHz)
*              [-x desired frequency offset (Hz)]
*              [-S number of seconds to skip at the beginning of the file]
*              [-d downsampling factor] (supports factor of 2 only)
*              [-r desired frequency resolution (Hz)] (only for buffer size)  
*              [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*       the -f argument specifies the data taking sampling frequency in MHz
*	the -d argument specifies optional downsampling by a factor of 2
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*  Jean-Luc Margot, Feb 2022
*******************************************************************************/

/* 
   $Log: $
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
  "$Id: $";

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[512];	/* command line assembled by processargs */

float *fftinbuf, *fftoutbuf;    /* arrays for transform products */

void processargs();
void open_file();
void copy_cmd_line();
void vector_power(float *data, int len);
void vector_window(float *data, int len);
void swap_freq(float *data, int len);
void swap_iandq(float *data, int len);
void zerofill(float *data, int len);
int  no_comma_in_string();	
void apply_linear_phase(float *data, double freq, double time, double timeint, int nsamples);
void apply_fast_phase(float *data, long n, int nsamples);

int main(int argc, char *argv[])
{
  long bufsize;		/* size of read buffer */
  int nsamples;		/* # of complex samples in each buffer */
  float freq;		/* frequency */
  double fsamp;		/* sampling frequency, MHz */
  double fquarter;	/* a quarter of the sampling frequency, Hz */
  double freqres;	/* frequency resolution, Hz */
  double foff;		/* frequency offset, Hz */
  double timeint;	/* sampling interval */ 
  double time;		/* time */ 
  int fftlen;		/* transform length, complex samples */
  int fftlen2;		/* twice the transform length */
  int downsample;	/* downsampling factor, dimensionless */
  int counter=0;	/* keeps track of number of transforms written */
  int open_flags;	/* flags required for open() call */
  int hanning;		/* apply Hanning window before fft routine */
  float nskipseconds;   /* optional number of seconds to skip at beginning of file */
  long nskipbytes;	/* number of bytes to skip at beginning of file */
  long n;               /* current sample count */               

  int i,j,k;
  short x;

  /* get the command line arguments */
  processargs(argc,argv,&infile,&outfile,&fsamp,&freqres,&downsample,&hanning,&nskipseconds,&foff);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open output file, stdout default */
  open_file(outfile,&fpoutput);

  /* open file input */
  open_flags = O_RDONLY;
  if((fdinput = open(infile, open_flags)) < 0 )
    {
      perror("open input file");
      exit(1);
    }

  /* compute transform parameters */
  fftlen = (int) rint(fsamp / freqres * 1e6);
  fftlen2 = 2 * fftlen;
  bufsize = fftlen * sizeof(float); 
  nsamples = (int) rint(bufsize / 4.0);

  /* describe what we are doing */
  fprintf(stderr,"\n%s\n\n",command_line);
  fprintf(stderr,"FFT length                     : %d\n",fftlen);
  fprintf(stderr,"Frequency resolution           : %e Hz\n",freqres);
  fprintf(stderr,"Processed bandwidth            : %e Hz\n",freqres*fftlen);

  fprintf(stderr,"Data required for one transform: %ld bytes\n",bufsize);
  fprintf(stderr,"Integration time for one sum   : %e s\n",1 / freqres);

  /* fsamp million samples per second during nskipseconds and 4 bytes per sample for real-sampling */
  nskipbytes = (long) rint(fsamp * 1e6 * nskipseconds * 4.0);
  if (nskipseconds != 0)
    {
      fprintf(stderr,"Skipping from BOF              : %f seconds\n",nskipseconds);
      fprintf(stderr,"Skipping from BOF              : %ld bytes\n",nskipbytes);
    }
  /* if xoff = fsamp / 4, can speed things up considerably */
  fquarter = fsamp * 1e6 / 4;
  if (downsample == 2)
    {
      fprintf(stderr,"Downsampling by a factor of 2. \n");  
    }
  if (foff == -fquarter)
    {
      fprintf(stderr,"Using accelerated downconversion due to fs=4fc. \n");  
      if (nsamples % 4 != 0)
        fprintf(stderr,"In this case make sure nsamples %d is a multiple of four.\n", nsamples);  
    }
  fprintf(stderr,"\n");
    
  /* skip unwanted bytes */
  if (nskipbytes != lseek(fdinput, nskipbytes, SEEK_SET))
    {
      fprintf(stderr,"Read error while skipping %ld bytes.  Check file size.\n",nskipbytes);
      exit(1);
    }

  /* allocate storage */
  fftinbuf  = (float *) malloc(fftlen2 * sizeof(float)); /* 2 * bufsize */
  fftoutbuf = (float *) malloc(fftlen2 * sizeof(float));
  if (!fftinbuf || !fftoutbuf)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* setup time counter */
  n = 0;
  time = 0;
  timeint = 1.0 / (fsamp * 1e6);

  /* infinite loop */
  while (1)
    {
      /* read one data buffer       */
      if (bufsize != read(fdinput, fftinbuf, bufsize))
        {
          fprintf(stderr,"Read error or EOF.\n");
          fprintf(stderr,"Wrote %d transforms\n",counter);
          exit(1);
        }

      /* put the real component in the correct locations in the complex array */
      /* and set the imaginary part to zero */
      for (i = 0, j = 0; i < fftlen; i++, j += 2) 
        {
          fftoutbuf[j] = fftinbuf[i];
          fftoutbuf[j+1] = 0;
        }

      /* optionally apply phase rotation and increment time */
      if (foff == -fquarter)
        {
          apply_fast_phase(fftoutbuf,n,nsamples);
          n += nsamples;
        }
      else if (foff != 0)
        {
          apply_linear_phase(fftoutbuf,foff,time,timeint,nsamples);
          time += timeint * nsamples;
        }

      /* write output */
      if (downsample == 2)
        {
          /* first copy outbuf back to inbuf */
          memcpy(fftinbuf,fftoutbuf,fftlen2 * sizeof(float));
          /* then average two consecutive samples */
          for (i = 0, k = 0; i < fftlen; i += 2, k += 2*downsample)
            {
              for (j = 0; j < 2*downsample; j += 2)
                {
                  fftoutbuf[i]   += fftinbuf[k+j];
                  fftoutbuf[i+1] += fftinbuf[k+j+1];
                }
            }
          if (fftlen != fwrite(fftoutbuf,sizeof(float),fftlen,fpoutput))
            fprintf(stderr,"Write error\n");
        }
      else
        if (fftlen2 != fwrite(fftoutbuf,sizeof(float),fftlen2,fpoutput))
          fprintf(stderr,"Write error\n");
      fflush(fpoutput);
      counter++;
    }

  free(fftinbuf);
  free(fftoutbuf);

  return 0;
}

/******************************************************************************/
/*	apply_fast_phase						      */
/******************************************************************************/
void apply_fast_phase(float *data, long n, int nsamples)
/* data		 data array */
/* n             current sample count */
/* nsamples	 number of complex samples in data array */
{
  /* apply a linear phase correction of freq Hz to the data using the
     an accelerated downconversion valid for fs=4fc (Considine 1983)
     See also Generating Complex Baseband and Analytic Bandpass Signals
     Rick Lyons November 2, 2011
     https://www.dsprelated.com/showarticle/153.php
  */
  int i,j;
  float  phase_r, phase_i;	/* real and imag phase components */
  float  data_r, data_i;	/* real and imag data components */

  for (i=0, j=1; i<2*nsamples; )
  {
    phase_r = 1;
    phase_i = 0;

    data_r = data[i];
    data_i = data[j];

    data[i] = data_r;
    data[j] = data_i;

    i+=2;
    j+=2;

    phase_r = 0;
    phase_i = -1;

    data_r = data[i];
    data_i = data[j];

    data[i] = data_i;
    data[j] = -data_r;

    i+=2;
    j+=2;

    phase_r = -1;
    phase_i = 0;

    data_r = data[i];
    data_i = data[j];

    data[i] = -data_r;
    data[j] = -data_i;

    i+=2;
    j+=2;

    phase_r = 0;
    phase_i = 1;

    data_r = data[i];
    data_i = data[j];

    data[i] = -data_i;
    data[j] = data_r;
    
    i+=2;
    j+=2;
    n+=4;
  }
  return;
}

/******************************************************************************/
/*	apply_linear_phase						      */
/******************************************************************************/
void apply_linear_phase(float *data, double freq, double time, double timeint, int nsamples)
/* data		 data array */
/* freq		 linear phase correction, Hz */
/* time		 start time of phase correction */
/* timeint	 sample period */
/* nsamples	 number of complex samples in data array */
{
  /* apply a linear phase correction of freq Hz to the data using the
     sample time spacing timeint
  */
  int i,j;
  float  phase_r, phase_i;	/* real and imag phase components */
  float  data_r, data_i;	/* real and imag data components */
  double phase;			/* phase correction in radians */
  double freq_rad;		/* frequency in rad/sec */
  double t;

  t        = time;
  freq_rad = 2 * M_PI * freq;

  for (i=0, j=1; i<2*nsamples; i+=2, j+=2)
  {
    phase   = freq_rad * t;
    phase_r = (float)cos( phase );
    phase_i = (float)sin( phase );

    data_r = data[i];
    data_i = data[j];

    data[i] = data_r*phase_r - data_i*phase_i;
    data[j] = data_r*phase_i + data_i*phase_r;

    t += timeint;
  }
  return;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,fsamp,freqres,downsample,hanning,nskipseconds,foff)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
double   *fsamp;
double   *freqres;
int     *downsample;
int     *hanning;
float   *nskipseconds;
double   *foff;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_fft program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "f:r:d:o:x:HS:"; /* options to search for :=> argument*/
  char *USAGE1="pfs_r2c -f sampling frequency (MHz) [-d downsampling factor] [-r desired frequency resolution (Hz)] [-w apply Hanning window before transform] [-S number of seconds to skip before applying first FFT] [-o outfile] [infile]";
  char *USAGE2="For phase rotation, also specify [-x desired frequency offset (Hz)] ";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *fsamp = 0;
  *freqres = 1;
  *downsample = 1;
  *hanning = 0;
  *nskipseconds = 0;    /* default is process entire file */

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
	sscanf(optarg,"%d",downsample);
	arg_count += 2;
	break;
	
      case 'S':
	sscanf(optarg,"%f",nskipseconds);
	arg_count += 2;
	break;
	
      case 'x':
	sscanf(optarg,"%lf",foff);
	arg_count += 2;
	break;
	
      case 'H':
	*hanning = 1;
	arg_count += 1;
	break;

      case '?':			 /*if not in myoptions, getopt rets ? */
	goto errout;
	break;
      }
  }
  
  if (arg_count < argc)		 /* 1st non-optioned param is infile */
    *infile = argv[arg_count];

  /* must specify a valid sampling frequency */
  if (*fsamp == 0) 
    {
      fprintf(stderr,"Must specify sampling frequency\n");
      goto errout;
    }

  /* must specify a valid downsampling factor */
  if (*downsample != 1 && *downsample != 2) 
    {
      fprintf(stderr,"Supports only downsampling by 1 or 2.  Use decimate.py for further downsampling.\n");
      goto errout;
    }

  return;

  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE1);
          fprintf(stderr,"%s\n",USAGE2);
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
  int	i;

  command_line[0] = '\n';

  for (i=0; i<argc; i++)
  {
    strcat (command_line, argv[i]);
    strcat (command_line, " ");
  }

  return;
}	

/******************************************************************************/
/*	vector_power							      */
/******************************************************************************/
void vector_power(float *data, int len)
{
  /* detects the complex array data of len complex samples placing the
     result in the bottom len samples of the input array
  */

  int i,j,k;

  for (i=0, j=0, k=1; i<len; i++, j+=2, k+=2)
    data[i] = data[j]*data[j] + data[k]*data[k];

  return;
}

/******************************************************************************/
/*	swap_freq							      */
/******************************************************************************/
void swap_freq(float *data, int len)
{
  /* swaps the location of the + and - frequencies for a continuous
     spectrum, the data array is 2*len samples long
  */
  int i,j;
  float temp;

  for (i=0, j=len; i<len; i++, j++)
    {
      temp    = data[i];
      data[i] = data[j];
      data[j] = temp;
    }

  return;
}

/******************************************************************************/
/*	swap_iandq							      */
/******************************************************************************/
void swap_iandq(float *data, int len)
{
  /* swaps the i and q components of each complex word to reverse the "phase"
     direction, the data array is len complex samples long or 2*len samples long
  */

  int i,j;
  float temp;

  for (i=0, j=1; i < 2*len; i+=2, j+=2)
    {
      temp    = data[i];
      data[i] = data[j];
      data[j] = temp;
    }

  return;
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
