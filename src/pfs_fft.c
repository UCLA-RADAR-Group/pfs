/*******************************************************************************
*  program pfs_fft
*  $Id: pfs_fft.c,v 4.2 2020/05/21 17:44:12 jlm Exp $
*  This programs performs spectral analysis on data acquired with
*  the portable fast sampler
*
*  usage:
*  	pfs_fft -m mode 
*               -f sampling frequency (MHz)
*              [-r desired frequency resolution (Hz)]  
*              [-d downsampling factor] 
*              [-n sum n transforms] 
*              [-l (dB output)]
*              [-b (binary output)]
*              [-t time series] 
*              [-x freqmin,freqmax (Hz)]
*              [-s scale to sigmas using smin,smax (Hz)]
*              [-c channel] 
*              [-i swap IQ before transform (invert freq axis)]
*              [-H apply Hanning window before transform]
*              [-C file of Chebyshev polynomial coefficients defining window to apply after transform] 
*              [-S number of seconds to skip before applying first FFT]
*              [-I dcoffi] [-Q dcoffq] 
*              [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*       the -f argument specifies the data taking sampling frequency in MHz
*	the -d argument specifies the factor by which to downsample the data
*	                (coherent sum before fft simulates lower sampling fr)
*	the -r argument specifies the desired frequency resolution in Hz
*       the -n argument specifies how many transforms to add
*			(incoherent sum after fft)
*       the -l argument specifies logarithmic (dB) output
*	the -t option indicates that (sums of) transforms ought to be written
*			one after the other until EOF
*       the -x option specifies an optional range of output frequencies
*       the -c argument specifies which channel (1 or 2) to process
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*  Jean-Luc Margot, Aug 2000
*******************************************************************************/

/* 
   $Log: pfs_fft.c,v $
   Revision 4.2  2020/05/21 17:44:12  jlm
   Added capability to read file of Chebyshev coefficients.

   Revision 4.0  2020/05/21 05:50:00  jlm
   Upgraded to FFTW3

   Revision 3.13  2017/03/05 05:35:27  jlm
   Added -C option to allow for application of Chebyshev polynomial after FFT.

   Revision 3.12  2016/12/24 22:21:06  jlm
   Improved computation of scaling to sigmas: (1) check bounds, (2) compute
   loop indices, (3) exclude points deviating by >3.5 sigmas from computation.

   Revision 3.11  2012/01/09 00:32:32  jlm
   Fixed uninitialized nskipbytes when nskipseconds = 0.

   Revision 3.10  2012/01/08 23:56:09  jlm
   Added -S option to allow skipping data at beginning of file.

   Revision 3.9  2009/11/16 19:09:36  jlm
   Added ifdef flag for Mac compilation

   Revision 3.8  2007/06/19 15:43:15  jao
   Increase dynamic range for parameters bufsize and sum

   Revision 3.7  2007/06/14 17:10:01  jlm
   Added option to apply Hanning window.

   Revision 3.6  2006/09/16 22:08:28  jlm
   Added -b option for output in binary format.  Previously the -t option
   was used to get a similar output, but with more restrictions.

   Revision 3.5  2003/11/19 01:33:24  cvs
   Improved output format.

   Revision 3.4  2003/11/09 03:28:12  cvs
   Added -i option for IQ swap.

   Revision 3.3  2003/09/18 19:51:43  cvs
   Abort if input file cannot be opened.

   Revision 3.2  2003/09/18 02:52:44  cvs
   Improved command-line option error handling.

   Revision 3.1  2003/02/25 22:45:34  cvs
   Sizes of fftinbuf and fftoutbuf must be identical.

   Revision 3.0  2003/02/22 07:48:08  cvs
   Switch to Joseph Jao's byte unpacking.

   Revision 1.12  2002/12/28 07:42:59  cvs
   Added capability to swap I and Q prior to FFT.

   Revision 1.11  2002/06/05 16:55:10  cvs
   Corrected one error message.

   Revision 1.10  2002/05/12 13:36:02  cvs
   Cosmetic changes.

   Revision 1.9  2002/05/02 05:48:59  cvs
   Mode added to usage line

   Revision 1.8  2002/04/27 20:22:47  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.7  2002/04/27 06:16:26  margot
   Added mode 16.

   Revision 1.6  2001/07/06 19:48:11  margot
   Added -s option for scaling to sigmas.

   Revision 1.5  2001/07/06 04:25:28  margot
   Added mode 32 for "unpacking" of 4 byte floating point quantities.

   Revision 1.4  2001/07/04 00:14:57  margot
   Added mode 8 for unpacking signed bytes.

   Revision 1.3  2000/10/30 05:38:44  margot
   Added variable nsamples.

   Revision 1.2  2000/09/15 22:10:38  margot
   Added -l and -x options.

   Revision 1.1  2000/09/15 19:36:37  margot
   Initial revision
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
#include <fftw3.h>

/* revision control variable */
static char const rcsid[] = 
"$Id: pfs_fft.c,v 4.2 2020/05/21 17:44:12 jlm Exp $";

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */
char   *chebfile;	        /* file of Chebyshev coefficients */

char	command_line[512];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();
void vector_power(float *data, int len);
void vector_window(float *data, int len);
void chebyshev_window(float *data, int len, double *chebcoeff, int degree);
void swap_freq(float *data, int len);
void swap_iandq(float *data, int len);
void zerofill(float *data, int len);
int  no_comma_in_string();	
double chebeval(double x, double c[], int degree);
int  read_cheb_coeffs(char *chebfile, double *chebcoeff);
void average(float *inbuf, int nsamples, double *i, double *q);

int main(int argc, char *argv[])
{
  int mode;
  long bufsize;		/* size of read buffer */
  char *buffer;		/* buffer for packed data */
  char *rcp;		/* buffer for unpacked data */
  float smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int levels;		/* # of levels for given quantization mode */
  int degree=0;         /* degree of Chebyshev polynomial, default none */

  float *fftinbuf, *fftoutbuf;
  float *total;

  double *chebcoeff;    /* array for polynomial coefficients */

  float freq;		/* frequency */
  float freqmin;	/* min frequency to output */
  float freqmax;	/* max frequency to output */
  float value;		/* value to output */
  float rmsmin;		/* min frequency for rms calculation */
  float rmsmax;		/* max frequency for rms calculation */
  double fsamp;		/* sampling frequency, MHz */
  double freqres;	/* frequency resolution, Hz */
  double mean,mean1;	/* needed for rms computation */
  double var,var1;	/* needed for rms computation */
  double sigma,sigma1;	/* needed for rms computation */
  int downsample;	/* downsampling factor, dimensionless */
  long long sum;	/* number of transforms to add, dimensionless */
  int timeseries;	/* process as time series, boolean */
  int dB;		/* write out results in dB */
  int fftlen;		/* transform length, complex samples */
  int chan;		/* channel to process (1 or 2) for dual pol data */
  int counter=0;	/* keeps track of number of transforms written */
  int open_flags;	/* flags required for open() call */
  int invert;		/* swap i and q before fft routine */
  int hanning;		/* apply Hanning window before fft routine */
  int swap = 1;		/* swap frequencies at output of fft routine */
  int binary;		/* write output as binary floating point quantities */
  float nskipseconds;   /* optional number of seconds to skip at beginning of file */
  long nskipbytes;	/* number of bytes to skip at beginning of file */
  int imin,imax;	/* indices for rms calculation */
  float	dcoffi,dcoffq;	/* user-provided dc offsets */
  int dcoffset=0;	/* compute and remove DC offset prior to FFT */
  
  fftwf_plan p;
  int i,j,k,l,n,n1;
  short x;

  /* get the command line arguments */
  processargs(argc,argv,&infile,&outfile,&mode,&fsamp,&freqres,&downsample,&sum,&binary,&timeseries,&chan,&freqmin,&freqmax,&rmsmin,&rmsmax,&dB,&invert,&hanning,&chebfile,&nskipseconds,&dcoffi,&dcoffq,&dcoffset);

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

  /* read Cheb coefficients, if requested */
  if (chebfile[0] != '-') 
    {
      chebcoeff = (double *) malloc(64 * sizeof(double));   /* allocate up to 64 coefficients */
      degree = read_cheb_coeffs(chebfile, chebcoeff);       /* read coeffs and return degree */
    }

  switch (mode)
    {
    case  -1: smpwd = 8; break;  
    case   1: smpwd = 8; break;
    case   2: smpwd = 4; break;
    case   3: smpwd = 2; break; 
    case   5: smpwd = 4; break;
    case   6: smpwd = 2; break;
    case   7: smpwd = 1; break;
    case   8: smpwd = 2; break; 
    case  16: smpwd = 1; break; 
    case  32: smpwd = 0.5; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* compute transform parameters */
  fftlen = (int) rint(fsamp / freqres * 1e6);
  bufsize = fftlen * 4 / smpwd; 
  fftlen = fftlen / downsample;

  /* describe what we are doing */
  fprintf(stderr,"\n%s\n\n",command_line);
  fprintf(stderr,"FFT length                     : %d\n",fftlen);
  fprintf(stderr,"Frequency resolution           : %e Hz\n",freqres);
  fprintf(stderr,"Processed bandwidth            : %e Hz\n",freqres*fftlen);
  if (rmsmin != 0 || rmsmax != 0)
    fprintf(stderr,"Scaling to rms power between   : [%e,%e] Hz\n\n",rmsmin,rmsmax);

  fprintf(stderr,"Data required for one transform: %ld bytes\n",bufsize);
  fprintf(stderr,"Number of transforms to add    : %qd\n",sum);
  fprintf(stderr,"Data required for one sum      : %qd bytes\n",sum * bufsize);
  fprintf(stderr,"Integration time for one sum   : %e s\n",sum / freqres);
  
  nskipbytes = (long) rint(fsamp * 1e6 * nskipseconds * 4.0 / smpwd);
  if (nskipseconds != 0)
    {
      fprintf(stderr,"Skipping from BOF              : %f seconds\n",nskipseconds);
      fprintf(stderr,"Skipping from BOF              : %ld bytes\n",nskipbytes);
    }
  if (chebfile[0] != '-')
    fprintf(stderr, "Degree of Chebyshev polynomial : %d\n",degree);    
  /* for (i = 0; i <= degree; i++) fprintf(stderr, "%d %lf\n", i, chebcoeff[i]); */
  fprintf(stderr,"\n");
    
  /* skip unwanted bytes */
  /* fsamp samples per second during nskipseconds, and 4/smpwd bytes per complex sample */
  if (nskipbytes != lseek(fdinput, nskipbytes, SEEK_SET))
    {
      fprintf(stderr,"Read error while skipping %ld bytes.  Check file size.\n",nskipbytes);
      exit(1);
    }

  /* verify that scaling request is sensible */
  if (rmsmin != 0 || rmsmax != 0)
    {
      if (rmsmin > rmsmax || rmsmin < -freqres*fftlen/2 || rmsmax > freqres*fftlen/2)
	{
	  fprintf(stderr,"Problem with -s parameters\n");
	  exit(1);
	}
    }

  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer    = (char *)  malloc(bufsize);
  fftinbuf  = (float *) malloc(2 * fftlen * sizeof(float));
  fftoutbuf = (float *) malloc(2 * fftlen * sizeof(float));
  total = (float *) malloc(fftlen * sizeof(float));
  rcp   = (char *)  malloc(2 * nsamples * sizeof(char));
  if (!buffer || !fftinbuf || !fftoutbuf || !total || !rcp)
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* compute fft plan */
  p = fftwf_plan_dft_1d(fftlen, (fftwf_complex *)fftinbuf, (fftwf_complex *)fftoutbuf, FFTW_FORWARD, FFTW_ESTIMATE);

  /* label used if time series is requested */
 loop:

  /* sum transforms */
  zerofill(total, fftlen);
  for (i = 0; i < sum; i++)
    {
      /* initialize fft array to zero */
      zerofill(fftinbuf, 2 * fftlen);
      
      /* read one data buffer       */
      if (bufsize != read(fdinput, buffer, bufsize))
	{
	  fprintf(stderr,"Read error or EOF.\n");
	  if (timeseries) fprintf(stderr,"Wrote %d transforms\n",counter);
	  exit(1);
	}

      /* unpack */
      switch (mode)
	{
	case 1:
	  unpack_pfs_2c2b(buffer, rcp, bufsize); 
	  break;
	case 2: 
	  unpack_pfs_2c4b(buffer, rcp, bufsize);
	  break;
	case 3: 
	  unpack_pfs_2c8b(buffer, rcp, bufsize);
	  break;
	case 5:
	  if (chan == 2) unpack_pfs_4c2b_lcp (buffer, rcp, bufsize);
	  else 		 unpack_pfs_4c2b_rcp (buffer, rcp, bufsize);
	  break;
	case 6: 
	  if (chan == 2) unpack_pfs_4c4b_lcp (buffer, rcp, bufsize);
	  else 		 unpack_pfs_4c4b_rcp (buffer, rcp, bufsize);
	  break;
	case 7:
	  if (chan == 2) unpack_pfs_4c8b_lcp (buffer, rcp, bufsize);
	  else 		 unpack_pfs_4c8b_rcp (buffer, rcp, bufsize);
	  break;
     	case 8: 
	  memcpy (rcp, buffer, bufsize);
	  break;
	case 16: 
	  for (i = 0, j = 0; i < bufsize; i+=sizeof(short), j++)
	    {
	      memcpy(&x,&buffer[i],sizeof(short));
	      fftinbuf[j] = (float) x;
	    }
	  break;
     	case 32: 
	  memcpy(fftinbuf,buffer,bufsize);
	  break;
	default: 
	  fprintf(stderr,"Mode not implemented yet\n"); 
	  exit(-1);
	}

      if (dcoffset)
	average(fftinbuf, fftlen, &dcoffi, &dcoffq);
	  
      /* deal with nonzero DC offsets if provided by user or if option -D was invoked */
      if (dcoffi != 0 || dcoffq != 0)
	for (k = 0; k < 2*fftlen; k += 2)
	  {
	    fftinbuf[k]   -= dcoffi;
	    fftinbuf[k+1] -= dcoffq; 
	  }
      
      /* downsample */
      if (mode != 16 && mode != 32)
	for (k = 0, l = 0; k < 2*fftlen; k += 2, l += 2*downsample)
	  {
	    for (j = 0; j < 2*downsample; j+=2)
	      {
		fftinbuf[k]   += (float) rcp[l+j];
		fftinbuf[k+1] += (float) rcp[l+j+1];
	      }
	  }

      /* transform, swap, and compute power */
      if (invert) swap_iandq(fftinbuf,fftlen); 
      if (hanning) vector_window(fftinbuf,fftlen);
      fftwf_execute(p); 
      if (swap) swap_freq(fftoutbuf,fftlen); 
      vector_power(fftoutbuf,fftlen);
      
      /* sum transforms */
      for (j = 0; j < fftlen; j++)
	total[j] += fftoutbuf[j];
    }
  
  /* set DC to average of neighboring values  */
  total[fftlen/2] = (total[fftlen/2-1]+total[fftlen/2+1]) / 2.0; 

  /* apply Chebyshev to detected power if needed */
  if (degree) chebyshev_window(total,fftlen,chebcoeff,degree);
  
  /* compute rms if needed */
  mean = 0;
  sigma = 1;
  if (rmsmin != 0 || rmsmax != 0)
    {
      /* identify relevant indices for rms power computation */
      imin = fftlen/2 + rmsmin/freqres; 
      imax = fftlen/2 + rmsmax/freqres; 
      mean1 = var1 = 0;
      n1 = 0;
      for (i = imin; i < imax; i++)
	{
	  mean1 += total[i];
	  var1  += total[i] * total[i];
	  n1++;
	}
      mean1  = mean1 / n1;
      var1   = var1 / n1;
      sigma1 = sqrt(var1 - mean1 * mean1);

      /* now redo calculation but exclude 3-sigma outliers */
      mean = var = 0;
      n = 0;
      for (i = imin; i < imax; i++)
	{
	  if (fabs((total[i] - mean1)/sigma1) > 3.5)
	    continue;
	  mean += total[i];
	  var  += total[i] * total[i];
	  n++;
	}
      mean  = mean / n;
      var   = var / n;
      sigma = sqrt(var - mean * mean);
      /*
      fprintf(stderr,"Computed mean,sigma with    %d outliers : %e +/- %e s\n",(n1-n),mean1,sigma1);
      fprintf(stderr,"Computed mean,sigma without %d outliers : %e +/- %e s\n",(n1-n),mean,sigma);
      */
    }
  
  /* write output */
  /* either time series */
  if (timeseries)
    {
      for (i = 0; i < fftlen; i++) total[i] = (total[i]-mean)/sigma;
      if (fftlen != fwrite(total,sizeof(float),fftlen,fpoutput))
	fprintf(stderr,"Write error\n");
      fflush(fpoutput);
      counter++;
      goto loop;
    }
  /* or standard output */
  /* or limited frequency range */
  else
    for (i = 0; i < fftlen; i++)
      {
	  freq = (i-fftlen/2)*freqres;
    
	  if ((freqmin == 0.0 && freqmax == 0.0) || (freq >= freqmin && freq <= freqmax)) 
	  {
	    value = (total[i]-mean)/sigma;
	    if (dB) value = 10*log10(value);

	    if (binary)
	      fwrite(&value,sizeof(float),1,fpoutput);
	    else
	      fprintf(fpoutput,"% .3f % .3e\n",freq,value);  
	  }
      }
  
  fftwf_destroy_plan(p);
  free(fftinbuf);
  free(fftoutbuf);
  
  return 0;
}

/******************************************************************************/
/*    average         							      */
/******************************************************************************/
void average(float *inbuf, int nsamples, double *i, double *q)
{
  int k;

  *i = 0;
  *q = 0;
  
  /* sum Is and Qs */
  for (k = 0; k < 2*nsamples; k += 2)
    {
      *i  += inbuf[k];
      *q  += inbuf[k+1];
    }

  /* divide by nsamples to get the average value */
  *i = *i / nsamples;
  *q = *q / nsamples;

  return;
}    
/******************************************************************************/
/*	vector_window							      */
/******************************************************************************/
void vector_window(float *data, int len)
{
  /* Hanning windows the data array of length 'len' (complex samples)
  */
  float  weight;		/* calculated weight */
  double n_minus_1;		/* weight calculation scale */
  int    i,j,k;

  n_minus_1 = 1.0/(double)(len - 1);

  for (i=0, j=0, k=1; i<len; i++, j+=2, k+=2)
  {
    weight = (float)(0.5 - 0.5 * cos( 2 * M_PI * (double)i * n_minus_1 ) );
    data[j] *= weight;
    data[k] *= weight;
  }
  return;
}

/******************************************************************************/
/*	chebyshev_window						      */
/******************************************************************************/
void chebyshev_window(float *data, int len, double *chebcoeff, int degree)
{
  /* Applies Chebyshev window the data array of length 'len' (floating point samples)
  */
  double x;
  double weight;		/* calculated weight */
  int    i;

  for (i=0; i<len; i++)
  {
    x = -0.5 + (double) i / (double) len;
    weight = chebeval(x, chebcoeff, degree);
    data[i] /= weight;
  }
  return;
}

/******************************************************************************/
/*	chebeval							      */
/******************************************************************************/
double chebeval(double x, double c[], int degree)
{
  /* 
     Evaluate a Chebyshev series at points x.
     Expects an array `c` of length at least degree + 1 = n + 1
     This function returns the value:
     p(x) = c_0 * T_0(x) + c_1 * T_1(x) + ... + c_n * T_n(x)
  */
  double x2;
  double c0, c1;
  double tmp;
  int i;
  
  x2 = 2*x;
  c0 = c[degree-1];
  c1 = c[degree];
  for (i = 2; i <= degree; i++)    
    {
      tmp = c0;
      c0 = c[degree-i] - c1;
      c1 = tmp + c1*x2;
    }
  return c0 + c1*x;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,fsamp,freqres,downsample,sum,binary,timeseries,chan,freqmin,freqmax,rmsmin,rmsmax,dB,invert,hanning,chebfile,nskipseconds,dcoffi,dcoffq,dcoffset)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
double   *fsamp;
double   *freqres;
int     *downsample;
long long     *sum;
int     *binary;
int     *timeseries;
int     *chan;
float   *freqmin;
float   *freqmax;
float   *rmsmin;
float   *rmsmax;
int     *dB;
int     *invert;
int     *hanning;
char    **chebfile;
float   *nskipseconds;
float   *dcoffi;
float   *dcoffq;
int     *dcoffset;
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

  char *myoptions = "m:f:d:r:n:tc:o:lbx:s:iHC:S:I:Q:D"; /* options to search for :=> argument*/
  char *USAGE1="pfs_fft -m mode -f sampling frequency (MHz) [-r desired frequency resolution (Hz)] [-d downsampling factor] [-n sum n transforms] [-l (dB output)] [-b (binary output)] [-t time series] [-x freqmin,freqmax (Hz)] [-s scale to sigmas using smin,smax (Hz)] [-c channel (1 or 2)] [-i swap IQ before transform (invert freq axis)] [-w apply Hanning window before transform] [-C file of Chebyshev polynomial coefficients defining window to apply after transform] [-S number of seconds to skip before applying first FFT] [-I dcoffi] [-Q dcoffq] [-D compute and remove DC offset prior to FFT] [-o outfile] [infile]";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t16: signed 16bit\n\t32: 32bit floats\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *fsamp = 0;
  *freqres = 1;
  *downsample = 1;
  *sum = 1;
  *binary = 0;
  *timeseries = 0;
  *chan = 1;
  *dB = 0;		/* default is linear output */
  *invert = 0;
  *hanning = 0;
  *chebfile = "-";
  *nskipseconds = 0;    /* default is process entire file */
  *freqmin = 0;		/* not set value */
  *freqmax = 0;		/* not set value */
  *rmsmin  = 0;		/* not set value */
  *rmsmax  = 0;		/* not set value */
  *dcoffi = 0;
  *dcoffq = 0;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
      {
      case 'I':
	sscanf(optarg,"%f",dcoffi);
	arg_count += 2;           /* two command line arguments */
	break;

      case 'Q':
	sscanf(optarg,"%f",dcoffq);
	arg_count += 2;           /* two command line arguments */
	break;

      case 'o':
	*outfile = optarg;	/* output file name */
	arg_count += 2;		/* two command line arguments */
	break;
	
      case 'm':
	sscanf(optarg,"%d",mode);
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
	
      case 'n':
	sscanf(optarg,"%qd",sum);
	arg_count += 2;
	break;

      case 'c':
	sscanf(optarg,"%d",chan);
	arg_count += 2;
	break;
	
      case 'S':
	sscanf(optarg,"%f",nskipseconds);
	arg_count += 2;
	break;
	
      case 'l':
	*dB = 1;
	arg_count += 1;
	break;

      case 'i':
	*invert = 1;
	arg_count += 1;
	break;

      case 'H':
	*hanning = 1;
	arg_count += 1;
	break;

      case 'C':
	*chebfile = optarg;	/* file name for Cheb coefficients */
	arg_count += 2;		/* two command line arguments */
	break;

      case 'D':
	*dcoffset = 1;		/* compute and remove DC offset prior to FFT */
	arg_count += 2;		/* two command line arguments */
	break;

      case 'b':
	*binary = 1;
	arg_count += 1;
	break;

      case 't':
	*timeseries = 1;
	arg_count += 1;
	break;

      case 'x':
	if ( no_comma_in_string(optarg) )
	  {
	    fprintf(stderr,"\nERROR: require comma between -x args\n");
	    goto errout;
	  }
	else
	  {
	    if (sscanf(optarg,"%f,%f",freqmin,freqmax) != 2)
	      goto errout;
	    arg_count += 2;          /* two command line arguments */
	  }
	break;
	
      case 's':
	if ( no_comma_in_string(optarg) )
	  {
	    fprintf(stderr,"\nERROR: require comma between -s args\n");
	    goto errout;
	  }
	else
	  {
	    if (sscanf(optarg,"%f,%f",rmsmin,rmsmax) != 2)
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

  /* must specify a valid mode */
  if (*mode == 0)
    {
      fprintf(stderr,"Must specify sampling mode\n");
      goto errout;
    }
  /* must specify a valid sampling frequency */
  if (*fsamp == 0) 
    {
      fprintf(stderr,"Must specify sampling frequency\n");
      goto errout;
    }
  /* must specify a valid channel */
  if (*chan != 1 && *chan != 2) goto errout;
  /* some combinations not implemented yet */
  if (*timeseries && *dB) 
    {
      fprintf(stderr,"Cannot have -t and -l simultaneously yet\n");
      goto errout;
    }
  if (*timeseries && (*freqmin != 0 || *freqmax !=0)) 
    {
      fprintf(stderr,"Cannot have -t and -x simultaneously yet\n");
      goto errout;
    }
  if (*downsample > 1 && (*mode == 16 || *mode == 32)) 
    {
      fprintf(stderr,"Cannot have -d with modes 16 or 32 yet\n");
      exit (1);
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
/*	read cheb coeffs						      */
/******************************************************************************/
int     read_cheb_coeffs(char *chebfile, double *chebcoeff)
{
  int degree = 0;
  FILE   *fpcheb;		/* pointer to file of Cheb coefficients */

  /* open the Cheb coeff file */
  fpcheb=fopen(chebfile,"r");
  if (fpcheb == NULL)
    {
      perror("read_cheb_coeff: cheb coefficients file open error");
      exit(1);
    }

  /* read coefficients */
  while (fscanf(fpcheb, "%lf", &chebcoeff[degree++]) == 1)
    continue;

  degree = degree - 1; /* EOF does not count */
  degree = degree - 1; /* number of coeffs = degree + 1 */

  return degree;
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
