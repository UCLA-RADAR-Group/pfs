/*******************************************************************************
*  program pfs_fft_2
*  $Id: pfs_fft_2.c,v 4.2 2020/05/21 17:47:53 jlm Exp $
*  This programs performs spectral analysis on data acquired with the portable
*  fast sampler (PFS), JPL clones of the PFS, and other data-taking devices.
*  It sums the powers obtained in two channels.
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
*              [-h fch1, write output in HDF5 format with starting frequency fch1 (MHz)]
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
*  Jean-Luc Margot, Mar 2017, based on pfs_fft 3.13
*******************************************************************************/

/* 
   $Log: pfs_fft_2.c,v $
   Revision 4.2  2020/05/21 17:47:53  jlm
   Added capability to read file of Chebyshev coefficients.

   Revision 4.0  2020/05/21 05:50:17  jlm
   Upgraded to FFTW3

   Revision 1.2  2017/05/26 04:33:07  jlm
   Allowed for non-integer values in nskipseconds.

   Revision 1.1  2017/03/05 06:34:27  jlm
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
#include <hdf5.h>

/* revision control variable */
static char const rcsid[] = 
"$Id: pfs_fft_2.c,v 4.2 2020/05/21 17:47:53 jlm Exp $";

#define UNDEFINED 0.987654321

FILE   *fpoutput;		/* pointer to output file */
int	fdinput1;		/* file descriptor for input file 1 */
int	fdinput2;		/* file descriptor for input file 2 */

char   *outfile;		/* output file name */
char   *infile1;	        /* input file name 1 */
char   *infile2;	        /* input file name 2 */
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
hid_t createHDF5File(const char* filename, size_t rows, size_t cols, double tsum, double freqres, double fch1);
void writeFloatLineToHDF5(hid_t dataset_id, const float* line_data, size_t current_row, size_t cols);

int main(int argc, char *argv[])
{
  int mode;
  long bufsize;		/* size of read buffer */
  char *buffer1;	/* buffer 1 for packed data */
  char *buffer2;	/* buffer 2 for packed data */
  char *rcp;		/* buffer for unpacked data */
  char *lcp;		/* buffer for unpacked data */
  float smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int levels;		/* # of levels for given quantization mode */
  int degree=0;         /* degree of Chebyshev polynomial, default none */

  float *fftinbuf1, *fftoutbuf1;
  float *fftinbuf2, *fftoutbuf2;
  float *total1,*total2;
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
  double tsum;		/* integration time for one sum */
  int timeseries;	/* process as time series, boolean */
  int dB;		/* write out results in dB */
  int fftlen;		/* transform length, complex samples */
  int fftout;		/* number of output (summed) transforms */
  int chan;		/* channel to process (1 or 2) for dual pol data */
  int counter=0;	/* keeps track of number of transforms written */
  int open_flags;	/* flags required for open() call */
  int invert;		/* swap i and q before fft routine */
  int hanning;		/* apply Hanning window before fft routine */
  double hdf5;		/* write output file in HDF5 format with starting frequency fch1 (Hz) */
  int swap = 1;		/* swap frequencies at output of fft routine */
  int binary;		/* write output as binary floating point quantities */
  float nskipseconds;     /* optional number of seconds to skip at beginning of file */
  long nskipbytes;	/* number of bytes to skip at beginning of file */
  long long inbytes;	/* size of input file in bytes */
  int imin,imax;	/* indices for rms calculation */
  
  fftwf_plan p1;
  fftwf_plan p2;
  int i,j,k,l,n,n1;
  short x;

  hid_t dataset_id;
  struct stat fileStat;

  /* get the command line arguments */
  processargs(argc,argv,&infile1,&infile2,&outfile,&mode,&fsamp,&freqres,&downsample,&sum,&binary,&timeseries,&chan,&freqmin,&freqmax,&rmsmin,&rmsmax,&dB,&invert,&hanning,&hdf5,&chebfile,&nskipseconds);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open file input */
  open_flags = O_RDONLY;
  if((fdinput1 = open(infile1, open_flags)) < 0 )
    {
      perror("open input file");
      exit(1);
    }
  if((fdinput2 = open(infile2, open_flags)) < 0 )
    {
      perror("open input file");
      exit(1);
    }
  if (fstat(fdinput1, &fileStat) == -1) {
      perror("fstat input file");
      exit(1);
    }
  inbytes = fileStat.st_size;

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
    case   8: smpwd = 2; break;
    case  16: smpwd = 1; break;
    case  32: smpwd = 0.5; break;
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* compute transform parameters */
  fftlen = (int) rint(fsamp / freqres * 1e6);
  bufsize = fftlen * 4 / smpwd;
  fftlen = fftlen / downsample;
  tsum = sum / freqres;

  /* compute number of output transforms */
  nskipbytes = (long) rint(fsamp * 1e6 * nskipseconds * 4.0 / smpwd);
  fftout = (inbytes - nskipbytes) / bufsize / sum;

  /* open output file, stdout default */
  if (hdf5 == UNDEFINED)
    open_file(outfile,&fpoutput);
  else
    dataset_id = createHDF5File(outfile, fftout, fftlen, tsum, freqres, hdf5);

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
  fprintf(stderr,"Integration time for one sum   : %e s\n",tsum);

  if (nskipseconds != 0)
    {
      fprintf(stderr,"Skipping from BOF              : %f seconds\n",nskipseconds);
      fprintf(stderr,"Skipping from BOF              : %ld bytes\n",nskipbytes);
    }
  fprintf(stderr,"Number of output (summed) ffts : %qd\n",fftout);

  if (chebfile[0] != '-')
    fprintf(stderr, "Degree of Chebyshev polynomial : %d\n",degree);
  /* for (i = 0; i <= degree; i++) fprintf(stderr, "%d %lf\n", i, chebcoeff[i]); */
  fprintf(stderr,"\n");
    
  /* skip unwanted bytes */
  /* fsamp samples per second during nskipseconds, and 4/smpwd bytes per complex sample */
  if (nskipbytes != lseek(fdinput1, nskipbytes, SEEK_SET))
    {
      fprintf(stderr,"Read error while skipping %ld bytes.  Check file size.\n",nskipbytes);
      exit(1);
    }
  if (nskipbytes != lseek(fdinput2, nskipbytes, SEEK_SET))
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
  buffer1    = (char *)  malloc(bufsize);
  buffer2    = (char *)  malloc(bufsize);
  fftinbuf1  = (float *) malloc(2 * fftlen * sizeof(float));
  fftinbuf2  = (float *) malloc(2 * fftlen * sizeof(float));
  fftoutbuf1 = (float *) malloc(2 * fftlen * sizeof(float));
  fftoutbuf2 = (float *) malloc(2 * fftlen * sizeof(float));
  total1 = (float *) malloc(fftlen * sizeof(float));
  total2 = (float *) malloc(fftlen * sizeof(float));
  total = (float *) malloc(fftlen * sizeof(float));
  rcp   = (char *)  malloc(2 * nsamples * sizeof(char));
  lcp   = (char *)  malloc(2 * nsamples * sizeof(char));
  if (!buffer2 || !fftinbuf2 || !fftoutbuf2 || !total || !rcp || !lcp)
    {
      fprintf(stderr,"Malloc error\n");
      exit(1);
    }

  /* compute fft plan */
  p1 = fftwf_plan_dft_1d(fftlen, (fftwf_complex *)fftinbuf1, (fftwf_complex *)fftoutbuf1, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftwf_plan_dft_1d(fftlen, (fftwf_complex *)fftinbuf2, (fftwf_complex *)fftoutbuf2, FFTW_FORWARD, FFTW_ESTIMATE);

  /* label used if time series is requested */
 loop:

  /* sum transforms */
  zerofill(total1, fftlen);
  zerofill(total2, fftlen);
  zerofill(total, fftlen);
  for (i = 0; i < sum; i++)
    {
      /* initialize fft array to zero */
      zerofill(fftinbuf1, 2 * fftlen);
      zerofill(fftinbuf2, 2 * fftlen);
      
      /* read one data buffer       */
      if (bufsize != read(fdinput1, buffer1, bufsize))
	{
	  fprintf(stderr,"Read error or EOF.\n");
	  if (timeseries) fprintf(stderr,"Wrote %d transforms\n",counter);
	  exit(1);
	}
      if (bufsize != read(fdinput2, buffer2, bufsize))
	{
	  fprintf(stderr,"Read error or EOF.\n");
	  if (timeseries) fprintf(stderr,"Wrote %d transforms\n",counter);
	  exit(1);
	}

      /* unpack */
      switch (mode)
	{
	case 1:
	  unpack_pfs_2c2b(buffer1, rcp, bufsize);
	  unpack_pfs_2c2b(buffer2, lcp, bufsize);
	  break;
	case 2: 
	  unpack_pfs_2c4b(buffer1, rcp, bufsize);
	  unpack_pfs_2c4b(buffer2, lcp, bufsize);
	  break;
	case 3: 
	  unpack_pfs_2c8b(buffer1, rcp, bufsize);
	  unpack_pfs_2c8b(buffer2, lcp, bufsize);
	  break;
	case 5:
	  unpack_pfs_4c2b_rcp (buffer1, rcp, bufsize);
	  unpack_pfs_4c2b_lcp (buffer2, lcp, bufsize);
	  break;
	case 6: 
	  unpack_pfs_4c4b_rcp (buffer1, rcp, bufsize);
	  unpack_pfs_4c4b_lcp (buffer2, lcp, bufsize);
	  break;
     	case 8: 
	  memcpy (rcp, buffer1, bufsize);
	  memcpy (lcp, buffer2, bufsize);
	  break;
	case 16: 
	  for (i = 0, j = 0; i < bufsize; i+=sizeof(short), j++)
	    {
	      memcpy(&x,&buffer1[i],sizeof(short));
	      fftinbuf1[j] = (float) x;
	      memcpy(&x,&buffer2[i],sizeof(short));
	      fftinbuf2[j] = (float) x;
	    }
	  break;
     	case 32: 
	  memcpy(fftinbuf1,buffer1,bufsize);
	  memcpy(fftinbuf2,buffer2,bufsize);
	  break;
	default: 
	  fprintf(stderr,"Mode not implemented yet\n");
	  exit(-1);
	}

      /* downsample */
      if (mode != 16 && mode != 32)
	for (k = 0, l = 0; k < 2*fftlen; k += 2, l += 2*downsample)
	  {
	    for (j = 0; j < 2*downsample; j+=2)
	      {
		fftinbuf1[k]   += (float) rcp[l+j];
		fftinbuf1[k+1] += (float) rcp[l+j+1];
		fftinbuf2[k]   += (float) lcp[l+j];
		fftinbuf2[k+1] += (float) lcp[l+j+1];
	      }
	  }

      /* transform, swap, and compute power */
      if (invert) swap_iandq(fftinbuf1,fftlen);
      if (invert) swap_iandq(fftinbuf2,fftlen);
      if (hanning) vector_window(fftinbuf1,fftlen);
      if (hanning) vector_window(fftinbuf2,fftlen);
      fftwf_execute(p1);
      fftwf_execute(p2);
      if (swap) swap_freq(fftoutbuf1,fftlen);
      if (swap) swap_freq(fftoutbuf2,fftlen);
      vector_power(fftoutbuf1,fftlen);
      vector_power(fftoutbuf2,fftlen);
      
      /* sum transforms */
      for (j = 0; j < fftlen; j++)
	{
	  total1[j] += fftoutbuf1[j];
	  total2[j] += fftoutbuf2[j];
	}
    }
  
  /* set DC to average of neighboring values  */
  total1[fftlen/2] = (total1[fftlen/2-1]+total1[fftlen/2+1]) / 2.0;
  total2[fftlen/2] = (total2[fftlen/2-1]+total2[fftlen/2+1]) / 2.0;

  /* sum powers */
  for (j = 0; j < fftlen; j++)
    total[j] = total1[j] + total2[j];

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
      if (hdf5 == UNDEFINED)
	{
	  if (fftlen != fwrite(total,sizeof(float),fftlen,fpoutput))
	    fprintf(stderr,"Write error\n");
	  fflush(fpoutput);
	}
      else
	{
	  if (counter < fftout)
	    writeFloatLineToHDF5(dataset_id, total, counter, fftlen);
	}
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
	      if (hdf5 == UNDEFINED)
		fwrite(&value,sizeof(float),1,fpoutput);
	      else
		/* writeFloatValueToHDF5(dataset_id, &value); */
		fprintf(stderr, "Writing a single transform to HDF5 is not implemented yet.  Try with -t -n 1.\n");
	    else
	      fprintf(fpoutput,"% .3f % .3e\n",freq,value);
	  }
      }

  /* close files */
  if (hdf5 == UNDEFINED)
    {
      fclose(fpoutput);
    }
  else
    {
      H5Dclose(dataset_id);
      H5Fclose(H5Iget_file_id(dataset_id));
    }
  
  fftwf_destroy_plan(p1);
  fftwf_destroy_plan(p2);
  free(fftinbuf1);
  free(fftoutbuf1);
  free(fftinbuf2);
  free(fftoutbuf2);
  
  return 0;
}

/******************************************************************************/
/*	writeFloatLineToHDF5						      */
/******************************************************************************/
void writeFloatLineToHDF5(hid_t dataset_id, const float* line_data, size_t current_row, size_t cols)
{
  hid_t dataspace_id, memspace_id;
  herr_t status;

  // Get the dataspace for the dataset
  dataspace_id = H5Dget_space(dataset_id);

  // Select the hyperslab for the new line
  hsize_t offset[3] = {current_row, 0, 0};
  hsize_t count[3] = {1, 1, cols};
  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
  if (status < 0)
    {
      fprintf(stderr,"HDF5 select hyperslab failed with error code %d\n", status);
      exit(1);
    }

  // Create the memory dataspace for the line of data
  hsize_t mem_dims[3] = {1, 1, cols};
  memspace_id = H5Screate_simple(3, mem_dims, NULL);

  // Write the line data to the dataset
  status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, line_data);
  if (status < 0)
    {
      fprintf(stderr,"HDF5 write failed with error code %d\n", status);
      exit(1);
    }

  // Close resources
  H5Sclose(dataspace_id);
}

/******************************************************************************/
/*	createHDF5File							      */
/******************************************************************************/
/* The structure and attributes of this HDF5 file are odd.
   They are meant to replicate Breakthrough Listen dynamic spectra. */
hid_t createHDF5File(const char* filename, size_t rows, size_t cols, double tsum, double freqres, double fch1)
{
  hid_t file_id, dataspace_id, dataset_id, attribute_id, attribute_type_id;
    int status;

    // Create a new HDF5 file
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create the dataspace for the dataset
    int rank = 3;
    hsize_t dims[3] = {rows, 1, cols};
    dataspace_id = H5Screate_simple(rank, dims, NULL);

    // Create the dataset
    dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_FLOAT, dataspace_id,
			    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create a scalar dataspace for the attributes
    dataspace_id = H5Screate(H5S_SCALAR);

    // Create the "nchans" attribute
    attribute_type_id = H5T_STD_I64LE;
    attribute_id = H5Acreate2(dataset_id, "nchans", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "nchans" attribute value
    int64_t nchans_attribute_value = cols;
    H5Awrite(attribute_id, attribute_type_id, &nchans_attribute_value);

    // Close the "nchans" attribute resources
    H5Aclose(attribute_id);

    // Create the "nbits" attribute
    attribute_type_id = H5T_STD_I64LE;
    attribute_id = H5Acreate2(dataset_id, "nbits", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "nbits" attribute value
    int64_t nbits_attribute_value = 32;
    H5Awrite(attribute_id, attribute_type_id, &nbits_attribute_value);

    // Close the "nbits" attribute resources
    H5Aclose(attribute_id);

    // Create the "data_type" attribute
    attribute_type_id = H5T_STD_I64LE;
    attribute_id = H5Acreate2(dataset_id, "data_type", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "data_type" attribute value
    int64_t data_type_attribute_value = 1;
    H5Awrite(attribute_id, attribute_type_id, &data_type_attribute_value);

    // Close the "data_type" attribute resources
    H5Aclose(attribute_id);

    // Create the "fch1" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "fch1", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "fch1" attribute value
    // Example values 1701.5625 for VEGAS band C channel 0, 1601.5625 for VEGAS band D channel 0
    double fch1_attribute_value = fch1;
    H5Awrite(attribute_id, attribute_type_id, &fch1_attribute_value);

    // Close the "fch1" attribute resources
    H5Aclose(attribute_id);

    // Create the "foff" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "foff", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "foff" attribute value
    double foff_attribute_value = -freqres / 1e6; /* must be in MHz to match fch1 units */
    H5Awrite(attribute_id, attribute_type_id, &foff_attribute_value);

    // Close the "foff" attribute resources
    H5Aclose(attribute_id);

    // Create the "tsamp" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "tsamp", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "tsamp" attribute value
    double tsamp_attribute_value = tsum;
    H5Awrite(attribute_id, attribute_type_id, &tsamp_attribute_value);

    // Close the "tsamp" attribute resources
    H5Aclose(attribute_id);

    // Create the "tstart" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "tstart", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "tstart" attribute value
    double tstart_attribute_value = 59332.933969907404;
    H5Awrite(attribute_id, attribute_type_id, &tstart_attribute_value);

    // Close the "tstart" attribute resources
    H5Aclose(attribute_id);

    // Create the "nifs" attribute
    attribute_type_id = H5T_STD_I64LE;
    attribute_id = H5Acreate2(dataset_id, "nifs", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "nifs" attribute value
    int64_t nifs_attribute_value = 1;
    H5Awrite(attribute_id, attribute_type_id, &nifs_attribute_value);

    // Close the "nifs" attribute resources
    H5Aclose(attribute_id);

    // Create the "telescope_id" attribute
    attribute_type_id = H5T_STD_I64LE;
    attribute_id = H5Acreate2(dataset_id, "telescope_id", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the "telescope_id" attribute value
    int64_t telescope_id_attribute_value = 6; /* for GBT */
    H5Awrite(attribute_id, attribute_type_id, &telescope_id_attribute_value);

    // Close the "telescope_id" attribute resources
    H5Aclose(attribute_id);

    // Create the "src_dej" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "src_dej", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    // Write the "src_dej" attribute value
    double src_dej_attribute_value = 0.0;
    H5Awrite(attribute_id, attribute_type_id, &src_dej_attribute_value);

    // Close the "src_dej" attribute resources
    H5Aclose(attribute_id);

    // Create the "src_raj" attribute
    attribute_type_id = H5T_IEEE_F64LE;
    attribute_id = H5Acreate2(dataset_id, "src_raj", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    // Write the "src_raj" attribute value
    double src_raj_attribute_value = 0.0;
    H5Awrite(attribute_id, attribute_type_id, &src_raj_attribute_value);

    // Close the "src_raj" attribute resources
    H5Aclose(attribute_id);

    // Create the string datatype for the attribute
    attribute_type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(attribute_type_id, H5T_VARIABLE);
    H5Tset_strpad(attribute_type_id, H5T_STR_NULLTERM);
    H5Tset_cset(attribute_type_id, H5T_CSET_ASCII);

    // Create the attribute
    attribute_id = H5Acreate2(file_id, "CLASS", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the attribute value
    const char* class_attribute_value = "FILTERBANK";
    H5Awrite(attribute_id, attribute_type_id, &class_attribute_value);

    // Close resources
    H5Aclose(attribute_id);

    // Create the attribute
    attribute_id = H5Acreate2(file_id, "VERSION", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the attribute value
    const char* version_attribute_value = "1.0";
    H5Awrite(attribute_id, attribute_type_id, &version_attribute_value);

    // Close resources
    H5Aclose(attribute_id);

    // Create the attribute
    attribute_id = H5Acreate2(dataset_id, "source_name", attribute_type_id,
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the attribute value
    const char* src_attribute_value = "test";
    H5Awrite(attribute_id, attribute_type_id, &src_attribute_value);

    // Close resources
    H5Aclose(attribute_id);
    H5Sclose(dataspace_id);

    // Create a string attribute for the dimension labels
    const char* col_labels[3] = {"frequency", "feed_id", "time"};

    // Create the dataspace for the dataset
    rank = 1;
    hsize_t dim[1] = {3};
    dataspace_id = H5Screate_simple(rank, dim, NULL);

    // Create the attribute
    attribute_id = H5Acreate2(dataset_id, "DIMENSION_LABELS", attribute_type_id,
                                dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write the attribute value
    H5Awrite(attribute_id, attribute_type_id, col_labels);

    // Close resources
    H5Aclose(attribute_id);
    H5Sclose(dataspace_id);
    H5Tclose(attribute_type_id);

    return dataset_id;
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
    /* fprintf(stderr,"%lf %lf %lf\n",weight,data[i],data[i]/weight); */
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
void	processargs(argc,argv,infile1,infile2,outfile,mode,fsamp,freqres,downsample,sum,binary,timeseries,chan,freqmin,freqmax,rmsmin,rmsmax,dB,invert,hanning,hdf5,chebfile,nskipseconds)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile1;		 /* input file name 1 */
char	**infile2;		 /* input file name 2 */
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
double  *hdf5;
char    **chebfile;
float     *nskipseconds;
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

  char *myoptions = "m:f:d:r:n:tc:h:o:lbx:s:iHC:S:"; /* options to search for :=> argument*/
  char *USAGE1="pfs_fft_2 -m mode -f sampling frequency (MHz) [-r desired frequency resolution (Hz)] [-d downsampling factor] [-n sum n transforms] [-l (dB output)] [-b (binary output)] [-t time series] [-x freqmin,freqmax (Hz)] [-s scale to sigmas using smin,smax (Hz)] [-c channel (1 or 2)] [-i swap IQ before transform (invert freq axis)] [-H apply Hanning window before transform] [-C file of Chebyshev polynomial coefficients defining window to apply after transform] [-S number of seconds to skip before applying first FFT] [-h fch1, write output in HDF5 format with starting frequency fch1 (MHz)] [-o outfile] infile1 infile2";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t16: signed 16bit\n\t32: 32bit floats\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile1  = "-";		 /* initialise to stdin, stdout */
  *infile2  = "-";		 /* initialise to stdin, stdout */
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
  *hdf5 = UNDEFINED;
  *chebfile = "-";
  *nskipseconds = 0;    /* default is process entire file */
  *freqmin = 0;		/* not set value */
  *freqmax = 0;		/* not set value */
  *rmsmin  = 0;		/* not set value */
  *rmsmax  = 0;		/* not set value */

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

      case 'h':
	sscanf(optarg,"%lf",hdf5);
	*binary = 1;
	arg_count += 2;
	break;

      case 'C':
 	*chebfile = optarg;	/* file name for Cheb coefficients */
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
  
  if (arg_count < argc)		 /* 1st non-optioned param is infile1 */
    *infile1 = argv[arg_count];
  arg_count += 1;

  if (arg_count < argc)		 /* 2nd non-optioned param is infile2 */
    *infile2 = argv[arg_count];
  arg_count += 1;
  
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
