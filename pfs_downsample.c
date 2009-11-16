/*******************************************************************************
*  program pfs_downsample
*  $Id$
*  This programs downsamples data from the portable fast sampler
*  by summing coherently a specified number of consecutive samples
*
*  usage:
*  	pfs_downsample -m mode -d downsampling factor 
*                      [-f scale fudge factor]
*                      [-b output signed bytes (default is floats)]
*                      [-a process all data files (default to 0)
*                      [-I dcoffi] [-Q dcoffq] 
*                      [-c channel] 
*                      [-i swap I/Q] 
*                      [-s number of complex samples to skip] 
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
*       4-byte floating point numbers or signed bytes if -b is used
*******************************************************************************/

/* 
   $Log$
   Revision 3.13  2007/06/14 21:33:27  jlm
   Clarified comment about end of loop condition with -a option.

   Revision 3.12  2007/06/14 20:42:05  jlm
   Introduced wordstoskip to continue straightening out wordstoskip situation.

   Revision 3.11  2007/06/14 19:08:51  jlm
   Starting to restore some order to bytestoskip situation.  More work ahead.

   Revision 3.10  2007/06/14 18:20:51  jlm
   Added quiet mode with -q option for those who desire to suppress
   informational messages.

   Revision 3.9  2007/06/14 18:11:59  jlm
   Restored printing of informational messages with verbose variable.

   Revision 3.8  2007/06/14 17:34:56  jlm
   Added clipping flag for detection and reporting of clipping situations.
   This frees up verbose variable that will be used for other informational
   messages.

   Revision 3.7  2005/02/22 19:53:31  jao
   Added '-i' option for I & Q swapping.

   Revision 3.6  2004/04/19 20:24:45  jao
   Added new option '-a' to process all data files.
   Take care half byte situation when skipping bytes in mode 1.

   Revision 3.5  2003/10/30 23:18:40  cvs
   Added large file support for output file.

   Revision 3.4  2003/05/30 20:52:59  cvs
   Fixed bug in summation with -m 32: j index should increment by 1, not 8.

   Revision 3.3  2003/05/29 22:49:43  cvs
   Fixed bug in bytestoskip: never initialized.
   Fixed bug in printout of buffers used: must be integers.

   Revision 3.2  2003/05/12 21:22:05  cvs
   Fixed bug in myoptions.

   Revision 3.1  2003/04/21 20:30:13  cvs
   Added option to skip a number of samples at the beginning of a file.
   Re-organized command-line options -b, -f, -s.
   Default output is now floating point numbers.

   Revision 3.0  2003/02/22 02:29:18  cvs
   Very extensive surgery by Joseph Jao for faster processing.
   Threads are used to read and process in parallel mode.
   Byte manipulation is used instead of float manipulation.

   Revision 1.10  2002/11/11 18:46:33  cvs
   Multiplied default buffer size by two.

   Revision 1.9  2002/06/05 16:54:38  cvs
   Added info on last buffer size.

   Revision 1.8  2002/05/26 05:07:11  cvs
   Better scheme for handling short buffers, including those at end of file.

   Revision 1.7  2002/05/26 00:43:13  cvs
   Added mode 32 for downsampling floats.

   Revision 1.6  2002/05/02 05:48:34  cvs
   Mode added to usage line

   Revision 1.5  2002/05/01 06:28:08  cvs
   Added downsampling of 8 bit data, modes 3 and 8.

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
#ifdef MAC
#include <fcntl.h>
#else
#include <asm/fcntl.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

int	fdinput;		/* file descriptor for input file */ 
int	fdoutput;		/* file descriptor to output file */

char	command_line[200];	/* command line assembled by processargs */
char	header[40];		/* data file name header */
int	ext;			/* data file extension */
int	open_rflags;	/* flags required for open() call for reading */
int	open_wflags;	/* flags required for open() call for writing */
struct  stat filestat;	/* input file status structure */
int	verbose = 1;    /* verbosity level */
int     clipping = 0;	/* flag for detection and reporting of clipping */
int	floats  = 1;    /* default output format is floating point */
int	allfiles = 0;   /* data file to be processed */
int	swapiq = 0;	/* swap I/Q */
int	downsample;	/* factor by which to downsample */
int     nsamples; 	/* # of complex samples in each buffer */
float	smpwd;		/* # of single pol complex samples in a 4 byte word */
float   wordstoskip=0.0;/* number of 4-byte words to skip */
float   bytestoskip=0.0;/* number of bytes to skip */
float   remainingbytestoskip=0.0;/* number of remaining bytes to skip after lseek call */

int	mode;		/* data acquisition mode */
int     chan;		/* channel to process (1 or 2) for dual pol data */
int	bufsize;	/* input buffer size */
float   scale; 		/* scaling factor to fit in a byte */
float	dcoffi,dcoffq;	/* dc offsets */

/* for thread ID */
pthread_t tid[3];

struct jdata {
    char    *bfrthr1;
    char    *bfrthr2;
    char   *chnthr1;
    char   *chnthr2;
    int	    bytesread;	/* number of bytes read from input file */
};


void *read_buf(void *rdata);
void *proc_buf(void *pdata);
void *iq_downsample (void *pdata);

void processargs();
void copy_cmd_line();


int main(int argc, char *argv[])
{
  unsigned char *buffer1, *buffer2;		/* buffer for packed data */
  char *channel1,*channel2;	/* buffer for unpacked data */
  float maxunpack;	/* maximum unpacked value from libunpack */
  float maxvalue;	/* maximum achievable value by downsampling */
  float fudge;		/* scale fudge factor */
  int samplestoskip;	/* number of complex samples to skip */
  int i;

  char   *outfile;	/* output file name */
  char   *infile;	/* input file name */

  struct jdata cntlbuf;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&downsample,&chan,&dcoffi,&dcoffq,&fudge,&samplestoskip);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* set mode */
  switch (mode)
    {
    case -1:  smpwd = 8; maxunpack =   +3; break;  
    case  1:  smpwd = 8; maxunpack =   +3; break;
    case  2:  smpwd = 4; maxunpack =  +15; break;
    case  3:  smpwd = 2; maxunpack = +255; break; 
    case  5:  smpwd = 4; maxunpack =   +3; break;
    case  6:  smpwd = 2; maxunpack =  +15; break;
    case  8:  smpwd = 2; maxunpack = +255; break;
    case 32:  smpwd = 0.5; maxunpack = +255; break;
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* open input file */
#ifdef LARGEFILE
  open_rflags = O_RDONLY|O_LARGEFILE;
#else
  open_rflags = O_RDONLY;
#endif

  if((fdinput = open(infile, open_rflags)) < 0 )
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
  
  /* test size compatibility */
  if (filestat.st_size % 4 != 0)
    if (verbose) fprintf(stderr,"Warning: file size %lld is not a multiple of 4\n", 
			 filestat.st_size);
  if (filestat.st_size % downsample != 0)
    if (verbose) fprintf(stderr,"Warning: file size %lld not a multiple of dwnsmplng factor\n",
			 filestat.st_size);

  /* skip samples if needed */
  /* but skip along 4-byte boundaries only */
  /* remaining bytes to skip, including half bytes in mode 1 (2c2b), are handled after unpacking */
  if (samplestoskip != 0)
  {
    wordstoskip = samplestoskip / smpwd;
    bytestoskip = wordstoskip * 4;
    if (verbose) fprintf(stderr, "Skipping %d complex samples, equivalent to %.1f words, equivalent to %.1f bytes\n", 
			 samplestoskip, wordstoskip, bytestoskip);
    
    /* skip desired amount of bytes */
    if ((int) bytestoskip != lseek(fdinput, (int) bytestoskip, SEEK_SET))
      {
        perror("lseek");
        fprintf(stderr, "Unable to skip %d bytes\n", (int) bytestoskip);
        exit(1);
      }

    /* compute the number of remaining samples to skip, if any.
       if this number is nonzero, it will be taken care of after unpacking */
    remainingbytestoskip = (int) ((wordstoskip - (int) wordstoskip) * 4);

    /* test new size compatibility */
    if ((filestat.st_size - (int) bytestoskip) % 4 != 0)
      if (verbose) fprintf(stderr,"Warning: file size %lld with %.1f byte skip not a multiple of 4\n",
			   filestat.st_size, bytestoskip);
    if ((filestat.st_size - (int) bytestoskip) % downsample != 0)
      if (verbose) fprintf(stderr,"Warning: file size %lld with %.1f byte skip not a multiple of dwnsmplng factor\n", 
			   filestat.st_size, bytestoskip);
  }

  /* open output file, stdout is default */
#ifdef LARGEFILE
  open_wflags = O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE;
#else
  open_wflags = O_RDWR | O_CREAT | O_TRUNC;
#endif

  if (outfile[0] == '-') {
     fdoutput=1;
  } else if ((fdoutput = open(outfile, open_wflags, 0660)) < 0 )
  {
     perror("open input file");
     exit(1);
  }

  /* compute dynamic range parameters */
  if (verbose) fprintf(stderr,"Downsampling file of size %d kB by %d\n", 
		       (int) (filestat.st_size / 1000), downsample);
  maxvalue = maxunpack * sqrt(downsample);
  scale = fudge * 0.25 * 128 / maxvalue;

  if (!floats && maxvalue > 255)
    {
      fprintf(stderr,"You may have a dynamic range problem\n");
      fprintf(stderr,"Turning clipping mode on so you can detect clipping instances\n");
      clipping = 1;
    }

  /* compute buffer size */
  /* we need a multiple of the downsampling factor, or order 1 MB */
  bufsize = (int) rint(1000000.0/downsample) * downsample;
  /* but the buffer size must be smaller than the file size */
  if (bufsize > (filestat.st_size - (int) bytestoskip)) 
    {
      bufsize = filestat.st_size - (int) bytestoskip;
      if (verbose) fprintf(stderr,"Reducing buffer size to file size minus bytes to skip: %d\n",
			   bufsize);
    }

  if (verbose) fprintf(stderr,"Using %d buffers of size %d\n", 
		       (int) floor((filestat.st_size - bytestoskip) / bufsize), bufsize);
  
  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer1 = (unsigned char *) malloc(bufsize);
  buffer2 = (unsigned char *) malloc(bufsize);

  /* for mode 32, data buffers are transferred as float. Others use char unpack */
  if (mode == 32) {
    channel1 = (char *) malloc(bufsize);
    channel2 = (char *) malloc(bufsize);
  } else {
    channel1 = (char *) malloc(2 * bufsize * smpwd / 4 * sizeof(char));
    channel2 = (char *) malloc(2 * bufsize * smpwd / 4 * sizeof(char));
  }

  if (!channel1 || !channel2 || !buffer1 || !buffer2) 
    fprintf(stderr,"Malloc error\n");


  if (nsamples % downsample != 0)
    fprintf(stderr,"Warning: # samples per buffer %d, downsampling factor %d\n",
	    nsamples,downsample);

    i = 1;

    /* read 1st buffer */
    cntlbuf.bfrthr1 = buffer1;
    cntlbuf.chnthr1 = channel1;
    pthread_create (&tid[0], NULL, read_buf, (void *)&cntlbuf);
    pthread_join (tid[0], NULL);

    /* read 2nd buffer & process 1st unpacking */
    cntlbuf.bfrthr1 = buffer2;
    cntlbuf.chnthr1 = channel2;
    cntlbuf.bfrthr2 = buffer1;
    cntlbuf.chnthr2 = channel1;
    pthread_create (&tid[1], NULL, proc_buf, (void *)&cntlbuf);
    pthread_create (&tid[0], NULL, read_buf, (void *)&cntlbuf);
    pthread_join (tid[0], NULL);
    pthread_join (tid[1], NULL);

    while (1) {
	/* buffer assignment for parallel processing */
	if (i++ % 2) {
	    cntlbuf.bfrthr1 = buffer1;
	    cntlbuf.chnthr1 = channel1;
	    cntlbuf.bfrthr2 = buffer2;
	    cntlbuf.chnthr2 = channel2;
	} else {
	    cntlbuf.bfrthr1 = buffer2;
	    cntlbuf.chnthr1 = channel2;
	    cntlbuf.bfrthr2 = buffer1;
	    cntlbuf.chnthr2 = channel1;
	}

        /* read buffer, process unpacking & downsample */
        pthread_create (&tid[1], NULL, proc_buf, (void *)&cntlbuf);
        pthread_create (&tid[0], NULL, read_buf, (void *)&cntlbuf);
        pthread_create (&tid[2], NULL, iq_downsample, (void *)&cntlbuf);

        pthread_join (tid[1], NULL);
        pthread_join (tid[0], NULL);
        pthread_join (tid[2], NULL);

	/* last buffer? */
	if (cntlbuf.bytesread != bufsize) {
	  break;
	}
    }

    /* last buffer == empty buffer? */
    if (cntlbuf.bytesread != 0) {
	if (i++ % 2) {
	    cntlbuf.bfrthr1 = buffer1;
	    cntlbuf.chnthr1 = channel1;
	    cntlbuf.bfrthr2 = buffer2;
	    cntlbuf.chnthr2 = channel2;
	} else {
	    cntlbuf.bfrthr1 = buffer2;
	    cntlbuf.chnthr1 = channel2;
	    cntlbuf.bfrthr2 = buffer1;
	    cntlbuf.chnthr2 = channel1;
	}

        /* process unpacking & downsample */
        pthread_create (&tid[1], NULL, proc_buf, (void *)&cntlbuf);
        pthread_create (&tid[2], NULL, iq_downsample, (void *)&cntlbuf);

        pthread_join (tid[1], NULL);
        pthread_join (tid[2], NULL);

        /* prepare to handle small size for the last buffer */
	bufsize = cntlbuf.bytesread;

	nsamples = (int) rint(bufsize * smpwd / 4.0);
	if (verbose) fprintf(stderr,"And one buffer of size %d\n", bufsize);
    }

    /* only process last buffer read size */
    if (i++ % 2) {
	cntlbuf.bfrthr1 = buffer1;
	cntlbuf.chnthr1 = channel1;
    } else {
	cntlbuf.bfrthr1 = buffer2;
	cntlbuf.chnthr1 = channel2;
    }
	
    /* process downsample only */
    pthread_create (&tid[2], NULL, iq_downsample, (void *)&cntlbuf);
    pthread_join (tid[2], NULL);


    /* clean up */
    free (buffer1);
    free (buffer2);
    free (channel1);
    free (channel2);

    close (fdinput);
    close (fdoutput);

  return 0;
}


void *read_buf (void *rdata) {
    struct jdata *rbuf = (struct jdata *)rdata;
    int datasz;
    char infile[80];

    /* read one buffer */
    if        ((rbuf->bytesread = read (fdinput, rbuf->bfrthr1, bufsize)) == -1) {
	perror("read");
	return;
    /* 03/05/04 SWJ - need more data from next data file? */
    } else if (rbuf->bytesread != bufsize && allfiles == 1) {
	close (fdinput);

	sprintf (infile, "%s.%3.3d", &header[0], ++ext);
	if (verbose) fprintf(stderr, "downsampling next file in sequence: \"%s\"\n", infile);

	/* the program must stop after reading the file with the highest extension */
	/* we assume that if the file with the next extension does not exist, 
	   we have reached the end */
	if ((fdinput = open(infile, open_rflags)) < 0) {
	    return;
	} else if (fstat (fdinput, &filestat) < 0) {
	    perror("fstat");
	    return;
	}

	/* make up the data buffer */
    	if ((datasz = read (fdinput, rbuf->bfrthr1 + rbuf->bytesread, bufsize - rbuf->bytesread )) == -1) {
	    perror ("read");
	    return;
	} else {
	    rbuf->bytesread = datasz + rbuf->bytesread;
	} 
    }
}


void *proc_buf (void *pdata) {
    struct jdata *pbuf = (struct jdata *)pdata;

    /* unpack and downsample */
    switch (mode)
      {
        case 1:
          unpack_pfs_2c2b (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          break;
        case 2:
          unpack_pfs_2c4b (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          break;
        case 3:
          unpack_pfs_2c8b (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          break;
        case 5:
          if (chan == 2) {
	    unpack_pfs_4c2b_lcp (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          } else {
            unpack_pfs_4c2b_rcp (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          }
          break;
        case 6:
          if (chan == 2) {
            unpack_pfs_4c4b_lcp (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
          } else {
            unpack_pfs_4c4b_rcp (pbuf->bfrthr2, pbuf->chnthr2, bufsize);
	  }
          break;
        case 8:
        case 32:
          memcpy (pbuf->chnthr2, pbuf->bfrthr2, bufsize);
          break;
        default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
      }
}


/******************************************************************************/
/*	iq_downsample							      */
/******************************************************************************/

void *iq_downsample (void *pdata) 
{
  struct jdata *pbuf = (struct jdata *)pdata;

  char	*inbuf  = (char *) pbuf->chnthr1; 
  float iq[2];

  /* accumulator larger enough to not cause overflow on all downsampled data */
  int	is  = 0,   qs  = 0;	/* is, qs  : char  accumulators for I & Q */
  float isf = 0.0, qsf = 0.0;	/* isf, qsf: float accumulators for I & Q */

  signed char *x;
  float *y;

  int j, k=0, l=0;
  int nbytes = 2 * nsamples / downsample;
  int nclipped = 0;
  int bcnt;

  float iscale = dcoffi * downsample * scale;
  float qscale = dcoffq * downsample * scale;

  if (floats) 
    y = (float *) malloc(4 * nbytes);
  else
    x = (signed char *) malloc(nbytes);

  bcnt = nsamples / downsample;

  /* 03/05/04 SWJ - need to skip 1st sample ?  */
  if (remainingbytestoskip > 0.0) {
    j = remainingbytestoskip;
    if (verbose) fprintf(stderr,"***** Skipping %d extra sample ***** \n", j);
    while (j > 0) {
      bcnt --;
      j --;

      /* skip I & Q */
      *inbuf++;
      *inbuf++;
    }

    /* byte skipping on begining of data segment only */
    remainingbytestoskip = 0.0;
  }

  for (; bcnt > 0; bcnt--)
  {
    if (mode == 32) {
	for (j = 0, isf = 0.0, qsf = 0.0; j < downsample; j += 1, k += 8) {
	  memcpy (&iq[0], &inbuf[k], 8);

	  /* sum Is and Qs */
	  isf  += iq[0];
	  qsf  += iq[1];
	}
    } else {
	for (j = 0, is = 0, qs = 0; j < downsample; j++) {
	  /* sum Is and Qs */
	  is  += *inbuf++;
	  qs  += *inbuf++;
	}

	isf = (float) is;
	qsf = (float) qs;
    }

    /* finished coherent sum */
    /* SWJ 12/07/04 added IQ swapping capability */
    if (floats)
    {
      /* scaling is unnecessary, but it makes */
      /* comparisons with other modes simpler */
      if (swapiq == 0) {
        y[l++] = scale * isf - iscale;
        y[l++] = scale * qsf - qscale;
      } else {
        y[l++] = scale * qsf - qscale;
        y[l++] = scale * isf - iscale;
      }
    }
    else
    {
      /* signed bytes have limited dynamic range */
      /* scale */
      is = scale * isf - iscale;
      qs = scale * qsf - qscale;

      /* compute clipped */
      if (is >  127) {is =  127; nclipped++;}
      if (is < -128) {is = -128; nclipped++;}
      if (qs >  127) {qs =  127; nclipped++;}
      if (qs < -128) {qs = -128; nclipped++;}
	      
      if (swapiq == 0) {
        x[l++] = (signed char) is;
        x[l++] = (signed char) qs;
      } else {
        x[l++] = (signed char) qs;
        x[l++] = (signed char) is;
      }
    }
  }

  /* 03/05/04 SWJ  if (l != nbytes) fprintf(stderr,"oops\n"); */
  
  /* print diagnostics */
  if (clipping && !floats) 
    fprintf(stderr,"this buffer: output samples %d nclipped %d\n",
	    nbytes,nclipped); 

  /* write it out */
  /* SWJ - replaced all nbytes with l */
  if (floats)
    {
      if (write(fdoutput, y, 4 * l) != 4 * l) perror ("Write floats");
      free(y);
    }
  else
    {
      if (write(fdoutput, x, l) != l) perror ("Write bytes");
      free(x);
    }

  return;
}    


/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,downsample,chan,dcoffi,dcoffq,fudge,samplestoskip)
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
int 	*samplestoskip;
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

  char *myoptions = "m:o:d:c:s:I:Q:b:f:axq"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_downsample -m mode -d downsampling factor [-s number of complex samples to skip] [-f scale fudge factor] [-b output byte quantities (default floats)] [-a downsample all data files] [-I dcoffi] [-Q dcoffq] [-c channel (1 or 2)] [-i (swap I/Q)] [-q (quiet mode)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t32: 32bit floats\n";
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
  *samplestoskip = 0;
  floats = 1;
  allfiles = 0;
  swapiq = 0;
  verbose = 1;

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
      *outfile = optarg;        /* output file name */
      arg_count += 2;           /* two command line arguments */
      break;

    case 'm':
      sscanf(optarg,"%d",mode);
      arg_count += 2;           /* two command line arguments */
      break;

    case 'd':
      sscanf(optarg,"%d",downsample);
      arg_count += 2;
      break;

    case 'c':
      sscanf(optarg,"%d",chan);
      arg_count += 2;
      break;

    case 'f':
      sscanf(optarg,"%f",fudge);
      arg_count += 2;
      break;

    case 'b':
      floats = 0;
      arg_count += 1;
      break;

    case 'a':
      allfiles = 1;
      arg_count += 1;
      break;

    case 'i':
      swapiq = 1;
      arg_count += 1;
      break;

    case 'q':
      verbose = 0;
      arg_count += 1;
      break;

    case 's':
      sscanf(optarg,"%d",samplestoskip);
      arg_count += 2;           /* two command line arguments */
      break;
  
    case '?':                    /*if not in myoptions, getopt rets ? */
      goto errout;
      break;
    }
  }   
      
  if (arg_count < argc)          /* 1st non-optioned param is infile */
    *infile = argv[arg_count];
    
  /* SWJ 03/05/04 */
  if (allfiles == 1) {          
        /* get data file information */
        /* get extension part of 'filename.ext' */
        strcpy (header, strrchr (*infile, '.'));
        ext = atoi (&header[1]);
      
        /* get 'name' part */
        strncpy (header, *infile, strlen(*infile) - strlen(header));
  } 

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

