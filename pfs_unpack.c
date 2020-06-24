/*******************************************************************************
*  program pfs_unpack
*  $Id: pfs_unpack.c,v 3.4 2009/11/16 19:07:49 jlm Exp $
*  This programs unpacks data from the portable fast sampler,
*  optionally applying a phase rotation to compensate for a frequency offset.
*  Default output is binary floating point quantities.
*
*  usage:
*  	pfs_unpack -m mode 
*                  [-a output ascii to stdout]
*                  [-d (detect and output magnitude)] 
*                  [-p (detect and output power)] 
*                  [-c channel] 
*                  [-o outfile] [infile]
*  for phase rotation, also specify
*                  [-f sampling frequency (MHz)]
*                  [-x desired frequency offset (Hz)]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m argument specifies the data acquisition mode
*       the -c argument specifies which channel (1 or 2) to process
*       the -a option allows text output instead of binary output
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log: pfs_unpack.c,v $
   Revision 3.4  2009/11/16 19:07:49  jlm
   Added ifdef flag for Mac compilation

   Revision 3.3  2007/06/18 20:50:06  jao
   Added option '-o -' & '-' for piping out/in to/from STDOUT & STDIN

   Revision 3.2  2003/09/18 02:44:54  cvs
   Now printing decimal places for mode 32.

   Revision 3.1  2003/05/29 22:51:40  cvs
   Fixed bug in allocation of rcp[]: should be 2*nsamples long
   Fixed bug in detect power: should operate on outbuf[], not rcp[]
   Fixed bug in transferring rcp to outbuf: loop should be 2*nsamples

   Revision 3.0  2003/02/22 02:09:08  cvs
   Revised by Joseph Jao to speed unpacking.

   Revision 2.8  2002/09/29 17:58:23  cvs
   Added -p option for square-law detection and power output

   Revision 2.7  2002/07/25 22:38:29  cvs
   Added mode 16 for 16-bit signed integers.

   Revision 2.6  2002/06/05 22:22:55  cvs
   Changed computation of outbufsize with -d option.

   Revision 2.5  2002/05/26 04:21:42  cvs
   Better scheme for handling short buffers, including those at end of file.

   Revision 2.4  2002/05/26 03:59:18  cvs
   Removed write statement unintentionally left over from earlier revision.

   Revision 2.3  2002/05/26 00:39:30  cvs
   Removed buggy adjustment to input buffer size

   Revision 2.2  2002/05/12 15:50:11  cvs
   Added detect option and better checks on buffer size.

   Revision 2.1  2002/05/12 13:42:26  cvs
   Added mode 32 for floats.

   Revision 2.0  2002/05/02 05:54:16  cvs
   Added capability to apply phase rotation with -f and -x options.
   Added apply_linear_phase() function.
   Now checking input file size and selecting a compatible buffer size.

   Revision 1.5  2002/05/02 03:38:25  cvs
   Added unpacking of mode 8, signed bytes.

   Revision 1.4  2002/04/27 20:21:44  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.3  2000/11/01 02:18:51  margot
   Added some fragile form of text output.

   Revision 1.2  2000/10/30 21:32:09  margot
   Added -c option

   Revision 1.1  2000/10/30 21:24:35  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <fcntl.h>
#else
#include <asm/fcntl.h>
#endif
#include <sys/stat.h>
#include <unistd.h>
#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id: pfs_unpack.c,v 3.4 2009/11/16 19:07:49 jlm Exp $";

int     fdoutput;		/* file descriptor for output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();
void apply_linear_phase(float *data, double freq, double time, double timeint, int nsamples);

int main(int argc, char *argv[])
{
  struct stat filestat;	/* input file status structure */
  int bufsize = 1000000;/* size of read buffer, default 1 MB, unless input file is smaller */
  int outbufsize;	/* output buffer size */
  int bytesread;	/* number of bytes read from input file */
  char *buffer;		/* buffer for packed data */
  char *rcp;		/* char buffer for unpacked data */
  float *outbuf;	/* float buffer for unpacked data */
  double fsamp;		/* sampling frequency, MHz */
  double foff;		/* frequency offset, Hz */
  double timeint;	/* sampling interval */ 
  double time;		/* time */ 
  int mode;
  float smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int chan;		/* channel to process (1 or 2) for dual pol data */
  int ascii;		/* text output */
  int mdetect;		/* magnitude output */
  int pdetect;		/* power output */
  char *format;		/* print format */
  int i,j;
  
  format = (char *) malloc(100);

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&chan,&ascii,&mdetect,&pdetect,&fsamp,&foff);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open file input */
#ifdef LARGEFILE
  if (infile[0] == '-')
    fdinput=0;
  else if((fdinput = open(infile, O_RDONLY|O_LARGEFILE)) < 0 )
    perror("open input file");

  if (outfile[0] == '-')
    fdoutput=1;
  else if((fdoutput = open(outfile, O_WRONLY|O_CREAT|O_LARGEFILE, 0644)) < 0 )
    perror("open output file");
#else
  if (infile[0] == '-')
    fdinput=1;
  else if((fdinput = open(infile, O_RDONLY)) < 0 )
    perror("open input file");

  if (outfile[0] == '-')
    fdoutput=1;
  else if((fdoutput = open(outfile, O_WRONLY|O_CREAT, 0644)) < 0 )
    perror("open output file");
#endif

  /* check file size */
  if (fstat (fdinput, &filestat) < 0)
    {
      perror("input file status");
      exit(1);
    }
  if (filestat.st_size % 4 != 0)
    fprintf(stderr,"Warning: file size %d is not a multiple of 4\n",
	    filestat.st_size);

  switch (mode)
    {
    case -1: smpwd = 8; break;  
    case  1: smpwd = 8; break;
    case  2: smpwd = 4; break;
    case  3: smpwd = 2; break; 
    case  5: smpwd = 4; break;
    case  6: smpwd = 2; break;
    case  8: smpwd = 2; break;
    case 16: smpwd = 1; break; 
    case 32: smpwd = 0.5; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = (int) rint(bufsize * smpwd / 4.0);
  outbufsize = 2 * nsamples * sizeof(float);
  outbuf = (float *) malloc(outbufsize);
  buffer = (char *) malloc(bufsize);
  rcp = (char *) malloc(2 * nsamples * sizeof(char));

  if (rcp == NULL || outbuf == NULL || buffer == NULL) 
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  /* setup time counter */
  time = 0;
  timeint = 1.0 / (fsamp * 1e6);

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
	  outbufsize = 2 * nsamples * sizeof(float);
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
	  if (chan == 2) 
	    unpack_pfs_4c2b_lcp(buffer, rcp, bufsize);
	  else 
	    unpack_pfs_4c2b_rcp(buffer, rcp, bufsize);
	  break;
	case 6:
	  if (chan == 2) 
	    unpack_pfs_4c4b_lcp(buffer, rcp, bufsize);
	  else 
	    unpack_pfs_4c4b_rcp(buffer, rcp, bufsize);
	  break;
     	case 8: 
	  memcpy (rcp, buffer, bufsize);
	  break;
     	case 16: 
	  unpack_pfs_signed16bits(buffer, outbuf, bufsize);
	  break;
	case 32:
	  memcpy (outbuf, buffer, bufsize);
	  break;
	default: 
	  fprintf(stderr,"mode not implemented yet\n"); 
	  exit(1);
	}
      if (mode != 16 && mode != 32)
	for (i = 0; i < 2 * nsamples; i++) 
	  outbuf[i] = (float) rcp[i];
      

      /* optionally apply phase rotation and increment time */
      if (foff != 0)
	{
	  apply_linear_phase(outbuf,foff,time,timeint,nsamples);
	  time += timeint * nsamples;
	}

      /* optionally compute magnitude */
      if (mdetect && !pdetect)
	{
	  outbufsize = nsamples * sizeof(float);
	  for (i = 0, j = 0; i < nsamples; i++, j+=2)
	    outbuf[i] = sqrt(outbuf[j]*outbuf[j]+outbuf[j+1]*outbuf[j+1]);
	}

      /* optionally compute power */
      if (pdetect && !mdetect)
	{
	  outbufsize = nsamples * sizeof(float);
	  for (i = 0, j = 0; i < nsamples; i++, j+=2)
	    outbuf[i] = outbuf[j]*outbuf[j]+outbuf[j+1]*outbuf[j+1];
	}
      
      /* write data to output file */
      if (ascii)
	{
	  if (mode == 32) 
	    sprintf(format, "%% .3f %% .3f\n");
	  else
	    sprintf(format, "%% .0f %% .0f\n");

	  if (mdetect || pdetect)
	    for (i = 0; i < nsamples; i++)
	      fprintf(stdout,"% .3f\n",outbuf[i]);
	  else
	    for (i = 0, j = 0; i < nsamples; i++, j+=2)
	      fprintf(stdout,format,outbuf[j],outbuf[j+1]);
	}
      else
	{
	  if (outbufsize != write(fdoutput,outbuf,outbufsize))
	    fprintf(stderr,"Write error\n");  
	}
    }

  return 0;
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
void	processargs(argc,argv,infile,outfile,mode,chan,ascii,mdetect,pdetect,fsamp,foff)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *chan;
int     *ascii;
int     *mdetect;
int     *pdetect;
double   *fsamp;
double   *foff;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_unpack program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:c:o:adpf:x:"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_unpack -m mode [-c channel (1 or 2)] [-d (detect and output magnitude)] [-p (detect and output power)] [-o outfile (- for stdout)] [infile (- for stdin)] ";
  char *USAGE2="For phase rotation, also specify [-f sampling frequency (MHz)] [-x desired frequency offset (Hz)] ";
  char *USAGE3="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t16: signed 16bit\n\t32: 32bit floats\n";

  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *chan  = 1;
  *ascii = 0;
  *mdetect = 0;
  *pdetect = 0;
  *foff  = 0;
  *fsamp = 0;

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
	
      case 'c':
	sscanf(optarg,"%d",chan);
	arg_count += 2;
	break;
	
      case 'f':
	sscanf(optarg,"%lf",fsamp);
	arg_count += 2;
	break;
	
      case 'x':
	sscanf(optarg,"%lf",foff);
	arg_count += 2;
	break;
	
      case 'a':
	*ascii = 1;
	arg_count += 1;
	break;

      case 'd':
	*mdetect = 1;
	arg_count += 1;
	break;
	
      case 'p':
	*pdetect = 1;
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
  if (*mode == 0 ) goto errout;

  /* must specify valid sampling frequency */
  if (*foff != 0 && *fsamp == 0) goto errout;
  
  return;

  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE1);
          fprintf(stderr,"%s\n",USAGE2);
          fprintf(stderr,"%s",USAGE3);
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

