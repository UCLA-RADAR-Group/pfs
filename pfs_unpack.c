/*******************************************************************************
*  program pfs_unpack
*  $Id$
*  This programs unpacks data from the portable fast sampler,
*  optionally applying a phase rotation to compensate for a frequency offset.
*  Default output is binary floating point quantities.
*
*  usage:
*  	pfs_unpack -m mode 
*                  [-a output ascii to stdout]
*                  [-d (detect and output magnitude)] 
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
   $Log$
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
#include <asm/fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include "unpack.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

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
  char *buffer;		/* buffer for packed data */
  float *rcp,*lcp;	/* buffer for unpacked data */
  double fsamp;		/* sampling frequency, MHz */
  double foff;		/* frequency offset, Hz */
  double timeint;	/* sampling interval */ 
  double time;		/* time */ 
  int mode;
  float smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int chan;		/* channel to process (1 or 2) for dual pol data */
  int ascii;		/* text output */
  int detect;		/* magnitude output */
  int i,j;

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&chan,&ascii,&detect,&fsamp,&foff);

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open file input */
#ifdef LARGEFILE
  if((fdinput = open(infile, O_RDONLY|O_LARGEFILE)) < 0 )
    perror("open input file");
  if((fdoutput = open(outfile, O_WRONLY|O_CREAT|O_LARGEFILE, 0644)) < 0 )
    perror("open output file");
#else
  if((fdinput = open(infile, O_RDONLY)) < 0 )
    perror("open input file");
  if((fdoutput = open(outfile, O_WRONLY|O_CREAT, 0644)) < 0 )
    perror("open output file");
#endif

  /* get file status and adjust buffer size if necessary */
  if (fstat (fdinput, &filestat) < 0)
    {
      perror("input file status");
      exit(1);
    }
  if (filestat.st_size < bufsize) bufsize = filestat.st_size;
  fprintf(stderr,"Unpacking file of size %d with %d byte buffers\n",
	  filestat.st_size,bufsize);
  if (filestat.st_size % bufsize != 0) 
    fprintf(stderr,"Warning: not an integer number of buffers\n");

  switch (mode)
    {
    case -1: smpwd = 8; break;  
    case  1: smpwd = 8; break;
    case  2: smpwd = 4; break;
    case  3: smpwd = 2; break; 
    case  5: smpwd = 4; break;
    case  6: smpwd = 2; break;
    case  8: smpwd = 2; break;
    case 32: smpwd = 0.5; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = (int) rint(bufsize * smpwd / 4.0);
  outbufsize = 2 * nsamples * sizeof(float);
  buffer = (char *) malloc(bufsize);
  rcp = (float *) malloc(outbufsize);
  lcp = (float *) malloc(outbufsize);
  if (lcp == NULL) 
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
      if (bufsize != read(fdinput, buffer, bufsize))
	{
	  perror("read");
	  fprintf(stderr,"Read error or EOF\n");
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
	  unpack_pfs_4c2b(buffer, rcp, lcp, bufsize);
	  if (chan == 2) memcpy(rcp, lcp, outbufsize);
	  break;
	case 6:
	  unpack_pfs_4c4b(buffer, rcp, lcp, bufsize);
	  if (chan == 2) memcpy(rcp, lcp, outbufsize);
	  break;
     	case 8: 
	  unpack_pfs_signedbytes(buffer, rcp, bufsize);
	  break;
     	case 32: 
	  memcpy(rcp,buffer,bufsize);
	  break;
	default: 
	  fprintf(stderr,"mode not implemented yet\n"); 
	  exit(1);
	}

      /* optionally apply phase rotation and increment time */
      if (foff != 0)
	{
	  apply_linear_phase(rcp,foff,time,timeint,nsamples);
	  time += timeint * nsamples;
	}

      /* optionally compute magnitude */
      if (detect)
	{
	  for (i = 0, j = 0; i < nsamples; i++, j+=2)
	    rcp[i] = sqrt(rcp[j]*rcp[j]+rcp[j+1]*rcp[j+1]);
	}
      
      /* write data to output file */
      if (ascii)
	{
	  if (detect)
	    for (i = 0; i < nsamples; i++)
	      fprintf(stdout,"% .3f\n",rcp[i]);
	  else
	    for (i = 0, j = 0; i < nsamples; i++, j+=2)
	      fprintf(stdout,"% .0f % .0f\n",rcp[j],rcp[j+1]);
	}
      else
	{
	  if (!detect)
	    {
	      if (outbufsize != write(fdoutput,rcp,outbufsize))
		fprintf(stderr,"Write error\n");  
	    }
	  else
	    {
	      if (outbufsize/2 != write(fdoutput,rcp,outbufsize/2))
		fprintf(stderr,"Write error\n");  
	    }
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
void	processargs(argc,argv,infile,outfile,mode,chan,ascii,detect,fsamp,foff)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *chan;
int     *ascii;
int     *detect;
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

  char *myoptions = "m:c:o:adf:x:"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_unpack -m mode [-c channel (1 or 2)] [-d (detect and output magnitude)] [-o outfile] [infile] ";
  char *USAGE2="For phase rotation, also specify [-f sampling frequency (MHz)] [-x desired frequency offset (Hz)] ";
  char *USAGE3="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n\t 8: signed bytes\n\t32: 32bit floats\n";

  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *chan  = 1;
  *ascii = 0;
  *detect = 0;
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
	*detect = 1;
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

