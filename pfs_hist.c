/*******************************************************************************
*  program pfs_hist
*  $Id$
*  This programs unpacks some data from the portable fast sampler
*  and prints and histogram of count values for all channels
*
*  usage:
*  	pfs_hist -m mode [-a (parse all data)] [-e (parse data at eof)] 
*               [-o outfile] [infile]
*
*  input:
*       the input parameters are typed in as command line arguments
*	the -m option specifies the data acquisition mode
*	the -e option specifies to parse data at the end of the file
*	the -a option specifies to parse all the data recorded
*                     (default is to parse the first megabyte)
*
*  output:
*	the -o option identifies the output file, stdout is default
*
*******************************************************************************/

/* 
   $Log$
   Revision 3.1  2007/06/19 17:00:21  jao
   Implemented '-a' option, mode 7 (4c8b) historgram, and 2's complemnt for mode 3 & 7

   Revision 3.0  2003/02/25 22:10:43  cvs
   Adapted to use Joseph Jao's byte unpacking.

   Revision 1.7  2002/05/01 05:29:33  cvs
   Fixed bug in handling of mode 3 (2c8b) histograms.

   Revision 1.6  2002/04/27 20:22:26  margot
   Removed obsolete routines specific to Golevka Sampling Box.

   Revision 1.5  2001/07/10 00:37:54  margot
   Adjusted input buffer size according to file size.

   Revision 1.4  2001/07/10 00:24:07  margot
   Added unpacking of signed bytes.

   Revision 1.3  2001/07/10 00:18:39  margot
   Added -e option for parsing at end of file.

   Revision 1.2  2000/10/30 05:21:41  margot
   Added variable nsamples.

   Revision 1.1  2000/10/30 04:46:22  margot
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

FILE   *fpoutput;		/* pointer to output file */
int	fdinput;		/* file descriptor for input file */

char   *outfile;		/* output file name */
char   *infile;		        /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_file();
void copy_cmd_line();


int main(int argc, char *argv[])
{
  struct stat filestat;	/* input file status structure */
  int mode;		/* data acquisition mode */
  int twoscmp = 0;	/* 2's complement (0 = FALSE)  */
  int bufsize = 1048576;/* size of read buffer, default 1 MB */
  char *buffer;		/* buffer for packed data */
  char *rcp,*lcp;	/* buffer for unpacked data */
  int smpwd;		/* # of single pol complex samples in a 4 byte word */
  int nsamples;		/* # of complex samples in each buffer */
  int levels;		/* # of levels for given quantization mode */
  int open_flags;	/* flags required for open() call */
  int parse_all;
  int parse_end;
  long long r_ihist[512], r_qhist[512];
  long long l_ihist[512], l_qhist[512];
  int i;

  /* initialization */
  for (i = 0; i < 512; i++) {
    r_ihist[i] = 0; r_qhist[i] = 0; l_ihist[i] = 0; l_qhist[i] = 0;
  }

  /* get the command line arguments and open the files */
  processargs(argc,argv,&infile,&outfile,&mode,&twoscmp,&parse_all,&parse_end);

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
  
  /* adjust buffer size if needed */
  if (filestat.st_size < bufsize)
    bufsize = filestat.st_size;

  switch (mode)
    {
    case -1: smpwd = 8; levels =   4; break;  
    case  1: smpwd = 8; levels =   4; break;
    case  2: smpwd = 4; levels =  16; break;
    case  3: smpwd = 2; levels = 256; break; 
    case  5: smpwd = 4; levels =   4; break;
    case  6: smpwd = 2; levels =  16; break;
    case  7: smpwd = 1; levels = 256; break;
    case  8: smpwd = 2; levels = 256; break; 
    default: fprintf(stderr,"Invalid mode\n"); exit(1);
    }

  /* allocate storage */
  nsamples = bufsize * smpwd / 4;
  buffer = (char *) malloc(bufsize);
  rcp = (char *) malloc(2 * nsamples * sizeof(char));
  lcp = (char *) malloc(2 * nsamples * sizeof(char));
  if (lcp == NULL) 
    {
      fprintf(stderr,"Malloc error\n"); 
      exit(1);
    }

  if (parse_end)
    lseek(fdinput, -bufsize, SEEK_END);

  /* read first buffer */
  do {
    if (bufsize != read(fdinput, buffer, bufsize)) {
      fprintf(stderr,"Read error\n");
      break;
    }  

    switch (mode) { 
      case 1:
        /* unpack & compute histogram */
        unpack_pfs_2c2b(buffer, rcp, bufsize);
      
        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels - 1] += 1; 
          r_qhist[(int)rcp[i+1] + levels - 1] += 1; 
        }

        break;
      case 2: 
        /* unpack & compute histogram */
        unpack_pfs_2c4b(buffer, rcp, bufsize);

        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels - 1] += 1; 
          r_qhist[(int)rcp[i+1] + levels - 1] += 1; 
        }

        break;
      case 3: 
        /* unpack & compute histogram */
        if (!twoscmp) {
          unpack_pfs_2c8b(buffer, rcp, bufsize);
        } else {
          unpack_pfs_2c8b_sb(buffer, rcp, bufsize);
        }

        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels/2] += 1; 
          r_qhist[(int)rcp[i+1] + levels/2] += 1; 
        }

        break;
      case 5:
        /* unpack & compute histogram */
        unpack_pfs_4c2b_rcp(buffer, rcp, bufsize);
        unpack_pfs_4c2b_lcp(buffer, lcp, bufsize);

        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels - 1] += 1; 
          r_qhist[(int)rcp[i+1] + levels - 1] += 1; 

          l_ihist[(int)lcp[i]   + levels - 1] += 1; 
          l_qhist[(int)lcp[i+1] + levels - 1] += 1; 
        }

        break;
      case 6:
        /* unpack & compute histogram */
        unpack_pfs_4c4b_rcp(buffer, rcp, bufsize);
        unpack_pfs_4c4b_lcp(buffer, lcp, bufsize);

        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels - 1] += 1; 
          r_qhist[(int)rcp[i+1] + levels - 1] += 1; 

          l_ihist[(int)lcp[i]   + levels - 1] += 1; 
          l_qhist[(int)lcp[i+1] + levels - 1] += 1; 
        }

        break;
      case 7: 
        /* unpack & compute histogram */
        if (!twoscmp) {
          unpack_pfs_4c8b_rcp(buffer, rcp, bufsize);
          unpack_pfs_4c8b_lcp(buffer, lcp, bufsize);
        } else {
          unpack_pfs_4c8b_rcp_sb(buffer, rcp, bufsize);
          unpack_pfs_4c8b_lcp_sb(buffer, lcp, bufsize);
        }

        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)rcp[i]   + levels/2] += 1; 
          r_qhist[(int)rcp[i+1] + levels/2] += 1; 

          l_ihist[(int)lcp[i]   + levels/2] += 1; 
          l_qhist[(int)lcp[i+1] + levels/2] += 1; 
        }

        break;
      case 8: 
        for (i = 0; i < 2*nsamples; i += 2) {
          r_ihist[(int)buffer[i]   + levels/2] += 1; 
          r_qhist[(int)buffer[i+1] + levels/2] += 1; 
        }
        
        break;
      default: fprintf(stderr,"mode not implemented yet\n"); exit(1);
    }
  } while (parse_all);


  /* print results */
  // mode 3 or 7 changes 256 -> 128 level for easy of display
  if (mode == 3 || mode == 8 || mode == 7) {  
      fprintf(fpoutput,"RCP hist\n");
 
      for (i = 0; i < levels; i ++) {
        fprintf(fpoutput,"%10d %15qd \t",i - levels/2,r_ihist[i]);
        fprintf(fpoutput,"%10d %15qd \n",i - levels/2,r_qhist[i]);
      }

      if (mode == 7) {     
        fprintf(fpoutput,"LCP hist\n");
 
        for (i = 0; i < levels; i ++) {
          fprintf(fpoutput,"%10d %15qd \t",i - levels/2,l_ihist[i]);
          fprintf(fpoutput,"%10d %15qd \n",i - levels/2,l_qhist[i]);
        }
      }
  } else {
      fprintf(fpoutput,"RCP hist\n");
 
      for (i = 0; i < 2 * levels; i += 2) {
        fprintf(fpoutput,"%10d %15qd \t",i - levels + 1,r_ihist[i]);
        fprintf(fpoutput,"%10d %15qd \n",i - levels + 1,r_qhist[i]);
      }

      if (mode > 4) {
        fprintf(fpoutput,"LCP hist\n");
 
        for (i = 0; i < 2 * levels; i += 2) {
          fprintf(fpoutput,"%10d %15qd \t",i - levels + 1,l_ihist[i]);
          fprintf(fpoutput,"%10d %15qd \n",i - levels + 1,l_qhist[i]);
        }
      }
  }

  return 0;
}


/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,mode,twoscmp,parse_all,parse_end)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* input file name */
char	**outfile;		 /* output file name */
int     *mode;
int     *twoscmp;
int     *parse_all;
int     *parse_end;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_hist program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "m:o:ae2"; 	 /* options to search for :=> argument*/
  char *USAGE1="pfs_hist -m mode [-2 (2's complement)] [-e (parse data at eof)] [-a (parse all data)] [-o outfile] [infile] ";
  char *USAGE2="Valid modes are\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n";
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *infile  = "-";		 /* initialise to stdin, stdout */
  *outfile = "-";

  *mode  = 0;                /* default value */
  *twoscmp  = 0;             /* default value */
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
	    
      case '2':
 	       *twoscmp = 1;
               arg_count ++;		/* one command line arguments */
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
  
  if (*twoscmp && *mode != 3 && *mode != 7) {
    fprintf(stderr,"2's complement is supported on mode 3 & 7 only\n"); 
    goto errout;
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

