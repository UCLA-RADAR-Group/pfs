/*******************************************************************************
*  program pfs_skipbytes
*  $Id$
*  This program reads nbyte-size records of data from the input file
*  optionally skipping over a user-specified number of bytes, and 
*  optionally writing out a user-specified window.
*
*  usage:
*  	pfs_skipbytes  -b nbytes [-n nreads] 
*                     [-s nskipbuffs,nskipbytes] 
*                     [-r startbyte,stopbyte]
*                     [-o outfile] infile
*
*  input:
*       the input parameters are typed in as command line arguments
*       -b number of bytes for each read   
*       -n number of reads to perform      (default 10)
*       -s number of buffers and bytes to skip (default 0)
*       -r desired start and stop bytes within a record.
*
*  output:
*	-o option identifies the output file (default is stdout)
*
*******************************************************************************/

/* 
   $Log$
   Revision 1.1  2002/04/27 06:26:37  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <asm/fcntl.h>

/* revision control variable */
static char const rcsid[] = 
"$Id$";

int	fdinput;		/* file descriptor for input file */
FILE   *fpinput;		/* pointer to input file */
FILE   *fpoutput;		/* pointer to output file */

char   *outfile;		/* output file name */
char   *infile;		       /* input file name */

char	command_line[200];	/* command line assembled by processargs */

void processargs();
void open_files();
void open_file();
void copy_cmd_line();

int main(int argc, char *argv[])
{
  int nbytes;			/* number of bytes per record */
  int nreads;			/* number of reads to perform */
  int nskips;
  int nskipbuffs;		
  int nskipbytes;
  int startbyte;
  int stopbyte;
  char *buffer;			/* buffer space */

  int i;
  int done;
  int buffcount;		/* buffer count */ 
  int open_flags;	/* flags required for open() call */

  /* get the command line arguments */
  processargs(argc,argv,&infile,&outfile,&nbytes,&nreads,&nskipbuffs,&nskipbytes,&startbyte,&stopbyte);

  /* check error conditions */
  if (nskipbytes > nbytes)
    exit(1);

  if (startbyte < 0 || startbyte > nbytes)
    exit(1);

  if (stopbyte < 0 || stopbyte > nbytes)
    exit(1);
  
  if (stopbyte < startbyte)
    exit(1);

  /* allocate data buffers */
  buffer = (char *) malloc(nbytes);
  if (buffer == NULL)
    {
      fprintf(stderr,"Unable to allocate %d-byte buffer.\n",nbytes);
      exit(1);
    }

  /* save the command line */
  copy_cmd_line(argc,argv,command_line);

  /* open input and output files */
  /* open_files(infile,outfile,&fpinput,&fpoutput); */

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

  /* Describe what we are doing to the user */
  fprintf(stderr,"Reading %d byte buffers\n",nbytes);
  fprintf(stderr,"Skipping first %d buffers\n",nskipbuffs);
  fprintf(stderr,"Plus an additional %d bytes\n",nskipbytes);
  fprintf(stderr,"Output bytes %d-%d from each record\n",startbyte,stopbyte);

  /* skip unwanted records */
  buffcount = 0;
  while (buffcount < nskipbuffs)
    {
      if (nbytes != read(fdinput, buffer, nbytes))
	fprintf(stderr,"Read error\n");
      buffcount++;
    }
  
  /* skip unwanted bytes */
  if (nskipbytes != 0)
    {
      if (nskipbytes != read(fdinput, buffer, nskipbytes))
	fprintf(stderr,"Read error\n");
    }

  /* read-write rest of data */
  while(buffcount < nreads + nskipbuffs)
    {
      if (nbytes != read(fdinput, buffer, nbytes))
	fprintf(stderr,"Read error\n");
      if (1 != fwrite(&buffer[startbyte],stopbyte-startbyte+1,1,fpoutput))
	fprintf(stderr,"Write error!\n"); 
      buffcount++;
    }

  return 0;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,infile,outfile,nbytes,nreads,nskipbuffs,nskipbytes,startbyte,stopbyte)
int	argc;
char	**argv;			 /* command line arguements */
char	**infile;		 /* array of input file names */
char	**outfile;		 /* output file name */
int     *nbytes;
int     *nreads;
int     *nskipbuffs;
int     *nskipbytes;
int     *startbyte;
int     *stopbyte;
{
  /* function to process a programs input command line.
     This is a template which has been customised for the skipbytes program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "b:n:s:r:o:"; 	 /* options to search for :=> argument*/
  char *USAGE="pfs_skipbytes -b nbytes [-n nreads] [-s nskipbuffs,nskipbytes] [-r startbyte,stopbyte] [-o outfile] infile";

  int  i;
  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */
  *outfile = "-";		 /* initialise to stdout */

  *nbytes  = -1;		 /* no default value */
  *nreads  = 10;
  *nskipbuffs = 0;
  *nskipbytes = 0;
  *startbyte = 1;
  *stopbyte = -1;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
    {
      case 'o':
 	       *outfile = optarg;	/* output file name */
               arg_count += 2;		/* two command line arguments */
	       break;

      case 'b':
 	       sscanf(optarg,"%d",nbytes);
               arg_count += 2;		/* two command line arguments */
	       break;
	    
      case 'n':
 	       sscanf(optarg,"%d",nreads);
               arg_count += 2;		/* two command line arguments */
	       break;

      case 's':
 	       sscanf(optarg,"%d,%d",nskipbuffs,nskipbytes);
               arg_count += 2;		/* two command line arguments */
	       break;

      case 'r':
 	       sscanf(optarg,"%d,%d",startbyte,stopbyte);
               arg_count += 2;		/* two command line arguments */
	       break;
	    
      case '?':			 /*if not in myoptions, getopt rets ? */
               goto errout;
               break;
    }
  }

  if (*nbytes < 0)
    goto errout;
  
  if (*stopbyte == -1)
    *stopbyte = *nbytes;

  if (arg_count < argc)	   /* non-optioned param is infile */
    *infile = argv[arg_count];

  return;
  
  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE);
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
/*	open files    							      */
/******************************************************************************/
void	open_files(infile,outfile,fpinput,fpoutput)
char	*infile;		/* array of input file names */
char	*outfile;		/* output file name */
FILE    **fpinput;		/* array of pointers to input files */
FILE    **fpoutput;		/* pointer to output file */
{

  int i;
  
  /* opens the input file */
  *fpinput=fopen(infile,"r");		/* open file input */
  if (*fpinput == NULL)
    {
      fprintf(stderr,"open_files: could not open input file\n");
      exit(1);
    }

  
  /* opens the output file, stdout is default */
  if (outfile[0] == '-')
    {
      *fpoutput=stdout;
    }
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

