/*******************************************************************************
*  program pfs_levels
*  $Id$
*  This program sets the programmable attenuators available on some PFS systems
*
*  usage:
*  	pfs_levels -a attenuation (1-15 dB) 
*
*  input:
*       the input parameters are typed in as command line arguments
*
*  output:
*
*******************************************************************************/

/* 
   $Log$
   Revision 1.1  2002/04/27 06:20:42  margot
   Initial revision

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#include "edt.h"
#include "edtinc.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

void processargs();

int main(int argc, char **argv)
{
    EdtDev *edt_p ;
    uint disable = 0x0;
    uint enable  = 0xf; 
    int attenuation;

    /* get the command line arguments */
    processargs(argc,argv,&attenuation);

    /* check sanity of attenuation */
    if (attenuation < 1 || attenuation > 15)
      {
	fprintf(stderr,"Attenuator levels must be 1-15 dB\n");
	exit(1);
      }

    if ((edt_p = edt_open("edt", 0)) == NULL)
    {
        perror("edt_open") ;
        exit(1) ;
    }
    else
      fprintf(stderr,"Device opened\n");
    
    /* enable level setting and turn off data acquisition */
    edt_reg_write(edt_p, PCD_FUNCT, enable);

    /* give it a few clock cycles to toggle flip-flops */
    usleep(1000);

    /* set levels according to user-specified argument */
    edt_reg_write(edt_p, PCD_FUNCT, attenuation);

    /* give it a few clock cycles to toggle flip-flops */
    usleep(1000);

    /* disable level setting and re-enable data acquisition */
    edt_reg_write(edt_p, PCD_FUNCT, disable);

    /* give it a few clock cycles to toggle flip-flops */
    usleep(1000);

    edt_close(edt_p);

    exit(0) ;
}

/******************************************************************************/
/*	processargs							      */
/******************************************************************************/
void	processargs(argc,argv,attenuation)
int	argc;
char	**argv;			 /* command line arguements */
int     *attenuation;		/* attenuator levels */
{
  /* function to process a programs input command line.
     This is a template which has been customised for the pfs_levels program:
	- the outfile name is set from the -o option
	- the infile name is set from the 1st unoptioned argument
  */

  int getopt();		/* c lib function returns next opt*/ 
  extern char *optarg; 	/* if arg with option, this pts to it*/
  extern int optind;	/* after call, ind into argv for next*/
  extern int opterr;    /* if 0, getopt won't output err mesg*/

  char *myoptions = "a:"; 	 /* options to search for :=> argument*/
  char *USAGE="pfs_levels -a attenuation (1-15 dB)";

  int  c;			 /* option letter returned by getopt  */
  int  arg_count = 1;		 /* optioned argument count */

  /* default parameters */
  opterr = 0;			 /* turn off there message */

  *attenuation = 15;

  /* loop over all the options in list */
  while ((c = getopt(argc,argv,myoptions)) != -1)
  { 
    switch (c) 
    {
      case 'a':
 	       sscanf(optarg,"%d",attenuation);
               arg_count += 2;		/* two command line arguments */
	       break;
	    
      case '?':			 /*if not in myoptions, getopt rets ? */
               goto errout;
               break;
    }
  }

return;
  
  /* here if illegal option or argument */
  errout: fprintf(stderr,"%s\n",rcsid);
          fprintf(stderr,"Usage: %s\n",USAGE);
	  exit(1);
}
