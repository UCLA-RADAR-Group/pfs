/*******************************************************************************
*  program pfs_trigger
*  $Id$
*  This program tests the pfs 1 pps and clock signals
*
*  usage:
*  	pfs_trigger  
*
*  input:
*       the input parameters are typed in as command line arguments
*
*  output:
*
*******************************************************************************/

/* 
   $Log$
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

int main(int argc, char **argv)
{
    EdtDev *edt_p ;
    unsigned char *buffer;
    unsigned char value;
    int bit[8];
    uint off = 0x0;
    uint on = 0x1; 

    if ((edt_p = edt_open("pcd", 0)) == NULL)
    {
        perror("edt_open") ;
        exit(1) ;
    }
    else
      fprintf(stderr,"Device opened\n");
    
    if (edt_configure_ring_buffers(edt_p, 1024*1024, 4, EDT_READ, NULL) == -1)
    {
        perror("edt_configure_ring_buffers") ;
        exit(1) ;
    }
    else
      fprintf(stderr,"Buffers configured\n");


    /* turn off trigger */
    edt_reg_write(edt_p, PCD_FUNCT, off);

    /* prepare buffers for free running mode */
    if (edt_start_buffers(edt_p, 0) == -1)
      {
        perror("edt_start_buffers") ;
        exit(1) ;
      }
    else
      fprintf(stderr,"Buffers started\n");

    printf("hit a key to start sampling on next second tick\n");
    getc(stdin);
    edt_reg_write(edt_p, PCD_FUNCT, on);

    for (;;)
      {
	fprintf(stderr,"Buffer %d\t",edt_done_count(edt_p));
	buffer = edt_wait_for_buffers(edt_p, 1) ;
	value  = buffer[0];
	bit[0] = (value & 0x1)? 1 : 0;
	bit[1] = (value & 0x2)? 1 : 0;
	bit[2] = (value & 0x4)? 1 : 0;
	bit[3] = (value & 0x8)? 1 : 0;
	/* fprintf(stdout,"%d\t",value); */
	fprintf(stdout,"chan A %2d %2d chan B %2d %2d\n",
		bit[0],bit[1],bit[2],bit[3]);
      }

    exit(0) ;
}
