/*******************************************************************************
*  program pfs_reset
*  $Id$
*  This program resets the PFS hardware
*
*  usage:
*  	pfs_reset  
*
*  input:
*       no input parameters
*
*  output:
*
*******************************************************************************/

/* 
   $Log$
   Revision 1.1  2002/05/12 23:16:24  cvs
   This program resets the PFS hardware.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "edtinc.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

int main(int argc, char **argv)
{
    EdtDev *edt_p ;
    uint off = 0x0;
    uint on = 0x1; 

    /* open device */
    if ((edt_p = edt_open("edt", 0)) == NULL)
    {
        perror("edt_open") ;
        exit(1) ;
    }
    else
      fprintf(stderr,"Device opened\n");
    
    /* turn off trigger */
    edt_reg_write(edt_p, PCD_FUNCT, off);
    fprintf(stderr,"Registers cleared\n");

    /* close device */
    edt_close(edt_p);
    fprintf(stderr,"Device closed\n");

    exit(0) ;
}
