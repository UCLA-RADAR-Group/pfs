#include "unpack.h"

/******************************************************************************/
/*	unpack_pfs_2c2b							      */
/******************************************************************************/
void unpack_pfs_2c2b(char *buf, float *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array contains 4*bufsize floats
    outbuf must have storage for at least 4*bufsize*sizeof(float)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  char value;
  float lookup[4] = {+3,+1,-1,-3}; 
  int i,k;
  
  k = 0;
  for (i = 0; i < bufsize; i += 4)
    {
      value = buf[i+1];
      outbuf[k+2] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+3] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+0] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+1] = lookup[value & 3];
      k += 4;

      value = buf[i+0];
      outbuf[k+2] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+3] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+0] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+1] = lookup[value & 3];
      k += 4;

      value = buf[i+3];
      outbuf[k+2] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+3] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+0] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+1] = lookup[value & 3];
      k += 4;

      value = buf[i+2];
      outbuf[k+2] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+3] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+0] = lookup[value & 3];
      value = value >> 2;
      outbuf[k+1] = lookup[value & 3];
      k += 4;
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_2c4b   						      */
/******************************************************************************/
void unpack_pfs_2c4b(char *buf, float *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 4-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains 2*bufsize floats
    output array must have been allocated for at least 2*bufsize*sizeof(float)
  */

  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  char value;
  float lookup[16] = {+15,+13,+11,+9,+7,+5,+3,+1,-1,-3,-5,-7,-9,-11,-13,-15}; 
  int i,k;
  
  k = 0;
  for (i = 0; i < bufsize; i += 4)
    {
      value = buf[i+1];
      outbuf[k++] = lookup[value & 15];
      value = value >> 4;
      outbuf[k++] = lookup[value & 15];

      value = buf[i+0];
      outbuf[k++] = lookup[value & 15];
      value = value >> 4;
      outbuf[k++] = lookup[value & 15];

      value = buf[i+3];
      outbuf[k++] = lookup[value & 15];
      value = value >> 4;
      outbuf[k++] = lookup[value & 15];

      value = buf[i+2];
      outbuf[k++] = lookup[value & 15];
      value = value >> 4;
      outbuf[k++] = lookup[value & 15];
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_2c8b   						      */
/******************************************************************************/
void unpack_pfs_2c8b(char *buf, float *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 8-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains bufsize floats
    output array must have been allocated for at least bufsize*sizeof(float)
  */

  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  int i,k;
  
  k = 0;
  for (i = 0; i < bufsize; i+=4)
    {
      outbuf[k++] = (float) (unsigned char) buf[i+0] - 128;
      outbuf[k++] = (float) (unsigned char) buf[i+1] - 128;
      outbuf[k++] = (float) (unsigned char) buf[i+2] - 128;
      outbuf[k++] = (float) (unsigned char) buf[i+3] - 128;
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_4c2b							      */
/******************************************************************************/
void unpack_pfs_4c2b(char *buf, float *rcp, float *lcp, int bufsize)
{
  /*
    unpacks 4-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp and lcp each contain 2*bufsize floats
    rcp and lcp must each have storage for at least 2*bufsize*sizeof(float)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /*          board 2 channel A, board 2 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I LCP-Q LCP-I */
  
  char value;
  float lookup[4] = {+3,+1,-1,-3};  
  int i,r,l;
  
  r = 0;
  l = 0;
  for (i = 0; i < bufsize; i += 4) 
    {
      value = buf[i+1];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+0];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+3];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+2];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_4c4b							      */
/******************************************************************************/
void unpack_pfs_4c4b(char *buf, float *rcp, float *lcp, int bufsize)
{
  /*
    unpacks 4-channel, 4-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp and lcp each contain bufsize floats
    rcp and lcp must each have storage for at least bufsize*sizeof(float)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /*          board 2 channel A, board 2 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I LCP-Q LCP-I */
  
  char value;
  float lookup[16] = {+15,+13,+11,+9,+7,+5,+3,+1,-1,-3,-5,-7,-9,-11,-13,-15}; 
  int i,r,l;
  
  r = 0;
  l = 0;
  for (i = 0; i < bufsize; i += 4) 
    {
      value = buf[i+1];
      lcp[l++] = lookup[value & 15];
      value = value >> 4;
      lcp[l++] = lookup[value & 15];

      value = buf[i+0];
      rcp[r++] = lookup[value & 15];
      value = value >> 4;
      rcp[r++] = lookup[value & 15];

      value = buf[i+3];
      lcp[l++] = lookup[value & 15];
      value = value >> 4;
      lcp[l++] = lookup[value & 15];

      value = buf[i+2];
      rcp[r++] = lookup[value & 15];
      value = value >> 4;
      rcp[r++] = lookup[value & 15];
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_signedbytes 						      */
/******************************************************************************/
void unpack_pfs_signedbytes(char *buf, float *outbuf, int bufsize)
{
  /*
    unpacks signed bytes, eg stream obtained by digital filtering
    of data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains bufsize floats
    output array must have been allocated for at least bufsize*sizeof(float)
  */

  int i;
  
  for (i = 0; i < bufsize; i++)
    outbuf[i] = (float) buf[i];
  
  return;
}

/******************************************************************************/
/*	unpack_pfs_signed16bits						      */
/******************************************************************************/
void unpack_pfs_signed16bits(char *buf, float *outbuf, int bufsize)
{
  /*
    unpacks signed 16 bit quantities, eg stream obtained from VME das
    input array buf is of size bufsize bytes
    output array outbuf contains 2*bufsize floats
    output array must have been allocated for at least 2*bufsize*sizeof(float)
  */

  int i,j;
  signed short int x;

  for (i = 0, j = 0; i < bufsize; i+=sizeof(signed short int), j++)
    {
      memcpy(&x,&buf[i],sizeof(signed short int));
      outbuf[j] = (float) x;
    }

  return;
}

/******************************************************************************/
/*unpack_pfs_4c2b_rcp      */
/******************************************************************************/
void unpack_pfs_4c2b_rcp(char *buf, float *rcp, int bufsize)
{
  /*
    unpacks 4-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp contains 2*bufsize floats
    rcp must each have storage for at least 2*bufsize*sizeof(float)
  */
  
  /* order is RCP-I, RCP-Q and LCP-I, LCP-Q */
  /* i.e. board 1 channel A, board 1 channel B */
  /* and  board 2 channel A, board 2 channel B */
  /* tested */
  
  char value;
  float lookup[4] = {+3,+1,-1,-3};  
  int i,r,l;
  
  r = 0;
  l = 0;
  for (i = 0; i < bufsize; i += 4) 
    {
      value = buf[i+1];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];

      value = buf[i+0];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];

      value = buf[i+3];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];

      value = buf[i+2];
      rcp[r++] = lookup[value & 3];
      value = value >> 2;
      rcp[r++] = lookup[value & 3];

    }

  return;
}


/******************************************************************************/
/*unpack_pfs_4c2b_lcp      */
/******************************************************************************/
void unpack_pfs_4c2b_lcp(char *buf, float *lcp, int bufsize)
{
  /*
    unpacks 4-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array lcp each contains 2*bufsize floats
    lcp must each have storage for at least 2*bufsize*sizeof(float)
  */
  
  char value;
  float lookup[4] = {+3,+1,-1,-3};  
  int i,l;
  
  l = 0;
  for (i = 0; i < bufsize; i += 4) 
    {
      value = buf[i+1];
      value = value >> 4;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+0];
      value = value >> 4;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+3];
      value = value >> 4;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];

      value = buf[i+2];
      value = value >> 4;
      lcp[l++] = lookup[value & 3];
      value = value >> 2;
      lcp[l++] = lookup[value & 3];
    }

  return;
}



/*****************************************************************************/
