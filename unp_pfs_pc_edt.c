#include "unpack.h"

#define DBG1


/******************************************************************************/
/*	unpack_pfs_2c2b							      */
/******************************************************************************/
void unpack_pfs_2c2b (unsigned char *buf, char *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array contains 4*bufsize bytes 
    outbuf must have storage for at least 4*bufsize*sizeof(char)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  unsigned char value, val2n;
  char lookup[13] = {+3,+1,-1,-3,+1,0,0,0,-1,0,0,0,-3}; 
  int i;
  
  for (i = 0; i < bufsize; i += 4)
  {
      value = buf[i+1];
      val2n = value >> 4;
      *outbuf++ = lookup[val2n & 3];
      *outbuf++ = lookup[val2n & 0x0C];
      *outbuf++ = lookup[value & 3];
      *outbuf++ = lookup[value & 0x0C];

      value = buf[i+0];
      val2n = value >> 4;
      *outbuf++ = lookup[val2n & 3];
      *outbuf++ = lookup[val2n & 0x0C];
      *outbuf++ = lookup[value & 3];
      *outbuf++ = lookup[value & 0x0C];

      value = buf[i+3];
      val2n = value >> 4;
      *outbuf++ = lookup[val2n & 3];
      *outbuf++ = lookup[val2n & 0x0C];
      *outbuf++ = lookup[value & 3];
      *outbuf++ = lookup[value & 0x0C];

      value = buf[i+2];
      val2n = value >> 4;
      *outbuf++ = lookup[val2n & 3];
      *outbuf++ = lookup[val2n & 0x0C];
      *outbuf++ = lookup[value & 3];
      *outbuf++ = lookup[value & 0x0C];
  }

  return;
}

/******************************************************************************/
/*	unpack_pfs_2c4b   						      */
/******************************************************************************/
void unpack_pfs_2c4b (unsigned char *buf, char *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 4-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains 2*bufsize bytess
    output array must have been allocated for at least 2*bufsize*sizeof(char)
  */

  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  unsigned char value;
  char lookup[16] = {+15,+13,+11,+9,+7,+5,+3,+1,-1,-3,-5,-7,-9,-11,-13,-15}; 
  int i;
  
  for (i = 0; i < bufsize; i += 4)
    {
      value = buf[i+1];
      *outbuf++ = lookup[value & 15];
      *outbuf++ = lookup[value >> 4];

      value = buf[i+0];
      *outbuf++ = lookup[value & 15];
      *outbuf++ = lookup[value >> 4];

      value = buf[i+3];
      *outbuf++ = lookup[value & 15];
      *outbuf++ = lookup[value >> 4];

      value = buf[i+2];
      *outbuf++ = lookup[value & 15];
      *outbuf++ = lookup[value >> 4];
    }

  return;
}

/******************************************************************************/
/*	unpack_pfs_2c8b   						      */
/******************************************************************************/
void unpack_pfs_2c8b (unsigned char *buf, char *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 8-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains bufsize bytes
    output array must have been allocated for at least bufsize*sizeof(char)
  */

  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  int i;

# ifdef DBG
  unsigned char value1;
  unsigned char value2;
  unsigned char value3;
  unsigned char value4;
  
  value1 = buf[0];
  value2 = buf[1];
  value3 = buf[2];
  value4 = buf[3];
  printf ("**** TEST 3  %x %x %x %x \n", value1, value2, value3, value4);
#endif

  for (i = 0; i < bufsize; i+=4)
  {
      *outbuf++ = (unsigned char) buf[i+0] - 128;
      *outbuf++ = (unsigned char) buf[i+1] - 128;
      *outbuf++ = (unsigned char) buf[i+2] - 128;
      *outbuf++ = (unsigned char) buf[i+3] - 128;
  }

  return;
}

/******************************************************************************/
/*	unpack_pfs_2c8b_sb - 2's compliment data format                           */
/******************************************************************************/
void unpack_pfs_2c8b_sb (char *buf, char *outbuf, int bufsize)
{
  /*
    unpacks 2-channel, 8-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array outbuf contains bufsize bytes
    output array must have been allocated for at least bufsize*sizeof(char)
  */

  /* order is board 1 channel A, board 1 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I */
  /* connect BBC sin and cos to I and Q, respectively, for positive freq */

  int i;
  
#ifdef DBG
  char value1;
  char value2;
  char value3;
  char value4;
  
  value1 = buf[0];
  value2 = buf[1];
  value3 = buf[2];
  value4 = buf[3];
  printf ("**** TEST 2  %d %d %d %d \n", value1, value2, value3, value4);
#endif

  for (i = 0; i < bufsize; i+=4)
  {
      *outbuf++ = buf[i+0];
      *outbuf++ = buf[i+1];
      *outbuf++ = buf[i+2];
      *outbuf++ = buf[i+3];
	  /* 
      *outbuf++ = buf[i+2];
      *outbuf++ = buf[i+3];
      *outbuf++ = buf[i+0];
      *outbuf++ = buf[i+1];
	  */
  }

  return;
}


/******************************************************************************/
/*	unpack_pfs_4c4b_rcp			      */
/******************************************************************************/
void unpack_pfs_4c4b_rcp (unsigned char *buf, char *rcp, int bufsize)
{
  /*
    unpacks 4-channel, 4-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp and lcp each contain bufsize bytes
    rcp and lcp must each have storage for at least bufsize*sizeof(char)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /*          board 2 channel A, board 2 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I LCP-Q LCP-I */
  
  unsigned char value;
  char lookup[16] = {+15,+13,+11,+9,+7,+5,+3,+1,-1,-3,-5,-7,-9,-11,-13,-15}; 
  int i;
  
  for (i = 0; i < bufsize; i += 4) 
  {
      value = buf[i+0];
      *rcp++ = lookup[value & 15];
      value = value >> 4;
      *rcp++ = lookup[value & 15];

      value = buf[i+2];
      *rcp++ = lookup[value & 15];
      value = value >> 4;
      *rcp++ = lookup[value & 15];
  }

  return;
}

/******************************************************************************/
/*	unpack_pfs_4c4b_lcp						      */
/******************************************************************************/
void unpack_pfs_4c4b_lcp (unsigned char *buf, char *lcp, int bufsize)
{
  /*
    unpacks 4-channel, 4-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp and lcp each contain bufsize bytes
    rcp and lcp must each have storage for at least bufsize*sizeof(char)
  */
  
  /* order is board 1 channel A, board 1 channel B */
  /*          board 2 channel A, board 2 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I LCP-Q LCP-I */
  
  unsigned char value;
  char lookup[16] = {+15,+13,+11,+9,+7,+5,+3,+1,-1,-3,-5,-7,-9,-11,-13,-15}; 
  int i;
  
  for (i = 0; i < bufsize; i += 4) 
  {
      value = buf[i+1];
      *lcp++ = lookup[value & 15];
      value = value >> 4;
      *lcp++ = lookup[value & 15];

      value = buf[i+3];
      *lcp++ = lookup[value & 15];
      value = value >> 4;
      *lcp++ = lookup[value & 15];
  }

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
void unpack_pfs_4c2b_rcp (unsigned char *buf, char *rcp, int bufsize)
{
  /*
    unpacks 4-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output arrays rcp contains 2*bufsize bytes
    rcp must each have storage for at least 2*bufsize*sizeof(char)
  */
  
  /* order is RCP-I, RCP-Q and LCP-I, LCP-Q */
  /* i.e. board 1 channel A, board 1 channel B */
  /* and  board 2 channel A, board 2 channel B */
  /* tested */
  
  unsigned char value;
  char lookup[13] = {3,1,-1,-3,1,0,0,0,-1,0,0,0,-3};
  int i;
  
  for (i = 0; i < bufsize; i += 4) 
  {
      value = buf[i+1];
      *rcp++ = lookup[value & 3];
      *rcp++ = lookup[value & 0x0C];

      value = buf[i+0];
      *rcp++ = lookup[value & 3];
      *rcp++ = lookup[value & 0x0C];

      value = buf[i+3];
      *rcp++ = lookup[value & 3];
      *rcp++ = lookup[value & 0x0C];

      value = buf[i+2];
      *rcp++ = lookup[value & 3];
      *rcp++ = lookup[value & 0x0C];
  }

  return;
}


/******************************************************************************/
/*unpack_pfs_4c2b_lcp      */
/******************************************************************************/
void unpack_pfs_4c2b_lcp (unsigned char *buf, char *lcp, int bufsize)
{
  /*
    unpacks 4-channel, 2-bit data from the portable fast sampler
    input array buf is of size bufsize bytes
    output array lcp each contains 2*bufsize bytes
    lcp must each have storage for at least 2*bufsize*sizeof(char)
  */
  
  unsigned char value;
  char lookup[13] = {3,1,-1,-3,1,0,0,0,-1,0,0,0,-3};
  int i;
  
  for (i = 0; i < bufsize; i += 4) 
  {
      value = buf[i+1] >> 4;
      *lcp++ = lookup[value & 3];
      *lcp++ = lookup[value & 0x0C];

      value = buf[i+0] >> 4;
      *lcp++ = lookup[value & 3];
      *lcp++ = lookup[value & 0x0C];

      value = buf[i+3] >> 4;
      *lcp++ = lookup[value & 3];
      *lcp++ = lookup[value & 0x0C];

      value = buf[i+2] >> 4;
      *lcp++ = lookup[value & 3];
      *lcp++ = lookup[value & 0x0C];
  }

  return;
}


/*
  unpacks 4-channel, 8-bit data from the portable fast sampler
  input array buf is of size bufsize bytes
  output arrays rcp and lcp each contain bufsize bytes
  rcp and lcp must each have storage for at least bufsize*sizeof(char)
*/
  
/******************************************************************************/
/*	unpack_pfs_4c8b_rcp			      */
/******************************************************************************/
void unpack_pfs_4c8b_rcp (unsigned char *buf, char *rcp, int bufsize)
{
  /* order is board 1 channel A, board 1 channel B */
  /*          board 2 channel A, board 2 channel B */
  /* with wiring of PFS, this corresponds to RCP-Q RCP-I LCP-Q LCP-I */

  int i;
  for (i = 0; i < bufsize; i += 4) {
      *rcp++ = (unsigned char)buf[i] - 128;
      *rcp++ = (unsigned char)buf[i+1] - 128;
  }
  return;
}

/******************************************************************************/
/*	unpack_pfs_4c8b_lcp						      */
/******************************************************************************/
void unpack_pfs_4c8b_lcp (unsigned char *buf, char *lcp, int bufsize)
{
  int i;
  for (i = 0; i < bufsize; i += 4) {
      *lcp++ = (unsigned char)buf[i+2] - 128;
      *lcp++ = (unsigned char)buf[i+3] - 128;
  }
  return;
}

/******************************************************************************/
/*	unpack_pfs_4c8b_rcp_sb			      */
/******************************************************************************/
void unpack_pfs_4c8b_rcp_sb (unsigned char *buf, char *rcp, int bufsize)
{
  int i;
  for (i = 0; i < bufsize; i += 4) {
      *rcp++ = buf[i];
      *rcp++ = buf[i+1];
  }
  return;
}

/******************************************************************************/
/*	unpack_pfs_4c8b_lcp_sb						      */
/******************************************************************************/
void unpack_pfs_4c8b_lcp_sb (unsigned char *buf, char *lcp, int bufsize)
{
  int i;
  for (i = 0; i < bufsize; i += 4) {
      *lcp++ = buf[i+2];
      *lcp++ = buf[i+3];
  }
  return;
}


