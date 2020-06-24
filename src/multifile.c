
#include <stdio.h>
#include <fcntl.h>
#include <sys/vfs.h>

#include "multifile.h"

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

/* a layer to supprt multiple file writes */

#define TERABYTE 500

#define TENGIGS 10000000000LL
#define TWOGIGS  (0x7fffffff)

static struct MULITCONST {
  long long maxfilesize;
  int largefile;
} multi = { TWOGIGS, O_LARGEFILE };

multi_config_maxfilesize( max )
long long max;
{
  if( max > TWOGIGS )
    multi.largefile = O_LARGEFILE;
  else
    multi.largefile = 0;
    
  if( max > TWOGIGS && multi.largefile == 0 )
    multi.maxfilesize = TWOGIGS;
  else
    multi.maxfilesize = max;
}

struct MULTIFILE *multi_open( prefix, openpar, mask, nfiles )
char *prefix;
unsigned int openpar;
int mask; 
int nfiles;
{
  struct MULTIFILE *m;
  char filename[256];
  struct statfs fs;
  long long avail, inc;
  int fd, n;

  m = (struct MULTIFILE *)malloc(sizeof(struct MULTIFILE));
  bzero( m, sizeof(struct MULTIFILE));

  for( n=0; n<nfiles; n++ ) {
    strncpy( m->name, prefix, sizeof(m->name) );
    sprintf( filename, "%s.%03d", prefix, n );
   
    if( (fd = open( filename, openpar|multi.largefile, mask ))<0 ) {
      perror("open"); 
      free(m);
      return(NULL);
    }
    m->fd[ m->max_file++ ] = fd;
  }
  m->max = multi.maxfilesize;
  return(m);
}

multi_close(m)
struct MULTIFILE *m;
{
  int i;
  char name[256];

  for( i=m->cur_file; i<m->max_file; i++ )
     close( m->fd[ i ] );

  for( i=m->cur_file+1; i<m->max_file; i++ ) {
      sprintf( name, "%s.%03d", m->name, i );
      if( unlink( name ) < 0 )
        perror("removing file");
      printf("removing unused file: %s\n", name );
  }
  
  free(m);
  return(0);
}

/*
  make the interface just like write
  cant have len > TWOGIGS
*/

int multi_write(m, buf, len )
struct MULTIFILE *m;
char *buf;
int len;
{
  long long off, wrt;
  int wlen, retlenlast, retlen;

   off = m->cur_off + len;
   retlen = 0;
   retlenlast = 0;

   /* printf("called multi_write %d\n", len );  */
   if( m->cur_file >= m->max_file )
     return(-1);

   if( off > m->max ) {
     wlen = m->max - m->cur_off; 


     /* printf( "big writing %d bytes at off %ld file %d, fd %d\n", wlen, off, 
	m->cur_file, m->fd[m->cur_file]); */


     if((retlenlast = write( m->fd[m->cur_file], buf, wlen )) != wlen ) {
       perror( "multi_write");
       return(-1);
     }
     m->cur_off = 0;
     wlen = len - wlen;
     /* printf(" closing file %d, fd %d\n", m->cur_file, m->fd[m->cur_file]); */
     /* close( m->fd[m->cur_file]); */
     if( ++m->cur_file >= m->max_file )
        return(retlenlast);
   } else
    wlen = len;

   /* printf( "writing %d bytes at file %d, off %ld, fd %d\n", wlen, off,
      m->cur_file, m->fd[m->cur_file]); */

  if((retlen = write( m->fd[m->cur_file], &buf[retlenlast], wlen )) != wlen) {
     perror( "multi_write");
     return(-1);
  }

  m->cur_off += wlen;
  return(retlen+retlenlast);
}

