
struct MULTIFILE {
  long long max;
  long long cur_off;
  int cur_file;
  int max_file;
  int fd[ 500 ]; /* biggest thing we could support is 2GIGS * 500 (TERABYTE) */
  char name[256];
};

int multi_config_maxfilesize( long long );
struct MULTIFILE *multi_open( char *, unsigned int, int, int ); 
int multi_close( struct MULTIFILE *);
int multi_write( struct MULTIFILE *, char *, int );

