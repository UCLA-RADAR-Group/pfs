/*******************************************************************************
*  program pfs_radar
*  $Id$
*  This program performs data acquitisition using the Portable Fast Sampler
*  and an EDT PCI CD-x0 interface card.
*
*  usage:
*       pfs_radar 
*       [-start yyyy,mm,dd,hh,mn,sc] 
*       [-secs sec] [-step sec] [-cycles c] 
*       [-files f] [-rings r] 
*       [-bytes b] [-log l]
*       -dir d [-dir d]... 
*
*  input:
*       the input parameters are typed in as command line arguments
*
*  output:
*       the output data is streamed to disk
*
*  Original program written by Jeff Hagen.
*  Written in the spirit of the Stewart Anderson wptape code,
*  this program reads a 16-bit wide EDT PCI CD-x0 card
*  and writes the data onto disk.  
*
*  compile with 
*  gcc -O2 -o pfs_radar -I/opt/EDTpcd pfs_radar.c /opt/EDTpcd/libedt.c -lthread
*
*******************************************************************************/

/* 
   $Log$
*/

#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include <termios.h>
#include <sys/time.h>
#include <sys/procfs.h>
#include <sys/priocntl.h>
#include <sys/rtpriocntl.h>
#include <sys/tspriocntl.h>

#include "edtinc.h"

/* revision control variable */
static char const rcsid[] = 
"$Id$";

void schedule_rt( int );

struct DISKWRITE { /* one of these for each diskbuffer allocated */
  char name[80];
  int fd;
  int len;
  char *out;
  thread_t proc;
};

struct RADAR { /* structure that holds the buffers and configuration */
  EdtDev *edt;
  unsigned int mode;
  int ameg;
  int twensec;
  int step;
  int cycles;
  int ringbufs;
  int diskbufs;
  struct DISKWRITE *dw;
  unsigned short **rings;
  time_t start;
  time_t startmone;
  time_t next;
  char timestr[80];
  char log[80];
  FILE *logfd;
  int pack;
} radar;

struct DIRLIST {
  struct DIRLIST *next;
  char name[80];
};

#define AMEG (1024*1024) /* default size of edt ring buffer */
#define RINGBUFS  8      /* default number of one meg edt ring buffers */
#define DISKBUFS  1      /* default number of pending disk writes */
#define TWENSECS  3600   /* default number of seconds to take */
#define AFEWSECS  3      /* interval bw key pressed and toggle EDT bit */

int ctlc_flag = 0;

int main(int argc, char *argv[])
{
  int i,cycle;
  unsigned char *p;
  struct DISKWRITE *w;
  int dcount, twensec;
  void *disk_write();
  struct RADAR *r;
  struct DIRLIST *dirhead, *dir;
  unsigned int off=0x00;
  unsigned int on=0x01;
  struct tm go;
  int time_set = 0;

#ifdef TIMER
  struct timeval   now;
  struct timeval   then;
  struct itimerval mytimer;
#endif

  r = &radar;
  bzero( r, sizeof(struct RADAR ));
  r->diskbufs = DISKBUFS;
  r->ringbufs = RINGBUFS;
  r->ameg = AMEG;
  r->pack = 1;
  r->twensec = 20 * TWENSECS;
  r->step = 0;
  r->cycles = 1;

  dirhead = NULL;
  for( i=1; i<argc; i++ ) {
    p = argv[i];
    if( strncasecmp( p, "-secs", strlen(p) ) == 0 ) {
      p = argv[++i];
      if( (r->twensec = atoi(p)*20)<=0 ) {
        fprintf(stderr, "bad value for -secs\n");
        pusage();
      }
    } else if( strncasecmp( p, "-step", strlen(p) ) == 0 ) {
      p = argv[++i];
      if((r->step = atoi(p))<=0 ) {
        fprintf(stderr, "bad value for -step\n");
        pusage();
      }
    } else if( strncasecmp( p, "-cycles", strlen(p) ) == 0 ) {
      p = argv[++i];
      if((r->cycles = atoi(p))<=0 ) {
        fprintf(stderr, "bad value for -cycles\n");
        pusage();
      }
    } else if( strncasecmp( p, "-files", strlen(p) ) == 0 ) {
      p = argv[++i];
      if((r->diskbufs = atoi(p))<=0 ) {
        fprintf(stderr, "bad value for -files\n");
        pusage();
      }
    } else if( strncasecmp( p, "-mode", strlen(p) ) == 0 ) {
      p = argv[++i];
      if((r->mode = atoi(p))<0 ) {
        fprintf(stderr, "bad value for -mode\n");
        pusage();
      }
      switch (r->mode)
	{
	case 1: break;
	case 2: break;
	case 3: break;
	case 5: break;
	case 6: break;
	default: fprintf(stderr,"invalid mode\n"); pusage();
	}
    } else if( strncasecmp( p, "-start", strlen(p) ) == 0 ) {
      p = argv[++i];
      if((sscanf(p,"%d,%d,%d,%d,%d,%d",
 		 &go.tm_year,&go.tm_mon,&go.tm_mday,
 		 &go.tm_hour,&go.tm_min,&go.tm_sec) != 6)) {
	fprintf(stderr, "bad value for -start\n");
	pusage();
      } else time_set = 1;
    } else if( strncasecmp( p, "-dir", strlen(p) ) == 0 ) {
      p = argv[++i];
      dir = (struct DIRLIST *)malloc( sizeof(struct DIRLIST));
      strcpy( dir->name, p );
      dir->next = dirhead;
      dirhead = dir;
      if( access( dir->name, W_OK|X_OK )) {
        fprintf(stderr, "unable to access directory %s\n", dir->name );
        pusage();
      }
    } else if( strncasecmp( p, "-rings", strlen(p) ) == 0 ) {
      p = argv[++i];
      if(( r->ringbufs = atoi(p))<=0 ) {
        fprintf(stderr, "bad value for -rings\n");
        pusage();
      }
    } else if( strncasecmp( p, "-bytes", strlen(p) ) == 0 ) {
      p = argv[++i];
      if(( r->ameg = atoi(p))<=0 ) {
        fprintf(stderr, "bad value for -bytes\n");
        pusage();
      }
    } else if( strncasecmp( p, "-log", strlen(p) ) == 0 ) {
      p = argv[++i];
      strcpy( r->log, p);
    } else if( strncasecmp( p, "-nopack", strlen(p) ) == 0 ) {
      r->pack = 0;
    }
  }

  if (r->cycles > 1 && r->step < r->twensec / 20 + AFEWSECS)
    {
      fprintf(stderr,"Step size must be bigger than duration of A/D\n");
      exit(1);
    }

  if( r->ringbufs <=0 )
    r->ringbufs = RINGBUFS;

  if( r->ameg <=0 )
    r->ameg = AMEG;

  if( r->diskbufs <=0 )
    r->diskbufs = DISKBUFS;

  if( !dirhead ) {
    fprintf(stderr, "At least one -dir switch is required\n");
    pusage();
  }
  
  /* set scheduling priority */
  schedule_rt(2); 
  /* start control-c interrupt handler */
  set_kb(1);

  printf("Starting the Portable Fast Sampler\n");

  /* open edt device */
  if ((r->edt = edt_open("pcd", 0)) == NULL)
    {
      perror("edt_open") ;
      exit(1) ;
    }

  /* block trigger immediately */
  edt_reg_write( r->edt, PCD_FUNCT, 0x00 | (r->mode << 1)); /* block trigger */

  /* allocate input buffers */
  allocate_ringbufs(r);

  /* allocate output buffers */
  allocate_writebufs(r);

  /* open log file */
  open_log(r);

  /* compute starting time */
  if (time_set)
    {
      printf("A/D will start on second tick at %04d %02d %02d %02d %02d %02d\n",
	     go.tm_year,go.tm_mon,go.tm_mday,
	     go.tm_hour,go.tm_min,go.tm_sec);
      /* take care of unusual definitions for month, day */
      go.tm_year = go.tm_year - 1900;
      go.tm_mon = go.tm_mon - 1;
      r->next = mktime(&go);
    }
  /* otherwise let user hit the go button */
  else
    {
      printf("Hit a key when ready to take data\n");
      printf("A/D will start on second tick following key_press event + %d seconds\n",AFEWSECS-1);
      getchar();

      /* get next start time */
      r->next = time(NULL) + AFEWSECS;
    }
  
  /* loop over desired number of cycles */
  for (cycle = 1; cycle <= r->cycles; cycle++)
    {
      /* define start times for this cycle and the next */
      r->start = r->next;
      r->startmone = r->start - 1;
      r->next = r->start + r->step;

      /* obtain time string and open files */
      get_tms(r->start,r->timestr); 
      open_files(r, dirhead);

      fprintf(stdout  , "\nCycle %d will start at %s\n" , cycle, r->timestr );
      fprintf(r->logfd, "\nCycle %d starting at %s\n" , cycle, r->timestr );
      fflush(r->logfd);

      edt_flush_fifo( r->edt );
      edt_start_buffers( r->edt, 0 );
      
      /* wait till .5 sec before expected pulse */
      wait_till_start(r->startmone); 
      /* arm trigger */
      edt_reg_write( r->edt, PCD_FUNCT, 0x01 | (r->mode << 1));  

#ifdef TIMER
      /* set current time */
      now.tv_sec = 100000;
      now.tv_usec = 0;
      /* set next timer expiration a long time from now */
      mytimer.it_value.tv_sec  = now.tv_sec;
      mytimer.it_value.tv_usec = now.tv_usec;
      /* reload timer with incremental values  */
      mytimer.it_interval.tv_sec = now.tv_sec;
      mytimer.it_interval.tv_usec = now.tv_usec;
      /* set timer going */
      if (!setitimer(ITIMER_REAL, &mytimer, NULL))
	perror("setitimer");
#endif

      /* main loop */
      twensec = r->twensec;
      dcount = 0;
      for( i=0; i<twensec; i++ ) {
	if( ctlc_flag )
	  break;
	if(!(p = edt_wait_for_buffers( r->edt, 1)))
	  printf("error \n");
	else {
	  if( i%100 == 0 )
	    printf("i = %d count = %d\n", i, edt_done_count(r->edt));
	  if( edt_ring_buffer_overrun(r->edt)) {
	    printf("overrun %d\n",i );
	    fprintf(r->logfd, "overrun %d\n" , i );
	  } else {
#ifdef TIMER
	    /* read current timer value */
	    getitimer(ITIMER_REAL, &mytimer);
	    now = mytimer.it_value;
	    fprintf(r->logfd,"%6d %06ld %06ld %20ld\n",
		    i,then.tv_usec,now.tv_usec,
		    1000000*(then.tv_sec-now.tv_sec)+then.tv_usec-now.tv_usec); 
	    then.tv_sec = now.tv_sec;
	    then.tv_usec = now.tv_usec;
#endif TIMER
	    w = &r->dw[dcount];
	    /* copy data from EDT to disk write output buffer */
	    memcpy(w->out,p,w->len);
	    if( w->proc )
	      thr_join( w->proc, NULL, NULL );
	    if( thr_create( NULL, 0, disk_write, w, THR_BOUND, &w->proc ))
	      perror("do_write");
	    if( ++dcount >=r->diskbufs )
	      dcount = 0;
	  }
	}
      }
      
      if( ctlc_flag )
	fprintf(r->logfd, "Stopped by user, read %d buffers\n\n\n", i );
      else
	fprintf(r->logfd, "Finished, read %d buffers\n\n\n", i );

#ifndef DEBUF      
      for( i=0; i< r->diskbufs; i++ ) {
	w = &r->dw[i];
	thr_join( w->proc, NULL, NULL );
      }
#endif      
      edt_stop_buffers( r->edt);
      /* clear trigger */
      edt_reg_write( r->edt, PCD_FUNCT, 0x00 | (r->mode << 1)); 
      edt_reset_ring_buffers(r->edt, 0);
      close_files(r);
    }
      
  edt_close(r->edt);
  fclose(r->logfd);
  set_kb(0);
}

/* disk_write thread */

void *disk_write( w )
struct DISKWRITE *w;
{
  int writ;

  if((writ = write( w->fd, w->out, w->len))!= w->len ) 
    printf(" disk write error: could only write %d bytes\n", writ );

  return(0);
}

/* courtesy of Stuart Anderson, swiped right out of proc_ut.c */

void
schedule_rt( int rt_priority )
{
    pcparms_t pcparms;
    rtparms_t *rtparmsp;
    pcinfo_t pcinfo;
    id_t rtID;
    short maxrtpri;
    rtinfo_t *rtinfop;


    /*
     * boost process to the RT scheduling class with rt_priority.
     */

    (void) strcpy (pcinfo.pc_clname, "RT");
    if (priocntl (0L, 0L, PC_GETCID, (caddr_t) &pcinfo) == -1L) {
	perror ("PC_GETCID failed for realtime class");
	exit (1);
    }
    
    rtID = pcinfo.pc_cid;
    rtinfop = (struct rtinfo *)pcinfo.pc_clinfo;
    maxrtpri = rtinfop->rt_maxpri;
    
    pcparms.pc_cid = rtID;
    rtparmsp = (struct rtparms *)pcparms.pc_clparms;
    rtparmsp->rt_pri = rt_priority;
    rtparmsp->rt_tqnsecs = RT_TQDEF;
    
    if (priocntl (P_PID, getpid(), PC_SETPARMS, (caddr_t) &pcparms) == -1) {
	perror ("PC_SETPARMS failed");
    }
}


/* 
  allocate output buffers using valloc and calling mlock
  r is the config structure
*/

allocate_writebufs(r)
struct RADAR *r;
{
  int pagesize, i, aout;
  struct DISKWRITE *w;

  pagesize = sysconf(_SC_PAGESIZE);
  r->ameg = r->ameg + r->ameg%pagesize;
  aout = (r->ameg) + (r->ameg)%pagesize;

  if(!(r->dw = (struct DISKWRITE *)valloc( 
     sizeof(struct DISKWRITE)*r->diskbufs))) {
      fprintf(stderr, "bad malloc allocating buffer\n");
      exit(1);
  }

  if( mlock( r->dw, sizeof(struct DISKWRITE)*r->diskbufs))
    perror("failed to mlock");

  for( i=0; i< r->diskbufs; i++ ) {
    w = &r->dw[i];
    if( !(w->out = ( char * ) valloc( aout ))) {
      fprintf(stderr, "bad valloc allocating buffer\n");
      exit(1);
    }
 
    if( mlock( w->out, aout ))
      perror("failed to mlock");

  }
}

/*
  open output files 
  r is the config structure
  head is the linked list of directory names to be used in opening files
*/

open_files(r, head)
struct RADAR *r;
struct DIRLIST *head;
{
  int i;
  struct DISKWRITE *w;
  struct DIRLIST *dir;
  char *tms;

  dir = head;

  for( i=0; i< r->diskbufs; i++ ) {
    w = &r->dw[i];

    if( !(dir = dir->next ))
      dir = head;
      
    sprintf(w->name, "%s/data%s.%02d", dir->name, r->timestr, i );
    if((w->fd = open(w->name, O_WRONLY|O_CREAT|O_LARGEFILE, 0664 )) < 0 )
      perror("write file open");
    w->len = r->ameg;
    w->proc = NULL;
  }
}


/*
  close output files 
  r is the config structure
*/

close_files(r)
struct RADAR *r;
{
  int i;
  struct DISKWRITE *w;

  for( i=0; i< r->diskbufs; i++ ) 
    {
      w = &r->dw[i];
      close(w->fd);
    }
}


/* build time string */

get_tms(time_t time, char *string)
{
  struct tm *ans;
 
  ans = gmtime(&time);
  /* correct for peculiar tm_mon : months since January - [0, 11] */
  /* by adding one to get usual [1,12] interval */
  ans->tm_mon += 1;
  sprintf( string, "%04d%02d%02d%02d%02d%02d",
    ans->tm_year+1900, ans->tm_mon, ans->tm_mday, 
    ans->tm_hour, 
    ans->tm_min, 
    ans->tm_sec );
}
  
  
/* 
  allocate input buffers using valloc and calling mlock
  r is the config structure
*/

allocate_ringbufs(r)
struct RADAR *r;
{
  int i;

  if(!(r->rings = (unsigned short **)valloc( 
     sizeof(unsigned short *)*r->ringbufs))) {
      fprintf(stderr, "bad malloc allocating ring buffer\n");
      exit(1);
  }

  if( mlock( r->rings, sizeof(unsigned short *)*r->ringbufs))
    perror("failed to mlock");

  for( i=0; i<r->ringbufs; i++ ) {
    if( !(r->rings[i] = ( unsigned short *)valloc(r->ameg)))
      printf("bad valloc\n");
    else
      if( mlock( r->rings[i], r->ameg ))
        perror("mlock data");
  }
  if( edt_configure_ring_buffers( r->edt, 
     r->ameg, r->ringbufs, EDT_READ, (void *)r->rings))
    edt_perror("configure edt card failure:");
}

open_log(r)
struct RADAR *r;
{
  if( r->log[0] == 0 )
    strcpy(r->log, "radar.log");

  if( (r->logfd = fopen( r->log, "a+") )== NULL ) {
    fprintf( stderr, "Failed to open log file %s\n", r->log );
    exit(1);
  }
  fprintf(r->logfd, "%s\n",rcsid);
  fprintf(r->logfd, "Input buffer size %d bytes\n", r->ameg );
  fprintf(r->logfd, "Input buffers, %d\n", r->ringbufs );
  fprintf(r->logfd, "Maximum time %d seconds\n", r->twensec/20 );
  fprintf(r->logfd, "Write buffers, %d\n", r->diskbufs );
  fprintf(r->logfd, "Data taking mode, %d\n", r->mode );
  fflush(r->logfd);
}


/* set keyboard:
   if you call it with a 1: set_kb(1);
   that is the same as -> stty cbreak,stty -echo.
   if you call it with a 0: set_kb(0);
   that is the same as -> stty -cbreak,stty echo.
 
   it operates on the standard input (0) which is
   assumed to be associated with a terminal.
*/

void do_ctlc()
{
  ctlc_flag = 1;
  set_kb(0);
}

set_kb(on)
int on;
{
  struct termios t;

  if(ioctl(0,TCGETS,&t))
    perror("ioctl:");

  if(on) {
    t.c_lflag &= ~ICANON;
    t.c_lflag &= ~ECHO;
    t.c_cc[VMIN] = 1;
    t.c_cc[VTIME] = 0;
    sigset( SIGINT, do_ctlc );

  } else {
    t.c_lflag |= ICANON;
    t.c_lflag |= ECHO;
    sigset( SIGINT, SIG_DFL );
  }

  if(ioctl(0,TCSETS,&t))
    perror("ioctl:");
}

/*
   Timing notes.

   Its important that we know the precise time the data was taken.
   A one second timing pulse is provided to trigger the hardware to start
   dumping data. From solaris we need to arm this trigger.  The start time 
   is encoded in the filenames. But it takes some time to open files, valloc 
   and mlock all the buffers.  

   Heres how we handle it.

   - We presume that the one second tick of the system clock is within +-0.5
     seconds of the timing pulse. This can be set and tested elsewhere. 
   - A start time in seconds is computed by adding an offset of AFEWSECS
     to the current one second time from time(NULL); 
   - Then the files are opened, memory is allocated and log file prepared.
   - Then the process polls the time() function till one second before the 
     start time by polling usleep(100), Suprisingly, this takes very little cpu.
   - Then we usleep(500000), a half a second.
   - Then enable data taking.
   - Then arm the trigger.

*/

/* wait till 0.5 sec before expected pulse */

wait_till_start( ttt )
time_t ttt;
{
  while( time(NULL) < ttt)
    usleep(100);

  usleep(500000);
}

pusage()
{
  fprintf( stderr, "Usage: radar -dir d [-dir d]... [-mode mode] [-secs sec] [-step sec] [-cycles c] [-files ff] [-rings r] [-nopack]\n");
  fprintf( stderr, "                                             (defaults)\n");
  /*  fprintf( stderr, "  -mode mode  0: 1c1b, 2: 1c2b, 4: 1c4b, 6: 1c8b, 8: 2c1b, 10: 2c2b, 12: 2c4b\n");*/
  fprintf( stderr, "  -mode mode\n\t 0: 2c1b (N/A)\n\t 1: 2c2b\n\t 2: 2c4b\n\t 3: 2c8b\n\t 4: 4c1b (N/A)\n\t 5: 4c2b\n\t 6: 4c4b\n\t 7: 4c8b (N/A)\n");
  fprintf( stderr, "  -secs sec   number of seconds of data to take (3600)\n");
  fprintf( stderr, "  -step sec   timestep between A/D cycles (0)\n");
  fprintf( stderr, "  -cycles c   number of repeat cycles (1)\n");
  fprintf( stderr, "  -files f    total number of files to open (1)\n");
  fprintf( stderr, "  -dir d      directory to use, multiples allowed\n");
  fprintf( stderr, "  -rings r    number of input buffers to use (8)\n");
  fprintf( stderr, "  -bytes b    size of input ring buffer (1048576 bytes)\n");
  fprintf( stderr, "  -log l      log file name \n");
  exit(1);
}

