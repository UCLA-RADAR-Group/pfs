#
#	Makefile to build pfs programs
#
#       The include file contains the definitions of environment variables
#       It is not included in the distribution so that you don't overwrite
#       your local settings when doing a software upgrade
#
include Makefile.inc
#
#
PROGRAMS=pfs_radar pfs_sample pfs_trigger pfs_reset pfs_levels pfs_hist pfs_stats pfs_unpack pfs_downsample pfs_fft pfs_dehop pfs_skipbytes
OBJECTS=pfs_radar.o pfs_sample.o pfs_trigger.o pfs_reset.o pfs_levels.o pfs_hist.o pfs_stats.o pfs_unpack.o pfs_downsample.o pfs_fft.o pfs_dehop.o pfs_skipbytes.o multifile.o libunpack.o
#
#
all: $(PROGRAMS)
#
#
# pfs_radar acquires data from the portable fast sampler
#
pfs_radar : pfs_radar.o multifile.o 
	$(CC) pfs_radar.o multifile.o  \
	-L/opt/EDTpcd -ledt \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_radar
#
# pfs_sample test samples some data from the portable fast sampler
#
pfs_sample : pfs_sample.o libunpack.o
	$(CC) pfs_sample.o libunpack.o \
	-L/opt/EDTpcd -ledt \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_sample
#
# pfs_trigger tests the 1 PPS and clock signals
#
pfs_trigger : pfs_trigger.o 
	$(CC) pfs_trigger.o \
	-L/opt/EDTpcd -ledt \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_trigger
#
#
# pfs_reset tests the 1 PPS and clock signals
#
pfs_reset : pfs_reset.o 
	$(CC) pfs_reset.o \
	-L/opt/EDTpcd -ledt \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_reset
#
# pfs_levels sets the PFS attenuator levels
#
pfs_levels : pfs_levels.o 
	$(CC) pfs_levels.o \
	-L/opt/EDTpcd -ledt \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_levels
#
# pfs_hist computes histograms of data from the portable fast sampler
#
pfs_hist : pfs_hist.o 
	$(CC) pfs_hist.o libunpack.o \
	$(LDFLAGS) \
	-o pfs_hist
#
# pfs_stats computes statistics of data from the portable fast sampler
#
pfs_stats : pfs_stats.o 
	$(CC) pfs_stats.o libunpack.o \
	$(LDFLAGS) \
	-o pfs_stats
#
# pfs_unpack unpacks data from the portable fast sampler
#
pfs_unpack : pfs_unpack.o 
	$(CC) pfs_unpack.o libunpack.o \
	$(LDFLAGS) \
	-o pfs_unpack
#
# pfs_downsample downsamples data from the portable fast sampler
#
pfs_downsample : pfs_downsample.o 
	$(CC) pfs_downsample.o libunpack.o \
	$(LDFLAGS) \
	-lpthread \
	-o pfs_downsample
#
# pfs_fft performs spectral analysis on data from the portable fast sampler
#
pfs_fft : pfs_fft.o 
	$(CC) pfs_fft.o libunpack.o \
	-lfftw \
	$(LDFLAGS) \
	-o pfs_fft
#
# pfs_dehop dehops fft spectra
#
pfs_dehop : pfs_dehop.o 
	$(CC) pfs_dehop.o \
	$(LDFLAGS) \
	-o pfs_dehop
#
# pfs_skipbytes skips over unwanted data
#
pfs_skipbytes : pfs_skipbytes.o 
	$(CC) pfs_skipbytes.o \
	$(LDFLAGS) \
	-o pfs_skipbytes
#
#
#
pfs_radar.o:	 pfs_radar.c ;	   $(CC) $(CFLAGS) -c pfs_radar.c -I/opt/EDTpcd
pfs_sample.o:	 pfs_sample.c ;	   $(CC) $(CFLAGS) -c pfs_sample.c -I/opt/EDTpcd
pfs_trigger.o:	 pfs_trigger.c ;   $(CC) $(CFLAGS) -c pfs_trigger.c -I/opt/EDTpcd
pfs_reset.o:	 pfs_reset.c ;   $(CC) $(CFLAGS) -c pfs_reset.c -I/opt/EDTpcd
pfs_levels.o:	 pfs_levels.c ;	   $(CC) $(CFLAGS) -c pfs_levels.c -I/opt/EDTpcd
pfs_hist.o:	 pfs_hist.c ;	   $(CC) $(CFLAGS) -c pfs_hist.c 
pfs_stats.o:	 pfs_stats.c ;	   $(CC) $(CFLAGS) -c pfs_stats.c 
pfs_unpack.o:    pfs_unpack.c ;    $(CC) $(CFLAGS) -c pfs_unpack.c 
pfs_downsample.o:pfs_downsample.c ;$(CC) $(CFLAGS) -c pfs_downsample.c 
pfs_fft.o:	 pfs_fft.c ;	   $(CC) $(CFLAGS) -c pfs_fft.c 
pfs_dehop.o:	 pfs_dehop.c ;     $(CC) $(CFLAGS) -c pfs_dehop.c 
pfs_skipbytes.o: pfs_skipbytes.c ; $(CC) $(CFLAGS) -c pfs_skipbytes.c 
multifile.o:	 multifile.c ;     $(CC) $(CFLAGS) -c multifile.c
libunpack.o:     unp_pfs_pc_edt.c; $(CC) $(CFLAGS) -c unp_pfs_pc_edt.c -o libunpack.o 
#
#
#
clean:
	/bin/rm -f a.out core $(OBJECTS) $(PROGRAMS) $(LIBOBJECTS) *~ \#*
#
install: $(PROGRAMS)
	@echo 'Installing programs : $(PROGRAMS)'
	@for f in $(PROGRAMS); \
	do install $$f $(DESTDIR)/$$f; \
	done;
#
distrib:
	tar cvf distrib.tar Makefile multifile.c multifile.h unpack.h unp_pfs_pc_edt.c pfs_radar.c pfs_sample.c pfs_trigger.c pfs_reset.c pfs_levels.c pfs_hist.c pfs_stats.c pfs_unpack.c pfs_downsample.c pfs_fft.c pfs_dehop.c pfs_skipbytes.c
