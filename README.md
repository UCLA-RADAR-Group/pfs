Portable Fast Sampler Software
====================================

These programs control the operation of the Portable Fast Sampler
(PFS) systems.  They also provide tools for initial data analysis 
(unpacking, digital filtering, spectral analysis, de-hopping, etc).
The code includes more than 8,000 lines of C code.  Most of this code 
has been incorporated in the JPL radar backend.  


If you use these programs, please cite the paper describing the instrument ([PDF](http://www.ursi.org/proceedings/procGA02/papers/p1949.pdf)).

To download and install, read the file [INSTALL](/INSTALL).

What you may need in addition to the software in this repository:
- [CVS](https://en.wikipedia.org/wiki/Concurrent_Versions_System), make, gcc.
- the EDT driver from www.edt.com, installed in its default /opt/EDTpcd location


If you want to take data with correct timestamps, please make sure
that the [Network Time Protocol (NTP)](https://en.wikipedia.org/wiki/Network_Time_Protocol) is installed and properly
configured.

Please send comments and bug reports to:

Jean-Luc Margot<br>
University of California, Los Angeles<br>
jlm@astro.ucla.edu
