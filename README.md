Portable Fast Sampler Software
====================================

These programs control the operation of the Portable Fast Sampler (PFS) systems that were in use at Arecibo (2000–2020), Goldstone (2001–2014), and Green Bank (2001–2017).  They also provide tools for initial data analysis (unpacking, digital filtering, spectral analysis, de-hopping, etc).  The code includes more than 8,000 lines of C code.  A substantial fraction of this code has been incorporated in the software that is used to operate and process data from the NASA JPL dual channel agilent receiver (DCAR) data-taking systems installed at Goldstone and Green Bank.

If you use this code, please cite the corresponding paper: J. L. Margot, [A Data-Taking System for Planetary Radar Applications](https://doi.org/10.1142/S225117172150001X), J. Astron. Instrum., <b>10</b>, 2021, ([PDF](https://ucla.box.com/s/t7mc5cflclor6si3kaqfqz46pxdlmjqc), [BibTeX](/pfs.bib)).

# Requirements

- Compilation and execution require bash, make, and gcc
  
- FFTW (version >= 3.3.8)

    - Install with your package manager:

       - RedHat/CentOS:   
       (include file at /usr/include/, library at /usr/lib64/)    
       ```sh
       sudo yum install fftw fftw-devel  
       ```
    
       - Debian/Ubuntu:  
       (include file at /usr/include/, library at /usr/lib/)  
       ```sh
       sudo apt-get install libfftw3-3 libfftw3-dev  
       ```
    
       - MacOS:  
       (include file at /opt/local/include/, library at /opt/local/lib/)  
       ```sh
       sudo port install fftw-3-single  
       ```
    
    - Install from source:  
       ```sh
       wget http://www.fftw.org/fftw-3.3.8.tar.gz  
       tar xvf fftw-3.3.8.tar.gz  
       ./configure --enable-float  
       make; make install
       ```

# Additional requirements for new PFS hardware installations

- EDT driver from [www.edt.com](www.edt.com), installed in its default /opt/EDTpcd location

- [Network Time Protocol (NTP)](https://en.wikipedia.org/wiki/Network_Time_Protocol) 

# Installation procedure

- Download repository:  

  With https:  
  ```sh
  git clone https://github.com/UCLA-RADAR-Group/pfs.git  
  ```
  
  With ssh:  
  ```sh
  git clone git@github.com:UCLA-RADAR-Group/pfs.git  
  ```

- Specify the target location for scripts and executables:  
  Edit the value of `GLOBDIR` in the top-level [Makefile](src/Makefile) if necessary (default is `$(HOME)/bin`)
  
- Compile:  
  ```sh
  cd src; make  
  ```

- Install:
  ```sh
  make install; make clean; cd ..
  ```

- Run tests:
  ```sh
  cd tests; make
  ```

- Successfully tested on CentOS 7.3.1611/7.5.1804, MacOS 10.11.6.

# Basic usage

Download a [users guide](https://seti.ucla.edu/jlm/research/pfs/pfs_usage.pdf).  

More information is available at the [PFS project webpage](https://seti.ucla.edu/jlm/research/pfs/).  

# Credits
Most of the software was written by Jean-Luc Margot, with software contributions from Jeff Hagen and Joseph Jao.  The design and construction of the first data-taking systems was funded by the National Astronomy and Ionosphere Center.

# Further reading  
Additional information about planetary radar astronomy can be found at [http://radarastronomy.org](http://radarastronomy.org).

# Contact  
Please send comments and bug reports to:  

Jean-Luc Margot  
University of California, Los Angeles  
jlm@astro.ucla.edu
