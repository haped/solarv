# About Solarv
SolaRV (solar radial velocity) is a precision ephemeris code to compute the relative velocity between a ground-based observatory on the earth and the Sun. Such information is, for example, required for absolute velocity calibration of spectroscopic observational data of the Sun. The code is based on the SPICE toolkit provided by NASA’s [Navigation and Ancillary Information Facility (NAIF)](http://naif.jpl.nasa.gov/naif/).

Solarv was initially developed for the data reduction code of LARS (Lars is an Absolute Reference Spectrograph) at the German Vacuum Tower Telescope (VTT), Tenerife (see [Doerr 2015)](http://adsabs.harvard.edu/abs/2015PhDT.......200D).

# Installation

## Requirements
- Linux (64 bit). Other systems might or might not work, but are not supported
- A csh shell (CSPICE' build script uses that)
- A C compiler that supports the C99 standard
- The cfitsio library

On Debian or Ubuntu-based Linux systems it should be sufficient to install the `libcfitsio-dev` and `tcsh` packages:

	sudo apt-get install libcfitsio-dev tcsh

## Preparations
Solarv depends on the CSPICE toolkit which needs to be installed prior to solarv. I recommend to consolidate all solarv related files and data in a single directory hierarchy, e.g. `~/local/solarv` for an installation in the users' home directory.

### Getting the source code
	# prepare directory tree
	mkdir ~/local
	cd ~/local
	# clone the repository from github
	git clone git@github.com:haped/solarv.git
	cd solarv
	# get the latest CSPICE source code distribution
	wget http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
	tar xzf cspice.tar.Z 

### Compiling CSPICE
	cd ~/local/solarv/cspice
	./makeall.csh

> **Note:**
> CSPICE is traditionally used in the same directory where it is compiled, there is no `make install`

### Compiling solarv
The solarv package doesn’t use an advanced build system, just a plain Makefile. Depending on your setup, you might have to apply some changes to `config.mk` or the `Makefile` itself.

	cd ~/local/solarv/src
	cp config.mk.default config.mk

Now, make sure you have the following definitions in your config.mk:

```
PKG_CONFIG = /usr/bin/pkg-config
CSPICE_PATH = $(HOME)/local/solarv/cspice
KERNEL_PATH = $(HOME)/local/solarv/kernels
BINARY_PATH = $(HOME)/local/bin
CC = c99
CFLAGS = -m64 -g -Wall -pedantic -O2
LDFLAGS = -m64
```

### Compiling solarv (continued)
	cd ~/local/solarv/src
	make
	make install
	make update-kernels
	export PATH=~/local/bin:$PATH

  > **Note:** You probably want to add ~/local/bin to your PATH environemnt varaible permanently by appending `export PATH=~/local/bin:$PATH` to your `~/.bash_profile`.

If these steps ran without any errors, solarv should be ready to use. Check it by running the following command. The output should look *exactly* the same.

	solarv -t -p 2016-01-01T12:00:00 hpc 0 0
	Date of observation..........  01 Jan 2016 12:00:00.00 UTC (JD 2457389.000000)
	Observer location............  VTT (28.30230 N, -16.51014 E, 2444 m)
	Terrestrial reference frame..  ITRF93
	Solar reference radius.......  695508000 m
	Apparent angular radius......  975.2690 arcsec
	Rotation model...............  fixed (no rotation w.r.t. interial frame)
	Sidereal rotation rate.......  0.0000 murad/s
	Position angle P.............  2.0880 deg
	Sub-observer Stonyhurst lat.. -3.0017 deg
	Sub-observer Stonyhurst lon..  0.0007 deg
	Sub-observer Carrington lon..  267.6156 deg
	Solar center distance........  147097208030 m / 0.983284103 AU
	Solar center radial velocity. -125.048 m/s
	Sun-Observer grav. redshift..  2.11231 ppm / 633.26 m/s
	Target HPC coordinates.......  0.0000, -0.0000 arcsec
	  Stonyhurst lon, lat........  0.0007, -3.0017 deg
	  Impact parameter...........  0 m / 0.00 arcsec
	  Heliocentric angle, mu.....  0.0000 deg, 1.0000
	  App. elevation, azimuth....  36.0186, 160.1698 deg
	  Zenith distance, airmass...  53.9814 deg, 1.6959
	  Distance...................  146401700030 m
	  Radial velocity............ -125.048 m/s

## Updating the SPICE kernels
TBD

# Preparing a site kernel
TBD

# Using solarv
A brief usage reference and some examples can be obtained when calling solarv  with the `-h` switch.

# License
The solarv source code is licensed under the [MIT](https://opensource.org/licenses/MIT) open source license.

Copyright © 2011-2020 Hans-Peter Doerr

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
