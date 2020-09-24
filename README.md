About Solarv
==========

This is the user manual for *solarv* (solar radial velocity), an ephemeris code to compute the radial velocity of the Sun for ground-based observatories. Solarv is based on the SPICE toolkit provided by NASA’s [Navigation and Ancillary Information Facility (NAIF)](http://naif.jpl.nasa.gov/naif/).

Solarv was initially developed for the data reduction code of LARS (Lars is an Absolute Reference Spectrograph) at the German Vacuum Tower Telescope (VTT), Tenerife (see [Doerr 2015)](http://adsabs.harvard.edu/abs/2015PhDT.......200D).

Installation
============

Requirements
------------

-   Linux (64 bit). Other systems might or might not work, but are not supported

-   A csh shell (CSPICE' build script uses that)

-   A C compiler that supports the C99 standard

-   The cfitsio library

Compiling the source
--------------------

The solarv source distribution includes the CSPICE distribution which needs to be installed prior to solarv following the instructions below.

I recommend to consolidate all solarv related files in a single directory hierarchy, e.g. `~/local/solarv` for a setup in the users' home directory. This is what we are going to use in the following instructions.

**Steps to compile CSPICE.**

    $ mkdir ~/local/solarv
    $ cd local/solarv
    # fetch cspice sources
    $ wget http://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
    $ tar xzf cspice.tar.Z
    $ cd cspice
    $ ./makeall.csh

> **Note**
>
> CSPICE is traditionally used in the same directory where it is compiled.

**Steps to compile and install solarv.**

The solarv package doesn’t use an advanced build system, just a plain Makefile. Depending on your setup, you might have to apply some changes to `config.mk` or the `Makefile` itself. Move the solarv source package to `~/local` and follow these steps:

    $ cd ~/local
    $ tar xjf solarv-XYZ.tar.bz2
    $ cd solarv/src
    $ cp config.mk.default config.mk

Edit `config.mk` to look like this:

    ## CSPICE INSTALL DIR
    CSPICE_PATH = $(HOME)/local/solarv/cspice
    ## where required SPICE kernels are
    KERNEL_PATH = $(HOME)/local/solarv/kernels
    ## where to put the 'solarv' binary
    BINARY_PATH = $(HOME)/local/bin
    ## the pinpoint utility is needed to update the
    ## stations.bsp binary kernel
    PINPOINT = $(CSPICE_PATH)/exe/pinpoint

    CC = c99
    CFLAGS = -m64 -g
    LDFLAGS = -m64

**Steps to compile and install solarv (continued).**

    $ cd ~/local/solarv/src
    $ mkdir -p ~/local/solarv/kernels
    $ make
    $ make install
    $ make update-kernels
    $ export PATH=~/local/bin:$PATH

If these steps ran without any errors, solarv should be ready to use. Check it by running the following command. The output should look *exactly* the same.

    $  solarv -t -p 2016-01-01T12:00:00 hpc 0 0
    Date of observation..........  01 Jan 2016 12:00:00.00 UTC (JD 2457389.000000)
    Observer location............  VTT (28.30239 N, -16.51005 E, 2413 m)
    Terrestrial reference frame..  ITRF93
    Solar reference radius.......  695508000 m
    Apparent angular radius......  975.2690 arcsec
    Rotation model...............  fixed (no rotation, fixed to inertial frame)
    Sidereal rotation rate.......  0.0000 murad/s
    Position angle P.............  2.0880 deg
    Sub-observer Stonyhurst lat.. -3.0017 deg
    Sub-observer Stonyhurst lon..  0.0007 deg
    Sub-observer Carrington lon..  267.6156 deg
    Solar center distance........  147097208053 m / 0.983284103 AU
    Solar center radial velocity. -125.047 m/s
    Sun-Observer grav. redshift..  2.11231 ppm / 633.26 m/s
    Target HPC coordinates....... -0.0000, 0.0000 arcsec
      Stonyhurst lon, lat........  0.0007, -3.0017 deg
      Impact parameter...........  0 m / 0.00 arcsec
      Heliocentric angle, mu.....  0.0000 deg, 1.0000
      App. elevation, azimuth....  36.0186, 160.1700 deg
      Zenith distance, airmass...  53.9814 deg, 1.6959
      Distance...................  146401700053 m
      Radial velocity............ -125.047 m/s

Updating the SPICE kernels
--------------------------
TBD

Preparing a site kernel
-----------------------
TBD

Using solarv
============
A  quick usage reference and some examples can be obtained when calling solarv  with the `-h`switch.


Implementation details
======================
TBD

Solar differential rotation
---------------------------

Gravitational redshift
----------------------

Airmass
-------

License
=======

The solarv source code is licensed under the [MIT](https://opensource.org/licenses/MIT) open source license.

Copyright © 2011-2020 Hans-Peter Doerr

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
