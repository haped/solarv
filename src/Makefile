## you may define these
CSPICE_PATH = /usr/local/cspice
#CFITSIO_PATH = /usr/local/pkg/cfitsio-3250
KERNEL_PATH = /usr/local/share/data/solarv/kernels
BINARY_PATH = /usr/local/bin
PINPOINT = pinpoint

## defining CFITSIO_PATH is optional, could also be pre-installed
## on the system
ifdef CFITSIO_PATH
	FITSIO_I = -I$(CFITSIO_PATH)/include
	FITSIO_L = -L$(CFITSIO_PATH)/lib -lcfitsio
else
	FITSIO_L = -lcfitsio
endif

## need to set SPICE_PATH in config.mk
ifdef CSPICE_PATH
	SPICE_I = -I$(CSPICE_PATH)/include
	SPICE_L = $(CSPICE_PATH)/lib/cspice.a
	PINPOINT = $(CSPICE_PATH)/exe/pinpoint
else
	SPICE_MISSING = 1
endif

ifndef KERNEL_PATH
	KERNEL_PATH = ../data/kernels
endif

CC = c99
CFLAGS = -Wall -pedantic -O2 -m64 $(SPICE_I) $(FITSIO_I) \
	-DKERNEL_PATH=\"$(KERNEL_PATH)\"
LD = $(CC) 
LDFLAGS = -m64 $(SPICE_L) $(FITSIO_L) -lm 
# solaris
#LDFLAGS = -m64 $(SPICE_L) $(FITSIO_L) -lm -lsocket -lnsl


all: solarv solarv.tm

install: solarv solarv.tm
	install -D solarv.tm $(KERNEL_PATH)/solarv.tm
	install -D solarv $(BINARY_PATH)/solarv

.PHONY: clean check

clean:
	rm -f solarv solarv.o solarv.tm stations.bsp

solarv.o: solarv.h solarv.c Makefile

solarv: solarv.o solarv.h
	$(LD) $< -o solarv $(LDFLAGS)


solarv.tm: solarv.tm.in Makefile
	@sed 's#KERNEL_SYSPATH#$(KERNEL_PATH)#g' solarv.tm.in > solarv.tm


stations.bsp: stations.defs
	$(PINPOINT) -def stations.defs \
	-pck $(KERNEL_PATH)/pck00010.tpc \
	-spk stations.bsp
