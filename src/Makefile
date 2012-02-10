## you may define these
CSPICE_PATH = /home/doerr/coding/lfc/solarv/cspice
#CFITSIO_PATH = 
#KERNEL_PATH = /usr/local/share/data/solarv/kernels


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
else
	SPICE_MISSING = 1
endif

ifndef KERNEL_PATH
	KERNEL_PATH = "../data/kernels"
endif

CC = c99
CFLAGS = -Wall -O0 -g $(SPICE_I) $(FITSIO_I) \
	-DKERNEL_PATH=\"$(KERNEL_PATH)\"
LD = $(CC)
LDFLAGS = $(SPICE_L) $(FITSIO_L) -lm


all: solarv solarv.tm


.PHONY: clean check

clean:
	rm -f solarv solarv.o solarv.tm

solarv.o: solarv.h solarv.c

solarv: solarv.o solarv.h
	$(LD) $< -o solarv $(LDFLAGS)


solarv.tm: solarv.tm.in Makefile
	@sed 's#KERNEL_SYSPATH#$(KERNEL_PATH)#g' solarv.tm.in > solarv.tm
