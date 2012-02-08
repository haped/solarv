## you may define these
CSPICE_PATH = /home/doerr/coding/lfc/solarv/cspice
#CFITSIO_PATH = 


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

CC = c99
CFLAGS = -Wall -O2 $(SPICE_I) $(FITSIO_I)
LD = $(CC)
LDFLAGS = $(SPICE_L) $(FITSIO_L) -lm -static


all: solarv


.PHONY: clean check

clean:
	rm -f solarv solarv.o

solarv.o: solarv.h solarv.c

solarv: solarv.o solarv.h
	$(LD) $< -o solarv $(LDFLAGS)
