include config.mk

FITSIO_I = `$(PKG_CONFIG) cfitsio --cflags`
FITSIO_L = `$(PKG_CONFIG) cfitsio --libs`

## need to set SPICE_PATH in config.mk
ifdef CSPICE_PATH
	SPICE_I = -I$(CSPICE_PATH)/include
	SPICE_L = $(CSPICE_PATH)/lib/cspice.a
else
	SPICE_MISSING = 1
endif

ifndef KERNEL_PATH
	KERNEL_PATH = ../data/kernels
endif

ifndef CC
	CC = c99
endif
CFLAGS += $(SPICE_I) $(FITSIO_I) \
	-DKERNEL_PATH=\"$(KERNEL_PATH)\"
LD = $(CC) 
LDFLAGS += $(SPICE_L) $(FITSIO_L) -lm 


all: solarv kernels/solarv.tm

install: solarv kernels/solarv.tm
	install -D -m 0644 kernels/* $(KERNEL_PATH)
	install -D solarv $(BINARY_PATH)/solarv

update-kernels:
	bash getkernels.sh || exit
	install -D -m 0644 kernels_tmp/* -t $(KERNEL_PATH)

.PHONY: clean check

clean:
	rm -f solarv solarv.o solarv.tm kernels/solarv.tm

solarv.o: solarv.h solarv.c Makefile

solarv: solarv.o solarv.h
	$(LD) $< -o solarv $(LDFLAGS)


kernels/solarv.tm: solarv.tm.in Makefile
	sed 's#KERNEL_SYSPATH#$(KERNEL_PATH)#g' solarv.tm.in > kernels/solarv.tm


kernels/stations.bsp: stations.defs
	rm -f kernels/stations.bsp
	$(PINPOINT) -def stations.defs \
	-pck $(KERNEL_PATH)/pck00010.tpc \
	-fk kernels/stations.tf \
	-spk kernels/stations.bsp
