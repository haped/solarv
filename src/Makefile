include config.mk

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

ifndef CC
	CC = c99
endif
CFLAGS += $(SPICE_I) $(FITSIO_I) \
	-DKERNEL_PATH=\"$(KERNEL_PATH)\"
LD = $(CC) 
LDFLAGS += $(SPICE_L) $(FITSIO_L) -lm 


all: solarv solarv.tm

install: solarv solarv.tm
	install -D -m 0644 solarv.tm $(KERNEL_PATH)/solarv.tm
	install -D -m 0644 stations.bsp $(KERNEL_PATH)/stations.bsp
	install -D -m 0644 stations.tf $(KERNEL_PATH)/stations.tf
	install -D -m 0644 earth_fixed_itrf93.tf $(KERNEL_PATH)/earth_fixed_itrf93.tf
	install -D -m 0644 earth_fixed.tf $(KERNEL_PATH)/earth_fixed.tf
	install -D -m 0644 heliospheric.tf $(KERNEL_PATH)/heliospheric.tf
	install -D solarv $(BINARY_PATH)/solarv

update-kernels:
	./getkernels.sh || exit
	install -D -m 0644 _tmp_/de421.bsp $(KERNEL_PATH)/de421.bsp
	install -D -m 0644 _tmp_/earth_latest_high_prec.bpc \
		$(KERNEL_PATH)/earth_latest_high_prec.bpc
	install -D -m 0644 _tmp_/naif0010.tls $(KERNEL_PATH)/naif0010.tls
	install -D -m 0644 _tmp_/pck00010.tpc $(KERNEL_PATH)/pck00010.tpc
	install -D -m 0644 _tmp_/earth_assoc_itrf93.tf $(KERNEL_PATH)/earth_assoc_itrf93.tf
	rm -rf _tmp_

.PHONY: clean check

clean:
	rm -f solarv solarv.o solarv.tm

solarv.o: solarv.h solarv.c Makefile

solarv: solarv.o solarv.h
	$(LD) $< -o solarv $(LDFLAGS)


solarv.tm: solarv.tm.in Makefile
	@sed 's#KERNEL_SYSPATH#$(KERNEL_PATH)#g' solarv.tm.in > solarv.tm


stations.bsp: stations.defs
	rm -f stations.bsp
	$(PINPOINT) -def stations.defs \
	-pck $(KERNEL_PATH)/pck00010.tpc \
	-spk stations.bsp
