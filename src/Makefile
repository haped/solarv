CSPICE_PATH = /home/doerr/coding/lfc/solarv/cspice
CSPICE_LIB = $(CSPICE_PATH)/lib/cspice.a
CSPICE_INC = -I$(CSPICE_PATH)/include

CC = c99
CFLAGS = -Wall -pedantic -O2 $(CSPICE_INC)
LD = $(CC)
LDFLAGS = $(CSPICE_LIB) -lm




all: solarv

solarv: solarv.o
	$(LD) $< -o solarv $(LDFLAGS)

.PHONY: clean


clean:
	rm -f solarv solarv.o


