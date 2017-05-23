# if you have non-standard search paths for package configs, add them to
# your PKG_CONFIG_PATH environment variable
PKG_CONFIG      = /bin/pkg-config

## CSPICE INSTALL DIR
CSPICE_PATH = $(HOME)/local/pkg/solarv/cspice
## where required SPICE kernels are
KERNEL_PATH = $(HOME)/local/share/solarv/kernels
## where to put the 'solarv' binary
BINARY_PATH = $(HOME)/local/bin
## the pinpoint utility is needed to update the 
## stations.bsp binary kernel
PINPOINT = $(CSPICE_PATH)/exe/pinpoint

CC = c99
CFLAGS = -m64 -g -Wall -pedantic -O2
LDFLAGS = -m64

