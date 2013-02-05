#!/bin/bash

die ()
{
    echo "$1"
    exit 1
}

kernels="http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0010.tls 
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc 
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de421.bsp
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_720101_070527.bpc"


mkdir _tmp_ || die "could not create tempdir"
( cd _tmp_ && wget $kernels ) || die "could not spice kernels"




