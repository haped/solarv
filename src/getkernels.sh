#!/bin/bash
owd=$PWD

die ()
{
    echo "$1"
    cd $owd
    exit 1
}

kernels="http://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0011.tls 
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc 
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_070425_370426_predict.bpc"

kern_tmp=kernels_tmp

test -d $kern_tmp || mkdir -p $kern_tmp

cd $kern_tmp || die "could not enter temp kernel dir"
rm -f *.tls *.tpc *.bsp *.tf *.bpc
wget $kernels || die "ERROR while getting the kernels"
cd $owd
