#!/bin/bash
owd=$PWD

die ()
{
    echo "$1"
    cd $owd
    exit 1
}

kernels="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls 
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc 
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc"

kern_tmp=kernels_tmp

test -d $kern_tmp || mkdir -p $kern_tmp

cd $kern_tmp || die "could not enter temp kernel dir"
rm -f *.tls *.tpc *.bsp *.tf *.bpc
wget $kernels || die "ERROR while getting the kernels"
cd $owd

