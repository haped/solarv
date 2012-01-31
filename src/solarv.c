/*
  Copyright (c) 2012 Hans-Peter Doerr <doerr@kis.uni-freiburg.de>,
    Kiepenheuer-Institut fuer Sonnenphysik, Freiburg, Germany.
 
  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:
 
  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.
 
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "solarv.h"

int main(int argc, char **argv)
{   
    SpiceDouble   et;
    char time_utc[80];
    SpiceDouble stepsize = 1; /* minutes */
    SpiceInt nsteps = 1;
    SpiceDouble
	lon = 0.0,
	lat = 0.0;


    /* commandline args */
    if (! (argc == 4 || argc == 6)) {
	fprintf (stderr, "usage: solarv <time utc> <lon> <lat> [<stepsize> <numsteps>]\n");
	return EXIT_FAILURE;
    }
    strncpy (time_utc, argv[1], 79);

    if (argc == 4) {
	lon = atof (argv[2]);
	lat = atof (argv[3]);
    }
    
    if (argc == 6) {
	stepsize = atof (argv[4]);
	nsteps = atoi (argv[5]);
    }

    /*
      load kernels: LSK, PCK, planet/satellite SPK and MGS spacecraft SPK
    */
    furnsh_c ("../data/kernels/naif0009.tls");
    furnsh_c ("../data/kernels/pck00010.tpc");
    furnsh_c ("../data/kernels/de421.bsp");
    //furnsh_c ("../data/kernels/de405.bsp");
    furnsh_c ("../data/kernels/earth_fixed.tf");
    furnsh_c ("../data/kernels/earth_assoc_itrf93.tf");
    furnsh_c ("../data/kernels/earth_000101_120411_120119.bpc");
    //furnsh_c ("../data/kernels/earthstns_fx_050714.bsp");
    furnsh_c ("../data/kernels/izana.bsp");
    furnsh_c ("../data/kernels/izana.txt");

    utc2et_c (time_utc, &et);

    /* das ist der offset zu den daten aus horizons ... wo kommt das her? */
    //et -= 66.1833;

    printf ("#jdate  B0(deg)  L0(deg)  c_vlos(m/s)  c_dist(m)\n");
    for (SpiceInt i  = 0; i < nsteps; ++i)
    {
	soleph_t eph;
	station_eph ("izana", et, lon, lat, &eph);

	printf ("%.10f  %.6f  %.6f  %.3f  %.3f  %.3f  %.3f\n",
		eph.jdate, eph.B0, eph.L0, eph.c_vlos, eph.c_dist, eph.vlos, eph.dist);

	printf ("\n"
		"julian date    : %f\n"
		"UTC date       : %s\n"
		"B0             : %f deg\n"
		"L0             : %f deg\n"
		"target posn    : %f, %f deg\n"
		"target velocity: %f m/s\n"
		"target distance: %f km\n"
		"center velocity: %f m/s\n"
		"center distance: %f km\n"
		"velocity diff  : %f m/s\n"
		"distance diff  : %f km\n",
		eph.jdate,
		eph.utcdate,
		eph.B0,
		eph.L0,
		eph.lon, eph.lat,
		eph.vlos, eph.dist / 1000,
		eph.c_vlos, eph.c_dist / 1000,
		(eph.vlos - eph.c_vlos),
		(eph.dist - eph.c_dist) / 1000);

	et += stepsize * 60;
    }
}


int station_eph (
    SpiceChar *station, /* Observer's NAIF BODY_NAME ("Izana")     */
    SpiceDouble et,     /* spice ephemeris time                    */
    SpiceDouble lon,    /* stonyhurst longitude   (deg)            */
    SpiceDouble lat,    /* stonyhurst latitude    (deg)            */
    soleph_t *eph       /* pt. to struct where to store ephem data */
    )
{
    SpiceDouble lt;
    SpiceDouble state_otc[6];     /* obs. to sun center state vector */
    SpiceDouble state_ctt_sun[6]; /* center to target state vector IAU_SUN  */
    SpiceDouble state_ctt_j2k[6]; /* center to target state vector J2000    */
    SpiceDouble state_ctt[6];     /* center to target state vector J2000    */
    SpiceDouble state_ott[6];     /* obs. to target state vector     */
    SpiceDouble los_otc[3];       /* los vector obs. to sun center   */
    SpiceDouble los_ott[3];       /* los vector obs. to target       */
    SpiceDouble xform[6][6];
    SpiceDouble subsc[3];
    SpiceDouble srfvec[3], losvec[3];
    SpiceDouble trgepc;
    SpiceDouble slat, slon, srad;

    /* remember julian date and ascii utc string of the event */
    eph->jdate = et / spd_c() + j2000_c();
    et2utc_c (et, "C", 2, 79, eph->utcdate);

    
    /* compute sub-observer point on the solar surface to derive B0, and L0 */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "IAU_SUN", "NONE", station,
	      subsc, &trgepc, srfvec);
    reclat_c (subsc, &srad, &slon, &slat);
    eph->B0 = slat / rpd_c();
    eph->L0 = slon / rpd_c();

    /* FIXME: implement calculation of P0 */
    /* P0 = .... */

    /* vlos and distance to sun's barycenter */
    /* FIXME: what about spkgeo_c()?? looks the same! */
    spkezr_c ("SUN", et, "J2000",  "NONE", station, state_otc, &lt);
    unorm_c (state_otc, los_otc, &(eph->c_dist));
    eph->c_dist *= 1000;
    eph->c_vlos = vdot_c (los_otc, &state_otc[3]) * 1000;
    
    /* for the following procedure, also see:
       http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html#
       Transforming%20a%20body-fixed%20state%20to%20the%20inertial%20J2000%20frame
    */

    /* state vector from sun center to target position on surface; velocity
     * of target is zero with repsect to the IAU_SUN frame. need to convert
     * it to the J2000 frame afterwards */
    SpiceDouble tlon = lon * rpd_c() + slon;
    SpiceDouble tlat = lat * rpd_c() + slat;
    latrec_c (SOLAR_RADIUS / 1000.0, tlon, tlat, state_ctt_sun);
    state_ctt_sun[3] = 0.0;
    state_ctt_sun[4] = 0.0;
    state_ctt_sun[5] = 0.0;

    /* rotate _ctt state vector to J2000 frame */
    sxform_c ("IAU_SUN", "J2000", et, xform);
    mxvg_c (xform, state_ctt_sun, 6, 6, state_ctt);

    /* finally, dist & velocity to target in meter, meter/second */
    vaddg_c (state_otc, state_ctt, 6, state_ott);
    unorm_c (state_ott, los_ott, &(eph->dist));
    eph->dist *= 1000;
    eph->vlos = vdot_c (los_ott, &state_ott[3]) * 1000;
	
    eph->lon = lon;
    eph->lat = lat;

    return RETURN_SUCCESS;
}
