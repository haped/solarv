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

    /* commandline args */
    if (! (argc == 2 || argc == 4)) {
	fprintf (stderr, "usage: solarv <time utc> [<stepsize> <numsteps>]\n");
	return EXIT_FAILURE;
    }
    strncpy (time_utc, argv[1], 79);
    if (argc == 4) {
	stepsize = atof (argv[2]);
	nsteps = atoi (argv[3]);
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
    et -= 66.1833;

    for (SpiceInt i  = 0; i < nsteps; ++i)
    {
	soleph_t eph;
	char utcstr[80];
	et2utc_c (et, "C", 2, 79, utcstr);

	SpiceDouble posx = 0, posy = 0;
	station_eph_xy ("izana", et, posx, posy, &eph);

	printf ("%.10f  %.6f  %.6f  %.3f  %.3f\n",
		eph.jdate, eph.B0, eph.L0, eph.s_vlos, eph.s_dist);
	et += stepsize * 60;
    }
}


int station_eph_xy (
    SpiceChar *station,
    SpiceDouble et,    
    SpiceDouble xpos,     
    SpiceDouble ypos,     
    soleph_t *eph      
    )
{
    SpiceDouble lt;
    SpiceDouble state [6];
    SpiceDouble subsc[3];
    SpiceDouble srfvec[3];
    SpiceDouble trgepc;
    SpiceDouble slat, slon, srad;
    SpiceDouble losvec[3];

    /* real julian date */
    eph->jdate = et / spd_c() + j2000_c();
    
    /* compute sub-observer point on the solar surface to derive B0, and L0 */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "IAU_SUN", "NONE", station,
	      subsc, &trgepc, srfvec);
    reclat_c (subsc, &srad, &slon, &slat);
    eph->B0 = slat * 180 / pi_c();
    eph->L0 = slon * 180 / pi_c();

    /* FIXME: implement calculation of P0 */
    /* P0 = .... */

    /* vlos and distance to sun's barycenter */
    spkezr_c ("SUN", et, "J2000",  "NONE", "IZANA", state, &lt);
    unorm_c (state, losvec, &(eph->s_dist));
    eph->s_dist *= 1000;
    eph->s_vlos = vdot_c (losvec, &state[3]) * 1000;


    /* compute target position epehems */


    return RETURN_SUCCESS;
}
