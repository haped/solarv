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
#include <strings.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>

#include "solarv.h"

void usage (FILE *stream)
{
    printf ("solarv v%s, %s by %s\n\n"
	    "compute precision radial velocities between earth-based "
	    "observatories\nand a given position on the sun.\n"
	    "\n"
	    "usage: solarv [options] <date> <lat lon> [<tstep> <nsteps]\n"
	    "options:\n"
	    "  -h          show this help\n"
	    "  -m model    set solar rotation model 'model'\n"
	    "              use 'list' to list available models\n"
	    "  -p          pretty-print ephemeris data\n"
	    "\n"
	    "'date' is a strftime() conforming string that might also include\n"
	    "a time zone. If no time zone is given, UTC is assumed.\n"
	    "Latitude and longitude are Stonyhurst coordinates given in degree.\n"
	    "The step width 'tstep' is in minutes.\n"
	    "\n"
	    "Examples\n"
	    "\n"
	    "Compute radial velocity at the western limb at 45deg latitude using\n"
	    "Snodgrass & Ulrich (1990) spectroscopy rotation model:\n"
	    "  solarv -p -m su90s 2010-01-01T00:00:00 90 45\n",
	    version, versiondate, author
	    
	);
}

int main(int argc, char **argv)
{   
    SpiceDouble et;
    char time_utc[80];
    SpiceDouble stepsize = 1; /* minutes */
    SpiceInt nsteps = 1;
    SpiceDouble lon = 0.0, lat = 0.0;
    bool fancy = false;
    int rotModel = fixed;

    int c; opterr = 0;
    while ((c = getopt (argc, argv, "+hm:pr:")) != -1)
    {
        switch (c)
        {
        case 'h': usage (stdout); return EXIT_SUCCESS;
        case 'm':
	    if (strcasecmp (optarg, "list") == 0) {
		printf ("available rotation models:\n");
		list_rotation_models (stdout);
		return EXIT_SUCCESS;
	    }
	    else if (strcasecmp (optarg, "rigid") == 0)
		rotModel = rigid;
	    else if (strcasecmp (optarg, "fixed") == 0)
		rotModel = fixed;
	    else if (strcasecmp (optarg, "su90s") == 0)
		rotModel = su90s;
	    else if (strcasecmp (optarg, "su90g") == 0)
		rotModel = su90g;
	    else if (strcasecmp (optarg, "su90m") == 0)
		rotModel = su90m;
	    else if (strcasecmp (optarg, "s84s") == 0)
		rotModel = s84s;
	    else {
		fprintf (stderr, "unknown rotation model: %s\n", optarg);
		list_rotation_models (stderr);
		return EXIT_FAILURE;
	    }					
	    break;
        case 'p': fancy = true; break;
        default: usage (stdout); return EXIT_FAILURE;
        }
    }

    /* commandline args */
    if (! (argc - optind == 3 || argc - optind == 5)) {
	usage (stdout);
	return EXIT_FAILURE;
    }
    strncpy (time_utc, argv[optind], 79);
    if (argc - optind == 3 || argc - optind == 5) {
	lon = atof (argv[optind + 1]);
	lat = atof (argv[optind + 2]);
    }
    if (argc - optind == 5) {
	stepsize = atof (argv[optind + 3]);
	nsteps = atoi (argv[optind + 4]);
    }
    
    /* required kernels are coded in solarv.tm  meta kernel */
    furnsh_c ("../data/kernels/solarv.tm");
    /* TODO: poss. to load a custom kernel here */

    str2et_c (time_utc, &et);

    if (! fancy)
	print_ephtable_head (stdout);
    for (SpiceInt i  = 0; i < nsteps; ++i)
    {
	soleph_t eph;
	soleph ("izana", et, lon, lat, &eph, rotModel);

	if (fancy) {
	    fancy_print_eph (stdout, &eph);
	    if (i < nsteps - 1) printf ("\n");
	}
	else
	    print_ephtable_row (stdout, &eph);

	et += stepsize * 60.0;
    }
    
    unload_c ("../data/kernels/solarv.tm");
    return EXIT_SUCCESS;
}



/*
 * this is the entry point to compute the solar ephemeris parameters for a
 * given position on the sun at a given time
 */
int soleph (
    SpiceChar *station, /* NAIF body name/code of the observer     */ 
    SpiceDouble et,     /* Spice ephemeris time of the observation */
    SpiceDouble lon,    /* stonyhurst longitude of target point    */
    SpiceDouble lat,    /* stonyhurst latitude of target point     */
    soleph_t *eph,
    int rotmodel)
{
    SpiceDouble subpoint[3];
    SpiceDouble srfvec[3];
    SpiceDouble trgepc;
    SpiceDouble subrad, sublon, sublat; /* lola cords of sub-observer point */
    SpiceDouble state_ots[6];
    SpiceDouble state_stt[6];
    SpiceDouble state_ott[6];
    SpiceDouble los_ott[3];
    SpiceDouble lonrd = lon * rpd_c();
    SpiceDouble latrd = lat * rpd_c();

    reset_soleph (eph);
    
    /* remember julian date and ascii utc string of the event */
    SpiceDouble deltaT;
    deltet_c (et, "ET", &deltaT);
    eph->jdate = unitim_c (et, "ET", "JDTDB") - deltaT / spd_c();
    et2utc_c (et, "C", 2, 79, eph->utcdate);
    strcpy (eph->station, station);

    eph->lon = lonrd;
    eph->lat = latrd;
    
    /* compute sub-observer point on the solar surface to derive B0, and L0 */
    /* FIXME: check if and which abberation correction is needed here. LT+S
     * should be correct because NONE would compute the geometric state at
     * 'et' while what we actually observe has happened ~8 mins before */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "IAU_SUN",
	      ABCORR, station, subpoint, &trgepc, srfvec);
    reclat_c (subpoint, &subrad, &sublon, &sublat);
    eph->B0 = sublat;
    eph->L0 = sublon; /* carrington longitude */

    /* compute solar barycenter state relative to observer */
    relstate_observer_sun (station, et, eph, state_ots);

    /* compute target state relative to sun barycenter */
    relstate_sun_target (station, et, lonrd, latrd, rotmodel, eph, state_stt);

    /* compute observer to target by adding ots + stt */
    vaddg_c (state_ots, state_stt, 6, state_ott);
    unorm_c (state_ott, los_ott, &(eph->dist));
    eph->vlos = vdot_c (los_ott, &state_ott[3]);
    
    /************************************************************************
     compute further ephemeris data from the parameters we gathered so far
    ************************************************************************/

    /* apparent diameter of the disc in arcsec */
    eph->rsun_as = RSUN / (eph->dist_sun * 1000) * dpr_c() * 3600.0;

    /* helicentric angle */
    /* FIXME: are we talkin about the same mu? */
    eph->mu = sqrt (1.0 - (eph->rho * 1000 / RSUN) * (eph->rho * 1000 / RSUN));

    return RETURN_SUCCESS;
}

/* compute relative state of observer to solar barycenter in the j2000
 * reference frame */
int relstate_observer_sun (
    SpiceChar *station,   /* NAIF body name of observer ("Izana") */
    SpiceDouble et,       /* ephemeris time of event              */
    soleph_t *eph,
    SpiceDouble *state_otc /* obs. to sun center state vector */
    )
{
    SpiceDouble lt;
    SpiceDouble los_otc[3];       /* los vector obs. to sun center   */
    
    /* vlos and distance to sun barycenter: km, km/s */
    spkezr_c ("SUN", et, "J2000",  "LT+S", station, state_otc, &lt);
    unorm_c (state_otc, los_otc, &(eph->dist_sun));
    eph->vlos_sun = vdot_c (los_otc, &state_otc[3]);

    return RETURN_SUCCESS;
}

/* compute relative state from solar center to target position on the
 * surface in IAU_SUN frame */
int relstate_sun_target (
    SpiceChar *station,
    SpiceDouble et,
    SpiceDouble lon,    /* stonyhurst lon in radian */
    SpiceDouble lat,    /* stonyhurst lat in radian */
    int rotmodel,
    soleph_t *eph,
    SpiceDouble *state_stt)
{
    SpiceDouble subpoint[3];
    SpiceDouble srfvec[3];
    SpiceDouble trgepc;
    SpiceDouble subrad, sublon, sublat; /* lola cords of sub-observer point   */
    SpiceDouble xform[6][6];            /* transf. matrix body-fixed -> j2000 */
    SpiceDouble state_stt_fixed[6];     /* sun-to-target state vector         */

    /* for the following procedure, also see:
       http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html#
       Transforming%20a%20body-fixed%20state%20to%20the%20inertial%20J2000%20frame
    */

    /* again we compute the sub-observer point, but this time in the fixed
     * Heliocentric Inertial frame. the sub-observer longitude computed in
     * such a way is then used as an offset to compute the rectangluar
     * target position on a sphere with solar radius. this coordinates are
     * then rotated in the j2000 frame and added to the j2000 state of the
     * solar barycenter to get the target position on the sun relative to
     * observer in j2000 coordinates */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "HCI", ABCORR, station,
	      subpoint, &trgepc, srfvec);
    reclat_c (subpoint, &subrad, &sublon, &sublat);
    
    /* rectangular, body fixed coorindates of target point given in
     * stonyhurst coordinates */
    SpiceDouble tlon = lon + sublon;
    SpiceDouble tlat = lat;
    latrec_c (RSUN / 1000.0, tlon, tlat, state_stt_fixed);

    /* radial distance from target to solar rotation axis */
    eph->rho = sqrt (state_stt_fixed[0] * state_stt_fixed[0] + 
		     state_stt_fixed[1] * state_stt_fixed[1]);

    /* additional veloctiy from (differential) rotation according to
     * rotation model specified */
    SpiceDouble omegas = omega_sun (lat, rotmodel);
    SpiceDouble omega[] = {0.0, 0.0, omegas};
    SpiceDouble radius[] = {state_stt_fixed[0], state_stt_fixed[1], 0.0};
    vcrss_c (omega, radius, &state_stt_fixed[3]);

    eph->omega = omegas * 1E6; /* we store omega in murad/s */
    eph->rotmodel = rotmodel;
    strcpy (eph->modelname, RotModels[rotmodel].name);
    strcpy (eph->modeldescr, RotModels[rotmodel].descr);

    /* now, rotate body-fixed state to j20000 */
    sxform_c ("HCI", "J2000", et, xform);
    mxvg_c (xform, state_stt_fixed, 6, 6, state_stt);

    /* helioprojective-cartesian coordinates (pointing coordinates). use the
     * stt_fixed vector again. this time we substract the sub-observer
     * latitude. sub-observer longitude is always zero in stondyhurst
     * coordinates per definition */
    latrec_c (RSUN / 1000.0, lon, lat - sublat, state_stt_fixed);

    /* Thomposn (2005) eq 4 */
    eph->x = state_stt_fixed[1] / eph->dist_sun * dpr_c () * 3600.0;
    eph->y = state_stt_fixed[2] / eph->dist_sun * dpr_c () * 3600.0;

    return RETURN_SUCCESS;
}

void reset_soleph (soleph_t *eph)
{
    eph->jdate = 0.0;
    strcpy (eph->utcdate, "na");
    strcpy (eph->station, "na");
    eph->B0 = 0.0;
    eph->L0 = 0.0;
    eph->P0 = 0.0;
    eph->dist_sun = 0.0;
    eph->vlos_sun = 0.0;
    eph->rsun_as = 0.0;
    eph->rotmodel = 0;
    strcpy (eph->modelname, "na");
    strcpy (eph->modeldescr, "na");

    eph->lon = 0.0;
    eph->lat = 0.0;
    eph->x = 0.0;
    eph->y = 0.0;
    eph->mu = 0.0;
    eph->dist = 0.0;
    eph->vlos = 0.0;
    eph->rho = 0.0;
    eph->omega = 0.0;
}


SpiceDouble omega_sun (SpiceDouble lat, int model)
{
    if ( (model > RotModel_END - 1 || model < 0)  && model != custom) {
	fprintf (stderr, "invalid rotation model: %i\n", model);
	exit (EXIT_FAILURE);
    }
    SpiceDouble
	A = RotModels[model].A * 1E-6,
	B = RotModels[model].B * 1E-6,
	C = RotModels[model].C * 1E-6;

    SpiceDouble sin2lat = sin (lat) * sin (lat);
    SpiceDouble sin4lat = sin2lat * sin2lat;
    SpiceDouble omega = A + B * sin2lat + C * sin4lat;

    return omega;
}

void print_ephtable_head (FILE *stream)
{
    fprintf (stream, "#jdate          B0(deg)  P0(deg)   L0(deg)   "
	     "vlos_t     d_t             vlos_sun  d_sun \n");
}

void print_ephtable_row (FILE *stream, soleph_t *eph)
{
    fprintf (stream, "%.6f %7.4f %9.4f %10.4f %9.4f  %14.3f  %9.4f  %13.3f\n",
	     eph->jdate, eph->B0 * dpr_c(), eph->P0 * dpr_c(), eph->L0 * dpr_c(),
	     eph->vlos * 1000, eph->dist, eph->vlos_sun * 1000, eph->dist_sun);
}

void fancy_print_eph (FILE *stream, soleph_t *eph)
{
	    printf ("Solar ephemeris for %s, %s, pos (%.3f, %.3f) deg\n"
		    "  julian date          : %f\n"
		    "  disk radius          : %.2f arcsec\n"
		    "  B0                   : %.4f deg\n"
		    "  P0                   : %.4f deg\n"
		    "  carrington L0        : %.4f deg\n"
		    "  center distance      : %.3f km\n"
		    "  center v_los         : %.3f m/s\n"
		    "  target position      : %.2f, %.2f arcsec\n"
		    "  heliocentric angle   : %.4f\n"
		    "  target distance      : %.3f km\n"
		    "  target v_los         : %.3f m/s\n"
		    "  solar rotation model : %s (%s)\n"
		    "  rotataion rate       : %.5f murad/s\n"
		    "  impact parameter     : %.3f km\n"
		    "  distance diff        : %.3f km\n",
		    eph->station, eph->utcdate,
		    eph->lon * dpr_c(), eph->lat * dpr_c(),
		    
		    eph->jdate,
		    eph->rsun_as,
		    eph->B0 * dpr_c(),
		    eph->P0 * dpr_c(),
		    eph->L0 * dpr_c(),
		    eph->dist_sun,
		    eph->vlos_sun * 1000,
		    eph->x, eph->y,
		    eph->mu,
		    eph->dist,
		    eph->vlos * 1000,
		    eph->modelname, eph->modeldescr,
		    eph->omega,
		    eph->rho,
		    (eph->dist - eph->dist_sun));
}


void list_rotation_models (FILE *stream)
{
    fprintf (stream, "name        A       B       C       description\n");
    fprintf (stream, "--------------------------------"
	    "--------------------------------\n");
    for (int i = 0; i < RotModel_END; ++i) {
	fprintf (stream, "%-10s % 6.4f % 6.4f % 6.4f  %s\n",
		 RotModels[i].name,
		 RotModels[i].A, RotModels[i].B, RotModels[i].C,
		 RotModels[i].descr);
    }
}
