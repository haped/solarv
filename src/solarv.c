/*
  Copyright (c) 2012, 2013 Hans-Peter Doerr <doerr@kis.uni-freiburg.de>,
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

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <libgen.h>
#include <ctype.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdarg.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <cfitsio/fitsio.h>

#include "solarv.h"

void usage (FILE *stream)
{
    fprintf (stream, 
	     "%s v%s (%s) -- Compute precise Sun-Observer radial velocity\n"
	     "\n"
	     "Usage:\n"
	     "  solarv [options] [request | file]\n"
	     "\n"
	     "Request specification:\n"
	     "  <timespec> <shg|hpc|hcr> <cord1> <cord2> [<nsteps> <tstep>]\n"
	     "\n"
	     "The request is either read from the commandline, from 'file' "
	     "or from STDIN if no\n"
	     "arguments are given. The position can be specified in "
	     "helio-projective cartesian\n"
	     "(hpc) coordinates in arcseconds from disk center, Stonyhurst "
	     "heliographic (shg) \n"
	     "coordinates in degree longitude and latitude or in helio-centric "
	     "radial (hcr)\n"
	     "coordinates given in cos(theta) from disk center and an angle "
	     "phi counter-\n"
	     "clockwise from heliographic north. The 'timespec' parameter "
	     "understands common\n"
	     "time strings like '2010-01-01 12:00:00 EST'. If no time zone is "
	     "given, UTC is\n"
	     "assumed. I highly recommend to use the ISO 8601 standard, "
	     "e.g.\n"
	     "'2012-01-01T00:00:00'. The optional 'nsteps' "
	     "and 'tstep' (in minutes) arguments\n"
	     "can be used to compute a set of values in a given time interval."
	     "\n\n"
	     "Options:\n"
	     "  -h            show this help\n"
	     "  -v            print program version\n"
	     "  -p            pretty-print output\n"
	     "  -m model      set solar rotation model 'model' [fixed]\n"
	     "                use 'list' to list available models\n"
	     "  -t            use the ITRF93 precision earth rotation kernel\n"
	     "                Note: this kernel has a limited time coverage and needs\n"
	     "                to be updated regularily\n"
	     "  -R radius     specifiy a different solar radius in km\n"
	     "  -O observer   set observer position. Can be any NAIF body code.\n"
	     "                pre-defined sites: 'VTT', 'SCHAUINSLAND', 'SST',\n"
	     "                'DST', 'MCMATH', 'BIGBEAR'\n"
	     "  -k kernel     load additional SPICE kernel 'kernel'\n"
	     "                this kernel will be loaded last in kernel pool\n"
	     "  -K kernel     load 'kernel' instead of the default meta kernel\n"
	     "  -i            Show information about loaded SPICE kernels\n"
	     "\n"
	     "Examples:\n"
	     "\n"
	     "Compute radial velocity at disk center, with no solar rotation. "
	     "Use the the high\n"
	     "precision ITRF93 earth rotation model, pretty-print output:\n"
	     "  solarv -t -p 2012-06-21T12:00:00 hpc 0 0"
	     "\n\n"
	     "Create a ephemeris table for the SST for one year, close to the "
	     "eastern limb\n"
	     "using the Snodgrass 1984 rotation model:\n"
	     "  solarv -t -O SST -m s84s 2012-01-01T12:00:00.0 hcr 0.003 90 365 1440\n"
	     "\n"
	     ,
	     _name, _version, _versiondate
	);
}

void program_info (FILE *stream)
{
    fprintf (stream, "%s v%s (%s)\n"
	     "SPICE toolkit version: %s\n"
	     "Copyright (C) %s\n"
	     "Licencse: MIT X11\n",
	     _name, _version, _versiondate, tkvrsn_c("toolkit"), _copyright);
}

void version_info (FILE *stream)
{
    fprintf (stream, "%s v%s (%s)\n",
	     _name, _version, _versiondate);
}

void dump_kernel_info (FILE *stream)
{
    const SpiceInt TYPLEN = 32;
    const SpiceInt SRCLEN = 128;
    SpiceInt  handle;
    SpiceInt nkernels;
    SpiceBoolean found;
    
    fprintf (stream, "Loaded SPICE kernels:\n");
    ktotal_c ("all", &nkernels);
    for (int i = 0; i < nkernels; ++i)
    {
	SpiceChar file[MAXPATH+1];
	SpiceChar filetype[TYPLEN+1];
	SpiceChar source[SRCLEN+1];
	
	kdata_c (i, "all", MAXPATH, TYPLEN, SRCLEN,
		 file, filetype, source, &handle, &found);
	kinfo_c (file, TYPLEN, SRCLEN, filetype, source, &handle, &found);
	fprintf (stream, "   %-5s: %s\n", filetype, file);
    }
}


/* get the body centered, body fixed coordinates of a location on earth */
void station_geopos (
    SpiceChar * station,
    SpiceDouble et,
    SpiceDouble *lon,
    SpiceDouble *lat,
    SpiceDouble *alt)
{
    SpiceDouble state_station[6];
    SpiceDouble lt;
    SpiceDouble abc[3], f, r_eq, r_pl;
    SpiceInt dim;

    spkezr_c (station, et, "EARTH_FIXED", "NONE", "EARTH", state_station, &lt);
    bodvrd_c( "EARTH", "RADII", 3, &dim, abc);
    r_eq = abc[0];
    r_pl = abc[2]; 
    f = (r_eq - r_pl) / r_eq;
    recgeo_c (state_station, r_eq, f, lon, lat, alt);
}

/* get inertial state of a target on the sun, given it's heliographic
 * coordinates an the inertial state of the solar barycenter */
int getstate_solar_target (
    SpiceDouble et,
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble omegas,
    SpiceDouble *sunstate,  /* in: inertial state of sun barycenter */
    SpiceDouble *tgtstate)  /* out: intertial state of target */
{
    SpiceDouble	relstate[6] = {0};
    SpiceDouble relstatei[6] = {0};
    SpiceDouble diffrot[6] = {0};
    SpiceDouble diffroti[6] = {0};

    latrec_c (RSUN, lon, lat, relstate);
    
    /* compute state velocity from given omega */
    SpiceDouble radius[] = {relstate[0], relstate[1], 0.0};
    SpiceDouble omegav[] = {0.0, 0.0, omegas};
    vcrss_c (omegav, radius, &diffrot[3]);

    /* 
     * transform relative state and the differential rotation component to
     * the inertial state. this procedure is a bit nasty, since HEEQ is not
     * an inertial state such that velocity is not preserved in the
     * transform. since the velocity of the target must be the same as that
     * of the sun center in the interital system (without rotation), we just
     * set it to zero. the velocity of the diffrot vector is preserved, as
     * it's spatial coordinates are zero.
    */
    SpiceDouble xform[6][6];
    sxform_c ("HEEQ", "J2000", et, xform);
    mxvg_c (xform, relstate, 6, 6, relstatei);
    mxvg_c (xform, diffrot, 6, 6, diffroti);

    /* make sure the velocity of the inertial relstate is zero */
    relstatei[3] = relstatei[4] = relstatei[5] = 0.0;

    /* get inertial state, add contribution from diff. rotation */
    vaddg_c (sunstate, relstatei, 6, tgtstate);
    vaddg_c (tgtstate, diffroti, 6, tgtstate);

    return SUCCESS;
}

/* get inertial state of an observer on any body */
int getstate_observer (
    SpiceChar *body,
    SpiceDouble et,
    SpiceDouble lon,  /* -pi .. pi     */
    SpiceDouble lat,  /* -pi/2 .. pi/2 */
    SpiceDouble alt,  /* kilometer     */
    SpiceDouble *state,
    SpiceDouble *bstate
    )
{
    SpiceDouble bodstate[6];   /* body inertial state             */
    SpiceDouble relstate[6];   /* observer state relative to body */
    SpiceDouble relstatei[6];  /* observer relative to inertial   */ 

    SpiceInt dim, frcode;
    SpiceDouble radi[3], equatr, polar, f;
    SpiceBoolean found;

    /* positon vector relative to body center */
    bodvrd_c (body, "RADII", 3, &dim, radi);
    equatr = radi[0];
    polar = radi[2];
    f = (equatr - polar) / equatr;
    georec_c (lon, lat, alt, equatr, f, relstate);

    /* find native reference frame of 'body' and use it to translate the
     * position vector to a inertial state */
    SpiceChar frname[MAXKEY + 1];
    cnmfrm_c (body, MAXKEY, &frcode, frname, &found);
    if (! found) {
	errmesg ("Could not find coordinate system for body %s\n", body);
	return FAILURE;
    }
    SpiceDouble xform[6][6];
    sxform_c (frname, "J2000", et, xform);
    mxvg_c (xform, relstate, 6, 6, relstatei);

    /* add observer state to body center inertial state */
    getstate_body (body, et, bodstate);
    vaddg_c (bodstate, relstatei, 6, state);

    /* let caller know the body barycenter state if he wants so */
    if (bstate) {
	memcpy (bstate, bodstate, 6 * sizeof(SpiceDouble));
    }

    return SUCCESS;
}

/* get J2000 inertial state of a body */
void getstate_body (
    SpiceChar *body,
    SpiceDouble et,
    SpiceDouble *state)
{
    SpiceDouble ltt;
    spkezr_c (body, et, "J2000", "NONE", "SSB",  state, &ltt);
}

/* get the relative state vector between two state vectors */
void relstate (
    SpiceDouble *sobs,
    SpiceDouble *stgt,
    SpiceDouble *srel,
    SpiceDouble *losv,  /* line of sight vector */
    SpiceDouble *dist,
    SpiceDouble *vrad,  /* radial velocity */
    SpiceDouble *lt)
{
    vsubg_c (stgt, sobs, 6, srel);
    unorm_c (srel, losv, dist);
    *vrad = vdot_c (losv, &srel[3]);
    *lt = *dist / clight_c();
}


int main (int argc, char **argv)
{   
    FILE *batchstream = stdin;
    bool batchmode = false;
    char addkernel[MAXPATH+1] = "na";
    char metakernel[MAXPATH+1] = "na";
    char earthkernel[MAXPATH+1] = "na";
    char observer[MAXKEY+1] = "VTT";
    char fitsfile[MAXPATH+1] = "na";
    bool fancy = false;
    int rotmodel = fixed;
    bool earth_itrf93 = false;
    int errorcode = SUCCESS;
    bool dumpinfo = false;
    
    int c; opterr = 0;
    while ((c = getopt (argc, argv, "+h:m:pr:O:vfK:k:tiR:")) != -1)
    {
        switch (c)
        {
	case 'i': dumpinfo = true; break;
        case 'h': usage (stdout); return EXIT_SUCCESS;
        case 'm':
	    if (strcasecmp (optarg, "list") == 0) {
		list_rotation_models (stdout);
		return EXIT_SUCCESS;
	    }
	    else if (strcasecmp (optarg, "rigid") == 0) rotmodel = rigid;
	    else if (strcasecmp (optarg, "crgt") == 0)	rotmodel = crgt;
	    else if (strcasecmp (optarg, "fixed") == 0) rotmodel = fixed;
	    else if (strcasecmp (optarg, "su90s") == 0)	rotmodel = su90s;
	    else if (strcasecmp (optarg, "su90g") == 0)	rotmodel = su90g;
	    else if (strcasecmp (optarg, "su90m") == 0)	rotmodel = su90m;
	    else if (strcasecmp (optarg, "s84s") == 0)	rotmodel = s84s;
	    else {
		errmesg ("Unknown Rotation Model: %s\n", optarg);
		return EXIT_FAILURE;
	    }					
	    break;
        case 'p': fancy = true; break;
	case 'O': strncpy (observer, optarg, MAXKEY); break;
	case 'v': program_info (stdout); return EXIT_SUCCESS; break;
	case 'k': strncpy (addkernel, optarg, MAXPATH); break;
	case 'K': strncpy (metakernel, optarg, MAXPATH); break;
	case 't': earth_itrf93 = true; break;
	case 'R': RSUN = atof (optarg);  break;
        default: usage (stdout); return EXIT_FAILURE;
        }
    }
    int posargs = argc - optind;
    
    if (posargs == 0) {
	batchmode = true;
	batchstream = stdin;
    }
    else if (posargs == 1) {
	batchmode = true;
	FILE *b = fopen (argv[optind], "r");
	if (NULL == b) {
	    errmesg ("Can't read input file '%s'\n", argv[optind]);
	    return EXIT_FAILURE;
	}
	batchstream = b;
    }
    else if (posargs != 0 && posargs != 1 && posargs != 4 && posargs != 6) {
	usage (stdout);
	return EXIT_FAILURE;
    }
    
    if (dumpinfo) {
	version_info (stdout);
	fprintf (stdout, "CSPICE toolkit version: %s\n", tkvrsn_c ("toolkit"));
    }

    if (earth_itrf93)
	snprintf (earthkernel, MAXPATH, "%s/earth_fixed_itrf93.tf", KERNEL_PATH);
    else
	snprintf (earthkernel, MAXPATH, "%s/earth_fixed.tf", KERNEL_PATH);


    /* load requested kernels */
    if (strcmp (metakernel, "na") == 0)
	strncpy (metakernel , KERNEL_PATH "/" METAKERNEL, MAXPATH);
    furnsh_c (metakernel);
    furnsh_c (earthkernel);

    if (strcmp (addkernel, "na") != 0) {
	printf ("Loading user-supplied additional kernel %s\n", addkernel);
	furnsh_c (addkernel);
    }

    if (dumpinfo)
	dump_kernel_info (stdout);
    
    /* the real work is done here */
    if (batchmode) {
	errorcode = mode_batch (observer, rotmodel,
				fancy,  fitsfile, batchstream, stdout);
    } else {
	if (! fancy)
	    print_ephtable_head (stdout, observer, rotmodel);
	errorcode = handle_request (observer, posargs, &argv[optind],
				    rotmodel, fancy, false, stdout);
    }
    
    unload_c (earthkernel);
    unload_c (metakernel);
    if (strcmp (addkernel, "na") != 0)
	unload_c (addkernel);
    
    if (FAILURE == errorcode) {
	fprintf (stderr, "There were Errors. Please Check!\n");
	return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

int anypos2lola (
    sunpos_t pos,
    SpiceDouble et,
    SpiceDouble *state_obs,
    SpiceDouble *state_sun,
    SpiceDouble *lon,
    SpiceDouble *lat)
{
    if (pos.type == cord_shg)
    {
	*lon = pos.x * rpd_c();
	*lat = pos.y * rpd_c();
    }
    else if (pos.type == cord_hpc)
    {
	SpiceDouble xp = pos.x * rpas();
	SpiceDouble yp = pos.y * rpas();
	SpiceBoolean on;
	pointing2lola (state_sun, et, xp, yp, state_obs, lon, lat, &on);
	if (! on) {
	    errmesg ("Coordinates xy=(%.3f, %.3f) not within "
		     "the visible disk\n", xp * aspr(), yp * aspr());
	    return FAILURE;
	}
    }
    else if (pos.type == cord_hcr)
    {
	if (pos.x < 0.0  || pos.x > 1.0) {
	    errmesg ("Value of mu not in range [0,1]\n");
	    return FAILURE;
	}
	SpiceDouble diff[3];
	vsub_c (state_sun, state_obs, diff);

	/* FIXME: this is approximate, do the real math */
	SpiceDouble phi = pos.y * rpd_c();
	SpiceDouble theta = acos (pos.x);
	SpiceDouble ang = RSUN / vnorm_c(diff);
	SpiceDouble rads = tan (ang);
	SpiceDouble xp = sin(phi) * sin (theta) * rads;
	SpiceDouble yp = cos(phi) * sin (theta) * rads;
	SpiceBoolean on;
	pointing2lola (state_sun, et, xp, yp, state_obs, lon, lat, &on);
	if (! on) {
	    errmesg ("Coordinates xy=(%.3f, %.3f) not within "
		     "the visible disk\n", xp * aspr(), yp * aspr());
	    return FAILURE;
	}
    }

    return SUCCESS;
}

/*
 * this is the entry point to compute the solar ephemeris parameters for a
 * given position on the sun at a given time
 */
int soleph (
    SpiceChar *observer, /* NAIF body name/code of the observer     */ 
    SpiceDouble et,      /* Spice ephemeris time of the observation */
    sunpos_t position,
    soleph_t *eph,
    int rotmodel)
{
    SpiceDouble subpoint[3];
    SpiceDouble srfvec[3];
    SpiceDouble trgepc;
    SpiceDouble subrad, sublon, sublat; /* lola cords of sub-observer point */
    SpiceDouble
	state_obs[6],       /* inertial state of the observer */
	state_sun[6],       /* inertial state of sun barycenter */
	state_tgt[6],       /* inertial state of the target */
	relstate_sun[6],    /* relative state observer -> sun barycenter */
	relstate_tgt[6];    /* relative state observer -> target */

    reset_soleph (eph);
    strncpy (eph->observer, observer, MAXKEY);
    /* we make heavy use of the observer-centric radial tangential normal
     * frame which is referenced to the observer station */
    pcpool_c ("FRAME_1803430_PRI_TARGET", 1, strlen (observer)-1,  observer);

    /* remember julian date and ascii utc string of the event */
    SpiceDouble deltaT;
    deltet_c (et, "ET", &deltaT);
    eph->et = et;
    /* we use the utc julain date here, so substract the delta */
    eph->jday = unitim_c (et, "ET", "JDTDB") - deltaT / spd_c();
    eph->mjd = eph->jday - 2400000.5;
    timout_c (et, "DD Mon YYYY HR:MN:SC.## UTC", 64, eph->utcdate);

    /* stonyhurst heliographic sub-observer coordinates */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "HEEQ",
    	      ABCORR, observer, subpoint, &trgepc, srfvec);
    reclat_c (subpoint, &subrad, &sublon, &sublat);
    eph->B0 = sublat;
    eph->L0hg = sublon;

    /* carrington sub-observer coordinates, carrington is 0...2pi */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "IAU_SUN",
    	      ABCORR, observer, subpoint, &trgepc, srfvec);
    reclat_c (subpoint, &subrad, &sublon, &sublat);
    eph->L0cr = (sublon < 0 ? 2 * pi_c() + sublon : sublon);

    /*
     * first, we compute the relative state between observer and sun
     * barycenter. we need this now to be able to translate between pointing
     * coordinates and heliographic coordinates
     */
    SpiceDouble los_tgt[3], los_sun[3], lt_sun, lt_tgt;
    //getstate_observer ("earth", et, 0, 0, 0, state_obs, NULL);
    getstate_body (observer, et, state_obs);
    getstate_body ("SUN", et, state_sun);
    memcpy (eph->state_sun, state_sun, 6 * sizeof(SpiceDouble));
    memcpy (eph->state_obs, state_obs, 6 * sizeof(SpiceDouble));

    relstate (state_obs, state_sun, relstate_sun,
	      los_sun, &eph->dist_sun, &eph->vlos_sun, &lt_sun);
    
    /* figure out target heliographic coordinates from use input */
    if (SUCCESS != (anypos2lola (position, et, state_obs, state_sun,
				 &eph->lon, &eph->lat))) {
	errmesg ("Could not process target coordinates\n");
	return FAILURE;
    }
    
    /* angular velocity from rotation model is needed to compute the full
     * target state */
    SpiceDouble omegas = omega_sun (eph->lat, rotmodel);
    eph->omega = omegas;
    eph->rotmodel = rotmodel;
    strncpy (eph->modelname, RotModels[rotmodel].name, MAXKEY);
    strncpy (eph->modeldescr, RotModels[rotmodel].descr, MAXKEY);

    /* TODO: handle light-time correction */
    getstate_solar_target
	(et, eph->lon, eph->lat, omegas, state_sun, state_tgt);
    relstate (state_obs, state_tgt, relstate_tgt,
	      los_tgt, &eph->dist, &eph->vlos, &lt_tgt);
    get_pointing (et, relstate_sun, relstate_tgt, &eph->x, &eph->y);


    /************************************************************************
     compute further ephemeris data from the parameters we gathered so far
    ************************************************************************/
    eph->rsun_ref = RSUN;
    eph->rsun_obs = asin (RSUN / eph->dist_sun);

    /* impact parameter (projected  distance to solar center) */
    SpiceDouble obsangle = vsep_c (relstate_tgt, relstate_sun);
    eph->rho = eph->dist * sin (obsangle);
    
    /* angle between sub-observer-point -> solar center, center -> target */
    SpiceDouble sunangle = asin (eph->rho / RSUN);

    /* heliocentric angle is the angle between the los and the surface
     * normal at the target */
    eph->mu = cos (obsangle + sunangle);

    /* P angle */
    SpiceDouble e[3] = {0, 0, 1}, sol[3];
    SpiceDouble xform[3][3];
    pxform_c ("EARTH_FIXED", "OBSCRTN", et, xform);
    mxv_c (xform, e, sol);
    eph->P0 = atan2 (sol[1], sol[2]);

    return SUCCESS;
}


void get_pointing (
    SpiceDouble et,
    SpiceDouble *relstate_sun,
    SpiceDouble *relstate_tgt,
    SpiceDouble *x,
    SpiceDouble *y)
{
    SpiceDouble sun[3], tgt[3];

    SpiceDouble xform[3][3];
    pxform_c ("J2000", "OBSCRTN", et, xform);
    mxv_c (xform, relstate_tgt, tgt);
    mxv_c (xform, relstate_sun, sun);

    SpiceDouble sun_x[3] = {sun[0], sun[1], 0};
    SpiceDouble sun_y[3] = {sun[0], 0, sun[2]};
    SpiceDouble tgt_x[3] = {tgt[0], tgt[1], 0};
    SpiceDouble tgt_y[3] = {tgt[0], 0, tgt[2]};

    /* FIXME: more elegant way to handle negative angles? */
    *x = vsep_c (sun_x, tgt_x);
    *y = vsep_c (sun_y, tgt_y);
    if (tgt_x[1] < sun_x[1]) *x = -*x;
    if (tgt_y[2] < sun_y[2]) *y = -*y;
}


/* pointing (angle from disk center) coordinatates to planetographic */
void pointing2lola (
    SpiceDouble *state_sun,
    SpiceDouble et,
    SpiceDouble x,           /* pointing angle from disk center */
    SpiceDouble y,           /* pointing angle from disk center */
    SpiceDouble *state_obs,  /* observer inertial state */
    SpiceDouble *lon,
    SpiceDouble *lat,
    SpiceBoolean *onbody)
{
    SpiceDouble obsbd[3], vsurf[3];
    vsub_c (state_obs, state_sun, obsbd);
    SpiceDouble sundist = vnorm_c (obsbd);

    SpiceDouble xform[3][3]; 
    pxform_c ("J2000", "OBSCRTN", et, xform);
    mxv_c (xform, obsbd, obsbd);

    SpiceDouble xdist = tan(x) * sundist;
    SpiceDouble ydist = tan(y) * sundist;
    SpiceDouble vpnt[3] = {-sundist, xdist, ydist};

    surfpt_c (obsbd, vpnt, RSUN, RSUN, RSUN, vsurf, onbody);
    pxform_c ("OBSCRTN", "HEEQ", et, xform);
    mxv_c (xform, vsurf, vsurf);

    SpiceDouble radius;
    reclat_c (vsurf, &radius, lon, lat);
}


void reset_soleph (soleph_t *eph)
{
    eph->jday = 0.0;
    eph->mjd = 0.0;
    strcpy (eph->utcdate, "na");
    strcpy (eph->observer, "na");
    eph->B0 = 0.0;
    eph->L0cr = 0.0;
    eph->L0hg = 0.0;
    eph->P0 = 0.0;
    eph->dist_sun = 0.0;
    eph->vlos_sun = 0.0;
    eph->rsun_ref = 0.0;
    eph->rsun_obs = 0.0;
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

    for (int i = 0; i < 6; ++i) {
	eph->state_obs[i] = 0.0;
	eph->state_sun[i] = 0.0;
	eph->state_rel[i] = 0.0;
    }
}

SpiceDouble omega_sun (SpiceDouble lat, int model)
{
    if ( (model > RotModel_END - 1 || model < 0)  && model != custom) {
	errmesg ("unknown rotation model %i, aborting\n", model);
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


void timestr_now (char *buf)
{
    time_t t = time (NULL);
    asctime_r (gmtime (&t), buf);
    buf[strlen (buf) - 1] = '\0'; // who had the idea with the trailing \n?
}

void print_ephtable_head (FILE *stream, SpiceChar *observer, SpiceInt rotmodel)
{
    char now[48];
    timestr_now (now);

    fprintf (stream,
	     "#*****************************************************"
	     "*************************\n"
	     "#  Ephemeris Data created on %s UTC by solarv v%s\n"
	     "#  SPICE Version      : %s\n"
	     "#  Reference Epoch    : J2000.0\n"
	     "#  Inert. Refr. Frame : ICRF/J2000, "
	     "Mean Equator and Equinox of Epoch\n"
	     "#  Abbr. Correction   : %s\n"
	     , now, _version, tkvrsn_c ("toolkit"), ABCORR);
    
    double lon, lat, alt;
    station_geopos (observer, 1.0, &lon, &lat, &alt);
    

    /* FIXME: really check if we are on earth or not */
    bool onEarth = true;
    if (0 != strcasecmp (observer, "VTT") &&
	0 != strcasecmp (observer, "SST") &&
	0 != strcasecmp (observer, "SCHAUINSLAND") &&
	0 != strcasecmp (observer, "DST") &&
	0 != strcasecmp (observer, "MCMATH") &&
	0 != strcasecmp (observer, "BIGBEAR")
	) {
	onEarth = false;
    }

    SpiceBoolean found = 0;
    SpiceInt frcode;
    SpiceChar frname[MAXKEY+1];
    if (onEarth)
	cnmfrm_c ("EARTH", MAXKEY, &frcode, frname, &found);
    else
	cnmfrm_c (observer, MAXKEY, &frcode, frname, &found);
    
    fprintf (stream,
	     "#*****************************************************"
	     "*************************\n"
	     "#  Observer Location  : %s", observer);
    if (onEarth) {
	fprintf (stream,
		 " (%.5f N, %.5f E, %.0f m)\n"
		 "#  Terr. Ref. Frame   : %s\n",
		 lat * dpr_c(), lon * dpr_c(), alt * 1000, frname);
    } else {
	fprintf (stream, "\n");
	/* fprintf (stream, "#  Obs. Ref. Frame   : %s\n", frname); */
    }
    SpiceDouble one_au;
    convrt_c(1.0, "AU", "KM", &one_au);
    fprintf (stream,
	     "#  Solar Radius       : %.0f km\n"
	     "#  Solar Rot. Model   : %s (%s)\n"
	     "#                       A, B, C = (%.4f, %.4f, %.4f) murad/s\n"
	     "#  One AU             : %.3f km\n",
	     RSUN, RotModels[rotmodel].name,
	     RotModels[rotmodel].descr,
	     RotModels[rotmodel].A,
	     RotModels[rotmodel].B,
	     RotModels[rotmodel].C,
	     one_au);
    
    fprintf (stream,
	     "#*****************************************************"
	     "*************************\n"
	     "#  Data Units        : km, km/s, rad, rad/s\n"
	     "#  Data Fields       : "
	     "1(utc string), 2(utc jd), 3(mjd), 4(P0), 5(L0), 6(B0),\n"
	     "#    7(rsun_obs), 8(x), "
	     "9(y), 10(lon), 11(lat), 12(rotrate), 13(rho), 14(mu),\n"
	     "#    15(dist), 16(vlos), 17(dist_sun), 18(vlos_sun)\n"
	     "#    19-24(sun inertial state), "
	     "25-30(observer intertial state)\n"
	     "#\n"
	);
}

void print_ephtable_row (FILE *stream, soleph_t *eph)
{
    SpiceChar utcstr[48];
    timout_c (eph->et, "YYYY-MM-DDTHR:MN:SC.### ::UTC", 48, utcstr);
    fprintf (stream,
    	     "%s % 16.7f % 10.6f "
    	     ,
    	     utcstr, eph->jday, eph->mjd);
    fprintf (stream,
	     "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E % 14.12E % 14.12E\n"
	     ,
	     eph->P0, eph->L0cr, eph->B0, eph->rsun_obs, eph->x, eph->y,
	     eph->lon, eph->lat, eph->omega, eph->rho, eph->mu,
	     eph->dist, eph->vlos,
	     eph->dist_sun, eph->vlos_sun,
	     eph->state_sun[0], eph->state_sun[1], eph->state_sun[2],
	     eph->state_sun[3], eph->state_sun[4], eph->state_sun[5],
	     eph->state_obs[0], eph->state_obs[1], eph->state_obs[2],
	     eph->state_obs[3], eph->state_obs[4], eph->state_obs[5]
    	);
}

void fancy_print_eph (FILE *stream, soleph_t *eph)
{
    double lon, lat, alt;
    station_geopos (eph->observer, eph->et, &lon, &lat, &alt);
    bool onEarth = true;
    if (0 != strcasecmp (eph->observer, "VTT") &&
	0 != strcasecmp (eph->observer, "SST") &&
	0 != strcasecmp (eph->observer, "SCHAUINSLAND") &&
	0 != strcasecmp (eph->observer, "DST") &&
	0 != strcasecmp (eph->observer, "MCMATH") &&
	0 != strcasecmp (eph->observer, "TEIDE") &&
	0 != strcasecmp (eph->observer, "BIGBEAR")
	) {
	onEarth = false;
    }
    SpiceInt frcode;
    SpiceBoolean found;
    SpiceChar frname[MAXKEY+1];
    if (onEarth)
	cnmfrm_c ("EARTH", MAXKEY, &frcode, frname, &found);
    else
	cnmfrm_c (eph->observer, MAXKEY, &frcode, frname, &found);

    fprintf (stream, 
	     "Date of Observation...........  %s (JD %.6f)\n",
	     eph->utcdate, eph->jday);

    if (onEarth)
	fprintf (stream,
		 "Observer Location.............  %s "
		 "(%3.5f N, %3.5f E, %.0f m)\n"
		 "Terrestr. Reference Frame.....  %s\n"
		 ,
		 eph->observer, lat * dpr_c(), lon * dpr_c(),
		 alt * 1000.0, 
		 frname);
    else
	fprintf (stream, "  Observer location...........  %s\n", eph->observer);
    
    SpiceDouble dist_au;
    convrt_c(eph->dist_sun, "KM", "AU", &dist_au);
    fprintf (stream, 
	     "Sun Reference Radius.......... % .0f m\n"
	     "Apparent Angular Radius....... % .4f arcsec\n"
	     "Rotation Model ...............  %s (%s)\n"
	     "Siderial Rotation Rate........  %.4f murad/s\n"
	     "Position Angle P.............. % -.4f deg\n"
	     "Sub-Obsrv. Stonyhurst lat..... % -.4f deg\n"
	     "Sub-Obsrv. Stonyhurst lon..... %  .4f deg\n"
	     "Sub-Obsrv. Carrington lon..... % -.4f deg\n"
	     "Solar Center Distance.........  %.0f m (%.9f AE)\n"
	     "Solar Center v_los............ % -.3f m/s\n"
	     "Helio-Projct. Cartesian x, y.. % -.4f, %.5f arcsec\n"
	     "Stonyhurst Heliogr. lon, lat.. % -.4f, %.5f deg\n"
	     "Impact parameter..............  %.0f m\n"
	     "Cos(theta) = mu............... % .4f\n"
	     "Distance......................  %.0f m\n"
	     "Line of Sight Velocity........ % -.3f m/s\n"
	     ,
	     eph->rsun_ref * 1000,
	     eph->rsun_obs * aspr() - 5E-5, /* round down to 4 digits */
	     eph->modelname, eph->modeldescr,
	     eph->omega * 1E6,
	     eph->P0 * dpr_c(),
	     eph->B0 * dpr_c(),
	     eph->L0hg * dpr_c(),
	     eph->L0cr * dpr_c(),
	     eph->dist_sun * 1000,
	     dist_au,
	     eph->vlos_sun * 1000,
	     eph->x * aspr(), eph->y * aspr(),
	     eph->lon * dpr_c(), eph->lat * dpr_c(),
	     eph->rho * 1000,
	     eph->mu,
	     eph->dist * 1000,
	     eph->vlos * 1000
	);
}


void list_rotation_models (FILE *stream)
{
    fprintf (stream, "Name        A       B       C       Description\n");
    fprintf (stream, "--------------------------------"
	    "------------------------------------------------\n");
    for (int i = 0; i < RotModel_END; ++i) {
	fprintf (stream, "%-10s % 6.4f % 6.4f % 6.4f  %s\n",
		 RotModels[i].name,
		 RotModels[i].A, RotModels[i].B, RotModels[i].C,
		 RotModels[i].descr);
    }
}

int handle_request (
    SpiceChar *observer,
    int argc,
    char **argv,
    int rotmodel,
    bool fancy,
    bool dofits,
    FILE *ostream)
{
    sunpos_t position;
    SpiceDouble et;
    size_t nsteps = 1;
    SpiceDouble stepsize = 60.0;
    fitsfile *fptr;
    int fstatus = 0;

    if (argc < 4 || argc > 6) {
	errmesg ("Insufficent input data in request\n");
	return FAILURE;
    }
    
    str2et_c (argv[0], &et);
    
    if (FAILURE == parse_sunpos (argv[1], argv[2], argv[3], &position)) {
	errmesg ("Could not parse target coordinates\n");
	return FAILURE;
    }
    if (argc == 6) {
	nsteps = atoi (argv[4]);
	stepsize = atof (argv[5]);
    }

    for (int i = 0; i < nsteps; ++i) {
	soleph_t eph;
	SpiceDouble ete = et + i * stepsize * 60.0;
	if (SUCCESS != soleph (observer, ete, position, &eph, rotmodel)) {
	    errmesg ("Could not compute ephemeris data\n");
	    return FAILURE;
	}
	if (fancy) {
	    fancy_print_eph (ostream, &eph);
	    if (i < nsteps - 1) fprintf (ostream, "\n");
	}
	else
	    print_ephtable_row (ostream, &eph);

	if (dofits) {
	    write_fits_ephtable_row (fptr, i+1, &eph, &fstatus);
	}
    }

    return SUCCESS;
}


int mode_batch (
    SpiceChar *observer,
    int rotmodel,
    bool fancy,
    char *fitsname,
    FILE *istream,
    FILE *ostream)
{
    char fitsfile[MAXPATH+1];

    /* we need a faked argv array to make handle_request happy */
    char **argv = malloc (6 * sizeof (char*));
    for (int i = 0; i < 6; ++i) {
	argv[i] = malloc (sizeof(char) * MAXPATH);
    }

    if (! fancy)
	print_ephtable_head (ostream, observer, rotmodel);

    /* bool dofits = false; */
    /* if (fitsname) { */
    /* 	dofits = true; */
    /* 	if (strcmp(fitsname, "auto") == 0) { */
    /* 	    char datestr[64+1]; */
    /* 	    timout_c (et, "YYYYMMDDTHRMNSC_UTC", 64, datestr); */
    /* 	    snprintf (fitsfile, MAXPATH, "solarv_%s_%s.fits", observer, datestr); */
    /* 	} else { */
    /* 	    strncpy (fitsfile, fitsname, MAXPATH); */
    /* 	} */
    /* } */
    /* if (dofits) { */
    /* 	fits_create_file (&fptr, fitsfile, &fstatus); */
    /* 	write_fits_ephtable_header (fptr, 1, &fstatus); */
    /* } */

    /* handle input stream line by line */
    ssize_t nread;
    char *line = 0;
    size_t len;
    size_t lineno = 1;
    while (-1 != (nread = getline(&line, &len, istream)))
    {
	char *tp = line;
	int argc = 0;
	while ( (tp = wordsep (tp, argv[argc++])) != 0) {
	    if (argc > 5) {
		errmesg ("Invalid request in input line %zu\n", lineno);
		return FAILURE;
	    }
	}
	argc--;
	
	if (SUCCESS != handle_request (observer, argc, argv,
				       rotmodel, fancy, false, ostream)) {
	    errmesg ("Invalid request in input line %zu\n", lineno);
	    return FAILURE;
	}
	lineno++;
    }

    for (int i = 0; i < 6; ++i) free (argv[i]);
    free (argv);
    
    /* if (dofits) { */
    /* 	fits_close_file (fptr, &fstatus); */
    /* } */

    return SUCCESS;
}

int write_fits_ephtable_header (fitsfile *fptr, long nrows, int *status)
{
    if (*status)
	return FAILURE;
    
    /* define table structure */
    int tfields = 28;
    char tname[] = "ephemeris table";
    char *ttype[] = { "date_obs", "jd_obs", "mjd_obs",
		      "observer", "obs_lon", "obs_lat", "obs_alt",
		      "solar_b0", "hgln_obs", "crln_obs", "solar_p0",
		      "rsun_obs", "rsun_ref",
		      "rotmodel", "modeldescr", "rotrate",
		      "state_sun", "dist_sun", "vlos_sun", 
		      "hg_lon", "hg_lat", "hpc_x", "hpc_y",
		      "mu", "impactparam", 
		      "state_obs", "dist_obs", "vlos_obs"
    };
    char *tform[] = { "32A", "D", "D",
		      "32A", "D", "D", "D",
		      "D", "D", "D", "D",
		      "D", "D",
		      "32A", "64A", "D",
		      "6D", "D", "D",
		      "D", "D", "D", "D", 
		      "D", "D",
		      "6D", "D", "D"
    };
    char *tunit[] = { "\0", "days", "days",
		      "\0", "rad", "rad", "rad",
		      "rad", "rad", "rad", "rad",
		      "rad", "km",
		      "\0", "\0", "rad/s",
		      "km, km/s", "km", "km/s",
		      "rad", "rad", "rad", "rad",
		      "\0", "km",
		      "km, km/s", "km", "km/s"
    };
    fits_create_tbl (fptr,
                     BINARY_TBL,       /* type              */
                     nrows,            /* nrows             */
                     tfields,          /* number of cols    */
                     ttype,            /* colnames          */
                     tform,            /* coltypes          */
                     tunit,            /* colunits          */
                     tname,            /* name of extension */
                     status);
    if (*status) {
        errmesg ("Error creating fitstable: %i\n", *status);
	fits_report_error (stderr, *status);
        return FAILURE;
    }
    
    return SUCCESS;
}


int write_fits_ephtable_row (
    fitsfile *fptr,
    long row,
    soleph_t *eph,
    int *status)
{
    long firstrow = row;
    long col = 1;

    SpiceDouble lon = 0, lat = 0, alt = 0; /* observer coordinates */


#define EPHTABLE_ADDCELL(type, ptr)					\
    fits_write_col (fptr, type, col++,					\
		    firstrow, 1, 1,					\
		    ptr, status)
    
    char *utcstr[] = {eph->utcdate};
    EPHTABLE_ADDCELL (TSTRING, utcstr);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->jday);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->mjd);

    char *observer[] = {eph->observer};
    EPHTABLE_ADDCELL (TSTRING, observer);
    EPHTABLE_ADDCELL (TDOUBLE, &lon);
    EPHTABLE_ADDCELL (TDOUBLE, &lat);
    EPHTABLE_ADDCELL (TDOUBLE, &alt);

    EPHTABLE_ADDCELL (TDOUBLE, &eph->B0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->L0hg);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->L0cr);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->P0);

    EPHTABLE_ADDCELL (TDOUBLE, &eph->rsun_obs);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rsun_ref);

    char *modelname[] = {eph->modelname};
    EPHTABLE_ADDCELL (TSTRING, modelname);
    char *modeldescr[] = {eph->modeldescr};
    EPHTABLE_ADDCELL (TSTRING, modeldescr);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->omega);
    
    fits_write_col (fptr, TDOUBLE, col++, row, 1, 6, eph->state_sun, status);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist_sun);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos_sun);

    EPHTABLE_ADDCELL (TDOUBLE, &eph->lon);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->lat);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->x);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->y);

    EPHTABLE_ADDCELL (TDOUBLE, &eph->mu);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rho);
    
    fits_write_col (fptr, TDOUBLE, col++, row, 1, 6, eph->state_obs, status);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos);
    
    if (*status) {
	errmesg ("couldn't add row to fitstable.");
	fits_report_error (stderr, *status);
        return FAILURE;
    }

    return SUCCESS;
}

int parse_sunpos (
    const char *type,
    const char *posx,
    const char *posy,
    sunpos_t *pos)
{
    if (strcasecmp (type, "shg") == 0)      pos->type = cord_shg;
    else if (strcasecmp (type, "hpc") == 0) pos->type = cord_hpc;
    else if (strcasecmp (type, "hcr") == 0) pos->type = cord_hcr;
    else {
	errmesg ("Unknown coordinate system: %s\n", type);
	return FAILURE;
    }

    pos->x = atof (posx);
    pos->y = atof (posy);

    return SUCCESS;
}


/* arcseconds per radian */
SpiceDouble aspr(void)
{
    return dpr_c() * 3600.0;
}

/* arcseconds per radian */
SpiceDouble rpas(void)
{
    return 1.0 / 3600.0 * rpd_c();
}

SpiceDouble dpas(void)
{
    return 1.0 / 3600.0;
}

/* helper to dump a full state vector */
void printstate (SpiceDouble *s)
{
    printf ("[% .6E, % .6E, % .6E, % .6E, % .6E, % .6E], absr=%.6E, absv=%.6E\n",
	    s[0], s[1], s[2],
	    s[3], s[4], s[5],
	    sqrt (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]),
	    sqrt (s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
    
}

/* helper to dump a vector */
void printvec (SpiceDouble *s)
{
    printf ("[% .6E, % .6E, % .6E]\n", s[0], s[1], s[2]);
}

/* word tokenizer to separate words at whitespace boundaries. substrings
   enclosed in double quotes are supported */
char* wordsep (char *str, char *token)
{
    enum {NONE, INWORD, INSTRING} state = NONE;
    char *token_start = str;
    char *w;
    size_t len;

    for (w = str; *w != '\0'; ++w)
    {
	char c = *w;

	switch (state)
	{
	case NONE:
	    if (isspace (c)) continue;

	    if (c == '"') {
		state = INSTRING;
		token_start = w + 1;
		continue;
	    }
	    state = INWORD;
	    token_start = w;
	    continue;

	case INSTRING:
	    if (c == '"') { /* closing quote */
		len = w - token_start;
		memcpy (token, token_start, len);
		token[len] = 0;
		return ++w; /* forward away from the quote */
	    }
	    continue;

	case INWORD:
	    if (isspace(c)) {
		len = w - token_start;
		memcpy (token, token_start, len);
		token[len] = 0;
		return w;
	    }
	    else if (*(w+1) == '\0') {
		len = w - token_start + 1;
		memcpy (token, token_start, len);
		token[len] = 0;
		return w;
	    }
	    
	    continue;
	}
    }

    /* we forgive a missing closing quote at the end of the sourc string */
    if (state == INSTRING) {
	len= w - token_start;
	memcpy (token, token_start, len);
	token[len] = 0;
	return w;
    }

    return NULL;
}

void errmesg (const char *mesg, ...)
{
    va_list ap;
    
    fprintf (stderr, "ERROR: ");
    va_start(ap, mesg);
    vfprintf (stderr, mesg, ap);
    va_end(ap);
}


