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
    printf ("%s v%s (%s) -- Compute precise Sun-Observer radial velocity\n"
	    "\n"
	    "Usage:\n"
	    "  solarv [options] [request | file]\n"
	    "\n"
            "Request specification:\n"
            "  <timespec> <xy|lola|muphi> <cord1> <cord2> [<nsteps> <tstep>]\n"
            "\n"
            "The request is either read from the commandline, from 'file' "
	    "or from STDIN if no arguments are given. The position can be "
	    "specified in helio-projective cartesian coordinates as "
	    "arcseconds from disk center (xy), Stonyhurst heliographic "
	    "coordinates in degree longitude and latitude (lola) or in a "
	    "polar coordinate system with the heliocentric parameter mu "
	    "and a position angle phi, counter-clockwise from solar north "
	    "(muphi). The 'timespec' parameter understands the common time "
	    "strings like '2010-01-01 12:00:00 EST'. If no time zone is "
	    "given, UTC is assumed. I highly recommend to use the ISO 8601 "
	    "standard, e.g. '2012-01-01T00:00:00'. The optional 'nsteps' "
	    "and 'tstep' (in minutes) arguments can be used to compute a "
	    "set of values in a given time interval."
	    "\n\n"
	    "Options:\n"
	    "  -h            show this help\n"
	    "  -m model      set solar rotation model 'model' [fixed]\n"
	    "                use 'list' to list available models\n"
	    "  -t            use the ITRF93 precision earth rotation kernel\n"
	    "                Note: this kernel has a limited time coverage and needs\n"
	    "                to be updated regularily\n"
	    "  -O observer   set observer position. Can be any NAIF body code.\n"
	    "                pre-defined sites: 'VTT', 'SCHAUINSLAND', 'SST',\n"
	    "                'DST', 'MCMATH'\n"
	    "  -p            pretty-print ephemeris data\n"
	    "  -v            print program version\n"
	    "  -k kernel     load additional SPICE kernel 'kernel'\n"
	    "                this kernel will be loaded last in kernel pool\n"
	    "  -K kernel     load 'kernel' instead of the default meta kernel\n"
	    "  -R radius     specifiy a different solar radius in meters\n"
	    "\n"
	    "Examples:\n"
	    "\n"
	    "Compute radial velocity at disk center, assuming a "
	    "carrington rotation model, pretty-print output:\n"
	    "  solarv -m carrington -p 2012-01-01T12:00:00 xy 0 0"
	    "\n\n"
	    "Compute radial velocity at the western limb using Snodgrass "
	    "& Ulrich (1990) spectroscopy rotation model for one year "
	    "with one sample per day:\n"
	    "  solarv -m su90s 2012-01-01T12:00:00.0 muphi 0.0 90 365 1440\n"
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
    SpiceDouble relstate[6], relstatei[6];
    
    /* state of target relative to solar center in j2000; sublon is the
     * stonyhurst longitude. by usng the HEEQ frame, we do not need the
     * sub-observer latitude, because HEEQ has is x-axis in the sub-
     * observer point*/
    latrec_c (RSUN, lon, lat, relstate);
    SpiceDouble xform[6][6];
    sxform_c ("HEEQ", "J2000", et, xform);
    mxvg_c (xform, relstate, 6, 6, relstatei);
    
    /* compute state velocity from given omega */
    SpiceDouble radius[] = {relstate[0], relstate[1], 0.0};
    SpiceDouble omegav[] = {0.0, 0.0, omegas};
    vcrss_c (omegav, radius, &relstatei[3]);
    vaddg_c (sunstate, relstatei, 6, tgtstate);

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
    char addkernel[MAXPATH + 1] = "na";
    char metakernel[MAXPATH + 1] = "na";
    char observer[MAXKEY + 1] = "VTT";
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
		printf ("available rotation models:\n");
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
		fprintf (stderr, "unknown rotation model: %s\n", optarg);
		list_rotation_models (stderr);
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
	    errmesg ("Can not open input file %s\n", argv[optind+posargs]);
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
    
    /* required kernels are coded in solarv.tm  meta kernel */
    if (strcmp (metakernel, "na") == 0) {
	strncpy (metakernel , KERNEL_PATH "/" METAKERNEL, MAXPATH);
	/* we only handle the -t flag, if we load the default meta kernel */
	if (earth_itrf93) {
	    furnsh_c (KERNEL_PATH "/earth_assoc_itrf93.tf");
	    furnsh_c (KERNEL_PATH "/earth_itrf.tf");
	} else {
	    furnsh_c (KERNEL_PATH "/earth_fixed.tf");
	}
    }
    furnsh_c (metakernel);
    if (strcmp (addkernel, "na") != 0) {
	printf ("Loading user-supplied additional kernel %s\n", addkernel);
	furnsh_c (addkernel);
    }

    if (dumpinfo) {
	dump_kernel_info (stdout);
    }

    /* the real work is done here */
    if (batchmode) {
	errorcode = mode_batch (observer, rotmodel,
				fancy,  batchstream, stdout);
    } else {
	if (! fancy)
	    print_ephtable_head (stdout, observer, rotmodel);
	errorcode = handle_request (observer, posargs, &argv[optind],
				    rotmodel, fancy, stdout);
    }
    
    unload_c (KERNEL_PATH "/" METAKERNEL);
    if (strcmp (addkernel, "na") != 0)
	unload_c (addkernel);
	
    if (FAILURE == errorcode) {
	fprintf (stderr, "Aborting, please check input!\n");
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
    if (pos.type == lola)
    {
	*lon = pos.x * rpd_c();
	*lat = pos.y * rpd_c();
    }
    else if (pos.type == xy)
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
    else if (pos.type == muphi)
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
    pcpool_c ( "FRAME_1803430_PRI_TARGET", 1, strlen (observer)-1,  observer);

    /* remember julian date and ascii utc string of the event */
    SpiceDouble deltaT;
    deltet_c (et, "ET", &deltaT);
    eph->et = et;
    /* we use the utc julain date here, so substract the delta */
    eph->jday = unitim_c (et, "ET", "JDTDB") - deltaT / spd_c();
    eph->mjd = eph->jday - 2400000.5; 
    et2utc_c (et, "C", 2, MAXKEY, eph->utcdate);

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
    eph->omega = omegas * 1E6;
    eph->rotmodel = rotmodel;
    strncpy (eph->modelname, RotModels[rotmodel].name, MAXKEY);
    strncpy (eph->modeldescr, RotModels[rotmodel].descr, MAXKEY);

    /* TODO: handle light-time correction */
    getstate_solar_target (et, eph->lon, eph->lat,
			   omegas, state_sun, state_tgt);
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
    SpiceDouble xform[3][3];
    SpiceDouble e[3] = {0, 0, 1}, sol[3];
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
    eph->P0 = -999999; /* we don't compute p0 yet */
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
	     "#  Ephemeris data created on %s UTC by solaRV v%s\n"
	     "#  SPICE version     : %s\n"
	     "#  Reference Epoch   : J2000.0\n"
	     "#  Reference Frame   : ICRF/J2000, "
	     "Mean Equator and Equinox of Epoch\n"
	     "#  Abbr. correction  : %s\n"
	     , now, _version, tkvrsn_c ("toolkit"), ABCORR);
    
    double lon, lat, alt;
    station_geopos (observer, 1.0, &lon, &lat, &alt);
    
    /* FIXME: really check if we are on earth or not */
    bool onEarth = true;
    if (0 != strcasecmp (observer, "VTT") &&
	0 != strcasecmp (observer, "SST") &&
	0 != strcasecmp (observer, "SCHAUINSLAND") &&
	0 != strcasecmp (observer, "DST") &&
	0 != strcasecmp (observer, "MCMATH")) {
	onEarth = false;
    }

    SpiceInt frcode;
    SpiceBoolean found;
    SpiceChar frname[MAXKEY+1];
    if (onEarth)
	cnmfrm_c ("EARTH", MAXKEY, &frcode, frname, &found);
    else
	cnmfrm_c (observer, MAXKEY, &frcode, frname, &found);
    
    fprintf (stream,
	     "#*****************************************************"
	     "*************************\n"
	     "#  Observer location : %s", observer);
    if (onEarth) {
	fprintf (stream,
		 " (%.5f N, %.5f E, %.0f m)\n"
		 "#  Terr. Ref. Frame  : %s\n",
		 lat * dpr_c(), lon * dpr_c(), alt * 1000, frname);
    } else {
	fprintf (stream, "\n");
	/* fprintf (stream, "#  Obs. Ref. Frame   : %s\n", frname); */
    }
    fprintf (stream,
	     "#  Solar radius      : %.0f km\n"
	     "#  Solar rot. Model  : %s (%s)\n",
	     RSUN, RotModels[rotmodel].name,
	     RotModels[rotmodel].descr);
    
    fprintf (stream,
	     "#*****************************************************"
	     "*************************\n"
	     "#  Data units        : km, km/s, rad\n"
	     "#  Data fields       : "
	     "1(date), 2(jd), 3(mjd), 4(P0), 5(L0), 6(B0),\n"
	     "#    7(rsun_obs), 8(x), "
	     "9(y), 10(lon), 11(lat), 12(rho), 13(mu), 14(dist),\n"
	     "#    15(vlos), 16(dist_sun), 17(vlos_sun)\n"
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
    fprintf (stream, "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E % 14.12E % 14.12E % 14.12E % 14.12E % 14.12E "
	     "% 14.12E %14.12E\n"
	     ,
	     eph->P0, eph->L0cr, eph->B0, eph->rsun_obs, eph->x, eph->y,
	     eph->lon, eph->lat, eph->rho, eph->mu, eph->dist, eph->vlos,
	     eph->dist_sun, eph->vlos_sun
    	);
}

void fancy_print_eph (FILE *stream, soleph_t *eph)
{
    double lon, lat, alt;

    station_geopos (eph->observer, eph->et, &lon, &lat, &alt);
    printf ("Solar ephemeris for %s\n"
	    "  Observer location..........  %s (%3.5f N, %3.5f E, %.0f m)\n"
	    "  TDB ephemeris time.........  %.3f\n"
	    "  UTC julian day.............  %f\n"
	    "  Modified julian day........  %f\n"
	    "  Abberation correction......  %s\n"
	    "  Sun Reference Radius....... % .0f m\n"
	    "  Angular Disk radius........ % .4f arcsec\n"
	    "  Solar Position angle P..... % -.4f deg\n"
	    "  Sub-obsrv. Latitude (B0)... % -.4f deg\n"
	    "  Sub-obsrv. Stonyhurst lon.. %  .9f deg\n"
	    "  Sub-obsrv. Carrington lon.. % -.4f deg\n"
	    "  Solar Barycenter distance..  %.0f m\n"
	    "  Solar Barycenter v_los..... % -.3f m/s\n"
	    "  Target Disk coordinates.... % -.5f, %.5f arcsec\n"
	    "  Target Lola coordinates.... % -.5f, %.5f deg\n"
	    "  Target Impact parameter....  %.0f m\n"
	    "  Target mu..................  %.6f\n"
	    "  Target distance............  %.0f m\n"
	    "  Target v_los............... % -.3f m/s\n"
	    "  Solar Rotation Model ......  %s (%s)\n"
	    "  Siderial rotation rate.....  %.5f murad/s\n"
	    ,
	    eph->utcdate,
	    eph->observer,
	    lat * dpr_c(), lon * dpr_c(), alt * 1000.0,
	    eph->et,
	    eph->jday,
	    eph->mjd,
	    ABCORR,
	    eph->rsun_ref * 1000,
	    eph->rsun_obs * aspr(), /* adopt this to displayed number of digits */
	    eph->P0 * dpr_c(),
	    eph->B0 * dpr_c(),
	    eph->L0hg * dpr_c(),
	    eph->L0cr * dpr_c(),
	    eph->dist_sun * 1000,
	    eph->vlos_sun * 1000,
	    eph->x * aspr(), eph->y * aspr(),
	    eph->lon * dpr_c(), eph->lat * dpr_c(),
	    eph->rho * 1000,
	    eph->mu,
	    eph->dist * 1000,
	    eph->vlos * 1000,
	    eph->modelname, eph->modeldescr,
	    eph->omega
	);
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

int handle_request (
    SpiceChar *observer,
    int argc,
    char **argv,
    int rotmodel,
    bool fancy,
    FILE *ostream)
{
    sunpos_t position;
    SpiceDouble et;
    size_t nsteps = 1;
    SpiceDouble stepsize = 60.0;

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
	if (SUCCESS != soleph (observer, et, position, &eph, rotmodel)) {
	    errmesg ("Could not compute ephemeris data\n");
	    return FAILURE;
	}
	if (fancy) {
	    fancy_print_eph (ostream, &eph);
	    if (i < nsteps - 1) fprintf (ostream, "\n");
	}
	else
	    print_ephtable_row (ostream, &eph);
	et += stepsize * 60.0;
    }

    return SUCCESS;
}


int mode_batch (
    SpiceChar *observer,
    int rotmodel,
    bool fancy,
    FILE *istream,
    FILE *ostream)
{
    /* we need a faked argv array to make handle_request happy */
    char **argv = malloc (6 * sizeof (char*));
    for (int i = 0; i < 6; ++i) {
	argv[i] = malloc (sizeof(char) * MAXPATH);
    }

    if (! fancy)
	print_ephtable_head (ostream, observer, rotmodel);

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
				       rotmodel, fancy, ostream)) {
	    errmesg ("Invalid request in input line %zu\n", lineno);
	    return FAILURE;
	}
	lineno++;
    }

    for (int i = 0; i < 6; ++i) free (argv[i]);
    free (argv);
    
    return SUCCESS;
}

int write_fits_ephtable_header (fitsfile *fptr, long nrows, int *status)
{
    if (*status)
	return FAILURE;
    
    /* define table structure */
    int tfields = 23;
    char tname[] = "ephemeris table";
    char *ttype[] = { "jd", "mjd",
		      "utc", "observer",
		      "B0", "L0hg", "L0cr", "P0",
		      "dist_sun", "vlos_sun", "rsun_as", "rsun_ref",
		      "modelname", "modeldescr",
		      "lon", "lat", "x", "y", "mu", "dist", "vlos",
		      "rho", "omega" };
    char *tform[] = { "D", "D",
		      "32A", "32A",
		      "D", "D", "D", "D",
		      "D", "D", "D", "D",
		      "32A", "32A",
		      "D", "D", "D", "D", "D", "D", "D",
		      "D", "D" };
    char *tunit[] = { "sec", "sec",
		      "\0", "\0",
		      "deg", "deg", "deg", "deg", "deg",
		      "km", "km/s", "arcsec", "km",
		      '\0', '\0',
		      "deg", "deg", "arcsec", "arcsec", '\0', "km", "km/s",
		      "km", "murad/s"};

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
        errmesg ("error creating fitstable: %i\n", *status);
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

#define EPHTABLE_ADDCELL(type, ptr)					\
    fits_write_col (fptr, type, col++,					\
		    firstrow, 1, 1,					\
		    (void*)ptr, status)
    

    //printf ("strlen (%s)=%zi\n", eph->observer, strlen (eph->observer));

    EPHTABLE_ADDCELL (TDOUBLE, &eph->jday);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->mjd);
    EPHTABLE_ADDCELL (TBYTE, eph->utcdate);
    EPHTABLE_ADDCELL (TSBYTE, eph->observer);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->B0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->L0hg);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->L0cr);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->P0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist_sun);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos_sun);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rsun_obs);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rsun_ref);
    EPHTABLE_ADDCELL (TBYTE, eph->modelname);
    EPHTABLE_ADDCELL (TBYTE, eph->modeldescr);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->lon);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->lat);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->x);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->y);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->mu);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rho);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->omega);
    
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
    if (strcasecmp (type, "lola") == 0)
	pos->type = lola;
    else if (strcasecmp (type, "xy") == 0)
	pos->type = xy;
    else if (strcasecmp (type, "muphi") == 0)
	pos->type = muphi;
    else {
	errmesg ("Unknown coordinate system: %s\n", type);
	return FAILURE;
    }

    pos->x = atof (posx);
    pos->y = atof (posy);

    return SUCCESS;
}


void errmesg (const char *mesg, ...)
{
    va_list ap;
    
    fprintf (stderr, "ERROR: ");
    va_start(ap, mesg);
    vfprintf (stderr, mesg, ap);
    va_end(ap);
}

int fitsframe_bcddate (
    fitsfile *fptr,
    long frameidx,
    long sy,
    char *utcstr,
    int *status)
{
    /* we assume 16 bit fits files here */
    const long bufsize = 32;
    unsigned short buf[bufsize];
    void *nulval = NULL;
    int anynul;
    long fpixel[] = {1, sy, frameidx};
    
    if (*status) return FAILURE;
    
    fits_read_pix (fptr, TSHORT, fpixel, bufsize/2,
		   nulval, (void*) buf, &anynul, status);
    if (*status) return FAILURE;
    
    /* bcd decoder borrowed from Kolja Klogowski */
    /*
      Format:
      buf[ 4] -> bcd[0]  ==  year * 100
      buf[ 5] -> bcd[1]  ==  year
      buf[ 6] -> bcd[2]  ==  month
      buf[ 7] -> bcd[3]  ==  day
      buf[ 8] -> bcd[4]  ==  hour
      buf[ 9] -> bcd[5]  ==  minute
      buf[10] -> bcd[6]  ==  second
      buf[11] -> bcd[7]  ==  mus * 10000
      buf[12] -> bcd[8]  ==  mus * 100
      buf[13] -> bcd[9]  ==  mus
    */
    unsigned char bcd[10];
    for (int i = 0; i < 10; ++i)
        bcd[i] = (unsigned char) (buf[4 + i] & 0xff);

    unsigned char val[10];
    for (int i = 0; i < 10; ++i)
        val[i] = (bcd[i] >> 4) * 10 + (bcd[i] & 0x0f);

    int year = 100 * (int)val[0] + (int)val[1];
    int mon = val[2];
    int day = val[3];
    int hour = val[4];
    int min = val[5];
    int sec = val[6];
    int msec = (int) (bcd[7]);

    /* some rather relaxed consistency checks. the used kernel probably
     * won't allow for larger time spans anyway */
    if (! ((1950 < year && 2100 > year) &&
	   (mon > 0 && mon <= 12)       &&
	   (day > 0 && day <= 31)       &&
	   (hour >= 0 && hour <= 24)    &&
	   (min >= 0 && min <= 59)      &&
	   (sec >= 0 && sec <= 59))) {
	errmesg ("timestamp of frame %li is invalid\n", frameidx);
	return FAILURE;
    }
	
    snprintf (utcstr, 63, "%4i-%02i-%02iT%02i:%02i:%02i.%i",
	      year, mon, day, hour, min, sec, msec);

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
    printf ("x=(%g, %g, %g), v=(%g, %g, %g), |x|=%g\n",
	    s[0], s[1], s[2],
	    s[3], s[4], s[5],
	    sqrt (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]));
}

/* helper to dump a vector */
void printvec (SpiceDouble *s)
{
    printf ("x=[%.9g, %.9g, %.9g]\n", s[0], s[1], s[2]);
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
