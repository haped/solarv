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
#include <libgen.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdarg.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <fitsio.h>

#include "solarv.h"

void usage (FILE *stream)
{
    printf ("%s v%s, %s\n"
	    "\n"
	    "compute precision radial velocities between earth-based "
	    "observatories\nand a given position on the sun.\n"
	    "\n"
	    "usage: solarv [options] <date> [<lola|xy> <lat> <lon>] [<tstep> <nsteps>]\n"
	    "options:\n"
	    "  -h            show this help\n"
	    "  -m model      set solar rotation model 'model'\n"
	    "                use 'list' to list available models\n"
	    /* "  -O observer   set observer location. can be any NAIF body code.\n" */
	    /* "                pre-defined sites: 'izana'\n" */
	    "  -p            pretty-print ephemeris data\n"
	    "  -v            print program version\n"
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
	    "  solarv -p -m su90s \"2010 Jan 1 12:00:00.0\" 90 45\n",
	    _name, _version, _versiondate
	    
	);
}

void version (void)
{
    printf ("%s v%s, %s\n"
	    "Copyright (c) 2012 %s\n",
	    _name, _version, _versiondate, _author);
}

int main (int argc, char **argv)
{   
    char fitsfile[MAXPATH+1] = "";
    char outdir[MAXPATH+1] = "";
    char time_utc[80];
    char observer[65] = "izana";
    SpiceDouble stepsize = 1; /* minutes */
    SpiceInt nsteps = 1;
    bool fancy = false;
    bool fitsmode = false;
    int rotmodel = fixed;
    int errorcode = RETURN_SUCCESS;

    sunpos_t position;

    int c; opterr = 0;
    while ((c = getopt (argc, argv, "+hm:pr:O:o:vf")) != -1)
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
		rotmodel = rigid;
	    else if (strcasecmp (optarg, "fixed") == 0)
		rotmodel = fixed;
	    else if (strcasecmp (optarg, "su90s") == 0)
		rotmodel = su90s;
	    else if (strcasecmp (optarg, "su90g") == 0)
		rotmodel = su90g;
	    else if (strcasecmp (optarg, "su90m") == 0)
		rotmodel = su90m;
	    else if (strcasecmp (optarg, "s84s") == 0)
		rotmodel = s84s;
	    else {
		fprintf (stderr, "unknown rotation model: %s\n", optarg);
		list_rotation_models (stderr);
		return EXIT_FAILURE;
	    }					
	    break;
        case 'p': fancy = true; break;
	case 'O': strncpy (observer, optarg, 64); break;
	case 'o': strncpy (outdir, optarg, MAXPATH); break;
	case 'v': version (); return EXIT_SUCCESS; break;
	case 'f': fitsmode = true; break;
        default: usage (stdout); return EXIT_FAILURE;
        }
    }

    /* commandline args */
    if (! (argc - optind == 4 || argc - optind == 6)) {
	usage (stdout);
	return EXIT_FAILURE;
    }
    if (fitsmode)
	strncpy (fitsfile, argv[optind], 1023);
    else
	strncpy (time_utc, argv[optind], 79);
    if (argc - optind == 4 || argc - optind == 6)
    {
	if (RETURN_FAILURE == parse_sunpos (argv[optind + 1], argv[optind + 2],
					    argv[optind + 3], &position))
	{
	    errmesg ("can't parse target coordinates, aborting\n");
	    return EXIT_FAILURE;
	}
    }
    if (argc - optind == 6 && ! fitsmode) {
	stepsize = atof (argv[optind + 4]);
	nsteps = atoi (argv[optind + 5]);
    }
    
    /* required kernels are coded in solarv.tm  meta kernel */
    furnsh_c ("../data/kernels/solarv.tm");
    /* TODO: poss. to load a custom kernel here */

    if (fitsmode) {
	errorcode = mode_fits (observer, fitsfile, outdir, position, rotmodel);
    }
    else {
	errorcode = mode_plain (observer, time_utc, position,
				stepsize, nsteps, fancy, stdout, rotmodel);
    }
    
    unload_c ("../data/kernels/solarv.tm");

    if (RETURN_FAILURE == errorcode) {
	errmesg ("there where errors, please check results\n");
	return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}



/*
 * this is the entry point to compute the solar ephemeris parameters for a
 * given position on the sun at a given time
 */
int soleph (
    SpiceChar *observer, /* NAIF body name/code of the observer     */ 
    SpiceDouble et,     /* Spice ephemeris time of the observation */
    sunpos_t position,
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

    reset_soleph (eph);

    /* remember julian date and ascii utc string of the event */
    SpiceDouble deltaT;
    deltet_c (et, "ET", &deltaT);
    /* we use the utc julain date here, so substract the delta */
    eph->jdate = unitim_c (et, "ET", "JDTDB") - deltaT / spd_c();
    et2utc_c (et, "C", 2, 79, eph->utcdate);
    strcpy (eph->observer, observer);

    /* compute sub-observer point on the solar surface to derive B0, and L0 */
    /* FIXME: check if and which abberation correction is needed here. LT+S
     * should be correct because NONE would compute the geometric state at
     * 'et' while what we actually observe has happened ~8 mins before */
    subpnt_c ("Near point: ellipsoid", "SUN", et, "IAU_SUN",
	      ABCORR, observer, subpoint, &trgepc, srfvec);
    reclat_c (subpoint, &subrad, &sublon, &sublat);
    eph->B0 = sublat * dpr_c();
    /* store carrington longitude as 0 ... 2 pi */
    eph->L0 = (sublon < 0 ? 2 * pi_c() + sublon : sublon) * dpr_c();

    /* compute solar barycenter state relative to observer */
    relstate_observer_sun (observer, et, eph, state_ots);

    /* compute target state relative to sun barycenter */
    relstate_sun_target (observer, et, position, rotmodel, eph, state_stt);

    /* compute observer to target by adding ots + stt */
    vaddg_c (state_ots, state_stt, 6, state_ott);
    unorm_c (state_ott, los_ott, &(eph->dist));
    eph->vlos = vdot_c (los_ott, &state_ott[3]);
    
    /************************************************************************
     compute further ephemeris data from the parameters we gathered so far
    ************************************************************************/

    /* apparent diameter of the disc in arcsec */
    eph->rsun_as = RSUN / (eph->dist_sun * 1000) * dpr_c() * 3600.0;

    /* helicentric impact paramter */
    eph->mu = sqrt (1.0 - (eph->rho_hc * 1000 / RSUN) *
		    (eph->rho_hc * 1000 / RSUN));

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

/* Thompson (2005) eq. (11) */
void lola2hcc (
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble L0,
    SpiceDouble B0,
    SpiceDouble *x,
    SpiceDouble *y,
    SpiceDouble *z)
{
    SpiceDouble costheta = cos (lat);
    SpiceDouble sintheta = sin (lat);
    SpiceDouble cosphimphi0 = cos (lon - L0);
    SpiceDouble cosB0 = cos (B0);
    SpiceDouble sinB0 = sin (B0);
    
    *x = RSUN * costheta * sin (lon - L0);
    *y = RSUN * (sintheta * cosB0 - costheta * cosphimphi0 * sinB0);
    *z = RSUN * (sintheta * sinB0 + costheta * cosphimphi0 * cosB0);
}

/* Thompson (2005) eq (12) */
void hcc2lola (
    SpiceDouble x,
    SpiceDouble y,
    SpiceDouble z,
    SpiceDouble B0,
    SpiceDouble *lon,
    SpiceDouble *lat)
{
    SpiceDouble cosB0 = cos (B0);
    SpiceDouble sinB0 = sin (B0);
    
    *lon = atan2 (x, z * cosB0 - y * sinB0);
    *lat = asin ((y * cosB0 + z * sinB0) / RSUN);
}

/* thomposn (2005) eq (15) */
void xy2hcc (
    SpiceDouble xp, /* helio-projective x */
    SpiceDouble yp, /* helio-projective y */
    SpiceDouble ds, /* dist obs-sun */
    SpiceDouble *x,
    SpiceDouble *y,
    SpiceDouble *z)
{
    SpiceDouble cosx = cos (xp);
    SpiceDouble sinx = sin (xp);
    SpiceDouble cosy = cos (yp);
    SpiceDouble siny = sin (yp);
    
    SpiceDouble q = ds * cosy * cosx;
    SpiceDouble dist = q * q - ds * ds + RSUN * RSUN;
    dist =  q - sqrt (dist);

    *x = dist * cosy * sinx;
    *y = dist * siny;
    *z = ds - dist * cosy * cosx;
}

/* Thompson (2005) eq (16) */
void hcc2xy (
    SpiceDouble x,
    SpiceDouble y,
    SpiceDouble z,
    SpiceDouble ds, /* dist obs-sun */
    SpiceDouble *xp,
    SpiceDouble *yp)
{
    SpiceDouble dt = sqrt (x*x + y*y + (ds - z) * (ds - z));
    *xp = atan2 (ds - z, x);
    *yp = asin (y / dt);
}

void lola2xy (
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble L0,
    SpiceDouble B0,
    SpiceDouble ds, /* dist obs - sun */
    SpiceDouble *xp,
    SpiceDouble *yp)
{
    SpiceDouble x, y, z;
    
    lola2hcc (lon, lat, L0, B0, &x, &y, &z);
    hcc2xy (x, y, z, ds, xp, yp);
}
 
void xy2lola (
    SpiceDouble xp,
    SpiceDouble yp,
    SpiceDouble ds,
    SpiceDouble B0,
    SpiceDouble *lon,
    SpiceDouble *lat)
{
    SpiceDouble x, y, z;
    
    xy2hcc (xp, yp, ds, &x, &y, &z);
    hcc2lola (x, y, z, B0, lon, lat);
}



/* compute relative state from solar center to target position on the
 * surface in IAU_SUN frame */
int relstate_sun_target (
    SpiceChar *station,
    SpiceDouble et,
    sunpos_t position,
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
    SpiceDouble lonrd, latrd;
    
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

    /* SpiceDouble subpnt_dist = vnorm_c (srfvec); */
    /* printf ("subobserver-distance=%f\n", subpnt_dist); */
    


    /* now that we know the heliographic sub-observer coordinates, we can
     * easily compute the helio-projective cartesian coordinates (x,y) from
     * the given position in (lon,lat) and vice versa */
    if (position.type == lola) {
	lonrd = position.x * rpd_c();
	latrd = position.y * rpd_c();
    }
    else if (position.type == xy) {
	SpiceDouble xp = position.x / 3600.0 * rpd_c();
	SpiceDouble yp = position.y / 3600.0 * rpd_c();
	xy2lola (xp, yp, eph->dist_sun * 1000, sublat, &lonrd, &latrd);
    }
    
    eph->lon = lonrd * dpr_c();
    eph->lat = latrd * dpr_c();

    
    /* rectangular, body fixed coorindates of target point given in
     * stonyhurst coordinates */
    SpiceDouble tlon = lonrd + sublon;
    SpiceDouble tlat = latrd;
    latrec_c (RSUN / 1000.0, tlon, tlat, state_stt_fixed);

    /* additional veloctiy from (differential) rotation according to
     * rotation model specified */
    SpiceDouble omegas = omega_sun (latrd, rotmodel);
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
    latrec_c (RSUN / 1000.0, lonrd, latrd - sublat, state_stt_fixed);

    /* Thompson (2005) eq 4 */
    eph->x = atan (state_stt_fixed[1] / eph->dist_sun) * dpr_c () * 3600.0;
    eph->y = atan (state_stt_fixed[2] / eph->dist_sun) * dpr_c () * 3600.0;

    /* radial distance from target to solar rotation axis */
    eph->rho_hc = sqrt (state_stt_fixed[1] * state_stt_fixed[1] + 
			state_stt_fixed[2] * state_stt_fixed[2]);

    return RETURN_SUCCESS;
}

void reset_soleph (soleph_t *eph)
{
    eph->jdate = 0.0;
    strcpy (eph->utcdate, "na");
    strcpy (eph->observer, "na");
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
    eph->rho_hc = 0.0;
    eph->omega = 0.0;
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

void print_ephtable_head (FILE *stream)
{
    fprintf (stream, "#jdate          B0(deg)  P0(deg)   L0(deg)   "
	     "vlos_t     d_t             vlos_sun  d_sun \n");
}

void print_ephtable_row (FILE *stream, soleph_t *eph)
{
    fprintf (stream, "%.6f %7.4f %9.4f %10.4f %9.4f  %14.3f  %9.4f  %13.3f\n",
	     eph->jdate, eph->B0, eph->P0, eph->L0,
	     eph->vlos * 1000, eph->dist, eph->vlos_sun * 1000, eph->dist_sun);
}

void fancy_print_eph (FILE *stream, soleph_t *eph)
{
	    printf ("Solar ephemeris for %s, %s\n"
		    "  UTC julian date      :% f\n"
		    "  disk radius          : %.2f arcsec\n"
		    "  B0                   :% -.4f deg\n"
		    //"  P0                   :% .4f deg\n"
		    "  carrington L0        :% .4f deg\n"
		    "  center distance      : %.3f km\n"
		    "  center v_los         :% .3f m/s\n"
		    "  disc coordinates     :% .2f, %.2f arcsec\n"
		    "  lola coordinates     :% .3f, %.3f deg\n"
		    "  cos(hel. angle) = mu : %.4f\n"
		    "  target distance      : %.3f km\n"
		    "  target v_los         :% .3f m/s\n"
		    "  solar rotation model : %s (%s)\n"
		    "  rotataion rate       : %.5f murad/s\n"
		    "  impact parameter     : %.3f km\n",
		    
		    eph->observer, eph->utcdate,
		    eph->jdate,
		    eph->rsun_as,
		    eph->B0,
		    //eph->P0,
		    eph->L0,
		    eph->dist_sun,
		    eph->vlos_sun * 1000,
		    eph->x, eph->y,
		    eph->lon, eph->lat,
		    eph->mu,
		    eph->dist,
		    eph->vlos * 1000,
		    eph->modelname, eph->modeldescr,
		    eph->omega,
		    eph->rho_hc);
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

int mode_plain (
    SpiceChar *observer, /* NAIF body name/code of the observer     */ 
    char *time_utc,      /* Spice ephemeris time of the observation */
    sunpos_t position,
    SpiceDouble stepsize,
    int nsteps,
    bool fancy,
    FILE * stream,
    int rotmodel)
{
    SpiceDouble et;
    str2et_c (time_utc, &et);

    if (! fancy)
	print_ephtable_head (stdout);
    for (SpiceInt i  = 0; i < nsteps; ++i)
    {
	soleph_t eph;
	soleph (observer, et, position, &eph, rotmodel);

	if (fancy) {
	    fancy_print_eph (stream, &eph);
	    if (i < nsteps - 1) fprintf (stream, "\n");
	}
	else
	    print_ephtable_row (stdout, &eph);

	et += stepsize * 60.0;
    }

    return RETURN_SUCCESS;
}

int mode_fits (
    SpiceChar *observer,
    char *infile,
    char *outdir,
    sunpos_t position,
    int rotmodel)
{
    char fitsfile2[MAXPATH+1];
    char fitsfile3[MAXPATH+1];
    char outfile[MAXPATH+1];
    char *basenam = '\0', *basedir = '\0';
    fitsfile *fptr, *fptrout;
    int status = 0;
    long naxes[3];
    int nfound;
    long nframes = 0;

    /* write the ephemeris table in the same directory and basename as the
     * input file but append '_ephem.fits' */
    strncpy (fitsfile2, infile, MAXPATH); 
    strncpy (fitsfile3, infile, MAXPATH); 
    basenam = basename (fitsfile2);
    char *lastp = rindex (basenam, '.'); /* remove file extension */
    if (lastp) *lastp = '\0';

    if (strcmp (outdir, "") == 0)
	basedir = dirname (fitsfile3);
    else
	basedir = outdir;
    if (strcmp (basedir, ".") == 0)
	snprintf (outfile, MAXPATH, "%s_ephem.fits", basenam);
    else
	snprintf (outfile, MAXPATH, "%s/%s_ephem.fits", basedir, basenam);

    printf ("reading  %s: ", infile);

    fits_open_file (&fptr, infile, READONLY, &status);

    /* check file compatibility */
    /* read the NAXIS1, NAXIS2 and NAXIS2 keyword to get image size */
    fits_read_keys_lng (fptr, "NAXIS", 1, 3, naxes, &nfound, &status);
    if (status)
    {
        errmesg ("error reading fits header from %s: %s\n", infile);
	fits_report_error (stderr, status);
        fits_close_file (fptr, &status);
        return RETURN_FAILURE;
    }
    if (nfound == 3) {
	printf ("%li frames\n", naxes[2]);
	nframes = naxes[2];
    }
    else if (nfound == 2) {
	printf ("1 frame\n");
	nframes = 1;
    }
    else {
	errmesg ("dimension=%i unsupported\n", nfound);
	return RETURN_FAILURE;
    }

    /* check if there should be a BCD header in the file; read DATE-OBS
     * keyword as a backup */
    bool haveBCDdate = false;
    bool haveExptime = false;
    char h_tsmode[81] = "na";
    fits_read_key (fptr, TSTRING, "TSMODE", h_tsmode, 0, &status);
    if (status == 202) status = 0; /* ignore unavailable keyword */
    if (strstr (h_tsmode, "Binary")) haveBCDdate = true;
    double h_exptime = -1;
    fits_read_key (fptr, TDOUBLE, "EXPTIME", &h_exptime, 0, &status);
    if (h_exptime != -1) haveExptime = true;
    if (status == 202) status = 0; /* ignore unavailable keyword */

    char h_dateobs[81] = "na";
    double h_deltime;
    if (! haveBCDdate)
    {
	fits_read_key (fptr, TSTRING, "DATE-OBS", h_dateobs, 0, &status);
	fits_read_key (fptr, TDOUBLE, "DELTIME", &h_deltime, 0, &status);
	if (status) {
	    errmesg ("no time of exposure information found in fits file, aborting\n");
	    fits_report_error (stderr, status);
	    return RETURN_FAILURE;
	}
	printf ("Warning: no binary timestamps available, using EXPTIME"
		" + DELTIME instead. This might be unreliable.\n");
    }
    
    fits_create_file (&fptrout, outfile, &status);
    fits_create_img (fptrout, BYTE_IMG, 0, 0, &status);
    fits_write_date (fptrout,  &status);
    write_fits_ephtable_header (fptrout, naxes[2], &status);
    if (status)  {
        errmesg ("couldn't write fitstable header", status);
	fits_report_error (stderr, status);
        fits_close_file (fptrout, &status);
        return RETURN_FAILURE;
    }
    
    /* in case the file does not provide the BCD timestamps, we have to
      extrapolate the time of exosure via T = DATE-OBS + i * (EXPTIME + DELTIME) */
    SpiceDouble et;
    if (! haveBCDdate)
	str2et_c (h_dateobs, &et);
    
    for (long i = 0; i < nframes; ++i)
    {
	soleph_t eph;

	if (haveBCDdate) {
	    char utcstr[64];
	    fitsframe_bcddate (fptr, i+1, naxes[1], utcstr, &status);
	    if (status) break;
	    str2et_c (utcstr, &et);
	}
	else
	    et += i * (h_exptime / 1000.0 + h_deltime / 1000.0);
	
	soleph (observer, et, position, &eph, rotmodel);
	printf ("   frame %ld: %s  exptime: %.2fs  vlos: %.2f m/s\n",
		i, eph.utcdate, h_exptime / 1000.0, eph.vlos * 1000); 
	
	write_fits_ephtable_row (fptrout, i+1, &eph, &status);
	if (status) break;
    }

    printf ("writing %s\n", outfile);
    fits_close_file (fptr, &status);
    fits_close_file (fptrout, &status);
    if (status) {
        errmesg ("couldn't write fitstable.");
	fits_report_error (stderr, status);
        return RETURN_FAILURE;
    }

    return RETURN_SUCCESS;
}


int write_fits_ephtable_header (fitsfile *fptr, long nrows, int *status)
{
    if (*status)
	return RETURN_FAILURE;
    
    /* define table structure */
    int tfields = 16;
    char tname[] = "ephemeris table";
    char *ttype[] = { "jdate",
		      //"utcdate", "observer",
		      "B0", "L0", "P0",
		      "dist_sun", "vlos_sun", "rsun_as",
		      //"modelname", "modeldescr",
		      "lon", "lat", "x", "y", "mu", "dist", "vlos",
		      "rho_hc", "omega" };
    char *tform[] = { "D",
		      //"32A", "32A",
		      "D", "D", "D",
		      "D", "D", "D",
		      //"32A", "32A",
		      "D", "D", "D", "D", "D", "D", "D",
		      "D", "D" };
    char *tunit[] = { "sec",
		      //'\0', '\0',
		      "deg", "deg", "deg",
		      "km", "km/s", "arcsec",
		      //'\0', '\0',
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
        return RETURN_FAILURE;
    }
    
    return RETURN_SUCCESS;
}


int write_fits_ephtable_row (
    fitsfile *fptr,
    long row,
    soleph_t *eph,
    int *status)
{
    long firstrow = row;
    long col = 1;

#define EPHTABLE_ADDCELL(type, ptr) fits_write_col (fptr, type, col++,	\
						    firstrow, 1, 1,	\
						    (void*)ptr, status)
    EPHTABLE_ADDCELL (TDOUBLE, &eph->jdate);
    /* EPHTABLE_ADDCELL (TSTRING, eph->utcdate); */
    /* EPHTABLE_ADDCELL (TSTRING, eph->observer); */
    EPHTABLE_ADDCELL (TDOUBLE, &eph->B0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->L0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->P0);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist_sun);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos_sun);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rsun_as);
    /* EPHTABLE_ADDCELL (TSTRING, eph->modelname); */
    /* EPHTABLE_ADDCELL (TSTRING, eph->modeldescr); */
    EPHTABLE_ADDCELL (TDOUBLE, &eph->lon);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->lat);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->x);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->y);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->mu);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->dist);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->vlos);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->rho_hc);
    EPHTABLE_ADDCELL (TDOUBLE, &eph->omega);
    
    if (*status) {
	errmesg ("couldn't add row to fitstable.");
	fits_report_error (stderr, *status);
        return RETURN_FAILURE;
    }

    return RETURN_SUCCESS;
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
    else {
	errmesg ("bad coordinate type: %s\n", type);
	return RETURN_FAILURE;
    }

    pos->x = atof (posx);
    pos->y = atof (posy);

    return RETURN_SUCCESS;
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
    
    if (*status) return RETURN_FAILURE;
    
    fits_read_pix (fptr, TSHORT, fpixel, bufsize/2,
		   nulval, (void*) buf, &anynul, status);
    if (*status) return RETURN_FAILURE;
    
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

    if (! ((1950 < year && 2100 > year) &&
	   (mon > 0 && mon <= 12)       &&
	   (day > 0 && day <= 31)       &&
	   (hour >= 0 && hour <= 24)    &&
	   (min >= 0 && min <= 59)      &&
	   (sec >= 0 && sec <= 59))) {
	errmesg ("timestamp of frame %li is invalid\n", frameidx);
	return RETURN_FAILURE;
    }
	
    snprintf (utcstr, 63, "%4i-%02i-%02iT%02i:%02i:%02i.%i",
	      year, mon, day, hour, min, sec, msec);

    return RETURN_SUCCESS;
}
