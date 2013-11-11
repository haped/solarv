/*
 Copyright (C) 2012, 2013 Hans-Peter Doerr <doerr@kis.uni-freiburg.de>,
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
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
*/

#ifndef _SOLARV_H_
#define _SOLARV_H_

#include <cfitsio/fitsio.h>

static const char _name[] = "Solarv";
static const char _author[] = "Hans-Peter Doerr";
static const char _version[] = "0.4.0a";
static const char _versiondate[] = "07 Nov 2013";
static const char _copyright[] = "2012, 2013 Hans-Peter Doerr";

#include "SpiceUsr.h"

#define METAKERNEL "solarv.tm"

#ifndef MAXPATH
#define MAXPATH (2048)
#endif

#ifdef MAXKEY
#undef MAXKEY
#endif
#define MAXKEY (64)

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif

#ifndef SUCCESS
#define SUCCESS (1)
#endif

#ifndef FAILURE
#define FAILURE (0)
#endif

/* solar photospheric radius in km from Brown and Dalsgaard (1998), APJ */
SpiceDouble RSUN = 6.95508E5;
SpiceChar ABCORR[8] = "None";

typedef struct
{
    double A;
    double B;
    double C;
    char name[33];
    char descr[129];
} rotmodel_t;

static const rotmodel_t RotModels[] =
{
    {0.000,    0.000,  0.000, "fixed", "no rotation, fixed to inertial frame"},
    {2.851,    0.000,  0.000, "rigid", "rigid body rotation, 2.851 murad/s"},
    {2.86532,  0.000,  0.000, "crgt", "rigid body rotation, Carrington rate"},
    {2.851,   -0.343, -0.474, "su90s", "Snodgrass & Ulrich (1990), spectroscopic"}, 
    {2.972,   -0.484, -0.361, "su90g", "Snodgrass & Ulrich (1990), supergranul."},
    {2.879,   -0.339, -0.485, "su90m", "Snodgrass & Ulrich (1990), magnetic"},
    {2.836,   -0.344, -0.504, "s84s", "Snodgrass (1984), spectrosc. MWO data"},
    {0.000,    0.000,  0.000, "custom", "Custom selected A, B, C coefficients"}
};
enum RotModel {fixed = 0,
	       rigid,
	       crgt,
	       su90s,
	       su90g,
	       su90m,
	       s84s,
	       RotModel_END,
	       custom /* this one can not be selcted */
};

typedef struct
{
    SpiceDouble x;
    SpiceDouble y;
    int type;
} sunpos_t;

enum PosType {
    cord_shg = 0,  /* Stonyhurst heliographic     */
    cord_hpc,      /* helio-projective cartesian  */
    cord_hcr       /* helio-centric radial        */
};


/* 
   ephemeris data structure. angles in radian, distances in km, velocities
   in km/s
*/
typedef struct
{
    /* common data; sun global parameters */
    SpiceDouble et;               /* SPICE ephemeris time                    */
    SpiceDouble jday;             /* julian day of the event                 */
    SpiceDouble mjd;              /* modified julian day                     */
    SpiceChar utcdate[MAXKEY+1];  /* ascii date in UTC                       */
    SpiceChar observer[MAXKEY+1]; /* NAIF station name                       */
    bool obs_on_earth;            /* true if observer is on earth            */
    SpiceDouble B0;               /* lat of sub-observer point               */
    SpiceDouble L0cr;             /* carrington lon of sub-observer point    */
    SpiceDouble L0hg;             /* stonyhurst lon of sub-observer point    */
    SpiceDouble P0;               /* Position angle of solar north           */
    SpiceDouble azimuth;          /* azimuth obs target as seen by observer  */
    SpiceDouble elev_app;
    SpiceDouble elev_true;
    SpiceDouble dist_sun;         /* distance obs. to solar center           */
    SpiceDouble vlos_sun;         /* radial velocity of obs to solar center  */
    SpiceDouble rsun_ref;         /* reference radius of the sun             */
    SpiceDouble rsun_obs;         /* apparent angular radius of the disk     */
    SpiceDouble gr_sun;           /* solar grav redshift at observer dist.   */
    SpiceDouble gr_obs;           /* grav redshift observer - sun            */
    SpiceDouble grs;
    int rotmodel;                 /* solar rotation model used               */
    char modelname[MAXKEY+1];     /* name of the rotation model              */
    char modeldescr[MAXKEY+1];    /* description of the rotation model       */
    								         
    /* target position parameters */				         
    SpiceDouble lon;              /* stonyhurst target longitude             */
    SpiceDouble lat;              /* stonyhurst target latitude              */
    SpiceDouble x;                /* helio-projective cartesian x            */
    SpiceDouble y;                /* helio-projective cartesian y            */
    SpiceDouble theta;            /* heliocentric angle                      */
    SpiceDouble mu;               /* heliocentric parameter of the target    */
    SpiceDouble dist;             /* distance to target                      */
    SpiceDouble vlos;             /* radial velocity of target               */
    SpiceDouble rho;              /* heliocentric impact parameter           */
    SpiceDouble omega;            /* angular velocity at target latitude     */
    SpiceDouble airmass;          /* observer target airmass index           */

    /* observer, target state vectors */
    SpiceDouble state_obs[6];
    SpiceDouble state_sun[6];
    SpiceDouble state_rel[6];
} soleph_t;

void station_geopos (
    SpiceChar * station,
    SpiceDouble et,
    SpiceDouble *lon,
    SpiceDouble *lat,
    SpiceDouble *alt);

void usage (FILE *stream);
void errmesg (const char *mesg, ...);

void dump_kernel_info (FILE *stream);
void list_rotation_models (FILE *stream);

char* wordsep (char *str, char *token);
int mode_batch (
    SpiceChar *observer,
    int rotmodel,
    bool fancy,
    char *fitsname,
    FILE *istream,
    FILE *ostream);

int parse_sunpos (
    const char *type,
    const char *posx,
    const char *posy,
    sunpos_t *pos);

void reset_soleph (soleph_t *eph);

int soleph (
    SpiceChar *station, /* NAIF body name/code of the observer     */ 
    SpiceDouble et,     /* Spice ephemeris time of the observation */
    sunpos_t position, 
    soleph_t *eph,
    int rotModel);

SpiceDouble omega_sun (SpiceDouble lat, int model);

int handle_request (
    SpiceChar *observer,
    int argc,
    char **argv,
    int rotmodel,
    bool fancy,
    fitsfile * fptr,
    FILE *ostream,
    size_t reqnum);

void getstate_body (
    SpiceChar *body,
    SpiceDouble et,
    SpiceDouble *state);

int getstate_observer (
    SpiceChar *body,
    SpiceDouble et,
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble alt,
    SpiceDouble *state,
    SpiceDouble *bodystate
    );

int getstate_solar_target (
    SpiceDouble et,
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble omegas,
    SpiceDouble *tstate,
    SpiceDouble *cstate);

void relstate (
    SpiceDouble *sobs,
    SpiceDouble *stgt,
    SpiceDouble *srel,
    SpiceDouble *losv,  /* line of sight vector */
    SpiceDouble *dist,
    SpiceDouble *vrad,  /* radial velocity */
    SpiceDouble *lt);

void get_pointing (
    SpiceDouble et,
    SpiceDouble *relstate_sun,
    SpiceDouble *relstate_tgt,
    SpiceDouble *x,
    SpiceDouble *y);

int getstate_pointing (
    SpiceChar *body,        /* target body name */
    SpiceDouble et,
    SpiceDouble x,           /* pointing angle from disk center */
    SpiceDouble y,           /* pointing angle from disk center */
    SpiceDouble *state_obs,  /* observer inertial state */
    SpiceDouble *state_tgt,
    SpiceBoolean *onbody);

void pointing2lola (
    SpiceDouble *state_sun,
    SpiceDouble et,
    SpiceDouble x,           /* pointing angle from disk center */
    SpiceDouble y,           /* pointing angle from disk center */
    SpiceDouble *state_obs,  /* observer inertial state */
    SpiceDouble *lon,
    SpiceDouble *lat,
    SpiceBoolean *onbody);

SpiceDouble grav_redshift (SpiceDouble GM, SpiceDouble r_emit, SpiceDouble r_obs);

void print_ephtable_head (
    FILE *stream,
    SpiceChar *observer,
    SpiceInt rotmodel);
void print_ephtable_row (FILE *stream, soleph_t *eph);
void fancy_print_eph (FILE *stream, soleph_t *eph);
int write_fits_ephtable_header (fitsfile *fptr, long nrows, int *status);
int write_fits_ephtable_row (
    fitsfile *fptr,
    long row,
    soleph_t *eph,
    int *status);
void printstate (SpiceDouble *s);
void printvec (SpiceDouble *s);

SpiceDouble aspr(void);
SpiceDouble dpas(void);
SpiceDouble rpas(void);


bool observer_on_earth(SpiceChar *observer);

#endif /* _SOLARV_H_ */
