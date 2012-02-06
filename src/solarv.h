/*
 Copyright (C) 2012 Hans-Peter Doerr <doerr@kis.uni-freiburg.de>,
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

#include "SpiceUsr.h"

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif

#ifndef RETURN_SUCCESS
#define RETURN_SUCCESS (1)
#endif

#ifndef RETURN_FAILURE
#define RETURN_FAILURE (0)
#endif


#define SOLAR_RADIUS (6.96E8) /* meters, Stix (2004) */

typedef struct
{
    double A;
    double B;
    double C;
    char name[12];
    char descr[128];
} rotmodel_t;

static const rotmodel_t RotModels[] =
{
    {2.851,  0.000,  0.000, "rigid", "rigid body rotation with 2.851 murad/s"},
    {0.000,  0.000,  0.000, "fixed", "no rotation, fixed to J2000 frame"},
    {2.851, -0.343, -0.474, "su90s", "Snodgrass & Ulrich (1990), spectroscopic"}, 
    {2.972, -0.484, -0.361, "su90g", "Snodgrass & Ulrich (1990), supergranulation"},
    {2.879, -0.339, -0.485, "su90m", "Snodgrass & Ulrich (1990), magnetic"},
    {2.836, -0.344, -0.504, "s84s", "Snodgrass (1984), spectroscopic MWO data"}
};
enum RotModel { rigid = 0,
		fixed,
		su90s,
		su90g,
		su90m,
		s84s,
		RotModel_END
};


/* 
   data structure to store epehemeris data of the target position on the
   solar surface (t_) and the sun barycenter (s_) with respect to the
   observer's position.

   angles in radian, distances in kilometers, velocity in meters/seconds
*/
typedef struct
{
    /* common data; sun global parameters */
    SpiceDouble jdate;     /* julian day of the event                    */
    SpiceChar utcdate[80]; /* ascii date in UTC                          */
    SpiceChar station[128];/* NAIF station name                          */
    SpiceDouble B0;        /* lat of sub-observer point                  */
    SpiceDouble L0;        /* lon of sub-observer point                  */
    SpiceDouble P0;        /* polar angle                                */
    SpiceDouble c_dist;    /* distance obs. to solar center              */
    SpiceDouble c_vlos;    /* radial velocity of obs to solar center     */
    SpiceDouble rsun_as;   /* apparent radius of the sun in arcsecs      */
    int rotmodel;          /* solar rotation model used                  */
    char modelname[12];    /* name of the rotation model                 */
    char modeldescr[128];  /* description of the rotation model          */
    								         
    /* target position parameters */				         
    SpiceDouble lon;       /* stonyhurst target longitude (deg)          */
    SpiceDouble lat;       /* stonyhurst target latitude  (deg)          */
    SpiceDouble x;         /* target x coordinate in as from disk center */
    SpiceDouble y;         /* target y coordinate in as from disk center */
    SpiceDouble mu;        /* heliocentric angle of the target           */
    SpiceDouble dist;      /* distance to target                         */
    SpiceDouble vlos;      /* radial velocity of target                  */
    SpiceDouble rho;       /* distance target to solar rotation axis     */
    SpiceDouble omega;     /* actually used omega value (murad/s)        */
} soleph_t;


void usage (FILE *stream);

int soleph (
    SpiceChar *station, /* NAIF body name/code of the observer     */ 
    SpiceDouble et,     /* Spice ephemeris time of the observation */
    SpiceDouble lon,    /* stonyhurst longitude of target point    */
    SpiceDouble lat,    /* stonyhurst latitude of target point     */
    soleph_t *eph);

int relstate_observer_sun (
    SpiceChar *station,
    SpiceDouble et,
    soleph_t *eph,
    SpiceDouble *state);

int relstate_sun_target (
    SpiceDouble et,
    SpiceDouble lon,
    SpiceDouble lat,
    soleph_t *eph,
    SpiceDouble *state);



int station_eph (
    SpiceChar *station, /* Observer's NAIF BODY_NAME ("Izana")     */
    SpiceDouble et,     /* spice ephemeris time                    */
    SpiceDouble lon,    /* stonyhurst longitude                    */
    SpiceDouble lat,    /* stonyhurst latitude                     */
    soleph_t *eph,       /* pt. to struct where to store ephem data */
    int rotModel
    );

int target_state (
    SpiceDouble lon,
    SpiceDouble lat,
    SpiceDouble lon0,
    SpiceDouble et,
    int model,
    SpiceDouble *state,
    soleph_t *eph
    );

SpiceDouble omega_sun (SpiceDouble lat, int model);

void print_ephtable_head (FILE *stream);
void print_ephtable_row (FILE *stream, soleph_t *eph);
void fancy_print_eph (FILE *stream, soleph_t *eph);
void list_rotation_models (FILE *stream);


#endif /* _SOLARV_H_ */
