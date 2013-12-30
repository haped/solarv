Observer-Centric Radial Tangential Normal (OBSCRTN) Frame

created: 2013 by H.-P. Doerr <doerr@kis.uni-freiburg.de>

This dynamic frame definition is an addition to, and was heavily inspired
by, the frame definitions in heliospheric.tf by W. Thompson.

The OBSCRTN frame simplifies the computation of several ephemeris data wrt
the solar disk (such as the p-angle). This file is solely intended to be
used with the solarv utility as it needs to be referenced to a PRI_TARGET
frame at runtime.

The FRAME_PRI_TARGET variable is set at runtime to match the observer
location.  Observer locations are defined in the stations.bf, stations.tf
files, which can be transfered in a .bsp kernel with the pinpoint tool.

     Definition of the Observer-Centric RTN Frame
 
              All vectors are geometric: no aberration corrections are used.
 
              The position of Observer relative to the Sun is the primary
              vector: the X axis points from the Sun center to Earth
 
              The solar rotation axis is the secondary vector: the Z axis is
	      the component of the solar north direction perpendicular to X.
 
              The Y axis is Z cross X, completing the right-handed reference
              frame.


\begindata

        FRAME_OBSCRTN                =  1803430
        FRAME_1803430_NAME           = 'OBSCRTN'
        FRAME_1803430_CLASS          =  5
        FRAME_1803430_CLASS_ID       =  18034301
        FRAME_1803430_CENTER         =  10
        FRAME_1803430_RELATIVE       = 'J2000'
        FRAME_1803430_DEF_STYLE      = 'PARAMETERIZED'
        FRAME_1803430_FAMILY         = 'TWO-VECTOR'
        FRAME_1803430_PRI_AXIS       = 'X'
        FRAME_1803430_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
        FRAME_1803430_PRI_OBSERVER   = 'SUN'
        FRAME_1803430_PRI_TARGET     = 'OBSERVER'
        FRAME_1803430_PRI_ABCORR     = 'NONE'
        FRAME_1803430_PRI_FRAME      = 'IAU_SUN'
        FRAME_1803430_SEC_AXIS       = 'Z'
        FRAME_1803430_SEC_VECTOR_DEF = 'CONSTANT'
        FRAME_1803430_SEC_FRAME      = 'IAU_SUN'
        FRAME_1803430_SEC_SPEC       = 'RECTANGULAR'
        FRAME_1803430_SEC_VECTOR      = ( 0, 0, 1 )

\begintext

