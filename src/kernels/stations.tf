KPL/FK
 
   FILE: kernels/stations.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.0.0 --- March 26, 2009
   PINPOINT RUN DATE/TIME:    2013-12-26T17:06:37
   PINPOINT DEFINITIONS FILE: stations.defs
   PINPOINT PCK FILE:         /usr/local/share/solarv/kernels/pck00010.tpc
   PINPOINT SPK FILE:         kernels/stations.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'VTT'
   NAIF_BODY_CODE                      += 399919
 
   NAIF_BODY_NAME                      += 'SST'
   NAIF_BODY_CODE                      += 399920
 
   NAIF_BODY_NAME                      += 'SCHAUINSLAND'
   NAIF_BODY_CODE                      += 399921
 
   NAIF_BODY_NAME                      += 'DST'
   NAIF_BODY_CODE                      += 399922
 
   NAIF_BODY_NAME                      += 'MCMATH'
   NAIF_BODY_CODE                      += 399923
 
   NAIF_BODY_NAME                      += 'BIGBEAR'
   NAIF_BODY_CODE                      += 399924
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame VTT_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame VTT_TOPO is centered at the site VTT
      which has Cartesian coordinates
 
         X (km):                  0.5390236406562E+04
         Y (km):                 -0.1597689356553E+04
         Z (km):                  0.3007196324301E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       -16.5100510000000
         Latitude  (deg):        28.3023900000000
         Altitude   (km):         0.2413000000002E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_VTT_TOPO                      =  1399919
   FRAME_1399919_NAME                  =  'VTT_TOPO'
   FRAME_1399919_CLASS                 =  4
   FRAME_1399919_CLASS_ID              =  1399919
   FRAME_1399919_CENTER                =  399919
 
   OBJECT_399919_FRAME                 =  'VTT_TOPO'
 
   TKFRAME_1399919_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399919_SPEC                =  'ANGLES'
   TKFRAME_1399919_UNITS               =  'DEGREES'
   TKFRAME_1399919_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399919_ANGLES              =  ( -343.4899490000000,
                                             -61.6976100000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame SST_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame SST_TOPO is centered at the site SST
      which has Cartesian coordinates
 
         X (km):                  0.5327383061976E+04
         Y (km):                 -0.1718720024125E+04
         Z (km):                  0.3051718618510E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       -17.8807360000000
         Latitude  (deg):        28.7597330000000
         Altitude   (km):         0.2360000000001E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_SST_TOPO                      =  1399920
   FRAME_1399920_NAME                  =  'SST_TOPO'
   FRAME_1399920_CLASS                 =  4
   FRAME_1399920_CLASS_ID              =  1399920
   FRAME_1399920_CENTER                =  399920
 
   OBJECT_399920_FRAME                 =  'SST_TOPO'
 
   TKFRAME_1399920_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399920_SPEC                =  'ANGLES'
   TKFRAME_1399920_UNITS               =  'DEGREES'
   TKFRAME_1399920_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399920_ANGLES              =  ( -342.1192640000000,
                                             -61.2402670000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame SCHAUINSLAND_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame SCHAUINSLAND_TOPO is centered at the site SCHAUINSLAND
      which has Cartesian coordinates
 
         X (km):                  0.4242984704318E+04
         Y (km):                  0.5891574709801E+03
         Z (km):                  0.4711352342574E+04
 
      and planetodetic coordinates
 
         Longitude (deg):         7.9052290000000
         Latitude  (deg):        47.9134710000000
         Altitude   (km):         0.1239000000001E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_SCHAUINSLAND_TOPO             =  1399921
   FRAME_1399921_NAME                  =  'SCHAUINSLAND_TOPO'
   FRAME_1399921_CLASS                 =  4
   FRAME_1399921_CLASS_ID              =  1399921
   FRAME_1399921_CENTER                =  399921
 
   OBJECT_399921_FRAME                 =  'SCHAUINSLAND_TOPO'
 
   TKFRAME_1399921_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399921_SPEC                =  'ANGLES'
   TKFRAME_1399921_UNITS               =  'DEGREES'
   TKFRAME_1399921_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399921_ANGLES              =  (   -7.9052290000000,
                                             -42.0865290000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame DST_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame DST_TOPO is centered at the site DST
      which has Cartesian coordinates
 
         X (km):                 -0.1463895180639E+04
         Y (km):                 -0.5166240553334E+04
         Z (km):                  0.3435666558701E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -105.8204980000000
         Latitude  (deg):        32.7872900000000
         Altitude   (km):         0.2799999999999E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_DST_TOPO                      =  1399922
   FRAME_1399922_NAME                  =  'DST_TOPO'
   FRAME_1399922_CLASS                 =  4
   FRAME_1399922_CLASS_ID              =  1399922
   FRAME_1399922_CENTER                =  399922
 
   OBJECT_399922_FRAME                 =  'DST_TOPO'
 
   TKFRAME_1399922_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399922_SPEC                =  'ANGLES'
   TKFRAME_1399922_UNITS               =  'DEGREES'
   TKFRAME_1399922_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399922_ANGLES              =  ( -254.1795020000000,
                                             -57.2127100000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame MCMATH_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame MCMATH_TOPO is centered at the site MCMATH
      which has Cartesian coordinates
 
         X (km):                 -0.1994182352967E+04
         Y (km):                 -0.5037951695800E+04
         Z (km):                  0.3357632130562E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -111.5952430000000
         Latitude  (deg):        31.9584460000000
         Altitude   (km):         0.2096000000000E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_MCMATH_TOPO                   =  1399923
   FRAME_1399923_NAME                  =  'MCMATH_TOPO'
   FRAME_1399923_CLASS                 =  4
   FRAME_1399923_CLASS_ID              =  1399923
   FRAME_1399923_CENTER                =  399923
 
   OBJECT_399923_FRAME                 =  'MCMATH_TOPO'
 
   TKFRAME_1399923_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399923_SPEC                =  'ANGLES'
   TKFRAME_1399923_UNITS               =  'DEGREES'
   TKFRAME_1399923_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399923_ANGLES              =  ( -248.4047570000000,
                                             -58.0415540000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame BIGBEAR_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame BIGBEAR_TOPO is centered at the site BIGBEAR
      which has Cartesian coordinates
 
         X (km):                 -0.2390085850365E+04
         Y (km):                 -0.4706785166442E+04
         Z (km):                  0.3571333642599E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -116.9212700000000
         Latitude  (deg):        34.2585180000000
         Altitude   (km):         0.2042999999999E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_BIGBEAR_TOPO                  =  1399924
   FRAME_1399924_NAME                  =  'BIGBEAR_TOPO'
   FRAME_1399924_CLASS                 =  4
   FRAME_1399924_CLASS_ID              =  1399924
   FRAME_1399924_CENTER                =  399924
 
   OBJECT_399924_FRAME                 =  'BIGBEAR_TOPO'
 
   TKFRAME_1399924_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399924_SPEC                =  'ANGLES'
   TKFRAME_1399924_UNITS               =  'DEGREES'
   TKFRAME_1399924_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399924_ANGLES              =  ( -243.0787300000000,
                                             -55.7414820000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file stations.defs
--------------------------------------------------------------------------------
 
begintext
 
        This is the definition for the solarv stations kernel. Stations are
        defined by their lat, lon, altitude co-ordinates.
 
        We also specify the *_UP, NORTH parameters so that pinpoint can
        add the corresponding topographic reference frames for each
        station.
 
        Data sources
          VTT:          google maps
          SST:          wikipedia
          SCHAUINSLAND: google maps
          DST:          google maps, site-info
          MCMATH:       google maps, site-info
 
 
        Last change: 2013, HP Doerr <doerr@kis.uni-freiburg.de>
 
 
begindata
         SITES      += 'VTT'
         VTT_CENTER = 399
         VTT_FRAME  = 'EARTH_FIXED'
         VTT_IDCODE = 399919
         VTT_LATLON = ( 28.30239, -16.510051, 2.413 )
         VTT_UP     = 'Z'
         VTT_NORTH  = 'X'
 
         SITES      +=  'SST'
         SST_CENTER = 399
         SST_FRAME  = 'EARTH_FIXED'
         SST_IDCODE = 399920
         SST_LATLON = ( 28.759733, -17.880736, 2.360 )
         SST_UP     = 'Z'
         SST_NORTH  = 'X'
 
         SITES      += , 'SCHAUINSLAND'
         SCHAUINSLAND_CENTER = 399
         SCHAUINSLAND_FRAME  = 'EARTH_FIXED'
         SCHAUINSLAND_IDCODE = 399921
         SCHAUINSLAND_LATLON = ( 47.913471, 7.905229, 1.239 )
         SCHAUINSLAND_UP = 'Z'
         SCHAUINSLAND_NORTH = 'X'
 
         SITES      += 'DST'
         DST_CENTER = 399
         DST_FRAME  = 'EARTH_FIXED'
         DST_IDCODE = 399922
         DST_LATLON = ( 32.78729, -105.820498, 2.800 )
         DST_UP     = 'Z'
         DST_NORTH  = 'X'
 
         SITES      += , 'MCMATH'
         MCMATH_CENTER = 399
         MCMATH_FRAME  = 'EARTH_FIXED'
         MCMATH_IDCODE = 399923
         MCMATH_LATLON = ( 31.958446, -111.595243, 2.096 )
         MCMATH_UP     = 'Z'
         MCMATH_NORTH  = 'X'
 
         SITES      += 'BIGBEAR'
         BIGBEAR_CENTER = 399
         BIGBEAR_FRAME  = 'EARTH_FIXED'
         BIGBEAR_IDCODE = 399924
         BIGBEAR_LATLON = ( 34.258518, -116.92127,  2.043)
         BIGBEAR_UP     = 'Z'
         BIGBEAR_NORTH  = 'X'
 
 
 
begintext
 
begintext
 
[End of definitions file]
 
