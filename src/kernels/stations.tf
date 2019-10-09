KPL/FK
 
   FILE: kernels/stations.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2020-02-13T17:37:17
   PINPOINT DEFINITIONS FILE: stations.defs
   PINPOINT PCK FILE:         /home/doerr/local/share/solarv/kernels/pck00010.tpc
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
 
   NAIF_BODY_NAME                      += 'GST'
   NAIF_BODY_CODE                      += 399924
 
   NAIF_BODY_NAME                      += 'DKIST'
   NAIF_BODY_CODE                      += 399925
 
   NAIF_BODY_NAME                      += 'IAG'
   NAIF_BODY_CODE                      += 399926
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame VTT_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame VTT_TOPO is centered at the
      site VTT, which has Cartesian coordinates
 
         X (km):                  0.5390264629397E+04
         Y (km):                 -0.1597706830465E+04
         Z (km):                  0.3007202236939E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       -16.5101400000000
         Latitude  (deg):        28.3023000000000
         Altitude   (km):         0.2444000000001E+01
 
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
   TKFRAME_1399919_ANGLES              =  ( -343.4898600000000,
                                             -61.6977000000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame SST_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame SST_TOPO is centered at the
      site SST, which has Cartesian coordinates
 
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
 
      Topocentric frame SCHAUINSLAND_TOPO is centered at the
      site SCHAUINSLAND, which has Cartesian coordinates
 
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
 
      Topocentric frame DST_TOPO is centered at the
      site DST, which has Cartesian coordinates
 
         X (km):                 -0.1463895180639E+04
         Y (km):                 -0.5166240553334E+04
         Z (km):                  0.3435666558701E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -105.8204980000000
         Latitude  (deg):        32.7872900000000
         Altitude   (km):         0.2800000000000E+01
 
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
 
      Topocentric frame MCMATH_TOPO is centered at the
      site MCMATH, which has Cartesian coordinates
 
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
 
   Topocentric frame GST_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame GST_TOPO is centered at the
      site GST, which has Cartesian coordinates
 
         X (km):                 -0.2390085850365E+04
         Y (km):                 -0.4706785166442E+04
         Z (km):                  0.3571333642599E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -116.9212700000000
         Latitude  (deg):        34.2585180000000
         Altitude   (km):         0.2043000000001E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_GST_TOPO                      =  1399924
   FRAME_1399924_NAME                  =  'GST_TOPO'
   FRAME_1399924_CLASS                 =  4
   FRAME_1399924_CLASS_ID              =  1399924
   FRAME_1399924_CENTER                =  399924
 
   OBJECT_399924_FRAME                 =  'GST_TOPO'
 
   TKFRAME_1399924_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399924_SPEC                =  'ANGLES'
   TKFRAME_1399924_UNITS               =  'DEGREES'
   TKFRAME_1399924_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399924_ANGLES              =  ( -243.0787300000000,
                                             -55.7414820000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame DKIST_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame DKIST_TOPO is centered at the
      site DKIST, which has Cartesian coordinates
 
         X (km):                 -0.5465999845833E+04
         Y (km):                 -0.2404405652428E+04
         Z (km):                  0.2242125713049E+04
 
      and planetodetic coordinates
 
         Longitude (deg):      -156.2560760000000
         Latitude  (deg):        20.7067420000000
         Altitude   (km):         0.3032000000000E+01
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_DKIST_TOPO                    =  1399925
   FRAME_1399925_NAME                  =  'DKIST_TOPO'
   FRAME_1399925_CLASS                 =  4
   FRAME_1399925_CLASS_ID              =  1399925
   FRAME_1399925_CENTER                =  399925
 
   OBJECT_399925_FRAME                 =  'DKIST_TOPO'
 
   TKFRAME_1399925_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399925_SPEC                =  'ANGLES'
   TKFRAME_1399925_UNITS               =  'DEGREES'
   TKFRAME_1399925_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399925_ANGLES              =  ( -203.7439240000000,
                                             -69.2932580000000,
                                             180.0000000000000 )
 
 
\begintext
 
   Topocentric frame IAG_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame IAG_TOPO is centered at the
      site IAG, which has Cartesian coordinates
 
         X (km):                  0.3913900299936E+04
         Y (km):                  0.6862529905339E+03
         Z (km):                  0.4972624006388E+04
 
      and planetodetic coordinates
 
         Longitude (deg):         9.9450000000000
         Latitude  (deg):        51.5593000000000
         Altitude   (km):         0.2010000000009E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame EARTH_FIXED.
 
 
\begindata
 
   FRAME_IAG_TOPO                      =  1399926
   FRAME_1399926_NAME                  =  'IAG_TOPO'
   FRAME_1399926_CLASS                 =  4
   FRAME_1399926_CLASS_ID              =  1399926
   FRAME_1399926_CENTER                =  399926
 
   OBJECT_399926_FRAME                 =  'IAG_TOPO'
 
   TKFRAME_1399926_RELATIVE            =  'EARTH_FIXED'
   TKFRAME_1399926_SPEC                =  'ANGLES'
   TKFRAME_1399926_UNITS               =  'DEGREES'
   TKFRAME_1399926_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399926_ANGLES              =  (   -9.9450000000000,
                                             -38.4407000000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file stations.defs
--------------------------------------------------------------------------------
 
begintext
 
        This is the definition for the solarv stations kernel. Stations are
        defined by their lat, lon, altitude coordinates.
 
        We also specify the *_UP, NORTH parameters so that pinpoint can
        add the corresponding topographic reference frames for each
        station.
 
        VTT: german vacuum tower telescope, tenerife
        SST: swedis solar telescope, la palma
        SCHAUINSLAND: old observatory of Kiepenheuer near Freiburg, Germany
        DST: Dunn solar telescope, NM, USA
        MCMATH: McMath Pierce solar facility, AZ, USA
        DKIST: Maui, Hawaii, USA
        GST: Goode solar telescope, Big Bear Lake, CA, USA
        IAG: Solar telscope of the Astronomical Institute, University of Goettingen, Germany
 
        Data sources
          VTT:          site survey with stationary GPS receiver
          SST:          wikipedia
          SCHAUINSLAND: google maps
          DST:          google maps, site-info
          MCMATH:       google maps, site-info
          DKIST:        google maps via wikipedia
          IAG:          google maps, altitude via gps
 
Last change: 2019, HP Doerr <doerr@mps.mpg.de>
 
 
begindata
         SITES      += 'VTT'
         VTT_CENTER = 399
         VTT_FRAME  = 'EARTH_FIXED'
         VTT_IDCODE = 399919
         VTT_LATLON = ( 28.30230, -16.51014, 2.444 )
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
 
         SITES      += 'GST'
         GST_CENTER = 399
         GST_FRAME  = 'EARTH_FIXED'
         GST_IDCODE = 399924
         GST_LATLON = ( 34.258518, -116.92127,  2.043)
         GST_UP     = 'Z'
         GST_NORTH  = 'X'
 
         SITES      += 'DKIST'
         DKIST_CENTER = 399
         DKIST_FRAME  = 'EARTH_FIXED'
         DKIST_IDCODE = 399925
         DKIST_LATLON = ( 20.706742, -156.256076,  3.032)
         DKIST_UP     = 'Z'
         DKIST_NORTH  = 'X'
 
         SITES      += 'IAG'
         IAG_CENTER = 399
         IAG_FRAME  = 'EARTH_FIXED'
         IAG_IDCODE = 399926
         IAG_LATLON = (51.55930, 9.94500, 0.201)
         IAG_UP     = 'Z'
         IAG_NORTH  = 'X'
 
 
begintext
 
begintext
 
[End of definitions file]
 
