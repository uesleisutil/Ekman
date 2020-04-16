"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  12 Apr 2019
Version:        1.7

Use this file to select what do you want Ekman Toolbox to do.
"""

####### Project section.
useDask              = True
projectName          = 'SC_2008'
projectAuthor        = 'Ueslei Adriano Sutil'

####### ROMS section. 
romsOriginalFilename = 'ocean_his.nc'
romsNewFilename      = 'ocean_his_new.nc'

selectRomsVars       = True
selectRomsBox        = True
selectRomsLevel      = True

romsZLevel           =-1 # Last sigma layer corresponds to ocean surface.
romsBox              = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]

romsTemp             = True
romsSalt             = False
romsZeta             = False
romsTKE              = False
romsRho              = False
romsLatent           = False
romsSensible         = False
romsLWRad            = False
romsSWRad            = False
romsEminusP          = False
romsEvaporation      = False
romsUwind            = False
romsVwind            = False
romsW                = False
romsOmega            = False
romsU                = False
romsV                = False
romsUbar             = False
romsVbar             = False

####### Sea-ice section.
iceOriginalFilename = 'ocean_his.nc'
iceNewFilename      = 'sea_ice_his.nc'

selectIceVars       = False
selectIceBox        = False

iceBox              = romsBox

iceAge              = True
iceA                = True
iceH                = True
iceV                = True
iceU                = True
iceSnowThick        = True
iceSurfaceTemp      = True
iceOceanMassFlux    = True
iceInteriorTemp     = True

####### WRF section.


####### WW3 section.