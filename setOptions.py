"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2020
Last modified:  03 Jan 2021
Version:        2.0

Use this file to select what do you want Ekman Toolbox to do.
"""
import numpy as np

####### Project section.
projectName          = 'SC_2008'
projectAuthor        = 'Ueslei Adriano Sutil'

###########################
####### ROMS section ######
###########################
romsOriginalFilename = 'ocean_his.nc'
romsNewFilename      = 'ocean_his_new.nc'

selectRomsBox        = True
romsBox              = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]

selectRomsLevel      = True
romsLevel            = np.arange(27, 28+1)

selectRomsTimeStep   = False
romsTimeStep         = np.arange(100, 255+1)

selectRomsVars       = False
romsTemp             = False
romsSalt             = False
romsTKE              = False
romsRho              = False
romsZeta             = False
romsW                = False
romsOmega            = False
romsLatent           = False
romsSensible         = False
romsLWRad            = False
romsSWRad            = False
romsEminusP          = False
romsEvaporation      = False
romsUwind            = False
romsVwind            = False
romsU                = False
romsV                = False
romsUbar             = False
romsVbar             = False

###########################
##### Sea-Ice section #####
###########################
iceOriginalFilename = 'ocean_his.nc'
iceNewFilename      = 'sea_ice_his.nc'

selectIceBox        = True
iceBox              = [-53, -40, -32, -23]

selectIceTimeStep   = True
iceTimeStep        = np.arange(200, 255+1)

selectIceVars       = False
iceAge              = False
iceA                = False
iceH                = False
iceV                = False
iceU                = False
iceSnowThick        = False
iceSurfaceTemp      = False
iceOceanMassFlux    = False
iceInteriorTemp     = False

###########################
####### WRF section #######
###########################
wrfOriginalFilename = 'wrf.nc'
wrfNewFilename      = 'wrf_new.nc'

selectWrfBox        = True
wrfBox              = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]

selectWrfLevel      = True
wrfLevel            = np.arange(27, 28+1)

selectWrfTimeStep   = True
wrfTimeStep         = np.arange(150, 162+1)

selectWrfVars       = True
wrfTemp             = False
wrfPotTemp          = False
wrfRh               = False
wrfTd               = False
wrfTwb              = False
wrfTv               = False
wrfPressure         = False
wrfAvo              = False
wrfPvo              = False
wrfDbz              = False
wrfGeopt            = False
wrfOmega            = False
wrfUnstaggeredU     = False
wrfUnstaggeredV     = False
wrfUnstaggeredW     = False
wrfUvmet            = False
wrfUvmet10m         = False
wrfLatent           = False
wrfSensible         = False
wrfSlp              = False
wrfRh2              = False
wrfTd2              = False
wrfTerrain          = False
wrfLandmask         = True

wrfU                = False
wrfV                = False

####### WW3 section.