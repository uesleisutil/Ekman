"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  03 Jan 2020
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

selectRomsVars       = True
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
romsUbar             = True
romsVbar             = False

###########################
##### Sea-Ice section #####
###########################
iceOriginalFilename = 'ocean_his.nc'
iceNewFilename      = 'sea_ice_his.nc'

selectIceBox        = True
iceBox              = [-53, -40, -32, -23]

selectIceTimeStep   = True
iceTimeStep        = np.arange(100, 255+1)

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

selectWrfVars       = False
selectWrfBox        = False
selectWrfLevel      = False

wrfLevel            = 1
wrfBox             = [-53, -40, -32, -23]  # [lon_min, lon_max, lat_min, lat_max]

wrfTemp             = True
wrfPotTemp          = True
wrfRH               = True
wrfU                = False
wrfV                = False

####### WW3 section.