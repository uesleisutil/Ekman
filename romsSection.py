"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  12 Apr 2019
Version:        1.5

Generates a new ROMS output file from scratch.
It is netCDF4 CF-compliant.
"""

from   netCDF4    import Dataset
from   setOptions import *
from   numpy      import dtype
import numpy      as np
import time

if romsTemp or romsSalt or romsZeta or romsTKE or romsLatent or romsSensible or romsLWRad or romsSWRad or romsEvaporation or romsEminusP or romsUwind or romsVwind or romsW or romsOmega == True:
    romsMassPoints = True
else:
    romsMassPoints = False   
if romsU or romsV or romsUbar or romsVbar== True:
    romsUVPoints = True
else:
    romsUVPoints = False

romsFillVal = 1.e+37

def romsVars(romsOriDir,romsNewDir):
    print('ROMS section activated.')
    # Original output file.
    romsRawFile             = Dataset(romsOriDir, mode='r')
    romsNewFile             = Dataset(romsNewDir, 'w', format='NETCDF4')
    romsNewFile.title       = "ROMS output file"
    romsNewFile.description = "Created with Ekman Toolbox at " + time.ctime(time.time())
    romsNewFile.link        = "https://github.com/uesleisutil/Ekman"

    # New ROMS output file.
    romsNewFile.createDimension('ocean_time', 0)
    romsNewOTdim           = romsNewFile.createVariable('ocean_time', dtype('double').char, ('ocean_time'))
    romsNewOTdim.long_name = 'seconds since 1900-01-01 00:00:00'
    romsNewOTdim.units     = 'seconds since 1900-01-01 00:00:00'

    # If a variable on mass point has been chosen.
    if romsMassPoints == True:
        eta_rho = romsRawFile.dimensions['eta_rho']
        xi_rho  = romsRawFile.dimensions['xi_rho']
        s_rho   = romsRawFile.dimensions['s_rho']
        s_w     = romsRawFile.dimensions['s_w']       
        lon_rho = romsRawFile.variables['lon_rho'][:,:]
        lat_rho = romsRawFile.variables['lat_rho'][:,:] 

        romsNewFile.createDimension('eta_rho', len(eta_rho))    
        romsNewFile.createDimension('xi_rho', len(xi_rho))
        romsNewFile.createDimension('s_rho', len(s_rho))
        romsNewFile.createDimension('s_w', len(s_w))

        romsNewLon               = romsNewFile.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
        romsNewLon.long_name     = 'Longitude on RHO-points'
        romsNewLon.units         = 'degree_east'
        romsNewLon.standard_name = 'longitude'
        romsNewLon[:, :]         = lon_rho

        romsNewLat               = romsNewFile.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=romsFillVal)
        romsNewLat.long_name     = 'Latitude on RHO-points'
        romsNewLat.units         = 'degree_north'
        romsNewLat.standard_name = 'latitude'
        romsNewLat[:, :]         = lat_rho

        # If ROMS potential temperature has been chosen.
        if romsTemp == True:
            print('Working on ROMS Temperature.')
            romsRawVar            = romsRawFile.variables['temp'][:,:,:,:]
            romsNewVar            = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name  = 'Potential Temperature'
            romsNewVar.units      = 'Celsius'
            romsNewVar._FillValue = '1.e+37f'
            romsNewVar[:,:,:,:]   = romsRawVar
            del romsRawVar, romsNewVar

        # If ROMS salinity has been chosen.
        if romsSalt == True:
            print('Working on ROMS Salinity.')
            romsRawVar            = romsRawFile.variables['salt'][:,:,:,:]
            romsNewVar            = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name  = 'Salinity'
            romsNewVar._FillValue = '1.e+37f'
            romsNewVar.units      = 'PSU'
            romsNewVar[:,:,:,:]   = romsRawVar
            del romsRawVar, romsNewVar
        
        if romsZeta == True:
            print('Working on ROMS Zeta.')
            romsRawVar            = romsRawFile.variables['zeta'][:,:,:]
            romsNewVar            = romsNewFile.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name  = 'Free-surface'
            romsNewVar.units      = 'meter'
            romsNewVar[:,:,:]     = romsRawVar
            del romsRawVar, romsNewVar            

        if romsTKE == True:
            print('Working on ROMS TKE.')
            romsRawVar            = romsRawFile.variables['tke'][:,:,:]
            romsNewVar            = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_w', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name  = 'Turbulent Kinectic Energy'
            romsNewVar.units      = 'm2 s-2'
            romsNewVar[:,:,:,:]   = romsRawVar
            del romsRawVar, romsNewVar               
          
        if romsLatent== True:
            print('Working on ROMS Latent Heat Flux.')
            romsRawVar                = romsRawFile.variables['latent'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('latent', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Latent Heat Flux'
            romsNewVar.negative_value = "upward flux = cooling"
            romsNewVar.positive_value = "downward flux = heating"
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsSensible== True:
            print('Working on ROMS Sensible Heat Flux.')
            romsRawVar                = romsRawFile.variables['sensible'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('sensible', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Sensible Heat Flux'
            romsNewVar.negative_value = "upward flux = cooling"
            romsNewVar.positive_value = "downward flux = heating"
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsLWRad== True:
            print('Working on ROMS Longwave Radiation Flux.')
            romsRawVar                = romsRawFile.variables['lwrad'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('lwrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Longwave Radiation Flux'
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsSWRad== True:
            print('Working on ROMS Shortwave Radiation Flux.')
            romsRawVar                = romsRawFile.variables['swrad'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('swrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Shortwave Radiation Flux'
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsEminusP== True:
            print('Working on ROMS Bulk Flux Surface Net Freshwater Flux.')
            romsRawVar                = romsRawFile.variables['EminusP'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('EminusP', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Bulk Flux Surface Net Freshwater Flux'
            romsNewVar.negative_value = "Upward = Freshening (Net Precipitation)"
            romsNewVar.positive_value = "Downward = Salting (Net Evaporation)"
            romsNewVar.units          = 'meter s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsEvaporation== True:
            print('Working on ROMS Evaporation Rate.')
            romsRawVar                = romsRawFile.variables['evaporation'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('evaporation', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Bulk Flux Surface Net Freshwater Flux'
            romsNewVar.negative_value = "Downward = Freshening (Condensation)"
            romsNewVar.positive_value = "Upward = Salting (Evaporation)"
            romsNewVar.units          = 'Kg m-2 s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsUwind== True:
            print('Working on ROMS U-wind Component.')
            romsRawVar                = romsRawFile.variables['Uwind'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('uwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Surface U-wind Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsVwind== True:
            print('Working on ROMS V-wind ;component.')
            romsRawVar                = romsRawFile.variables['Vwind'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('vwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Surface V-wind Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar                   

        if romsW== True:
            print('Working on ROMS Vertical Momentum Component.')
            romsRawVar                = romsRawFile.variables['w'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_w','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Vertical Momentum Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar     

        if romsOmega== True:
            print('Working on ROMS S-coordinate Vertical Momentum Component.')
            romsRawVar                = romsRawFile.variables['omega'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_w', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'S-coordinate Vertical Momentum Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:,:]       = romsRawVar
            del romsRawVar, romsNewVar     


    if romsMassPoints == False:
        print('No variables on mass-points has been chosen. Continuing.')
        pass

    # If a variable on U and/or V point has been chosen.
    if romsUVPoints == True:
        eta_u = romsRawFile.dimensions['eta_u']
        xi_u  = romsRawFile.dimensions['xi_u']
        lon_u = romsRawFile.variables['lon_u'][:,:]
        lat_u = romsRawFile.variables['lat_u'][:,:] 

        eta_v = romsRawFile.dimensions['eta_v']
        xi_v  = romsRawFile.dimensions['xi_v']
        lon_v = romsRawFile.variables['lon_v'][:,:]
        lat_v = romsRawFile.variables['lat_v'][:,:] 

        s_rho = romsRawFile.dimensions['s_rho']

        romsNewFile.createDimension('eta_u', len(eta_u))    
        romsNewFile.createDimension('xi_u', len(xi_u))

        romsNewFile.createDimension('eta_v', len(eta_v))    
        romsNewFile.createDimension('xi_v', len(xi_v))

        if romsMassPoints == True:
            pass
        else:
            romsNewFile.createDimension('s_rho', len(s_rho))

        romsNewLonU               = romsNewFile.createVariable('lon_u', 'd', ('eta_u', 'xi_u'), zlib=True, fill_value=romsFillVal)
        romsNewLonU.long_name     = 'Longitude on U-points'
        romsNewLonU.units         = 'degree_east'
        romsNewLonU.standard_name = 'longitude'
        romsNewLonU[:, :]         = lon_u

        romsNewLonV               = romsNewFile.createVariable('lon_v', 'd', ('eta_v', 'xi_v'), zlib=True, fill_value=romsFillVal)
        romsNewLonV.long_name     = 'Longitude on V-points'
        romsNewLonV.units         = 'degree_east'
        romsNewLonV.standard_name = 'longitude'
        romsNewLonV[:, :]         = lon_v

        romsNewLatU               = romsNewFile.createVariable('lat_u', 'd', ('eta_u', 'xi_u'), fill_value=romsFillVal)
        romsNewLatU.long_name     = 'Latitude on U-points'
        romsNewLatU.units         = 'degree_north'
        romsNewLatU.standard_name = 'latitude'
        romsNewLatU[:, :]         = lat_u

        romsNewLatV               = romsNewFile.createVariable('lat_v', 'd', ('eta_v', 'xi_v'), fill_value=romsFillVal)
        romsNewLatV.long_name     = 'Latitude on U-points'
        romsNewLatV.units         = 'degree_north'
        romsNewLatV.standard_name = 'latitude'
        romsNewLatV[:, :]         = lat_v

        # If ROMS V-velocity has been selected.
        if romsV == True:
            print('Working on ROMS V-velocity')
            romsRawVar           = romsRawFile.variables['v'][:,:,:,:]
            romsNewVar           = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'V-Velocity'
            romsNewVar.units     = 'm/s'
            romsNewVar[:,:,:,:]  = romsRawVar
            del romsRawVar, romsNewVar

        if romsU == True:
            print('Working on ROMS U-Velocity.')
            romsRawVar           = romsRawFile.variables['u'][:,:,:,:]
            romsNewVar           = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'U-Velocity'
            romsNewVar.units     = 'm/s'
            romsNewVar[:,:,:,:]  = romsRawVar
            del romsRawVar, romsNewVar

        if romsUbar == True:
            print('Working on ROMS Ubar.')
            romsRawVar           = romsRawFile.variables['ubar'][:,:,:]
            romsNewVar           = romsNewFile.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'Vertically Integrated U-momentum Component'
            romsNewVar.units     = 'm s-1'
            romsNewVar[:,:,:]   = romsRawVar
            del romsRawVar, romsNewVar

        if romsVbar == True:
            print('Working on ROMS Vbar.')
            romsRawVar           = romsRawFile.variables['vbar'][:,:,:]
            romsNewVar           = romsNewFile.createVariable('vbar', 'f', ('ocean_time', 'eta_v', 'xi_v'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'Vertically Integrated V-momentum Component'
            romsNewVar.units     = 'm s-1'
            romsNewVar[:,:,:]   = romsRawVar
            del romsRawVar, romsNewVar
    if romsUVPoints == False:
        print('No variable on U or V-points has been chosen. Continuing.')
