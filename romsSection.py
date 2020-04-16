"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  13 Apr 2019
Version:        1.8

This file generates a new ROMS output file from scratch.
It is netCDF4 CF-compliant.

WARNING: Do not change anything in this file.
"""

from   netCDF4    import Dataset
from   setOptions import *
from   numpy      import dtype
from   matplotlib import path 
import numpy      as np
import time
import netCDF4

if romsTemp or romsSalt or romsZeta or romsTKE or romsLatent or romsSensible or romsLWRad or romsSWRad or romsEvaporation or romsEminusP or romsUwind or romsVwind or romsW or romsOmega or romsRho == True:
    romsMassPoints = True
else:
    romsMassPoints = False   
if romsU or romsV or romsUbar or romsVbar== True:
    romsUVPoints = True
else:
    romsUVPoints = False

romsFillVal = 1.e+37

def bbox2ij(lon,lat,romsBox=[-160., -155., 18., 23.]):
    """Return indices for i,j that will completely cover the specified bounding box.

    i0,i1,j0,j1 = bbox2ij(lon,lat,romsBox)
    
    lon,lat = 2D arrays that are the target of the subset
    romsBox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

    Example
    -------  
    >>> i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
    >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
    """
    romsBox=np.array(romsBox)
    mypath=np.array([romsBox[[0,1,1,0]],romsBox[[2,2,3,3]]]).T
    p = path.Path(mypath)
    points = np.vstack((lon.flatten(),lat.flatten())).T
    n,m = np.shape(lon)
    inside = p.contains_points(points).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside])

def romsVars(romsOriDir,romsNewDir):
    """
    Generates a new ROMS output file from scratch.
    """
    # Original output file.
    romsRawFile             = Dataset(romsOriDir, mode='r')
    romsNewFile             = Dataset(romsNewDir, 'w', format='NETCDF4')   
    romsNewFile.title       = "ROMS output file made by "+projectAuthor
    romsNewFile.description = "Created with Ekman Toolbox at " + time.ctime(time.time())
    romsNewFile.link        = "https://github.com/uesleisutil/Ekman"

    # New ROMS output file.
    romsNewFile.createDimension('ocean_time', 0)
    roms_time              = romsRawFile.variables['ocean_time']
    romsNewOTdim           = romsNewFile.createVariable('ocean_time', dtype('double').char, ('ocean_time'))
    romsNewOTdim.long_name = roms_time.units
    romsNewOTdim.units     = roms_time.units

    # If a variable on mass point has been chosen.
    if romsMassPoints == True:
        s_rho   = romsRawFile.dimensions['s_rho']
        s_w     = romsRawFile.dimensions['s_w']       
        if selectRomsBox == True:
            lon_rho     = romsRawFile.variables['lon_rho'][:,:]
            lat_rho     = romsRawFile.variables['lat_rho'][:,:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,romsBox)
            lon_rho     = romsRawFile.variables['lon_rho'][j0:j1, i0:i1]
            lat_rho     = romsRawFile.variables['lat_rho'][j0:j1, i0:i1]  
            romsNewFile.createDimension('eta_rho', len(lon_rho[:,0]))    
            romsNewFile.createDimension('xi_rho', len(lon_rho[0,:]))           
        else:            
            lon_rho = romsRawFile.variables['lon_rho'][:,:]
            lat_rho = romsRawFile.variables['lat_rho'][:,:] 
            eta_rho = romsRawFile.dimensions['eta_rho']
            xi_rho  = romsRawFile.dimensions['xi_rho']
            romsNewFile.createDimension('eta_rho', len(eta_rho))    
            romsNewFile.createDimension('xi_rho', len(xi_rho))   

        romsNewFile.createDimension('s_rho', len(s_rho))
        romsNewFile.createDimension('s_w', len(s_w))

        romsNewLon               = romsNewFile.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
        romsNewLon.long_name     = 'Longitude on RHO-points'
        romsNewLon.units         = 'degree_east'
        romsNewLon.standard_name = 'longitude'
        romsNewLon[:,:]          = lon_rho

        romsNewLat               = romsNewFile.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=romsFillVal)
        romsNewLat.long_name     = 'Latitude on RHO-points'
        romsNewLat.units         = 'degree_north'
        romsNewLat.standard_name = 'latitude'
        romsNewLat[:, :]         = lat_rho

        # If ROMS potential temperature has been chosen.
        if romsTemp == True:
            print('Working on ROMS Temperature.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['temp'][:,romsZLevel,j0:j1, i0:i1]  
                romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['temp'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['temp'][:,:,j0:j1, i0:i1]
                romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['temp'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'Potential Temperature'
            romsNewVar.units     = 'Celsius'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar

        # If ROMS salinity has been chosen.
        if romsSalt == True:
            print('Working on ROMS Salinity.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['salt'][:,romsZLevel,j0:j1, i0:i1]  
                romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['salt'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['salt'][:,:,j0:j1, i0:i1]
                romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['salt'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'Salinity'
            romsNewVar.units     = 'PSU'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar

        if romsTKE == True:
            print('Working on ROMS Turbulent Kinetic Energy.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['tke'][:,romsZLevel,j0:j1, i0:i1]  
                romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['tke'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['tke'][:,:,j0:j1, i0:i1]
                romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['tke'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'Turbulent Kinectic Energy'
            romsNewVar.units     = 'm2 s-2'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar         

        if romsRho == True:
            print('Working on ROMS Density Anomaly.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['rho'][:,romsZLevel,j0:j1, i0:i1]  
                romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['rho'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['rho'][:,:,j0:j1, i0:i1]
                romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['rho'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'Density Anomaly'
            romsNewVar.units     = 'kilogram meter-3'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar

        if romsZeta == True:
            print('Working on ROMS Zeta.')
            if selectRomsBox == True:
                romsRawVar        = romsRawFile.variables['zeta'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar        = romsRawFile.variables['zeta'][:,:,:]
            romsNewVar            = romsNewFile.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name  = 'Free-surface'
            romsNewVar.units      = 'meter'
            romsNewVar[:,:,:]     = romsRawVar
            del romsRawVar, romsNewVar            

        if romsLatent== True:
            print('Working on ROMS Latent Heat Flux.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['latent'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['latent'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('latent', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Latent Heat Flux'
            romsNewVar.negative_value = "upward flux = cooling"
            romsNewVar.positive_value = "downward flux = heating"
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsSensible== True:
            print('Working on ROMS Sensible Heat Flux.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['sensible'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['sensible'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('sensible', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Sensible Heat Flux'
            romsNewVar.negative_value = "upward flux = cooling"
            romsNewVar.positive_value = "downward flux = heating"
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsLWRad== True:
            print('Working on ROMS Longwave Radiation Flux.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['lwrad'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['lwrad'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('lwrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Longwave Radiation Flux'
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsSWRad== True:
            print('Working on ROMS Shortwave Radiation Flux.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['swrad'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['swrad'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('swrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Net Shortwave Radiation Flux'
            romsNewVar.units          = 'Watts m-2'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsEminusP== True:
            print('Working on ROMS Bulk Flux Surface Net Freshwater Flux.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['EminusP'][:,j0:j1, i0:i1]  
            else:                       
                romsRawVar            = romsRawFile.variables['EminusP'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('EminusP', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Bulk Flux Surface Net Freshwater Flux'
            romsNewVar.negative_value = "Upward = Freshening (Net Precipitation)"
            romsNewVar.positive_value = "Downward = Salting (Net Evaporation)"
            romsNewVar.units          = 'meter s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsEvaporation== True:
            print('Working on ROMS Evaporation Rate.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['evaporation'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['evaporation'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('evaporation', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Bulk Flux Surface Net Freshwater Flux'
            romsNewVar.negative_value = "Downward = Freshening (Condensation)"
            romsNewVar.positive_value = "Upward = Salting (Evaporation)"
            romsNewVar.units          = 'Kg m-2 s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsUwind== True:
            print('Working on ROMS U-wind Component.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['Uwind'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['Uwind'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('Uwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Surface U-wind Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar  

        if romsVwind== True:
            print('Working on ROMS V-wind ;component.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['Vwind'][:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['Vwind'][:,:,:]
            romsNewVar                = romsNewFile.createVariable('Vwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Surface V-wind Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:]         = romsRawVar
            del romsRawVar, romsNewVar                   

        if romsW== True:
            print('Working on ROMS Vertical Momentum Component.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['w'][:,:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['w'][:,:,:,:]
            romsNewVar                = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_w','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'Vertical Momentum Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:,:]       = romsRawVar
            del romsRawVar, romsNewVar     

        if romsOmega== True:
            print('Working on ROMS S-coordinate Vertical Momentum Component.')
            if selectRomsBox == True:
                romsRawVar            = romsRawFile.variables['omega'][:,:,j0:j1, i0:i1]  
            else:              
                romsRawVar            = romsRawFile.variables['omega'][:,:,:,:]
            romsNewVar                = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_w', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name      = 'S-coordinate Vertical Momentum Component'
            romsNewVar.units          = 'm s-1'
            romsNewVar[:,:,:,:]       = romsRawVar
            del romsRawVar, romsNewVar     

    if romsMassPoints == False:
        print('No ROMS variables on mass-points has been chosen. Continuing.')
        pass

    # If a variable on U and/or V point has been chosen.
    if romsUVPoints == True:
        if selectRomsBox == True:
            lon_u = romsRawFile.variables['lon_u'][:, :]
            lat_u = romsRawFile.variables['lat_u'][:, :]
            i0_u,i1_u,j0_u,j1_u = bbox2ij(lon_u,lat_u,romsBox)
            lon_u = romsRawFile.variables['lon_u'][j0_u:j1_u, i0_u:i1_u]
            lat_u = romsRawFile.variables['lat_u'][j0_u:j1_u, i0_u:i1_u]
            romsNewFile.createDimension('eta_u', len(lon_u[:,0]))    
            romsNewFile.createDimension('xi_u', len(lon_u[0,:])) 

            lon_v = romsRawFile.variables['lon_v'][:, :]
            lat_v = romsRawFile.variables['lat_v'][:, :]
            i0_v,i1_v,j0_v,j1_v = bbox2ij(lon_v,lat_v,romsBox)
            lon_v = romsRawFile.variables['lon_v'][j0_v:j1_v, i0_v:i1_v]
            lat_v = romsRawFile.variables['lat_v'][j0_v:j1_v, i0_v:i1_v]
            romsNewFile.createDimension('eta_v', len(lon_v[:,0]))    
            romsNewFile.createDimension('xi_v', len(lon_v[0,:])) 

        else:
            eta_u = romsRawFile.dimensions['eta_u']
            xi_u  = romsRawFile.dimensions['xi_u']
            lon_u = romsRawFile.variables['lon_u'][:,:]
            lat_u = romsRawFile.variables['lat_u'][:,:] 
            romsNewFile.createDimension('eta_u', len(eta_u))    
            romsNewFile.createDimension('xi_u', len(xi_u)) 
         
            eta_v = romsRawFile.dimensions['eta_v']
            xi_v  = romsRawFile.dimensions['xi_v']
            lon_v = romsRawFile.variables['lon_v'][:,:]
            lat_v = romsRawFile.variables['lat_v'][:,:] 
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
            print('Working on ROMS V-velocity.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['v'][:,romsZLevel,j0_v:j1_v, i0_v:i1_v] 
                romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['v'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['v'][:,:,j0_v:j1_v, i0_v:i1_v] 
                romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['v'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'V-Velocity'
            romsNewVar.units     = 'm s'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar    

        if romsU == True:
            print('Working on ROMS U-Velocity.')
            if selectRomsBox == True and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['u'][:,romsZLevel,j0_v:j1_v, i0_v:i1_v] 
                romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['v'][:,:,:,:]
                romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == True and selectRomsLevel == False:           
                romsRawVar = romsRawFile.variables['u'][:,:,j0_v:j1_v, i0_v:i1_v] 
                romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
            if selectRomsBox == False and selectRomsLevel == True:
                romsRawVar = romsRawFile.variables['u'][:,romsZLevel,:,:]  
                romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho','eta_rho', 'xi_rho'), zlib=True, fill_value=romsFillVal)
           
            romsNewVar.long_name = 'U-Velocity'
            romsNewVar.units     = 'm s'
            if selectRomsLevel == False:
                romsNewVar[:,:,:,:] = romsRawVar
            if selectRomsLevel == True:
                romsNewVar[:,:,:] = romsRawVar                
            del romsRawVar, romsNewVar    

        if romsUbar == True:
            print('Working on ROMS Ubar.')
            if selectRomsBox == True:
                romsRawVar       = romsRawFile.variables['ubar'][:,j0_u:j1_u, i0_u:i1_u]    
            else:              
                romsRawVar       = romsRawFile.variables['ubar'][:,:,:]
            romsNewVar           = romsNewFile.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'Vertically Integrated U-momentum Component'
            romsNewVar.units     = 'm s-1'
            romsNewVar[:,:,:]  = romsRawVar
            del romsRawVar, romsNewVar

        if romsVbar == True:
            print('Working on ROMS Vbar.')
            if selectRomsBox == True:
                romsRawVar       = romsRawFile.variables['vbar'][:,j0_v:j1_v, i0_v:i1_v]    
            else:              
                romsRawVar       = romsRawFile.variables['vbar'][:,:,:]
            romsNewVar           = romsNewFile.createVariable('vbar', 'f', ('ocean_time', 'eta_v', 'xi_v'), zlib=True, fill_value=romsFillVal)
            romsNewVar.long_name = 'Vertically Integrated V-momentum Component'
            romsNewVar.units     = 'm s-1'
            romsNewVar[:,:,:]  = romsRawVar
            del romsRawVar, romsNewVar
    if romsUVPoints == False:
        print('No ROMS variable on U or V-points has been chosen. Continuing.')
