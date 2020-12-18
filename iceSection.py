"""
Author: Ueslei Adriano Sutil
Created: 12 Apr 2019
Last modified: 12 Apr 2019
Version: 1.1

This file generates a new Sea Ice output file from scratch.
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

iceFillVal = 1.e+37

def bbox2ij(lon,lat,iceBox=[-160., -155., 18., 23.]):
    """Return indices for i,j that will completely cover the specified bounding box.

    i0,i1,j0,j1 = bbox2ij(lon,lat,iceBox)
    
    lon,lat = 2D arrays that are the target of the subset
    iceBox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

    Example
    -------  
    >>> i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
    >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
    """
    iceBox=np.array(iceBox)
    mypath=np.array([iceBox[[0,1,1,0]],iceBox[[2,2,3,3]]]).T
    p = path.Path(mypath)
    points = np.vstack((lon.flatten(),lat.flatten())).T
    n,m = np.shape(lon)
    inside = p.contains_points(points).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside])

def iceVars(iceOriDir,iceNewDir):
    """
    Generates a new Sea-ice model output file from scratch.
    """
    print('Sea-ice section activated.')
    # Original output file.
    iceRawFile             = Dataset(iceOriDir, mode='r')
    iceNewFile             = Dataset(iceNewDir, 'w', format='NETCDF4')
    iceNewFile.title       = "Budgell Sea-ice output file"
    iceNewFile.description = "Created with Ekman Toolbox at " + time.ctime(time.time())
    iceNewFile.link        = "https://github.com/uesleisutil/Ekman"

    # New Sea-ice output file.
    iceNewFile.createDimension('ocean_time', 0)
    ice_time              = iceRawFile.variables['ocean_time']
    iceNewOTdim           = iceNewFile.createVariable('ocean_time', dtype('double').char, ('ocean_time'))
    iceNewOTdim.long_name = ice_time.units
    iceNewOTdim.units     = ice_time.units

    s_rho   = iceRawFile.dimensions['s_rho']
    s_w     = iceRawFile.dimensions['s_w'] 
          
    if selectIceBox == True:
        lon_rho     = iceRawFile.variables['lon_rho'][:,:]
        lat_rho     = iceRawFile.variables['lat_rho'][:,:]
        i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,iceBox)
        lon_rho     = iceRawFile.variables['lon_rho'][j0:j1, i0:i1]
        lat_rho     = iceRawFile.variables['lat_rho'][j0:j1, i0:i1]  
        iceNewFile.createDimension('eta_rho', len(lon_rho[:,0]))    
        iceNewFile.createDimension('xi_rho', len(lon_rho[0,:]))           
    else:            
        lon_rho = iceNewFile.variables['lon_rho'][:,:]
        lat_rho = iceNewFile.variables['lat_rho'][:,:] 
        eta_rho = iceNewFile.dimensions['eta_rho']
        xi_rho  = iceNewFile.dimensions['xi_rho']
        iceNewFile.createDimension('eta_rho', len(eta_rho))    
        iceNewFile.createDimension('xi_rho', len(xi_rho))  
    
    iceNewLon               = iceNewFile.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
    iceNewLon.long_name     = 'Longitude on RHO-points'
    iceNewLon.units         = 'degree_east'
    iceNewLon.standard_name = 'longitude'
    iceNewLon[:,:]          = lon_rho

    iceNewLat               = iceNewFile.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=iceFillVal)
    iceNewLat.long_name     = 'Latitude on RHO-points'
    iceNewLat.units         = 'degree_north'
    iceNewLat.standard_name = 'latitude'
    iceNewLat[:, :]         = lat_rho    

    if iceAge == True:
        print('Working on Sea-ice Age.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['ageice'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['ageice'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('ageice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Sea-ice age'
        iceNewVar.units          = 's'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceA == True:
        print('Working on Sea-ice fraction of cell covered by ice.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['aice'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['aice'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('aice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Fraction of Cell Covered by Ice'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceH == True:
        print('Working on Sea-ice average ice thickness in cell.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['hice'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['hice'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('hice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Average Ice Thickness in Cell'
        iceNewVar.units          = 'meter'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceV == True:
        print('Working on Sea-ice V-velocity.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['vice'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['vice'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('vice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'V-component of Ice Velocity'
        iceNewVar.units          = 'm s-1'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceU == True:
        print('Working on Sea-ice U-velocity.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['uice'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['uice'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('vice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'U-component of Ice Velocity'
        iceNewVar.units          = 'm s-1'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceSnowThick == True:
        print('Working on Sea-ice U-velocity.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['snow_thick'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['snow_thick'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('snow_thick', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Sea-cover Thickness'
        iceNewVar.units          = 'meter'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceSurfaceTemp == True:
        print('Working on Sea-ice Surface Temperature.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['tisrf'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['tisrf'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('tisrf', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Sea-ice Surface Temperature'
        iceNewVar.units          = 'Degree Celsius'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceOceanMassFlux == True:
        print('Working on Ice-Ocean Mass Flux.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['iomflx'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['iomflx'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('iomflx', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Ice-Ocean Mass Flux'
        iceNewVar.units          = 'm s-1'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  

    if iceInteriorTemp == True:
        print('Working on Interior Ice Temperature.')
        if selectIceBox == True:
            iceRawVar            = iceRawFile.variables['ti'][:,j0:j1, i0:i1]  
        else:              
            iceRawVar            = iceRawFile.variables['ti'][:,:,:]
        iceNewVar                = iceNewFile.createVariable('ti', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True, fill_value=iceFillVal)
        iceNewVar.long_name      = 'Interior Ice Temperature'
        iceNewVar.units          = 'Degree Celcius'
        iceNewVar[:,:,:]         = iceRawVar
        del iceRawVar, iceNewVar  




