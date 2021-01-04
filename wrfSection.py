"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  03 Jan 2020
Version:        2.0

This file generates a new WRF output file from scratch.
It is netCDF4 CF-compliant.

WARNING: Do not change anything in this file.
"""

from   netCDF4    import Dataset
from   setOptions import *
from   numpy      import dtype
from   matplotlib import path 
from   wrf        import getvar
from progress.bar import IncrementalBar
import numpy      as np
import time

if wrfTemp or wrfPotTemp or wrfRH == True:
    wrfMassPoints = True
else:
    wrfMassPoints = False   
if wrfU or wrfV == True:
    wrfUVPoints = True
else:
    wrfUVPoints = False

wrfFillVal = 1.e+37

def bbox2ij(lon,lat,wrfBox=[-160., -155., 18., 23.]):
    """Return indices for i,j that will completely cover the specified bounding box.

    i0,i1,j0,j1 = bbox2ij(lon,lat,wrfBox)
    
    lon,lat = 2D arrays that are the target of the subset
    wrfBox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

    Example
    -------  
    >>> i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
    >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
    """
    wrfBox = np.array(wrfBox)
    mypath = np.array([wrfBox[[0,1,1,0]],wrfBox[[2,2,3,3]]]).T
    p = path.Path(mypath)
    points = np.vstack((lon.flatten(),lat.flatten())).T
    n,m = np.shape(lon)
    inside = p.contains_points(points).reshape((n,m))
    ii,jj = np.meshgrid(range(m),range(n))
    return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside])

def wrfVars(wrfOriDir,wrfNewDir):
    """
    Generates a new WRF output file from scratch.
    """
    # Original output file.
    wrfRawFile             = Dataset(wrfOriDir, mode='r')
    wrfNewFile             = Dataset(wrfNewDir, 'w', format='NETCDF4')   
    wrfNewFile.title       = "WRF output file made by "+projectAuthor
    wrfNewFile.description = "Created with Ekman Toolbox in " + time.ctime(time.time())
    wrfNewFile.link        = "https://github.com/uesleisutil/Ekman"

    # New wrf output file.
    wrfNewFile.createDimension('Times', 0)
    wrf_time              = wrfRawFile.variables['Times']
    wrfNewOTdim           = wrfNewFile.createVariable('Times', dtype('double').char, ('Times'))


    if wrfMassPoints == True:
        if selectWrfBox == True:
            print("Bounding box selected. Creating new XLAT and XLONG variables.")
            lon_wrf    = wrfRawFile.variables['XLONG'][:,:,:]
            lat_wrf    = wrfRawFile.variables['XLAT'][:,:]
            i0,i1,j0,j1 = bbox2ij(lon_wrf[0,:,:],lat_wrf[0,:,:],wrfBox)
            lon_wrf    = wrfRawFile.variables['XLONG'][:,j0:j1, i0:i1]
            lat_wrf    = wrfRawFile.variables['XLAT'][:,j0:j1, i0:i1]  
            wrfNewFile.createDimension('Time', len(lon_wrf[:,0,0]))    
            wrfNewFile.createDimension('south_north', len(lat_wrf[0,:,0]))      
            wrfNewFile.createDimension('west_east', len(lon_wrf[0,0,:]))               
        else:
            print("Copying XLAT and XLONG variables.")
            lon_wrf = wrfRawFile.variables['XLONG'][:,:,:]
            lat_wrf = wrfRawFile.variables['XLAT'][:,:,:] 
            wrfNewFile.createDimension('Time', len(lon_wrf[:,0,0]))    
            wrfNewFile.createDimension('south_north', len(lat_wrf[0,:,0]))      
            wrfNewFile.createDimension('west_east', len(lon_wrf[0,0,:]))   

        wrfNewLon               = wrfNewFile.createVariable('XLONG', 'd', ('Time','south_north', 'west_east'), fill_value=wrfFillVal)
        wrfNewLon.long_name     = 'Longitude on RHO-points'
        wrfNewLon.units         = 'degree_east'
        wrfNewLon.standard_name = 'longitude'
        wrfNewLon[:,:,:]        = lon_wrf

        wrfNewLat               = wrfNewFile.createVariable('XLAT', 'd', ('Time','south_north', 'west_east'), fill_value=wrfFillVal)
        wrfNewLat.long_name     = 'LATITUDE, SOUTH IS NEGATIVE'
        wrfNewLat.units         = 'degree_north'
        wrfNewLat.standard_name = 'latitude'
        wrfNewLat.stagger       = ' '
        wrfNewLat[:,:, :]       = lat_wrf

        # If WRF temperature has been chosen.
        if wrfTemp == True:
            print('Working on WRF Temperature.')
            foo        = wrfRawFile.variables['T'][:,:,0,0]
            ntimes     = foo[:,0]
            levels     = foo[:,0]
            del foo
            bar        = IncrementalBar('WRF Temperature:', max=len(ntimes))
            for i in range(0,len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1]                     
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Celsius'                   
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar   
                   
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])   
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[:,:,:]
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Celsius'    
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[:,:,:]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                          
                elif selectWrfBox == True and selectWrfLevel == False: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[:,j0:j1, i0:i1]                        
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Celsius'   
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[:,j0:j1, i0:i1]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                              
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar  = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Celsius'    
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar[i,:,:] = wrfRawVar  
                bar.next()
            bar.finish()

        # If WRF potential temperature has been chosen.
        if wrfPotTemp == True:
            print('Working on WRF Potentital Temperature.')
            foo        = wrfRawFile.variables['T'][:,:,0,0]
            ntimes     = foo[:,0]
            levels     = foo[:,0]
            del foo
            bar        = IncrementalBar('WRF Potential Temperature:', max=len(ntimes))
            for i in range(0,len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1]                     
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Celsius'                   
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar   
                   
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])   
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[:,:,:]
                        wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Celsius'    
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[:,:,:]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                          
                elif selectWrfBox == True and selectWrfLevel == False: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[:,j0:j1, i0:i1]                        
                        wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Celsius'   
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[:,j0:j1, i0:i1]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                              
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar  = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Celsius'    
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar[i,:,:] = wrfRawVar  
                bar.next()
            bar.finish()

        # If WRF Relative Humidity has been chosen.
        if wrfRH == True:
            print('Working on WRF Relative Humidity.')
            foo        = wrfRawFile.variables['T'][:,:,0,0]
            ntimes     = foo[:,0]
            levels     = foo[:,0]
            del foo
            bar        = IncrementalBar('WRF Relative Humidity:', max=len(ntimes))
            for i in range(0,len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1]                     
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'                   
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[wrfLevel,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar   
                   
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])   
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[:,:,:]
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'    
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[:,:,:]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                          
                elif selectWrfBox == True and selectWrfLevel == False: 
                    if i == 0:
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[:,j0:j1, i0:i1]                        
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time','bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'   
                    else:                
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[:,j0:j1, i0:i1]
                        wrfNewVar[i,:,:,:] = wrfRawVar                                              
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if i == 0:
                        wrfNewVar  = np.zeros([len(ntimes),len(lat_wrf), len(lon_wrf)])  
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar[i,:,:] = wrfRawVar 
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'    
                    else:               
                        wrfRawVar = getvar(wrfRawFile,'rh',meta=False, timeidx=i)[wrfLevel,:,:]                 
                        wrfNewVar[i,:,:] = wrfRawVar  
                bar.next()
            bar.finish()            
    if wrfMassPoints == False:
        print('No WRF variables on mass-points has been chosen. Continuing.')
        pass

