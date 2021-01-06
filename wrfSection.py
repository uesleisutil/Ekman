"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2019
Last modified:  05 Jan 2021
Version:        2.12

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

if wrfTemp or wrfPotTemp or wrfRh or wrfTd or wrfTwb or wrfTv or wrfSST or wrfPressure or wrfUnstaggeredU or wrfUnstaggeredV or wrfUnstaggeredW or wrfUvmet or wrfUvmet10m or wrfLatent or wrfSensible or wrfSlp or wrfAvo or wrfDbz or wrfGeopt or wrfOmega or wrfPvo or wrfTerrain or wrfRh2 or wrfTd2 or wrfLandmask== True:
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
    wrfRawFile             = Dataset(wrfOriDir, mode='r')
    wrfNewFile             = Dataset(wrfNewDir, 'w', format='NETCDF4')   
    wrfNewFile.title       = "WRF output file made by "+projectAuthor
    wrfNewFile.description = "Created with Ekman Toolbox in " + time.ctime(time.time())
    wrfNewFile.link        = "https://github.com/uesleisutil/Ekman"
    
    # If a variable on mass point has been chosen.
    if wrfMassPoints == True:
        # Define Lat and Lon variables.
        if selectWrfBox == True:
            lon_wrf    = wrfRawFile.variables['XLONG'][0,:,:]
            lat_wrf    = wrfRawFile.variables['XLAT'][0,:,:]
            i0,i1,j0,j1 = bbox2ij(lon_wrf[:,:],lat_wrf[:,:],wrfBox)
            lon_wrf    = wrfRawFile.variables['XLONG'][0,j0:j1, i0:i1]
            lat_wrf    = wrfRawFile.variables['XLAT'][0,j0:j1, i0:i1]   
            wrfNewFile.createDimension('south_north', len(lat_wrf[:,0]))      
            wrfNewFile.createDimension('west_east', len(lon_wrf[0,:]))
            print("Bounding box selected. New domain limits are: Longitude "+str(wrfBox[0])+"/"+str(wrfBox[1])+" and Latitude "+str(wrfBox[2])+"/"+str(wrfBox[3])+".")    
        else:
            print("No bounding box selected: Using XLAT and XLONG variables from input file.")
            lon_wrf = wrfRawFile.variables['XLONG'][0,:,:]
            lat_wrf = wrfRawFile.variables['XLAT'][0,:,:]   
            wrfNewFile.createDimension('south_north', len(lat_wrf[:,0]))      
            wrfNewFile.createDimension('west_east', len(lon_wrf[0,:]))   

        wrfNewLon               = wrfNewFile.createVariable('XLONG', 'd', ('south_north', 'west_east'), fill_value=wrfFillVal)
        wrfNewLon.long_name     = 'Longitude on RHO-points'
        wrfNewLon.units         = 'degree_east'
        wrfNewLon.standard_name = 'longitude'
        wrfNewLon[:,:]        = lon_wrf

        wrfNewLat               = wrfNewFile.createVariable('XLAT', 'd', ('south_north', 'west_east'), fill_value=wrfFillVal)
        wrfNewLat.units         = 'degree_north'
        wrfNewLat.standard_name = 'latitude'
        wrfNewLat[:, :]       = lat_wrf

        # Define vertical levels and time-steps. 
        if selectWrfLevel == True and len(wrfLevel) == 1:
            print("One vertical level selected: Working on vertical level "+str(wrfLevel)+".")
            wrfNewFile.createDimension('bottom_top', len(wrfLevel))
        if selectWrfLevel == True and len(wrfLevel) > 1:
            print("Multiple vertical levels selected: Working from level "+str(wrfLevel[0])+" to "+str(wrfLevel[-1])+".")
            wrfNewFile.createDimension('bottom_top', len(wrfLevel))
        if selectWrfLevel == False:
            print("No selected vertical levels specified: Using entire vertical level from input file.")
            levels = getvar(wrfRawFile,'z',meta=False)[:,0,0]
            wrfNewFile.createDimension('bottom_top', len(levels))
        if selectWrfTimeStep == True:
            ntimes = wrfTimeStep
            print("Time-step selected: Working from time-step "+str(ntimes[0])+" to "+str(ntimes[-1])+".")
            wrfNewFile.createDimension('Time', len(ntimes))
        else:
            print("No time-step selected. Working with entire time-step.")
            ntimes = getvar(wrfRawFile,'LH',meta=False, timeidx=None)[:,0,0]
            ntimes = np.arange(np.argmin(ntimes), len(ntimes)) 
            wrfNewFile.createDimension('Time', len(ntimes))
            
        # If WRF Temperature has been chosen.
        if wrfTemp == True:
            print('Working on WRF Temperature.')
            bar = IncrementalBar('WRF Temperature:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'temp',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = wrfRawFile.variables['temp'][i,wrfStart:wrfStop,:, :] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('temp', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = wrfRawFile.variables['temp'][i,wrfStart:wrfStop,:, :]  
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish()  

        # If WRF Potential Temperature has been chosen.
        if wrfPotTemp == True:
            print('Working on WRF Potential Temperature.')
            bar = IncrementalBar('WRF Potential Temperature:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Potential Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'th',units="degC",meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = wrfRawFile.variables['th'][i,wrfStart:wrfStop,:, :] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('th', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = wrfRawFile.variables['th'][i,wrfStart:wrfStop,:, :]  
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish()  

        # If WRF Relative Humidity has been chosen.
        if wrfRh == True:
            print('Working on WRF Relative Humidity.')
            bar = IncrementalBar('WRF Relative Humidity:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Relative Humidity'
                            wrfNewVar.units     = '%'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Relative Humidity'
                            wrfNewVar.units     = '%'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Relative Humidity'
                        wrfNewVar.units     = '%'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Relative Humidity'
                            wrfNewVar.units     = '%'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'rh', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = wrfRawFile.variables['rh'][i,wrfStart:wrfStop,:, :] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('rh', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Relative Humidity'
                            wrfNewVar.units     = '%'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = wrfRawFile.variables['rh'][i,wrfStart:wrfStop,:, :]  
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Dew Point Temperature has been chosen.
        if wrfTd == True:
            print('Working on WRF Dew Point Temperature.')
            bar = IncrementalBar('WRF Dew Point Temperature:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Dew Point Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Dew Point Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Dew Point Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Dew Point Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Dew Point Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'td', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('td', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Dew Point Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'td',units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Wet Bulb Temperature has been chosen.
        if wrfTwb == True:
            print('Working on WRF Wet Bulb Temperature.')
            bar = IncrementalBar('WRF Wet Bulb Temperature:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Wet Bulb Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Wet Bulb Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Wet Bulb Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Wet Bulb Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Wet Bulb Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'twb', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('twb', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Wet Bulb Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'twb',units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 


        # If WRF Virtual Temperature has been chosen.
        if wrfTv == True:
            print('Working on WRF Virtual Temperature.')
            bar = IncrementalBar('WRF Virtual Temperature:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Virtual Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Virtual Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Virtual Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Virtual Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Virtual Temperature'
                            wrfNewVar.units     = 'Degree Celsius'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'tv', units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('tv', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Virtual Temperature'
                            wrfNewVar.units     = 'Degree Celsius'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'tv',units="degC", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Pressure has been chosen.
        if wrfPressure == True:
            print('Working on WRF Full Model Pressure.')
            bar        = IncrementalBar('WRF Full Model Pressure:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Full Model Pressure'
                            wrfNewVar.units     = 'hPa'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Full Model Pressure'
                            wrfNewVar.units     = 'hPa'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Full Model Pressure'
                        wrfNewVar.units     = 'hPa'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Full Model Pressure'
                        wrfNewVar.units     = 'hPa'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Full Model Pressure'
                            wrfNewVar.units     = 'hPa'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'p', units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('p', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Full Model Pressure'
                            wrfNewVar.units     = 'hPa'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'p',units="hPa", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Absolute Vorticity has been chosen.
        if wrfAvo == True:
            print('Working on WRF Absolute Vorticity.')
            bar        = IncrementalBar('WRF Absolute Vorticity:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'avo',  meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Absolute Vorticity'
                            wrfNewVar.units     = '10-5 s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Absolute Vorticity'
                            wrfNewVar.units     = '10-5 s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'avo',meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Absolute Vorticity'
                        wrfNewVar.units     = '10-5 s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Absolute Vorticity'
                        wrfNewVar.units     = '10-5 s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Absolute Vorticity'
                            wrfNewVar.units     = '10-5 s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('avo', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Absolute Vorticity'
                            wrfNewVar.units     = '10-5 s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'avo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Potential Vorticity has been chosen.
        if wrfPvo == True:
            print('Working on WRF Potential Vorticity.')
            bar        = IncrementalBar('WRF Potential Vorticity:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'pvo',  meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Vorticity'
                            wrfNewVar.units     = '10-5 s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Vorticity'
                            wrfNewVar.units     = '10-5 s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'pvo',meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Potential Vorticity'
                        wrfNewVar.units     = '10-5 s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Potential Vorticity'
                        wrfNewVar.units     = '10-5 s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Vorticity'
                            wrfNewVar.units     = '10-5 s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('pvo', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Potential Vorticity'
                            wrfNewVar.units     = '10-5 s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'pvo', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 


        # If WRF Reflectivity has been chosen.
        if wrfDbz == True:
            print('Working on WRF Reflectivity.')
            bar        = IncrementalBar('WRF Reflectivity:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'dbz',  meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Reflectivity'
                            wrfNewVar.units     = 'dBZ'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Reflectivity'
                            wrfNewVar.units     = 'dBZ'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'dbz',meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Reflectivity'
                        wrfNewVar.units     = 'dBZ'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Reflectivity'
                        wrfNewVar.units     = 'dBZ'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Reflectivity'
                            wrfNewVar.units     = 'dBZ'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('dbz', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Reflectivity'
                            wrfNewVar.units     = 'dBZ'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'dbz', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Geopotential Height for the Mass Grid has been chosen.
        if wrfGeopt == True:
            print('Working on WRF Geopotential Height for the Mass Grid.')
            bar        = IncrementalBar('WRF Geopotential Height for the Mass Grid:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'geopt',  meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                            wrfNewVar.units     = 'm2 s-2'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                            wrfNewVar.units     = 'm2 s-2'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'geopt',meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                        wrfNewVar.units     = 'm2 s-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                        wrfNewVar.units     = 'm2 s-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                            wrfNewVar.units     = 'm2 s-2'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('geopt', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Geopotential Height for the Mass Grid'
                            wrfNewVar.units     = 'm2 s-2'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'geopt', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF Omega has been chosen.
        if wrfOmega == True:
            print('Working on WRF Omega.')
            bar        = IncrementalBar('WRF Omega:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'omega',  meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Omega'
                            wrfNewVar.units     = 'Pa s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Omega'
                            wrfNewVar.units     = 'Pa s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'omega',meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Omega'
                        wrfNewVar.units     = 'Pa s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Omega'
                        wrfNewVar.units     = 'Pa s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Omega'
                            wrfNewVar.units     = 'Pa s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('omega', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'Omega'
                            wrfNewVar.units     = 'Pa s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'omega', meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 


        # If WRF U-component of Wind on Mass Points has been chosen.
        if wrfUnstaggeredU == True:
            print('Working on WRF U-component of Wind on Mass Points.')
            bar        = IncrementalBar('WRF U-component of Wind on Mass Points:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile, 'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'ua', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('ua', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'ua',units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF V-component of Wind on Mass Points has been chosen.
        if wrfUnstaggeredV == True:
            print('Working on WRF V-component of Wind on Mass Points.')
            bar        = IncrementalBar('WRF V-component of Wind on Mass Points:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile, 'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'va', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('va', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'V-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'va',units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

       # If WRF W-component of Wind on Mass Points has been chosen.
        if wrfUnstaggeredW == True:
            print('Working on WRF W-component of Wind on Mass Points.')
            bar        = IncrementalBar('WRF W-component of Wind on Mass Points:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[i,:,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile, 'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfLevel,:,:]   
                            wrfNewVar[i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'wa', units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewVar = wrfNewFile.createVariable('wa', 'f', ('Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'W-component of Wind on Mass Points'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'wa',units="m s-1", meta=False, timeidx=ntimes[0]+i)[wrfStart:wrfStop,:, :]
                            wrfNewVar[i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

       # If WRF U and V Components of Wind Rotated to Earth Coordinates has been chosen.
        if wrfUvmet == True:
            print('Working on WRF U and V Components of Wind Rotated to Earth Coordinates.')
            bar        = IncrementalBar('WRF U and V Components of Wind Rotated to Earth Coordinates:', max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True and selectWrfLevel == True:
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,wrfLevel,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([2,len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewFile.createDimension('UV', 2)
                            wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[:,i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,wrfLevel,j0:j1, i0:i1]  
                            wrfNewVar[:,i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar = np.zeros([2,len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewFile.createDimension('UV', 2)
                            wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[:,i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfStart:wrfStop,j0:j1, i0:i1] 
                            wrfNewVar[:,i,:,:] = wrfRawVar                              
                elif selectWrfBox == False and selectWrfLevel == False:    
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, :,:,:] 
                        wrfNewVar = np.zeros([2, len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewFile.createDimension('UV', 2)
                        wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[:,i,:,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, :,:,:] 
                        wrfNewVar[:,i,:,:,:] = wrfRawVar                                                                                
                elif selectWrfBox == True and selectWrfLevel == False:  
                    if i == np.argmin(ntimes):             
                        wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, :,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([2,len(ntimes),len(levels),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewFile.createDimension('UV', 2)
                        wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time', 'bottom_top', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[:,i,:,:,:] = wrfRawVar  
                    else: 
                        wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, :,j0:j1, i0:i1] 
                        wrfNewVar[i,:,:,:] = wrfRawVar  
                elif selectWrfBox == False and selectWrfLevel == True: 
                    if len(wrfLevel) == 1:
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile, 'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfLevel,:,:]   
                            wrfNewVar = np.zeros([2,len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewFile.createDimension('UV', 2)
                            wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                            wrfNewVar.units     = 'm s-1'
                            wrfNewVar[:,i,:,:] = wrfRawVar  
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfLevel,:,:]   
                            wrfNewVar[:,i,:,:] = wrfRawVar                                                             
                    else:
                        wrfStart  = slice(min(wrfLevel),max(wrfLevel)+1).start
                        wrfStop   = slice(min(wrfLevel),max(wrfLevel)+1).stop
                        if i == np.argmin(ntimes):
                            wrfRawVar = getvar(wrfRawFile,'uvmet', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfStart:wrfStop,:, :]
                            wrfNewVar = np.zeros([2,len(ntimes),len(wrfLevel),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                            wrfNewFile.createDimension('UV', 2)
                            wrfNewVar = wrfNewFile.createVariable('uvmet', 'f', ('UV','Time','bottom_top','south_north', 'west_east'), fill_value=wrfFillVal)
                            wrfNewVar.long_name = 'U and V Components of Wind Rotated to Earth Coordinates'
                            wrfNewVar.units     = 'm s-1'                     
                            wrfNewVar[:,i,:,:,:] = wrfRawVar
                        else: 
                            wrfRawVar = getvar(wrfRawFile,'uvmet',units="m s-1", meta=False, timeidx=ntimes[0]+i)[:, wrfStart:wrfStop,:, :]
                            wrfNewVar[:,i,:,:,:] = wrfRawVar  
                bar.next()
            bar.finish() 

        # If WRF 10 m U and V Components of Wind Rotated to Earth Coordinates been chosen.                                               
        if wrfUvmet10m == True:
            print('Working on  10 m U and V Components of Wind Rotated to Earth Coordinates .')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'uvmet10', units="m s-1",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([2,len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewFile.createDimension('UV', 2)
                        wrfNewVar = wrfNewFile.createVariable('uvmet10', 'f', ('UV','Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '10 m U and V Components of Wind Rotated to Earth Coordinates '
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[:,i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'uvmet10', units="m s-1",meta=False, timeidx=ntimes[0]+i)[:,j0:j1, i0:i1] 
                        wrfNewVar[:,i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'uvmet10', units="m s-1",meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar = np.zeros([2, len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewFile.createDimension('UV', 2)
                        wrfNewVar = wrfNewFile.createVariable('uvmet10', 'f', ('UV', 'Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '10 m U and V Components of Wind Rotated to Earth Coordinates '
                        wrfNewVar.units     = 'm s-1'
                        wrfNewVar[:,i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'uvmet10', units="m s-1", meta=False, timeidx=ntimes[0]+i)[:,:,:] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF Latent Heat Flux has been chosen.                                               
        if wrfLatent == True:
            print('Working on WRF Latent Heat Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'LH', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('LH', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Latent Heat Flux'
                        wrfNewVar.units     = 'W m-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'LH', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'LH', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('LH', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Latent Heat Flux'
                        wrfNewVar.units     = 'W m-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'LH', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF Sensible Heat Flux has been chosen.                                               
        if wrfSensible == True:
            print('Working on WRF Sensible Heat Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'HFX', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('HFX', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sensible Heat Flux'
                        wrfNewVar.units     = 'W m-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'HFX', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'HFX', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('HFX', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sensible Heat Flux'
                        wrfNewVar.units     = 'W m-2'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'HFX', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF Sea Level Pressure has been chosen.                                               
        if wrfSlp == True:
            print('Working on WRF Sea Level Pressure.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'slp', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('slp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sea Level Pressure'
                        wrfNewVar.units     = 'hPa'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'slp', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'slp', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('slp', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sea Level Pressure'
                        wrfNewVar.units     = 'hPa'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'slp', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF 2m Relative Humidity has been chosen. 
        if wrfRh2 == True:
            print('Working on WRF 2m Relative Humidity.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'rh2', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('rh2', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '2m Relative Humidity'
                        wrfNewVar.units     = '%'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'rh2', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'rh2', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('rh2', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '2m Relative Humidity'
                        wrfNewVar.units     = '%'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'rh2', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF Model Terrain Height has been chosen. 
        if wrfTerrain == True:
            print('Working on WRF Model Terrain Height.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'ter', units='m', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('ter', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Model Terrain Height'
                        wrfNewVar.units     = 'm'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'ter', units='m', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'ter', units='m', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('ter', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Model Terrain Height'
                        wrfNewVar.units     = 'm'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'ter', units='m',meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF 2m Dew Point Temperature has been chosen. 
        if wrfTd2 == True:
            print('Working on WRF 2m Dew Point Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'td2', units='degC', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('td2', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '2m Dew Point Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'td2', units='degC', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'td2', units='degC', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('td2', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = '2m Dew Point Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'td2', units='degC', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 

        # If WRF Landmask has been chosen. 
        if wrfLandmask == True:
            print('Working on WRF Landmask.')
            if selectWrfBox == True:
                    wrfRawVar = getvar(wrfRawFile,'LANDMASK', meta=False)[j0:j1, i0:i1] 
                    wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                    wrfNewVar = wrfNewFile.createVariable('LANDMASK', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                    wrfNewVar.long_name = 'Land mask'
                    wrfNewVar.units     = '1=land, 0=water'
                    wrfNewVar[:,:] = wrfRawVar                         
            else: 
                    wrfRawVar = getvar(wrfRawFile,'tLANDMASKd2', meta=False)[:, :] 
                    wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                    wrfNewVar = wrfNewFile.createVariable('LANDMASK', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                    wrfNewVar.long_name = 'Land mask'
                    wrfNewVar.units     = '1=land, 0=water'
                    wrfNewVar[:,:] = wrfRawVar  

        # If WRF Sea Surface Temperature has been chosen.                                               
        if wrfSST == True:
            print('Working on WRF Sea Surface Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectWrfBox == True:
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'SST', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('SST', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sea Surface Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'SST', meta=False, timeidx=ntimes[0]+i)[j0:j1, i0:i1] 
                        wrfNewVar[i,:,:] = wrfRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        wrfRawVar = getvar(wrfRawFile,'SST', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar = np.zeros([len(ntimes),len(lat_wrf[:,0]), len(lon_wrf[0,:])])
                        wrfNewVar = wrfNewFile.createVariable('SST', 'f', ('Time', 'south_north', 'west_east'), fill_value=wrfFillVal)
                        wrfNewVar.long_name = 'Sea Surface Temperature'
                        wrfNewVar.units     = 'Degree Celsius'
                        wrfNewVar[i,:,:] = wrfRawVar  
                    else:
                        wrfRawVar = getvar(wrfRawFile,'SST', meta=False, timeidx=ntimes[0]+i)[:, :] 
                        wrfNewVar[i,:,:] = wrfRawVar                     
                bar.next()
            bar.finish() 