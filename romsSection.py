"""
Author:         Ueslei Adriano Sutil
Created:        08 Apr 2020
Last modified:  06 Jan 2021
Version:        2.12

This file generates a new ROMS output file from scratch.
It is netCDF4 CF-compliant.

WARNING: Do not change anything in this file.
"""

from   netCDF4      import Dataset
from   setOptions   import *
from   matplotlib   import path 
from   progress.bar import IncrementalBar
import numpy        as     np
import time

if romsSST or romsTemp or romsSalt or romsZeta or romsTKE or romsLatent or romsSensible or romsLWRad or romsSWRad or romsEvaporation or romsEminusP or romsUwind or romsVwind or romsW or romsOmega or romsRho == True:
    romsMassPoints = True
else:
    romsMassPoints = False   
if romsU or romsV or romsUbar or romsVbar == True:
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
    romsNewFile.description = "Created with Ekman Toolbox in " + time.ctime(time.time())
    romsNewFile.link        = "https://github.com/uesleisutil/Ekman"

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
            print("Bounding box selected. New domain limits are: Longitude "+str(romsBox[0])+"/"+str(romsBox[1])+" and Latitude "+str(romsBox[2])+"/"+str(romsBox[3])+".")
        else:            
            print("No bounding box selected: Using XLAT and XLONG variables from input file.")
            lon_rho = romsRawFile.variables['lon_rho'][:,:]
            lat_rho = romsRawFile.variables['lat_rho'][:,:] 
            eta_rho = romsRawFile.dimensions['eta_rho']
            xi_rho  = romsRawFile.dimensions['xi_rho']
            romsNewFile.createDimension('eta_rho', len(eta_rho))    
            romsNewFile.createDimension('xi_rho', len(xi_rho))  
        if selectRomsLevel == True:
            romsNewFile.createDimension('s_rho', len(romsLevel))
        else:
            romsNewFile.createDimension('s_rho', len(s_rho))

        romsNewFile.createDimension('s_w', len(s_w))

        romsNewLon               = romsNewFile.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=romsFillVal)
        romsNewLon.long_name     = 'Longitude on RHO-points'
        romsNewLon.units         = 'degree_east'
        romsNewLon.standard_name = 'longitude'
        romsNewLon[:,:]          = lon_rho

        romsNewLat               = romsNewFile.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=romsFillVal)
        romsNewLat.long_name     = 'Latitude on RHO-points'
        romsNewLat.units         = 'degree_north'
        romsNewLat.standard_name = 'latitude'
        romsNewLat[:, :]         = lat_rho

        # Define vertical levels and time-steps. 
        levels = len(romsRawFile.variables['s_rho'][:]) 
        if selectRomsLevel == True and len(romsLevel) == 1:
            print("One vertical level selected: Working on vertical level "+str(romsLevel)+".")
        if selectRomsLevel == True and len(romsLevel) > 1:
            print("Multiple vertical levels selected: Working from level "+str(romsLevel[0])+" to "+str(romsLevel[-1])+".")
        if selectRomsLevel == False:
            print("No selected vertical levels specified: Using entire vertical level from input file.")
        if selectRomsTimeStep == True:
            ntimes = romsTimeStep
            romsNewFile.createDimension('ocean_time', 0)
            print("Time-step selected: Working from time-step "+str(ntimes[0])+" to "+str(ntimes[-1])+".")
        else:
            ntimes = romsRawFile.variables['ocean_time'][:]
            ntimes = np.arange(np.argmin(ntimes), len(ntimes)) 
            romsNewFile.createDimension('ocean_time', 0)
            print("No time-step selected. Working with entire time-step.")

        # If ROMS Sea Surface Temperature has been chosen.                                               
        if romsSST == True:
            print('Working on ROMS Sea Surface Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,-1,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('sst', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Sea Surface Temperature'
                        romsNewVar.units     = 'Degree Celsius'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,-1,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,-1,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('sst', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Sea Surface Temperature'
                        romsNewVar.units     = 'Degree Celsius'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,-1,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 


        # If ROMS Potential Temperature has been chosen.
        if romsTemp == True:
            print('Working on ROMS Potential Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Potential Temperature'
                            romsNewVar.units     = 'Degree Celsius'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['temp'][i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Potential Temperature'
                            romsNewVar.units     = 'Degree Celsius'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Potential Temperature'
                        romsNewVar.units     = 'Degree Celsius'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Potential Temperature'
                        romsNewVar.units     = 'Degree Celsius'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Potential Temperature'
                            romsNewVar.units     = 'Degree Celsius'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('temp', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Potential Temperature'
                            romsNewVar.units     = 'Degree Celsius'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['temp'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  

        # If ROMS Salinity has been chosen.
        if romsSalt == True:
            print('Working on ROMS Salinity.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Salinity'
                            romsNewVar.units     = 'PSU'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Salinity'
                            romsNewVar.units     = 'PSU'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Salinity'
                        romsNewVar.units     = 'PSU'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Salinity'
                        romsNewVar.units     = 'PSU'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Salinity'
                            romsNewVar.units     = 'PSU'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('salt', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Salinity'
                            romsNewVar.units     = 'PSU'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['salt'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()    

        # If ROMS Turbulent Kinectic Energy has been chosen.
        if romsTKE == True:
            print('Working on ROMS Turbulent Kinectic Energy.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Turbulent Kinectic Energy'
                            romsNewVar.units     = 'm2 s-2'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Turbulent Kinectic Energy'
                            romsNewVar.units     = 'm2 s-2'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Turbulent Kinectic Energy'
                        romsNewVar.units     = 'm2 s-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Turbulent Kinectic Energy'
                        romsNewVar.units     = 'm2 s-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Turbulent Kinectic Energy'
                            romsNewVar.units     = 'm2 s-2'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('tke', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Turbulent Kinectic Energy'
                            romsNewVar.units     = 'm2 s-2'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['tke'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish() 

        # If ROMS Density Anomaly has been chosen.
        if romsRho == True:
            print('Working on ROMS Density Anomaly.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Density Anomaly'
                            romsNewVar.units     = 'kilogram meter-3'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Density Anomaly'
                            romsNewVar.units     = 'kilogram meter-3'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Density Anomaly'
                        romsNewVar.units     = 'kilogram meter-3'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Density Anomaly'
                        romsNewVar.units     = 'kilogram meter-3'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Density Anomaly'
                            romsNewVar.units     = 'kilogram meter-3'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('rho', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Density Anomaly'
                            romsNewVar.units     = 'kilogram meter-3'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['rho'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish() 

        # If ROMS Vertical Momentum Component has been chosen.
        if romsW == True:
            print('Working on ROMS Vertical Momentum Component.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertical Momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertical Momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['w'][i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('w', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['w'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  

        # If ROMS S-coordinate Vertical Momentum Component has been chosen.
        if romsOmega == True:
            print('Working on ROMS S-coordinate Vertical Momentum Component.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsLevel,j0:j1, i0:i1]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsStart:romsStop,j0:j1, i0:i1] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['omega'][i,romsStart:romsStop,j0:j1, i0:i1]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,:,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,:,j0:j1, i0:i1] 
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                            romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_rho[:,0]), len(lon_rho[0,:])])
                            romsNewVar = romsNewFile.createVariable('omega', 'f', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'S-coordinate Vertical Momentum Component'
                            romsNewVar.units     = 'm s-1'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['omega'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  

        # If ROMS Free-surface has been chosen.                                               
        if romsZeta == True:
            print('Working on ROMS Free-surface.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['zeta'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Free-surface'
                        romsNewVar.units     = 'meters'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['zeta'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['zeta'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('zeta', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Free-surface'
                        romsNewVar.units     = 'meters'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['zeta'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                  
                bar.next()
            bar.finish()  

        # If ROMS Latent Heat Flux has been chosen.                                               
        if romsLatent == True:
            print('Working on ROMS Latent Heat Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['latent'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('latent', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Latent Heat Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['latent'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['latent'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('latent', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Latent Heat Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['latent'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                
                bar.next()
            bar.finish() 

        # If ROMS Sensible Heat Flux has been chosen.                                               
        if romsSensible == True:
            print('Working on ROMS Sensible Heat Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['sensible'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('sensible', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Sensible Heat Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar.negative_value = "Upward flux = Cooling"
                        romsNewVar.positive_value = "Fownward flux = Heating"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['sensible'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['sensible'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('sensible', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Sensible Heat Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar.negative_value = "Upward flux = Cooling"
                        romsNewVar.positive_value = "Downward flux = Heating"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['sensible'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                  
                bar.next()
            bar.finish() 

        # If ROMS Net Longwave Radiation Flux has been chosen.                                               
        if romsLWRad == True:
            print('Working on ROMS Net Longwave Radiation Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['lwrad'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('lwrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Net Longwave Radiation Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['lwrad'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['lwrad'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('lwrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Net Longwave Radiation Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['lwrad'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                   
                bar.next()
            bar.finish() 

        # If ROMS Net Shortwave Radiation Flux has been chosen.                                               
        if romsSWRad == True:
            print('Working on ROMS Net Shortwave Radiation Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['swrad'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('swrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Net Shortwave Radiation Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['swrad'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['swrad'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('swrad', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Net Shortwave Radiation Flux'
                        romsNewVar.units     = 'W m-2'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['swrad'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 

        # If ROMS Bulk Flux Surface Net Freshwater Flux has been chosen.                                               
        if romsEminusP == True:
            print('Working on ROMS Bulk Flux Surface Net Freshwater Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['EminusP'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('EminusP', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Bulk Flux Surface Net Freshwater Flux'
                        romsNewVar.units     = 'meter s-1'
                        romsNewVar.negative_value = "Upward = Freshening (Net Precipitation)"
                        romsNewVar.positive_value = "Downward = Salting (Net Evaporation)"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['EminusP'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['EminusP'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('EminusP', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Bulk Flux Surface Net Freshwater Flux'
                        romsNewVar.units     = 'meter s-1'
                        romsNewVar.negative_value = "Upward = Freshening (Net Precipitation)"
                        romsNewVar.positive_value = "Downward = Salting (Net Evaporation)"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['EminusP'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 

        # If ROMS Evaporation Rate has been chosen.                                               
        if romsEvaporation == True:
            print('Working on ROMS Evaporation Rate.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['evaporation'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('evaporation', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Evaporation Rate'
                        romsNewVar.units     = 'Kg m-2 s-1'
                        romsNewVar.negative_value = "Downward = Freshening (Condensation)"
                        romsNewVar.positive_value = "Upward = Salting (Evaporation)"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['evaporation'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['evaporation'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('evaporation', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Evaporation Rate'
                        romsNewVar.units     = 'Kg m-2 s-1'
                        romsNewVar.negative_value = "Downward = Freshening (Condensation)"
                        romsNewVar.positive_value = "Upward = Salting (Evaporation)"
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['evaporation'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 

        # If ROMS U-wind Component has been chosen.                                               
        if romsUwind == True:
            print('Working on ROMS U-wind Component.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['Uwind'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('Uwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Surface U-wind Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['Uwind'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['Uwind'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('Uwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Surface U-wind Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['Uwind'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 

        # If ROMS V-wind Component has been chosen.                                               
        if romsVwind == True:
            print('Working on ROMS V-wind Component.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['Vwind'][ntimes[0]+i,j0:j1, i0:i1]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('Vwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Surface V-wind Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['Vwind'][ntimes[0]+i,j0:j1,i0:i1]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['Vwind'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        romsNewVar = romsNewFile.createVariable('Vwind', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Surface V-wind Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['Vwind'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 
   
    elif romsUVPoints == True:     
        if selectRomsBox == True:
            if romsU or romsUbar == True:
                lon_u = romsRawFile.variables['lon_u'][:, :]
                lat_u = romsRawFile.variables['lat_u'][:, :]
                i0_u,i1_u,j0_u,j1_u = bbox2ij(lon_u,lat_u,romsBox)
                lon_u = romsRawFile.variables['lon_u'][j0_u:j1_u, i0_u:i1_u]
                lat_u = romsRawFile.variables['lat_u'][j0_u:j1_u, i0_u:i1_u]
                romsNewFile.createDimension('eta_u', len(lon_u[:,0]))    
                romsNewFile.createDimension('xi_u', len(lon_u[0,:])) 
            if romsV or romsVbar == True:
                lon_v = romsRawFile.variables['lon_v'][:, :]
                lat_v = romsRawFile.variables['lat_v'][:, :]
                i0_v,i1_v,j0_v,j1_v = bbox2ij(lon_v,lat_v,romsBox)
                lon_v = romsRawFile.variables['lon_v'][j0_v:j1_v, i0_v:i1_v]
                lat_v = romsRawFile.variables['lat_v'][j0_v:j1_v, i0_v:i1_v]
                romsNewFile.createDimension('eta_v', len(lon_v[:,0]))    
                romsNewFile.createDimension('xi_v', len(lon_v[0,:])) 
            print("Bounding box selected. New domain limits are: Longitude "+str(romsBox[0])+"/"+str(romsBox[1])+" and Latitude "+str(romsBox[2])+"/"+str(romsBox[3])+".")
        else: 
            print("No bounding box selected: Using XLAT and XLONG variables from input file.")
            if romsU or romsUbar == True:           
                eta_u = romsRawFile.dimensions['eta_u']
                xi_u  = romsRawFile.dimensions['xi_u']
                lon_u = romsRawFile.variables['lon_u'][:,:]
                lat_u = romsRawFile.variables['lat_u'][:,:] 
                romsNewFile.createDimension('eta_u', len(eta_u))    
                romsNewFile.createDimension('xi_u', len(xi_u)) 
            if romsV or romsVbar == True:
                eta_v = romsRawFile.dimensions['eta_v']
                xi_v  = romsRawFile.dimensions['xi_v']
                lon_v = romsRawFile.variables['lon_v'][:,:]
                lat_v = romsRawFile.variables['lat_v'][:,:] 
                romsNewFile.createDimension('eta_v', len(eta_v))    
                romsNewFile.createDimension('xi_v', len(xi_v)) 

        if selectRomsLevel == True:
            romsNewFile.createDimension('s_rho', len(romsLevel))
        else:
            s_rho = romsRawFile.dimensions['s_rho']
            romsNewFile.createDimension('s_rho', len(s_rho))

        # Define vertical levels and time-steps. 
        levels = len(romsRawFile.variables['s_rho'][:])        
        if selectRomsTimeStep == True:
            ntimes = romsTimeStep
            print("Time-step selected: Working from time-step "+str(np.argmin(ntimes))+" to "+str(np.argmax(ntimes))+".")
        else:
            ntimes = romsRawFile.variables['ocean_time'][:]
            print("No time-step selected. Working with entire time-step.")      
        if selectRomsLevel and len(romsLevel) == 1 and romsU or romsV == True:
            print("One vertical level selected: Working on level "+str(romsLevel)+".")
        if selectRomsLevel and len(romsLevel) > 1 and romsU or romsV == True:
            print("Multiple vertical levels selected: Working from level "+str(romsLevel[0])+" to "+str(romsLevel[-1])+".")
        if selectRomsLevel == False and romsU or romsV == True:
            print("No selected vertical levels specified: Using entire vertical level from input file.")

        s_w     = romsRawFile.dimensions['s_w']      
        romsNewFile.createDimension('s_w', len(s_w))

        # Create lat and lon variables.
        if romsU or romsUbar == True:
            romsNewLonU               = romsNewFile.createVariable('lon_u', 'd', ('eta_u', 'xi_u'), fill_value=romsFillVal)
            romsNewLonU.long_name     = 'Longitude on U-points'
            romsNewLonU.units         = 'degree_east'
            romsNewLonU.standard_name = 'longitude'
            romsNewLonU[:, :]         = lon_u
            romsNewLatU               = romsNewFile.createVariable('lat_u', 'd', ('eta_u', 'xi_u'), fill_value=romsFillVal)
            romsNewLatU.long_name     = 'Latitude on U-points'
            romsNewLatU.units         = 'degree_north'
            romsNewLatU.standard_name = 'latitude'
            romsNewLatU[:, :]         = lat_u
        if romsV or romsVbar == True:
            romsNewLonV               = romsNewFile.createVariable('lon_v', 'd', ('eta_v', 'xi_v'), fill_value=romsFillVal)
            romsNewLonV.long_name     = 'Longitude on V-points'
            romsNewLonV.units         = 'degree_east'
            romsNewLonV.standard_name = 'longitude'
            romsNewLonV[:, :]         = lon_v
            romsNewLatV               = romsNewFile.createVariable('lat_v', 'd', ('eta_v', 'xi_v'), fill_value=romsFillVal)
            romsNewLatV.long_name     = 'Latitude on U-points'
            romsNewLatV.units         = 'degree_north'
            romsNewLatV.standard_name = 'latitude'
            romsNewLatV[:, :]         = lat_v

        # If ROMS V-wind Component has been chosen.
        if romsV == True:
            print('Working on ROMS V-wind Component.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,j0_v:j1_v, i0_v:i1_v]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,j0_v:j1_v, i0_v:i1_v]   
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,j0_v:j1_v, i0_v:i1_v]  
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_v[:,0]), len(lon_v[0,:])])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,j0_v:j1_v, i0_v:i1_v]   
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'V-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,j0_v:j1_v, i0_v:i1_v]  
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'V-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,j0_v:j1_v, i0_v:i1_v]  
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_v[:,0]), len(lon_v[0,:])])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  


        # If ROMS V-wind Component has been chosen.
        if romsU == True:
            print('Working on ROMS V-wind Component.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsLevel,j0_u:j1_u, i0_u:i1_u]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_u), len(lon_u)])
                            romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsLevel,j0_u:j1_u, i0_u:i1_u]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsStart:romsStop,j0_u:j1_u, i0_u:i1_u]  
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_u[:,0]), len(lon_u[0,:])])
                            romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsStart:romsStop,j0_u:j1_u, i0_u:i1_u]  
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_u), len(lon_u)])
                        romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'V-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,:,j0_u:j1_u, i0_u:i1_u]  
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_u), len(lon_u)])
                        romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'V-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,:,j0_u:j1_u, i0_u:i1_u]  
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_u), len(lon_u)])
                            romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_u[:,0]), len(lon_u[0,:])])
                            romsNewVar = romsNewFile.createVariable('u', 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'V-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['u'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  

        # If ROMS U-wind Component has been chosen.
        if romsV == True:
            print('Working on ROMS U-wind Component.')
            bar = IncrementalBar(max=len(ntimes))
            
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True and selectRomsLevel == True:
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,j0_v:j1_v, i0_v:i1_v]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'U-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,j0_v:j1_v, i0_v:i1_v]   
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,j0_v:j1_v, i0_v:i1_v]  
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_v[:,0]), len(lon_v[0,:])])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'U-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,j0_v:j1_v, i0_v:i1_v]   
                            romsNewVar[i,:,:,:] = romsRawVar                              
                elif selectRomsBox == False and selectRomsLevel == False:    
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,:,:]
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'U-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,:,:] 
                        romsNewVar[i,:,:] = romsRawVar                                                                                
                elif selectRomsBox == True and selectRomsLevel == False:  
                    if i == np.argmin(ntimes):             
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,j0_v:j1_v, i0_v:i1_v]  
                        romsNewVar = np.zeros([len(ntimes),levels,len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'U-wind Component'
                        romsNewVar.units     = 'm s'
                        romsNewVar[i,:,:] = romsRawVar  
                    else: 
                        romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,:,j0_v:j1_v, i0_v:i1_v]  
                        romsNewVar[i,:,:] = romsRawVar  
                elif selectRomsBox == False and selectRomsLevel == True: 
                    if len(romsLevel) == 1:
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'U-wind Component'
                            romsNewVar.units     = 'm s'
                            romsNewVar[i,:,:] = romsRawVar  
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsLevel,:, :]  
                            romsNewVar[i,:,:] = romsRawVar                                                             
                    else:
                        romsStart  = slice(min(romsLevel),max(romsLevel)+1).start
                        romsStop   = slice(min(romsLevel),max(romsLevel)+1).stop
                        if i == np.argmin(ntimes):
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,:, :] 
                            romsNewVar = np.zeros([len(ntimes),len(romsLevel),len(lat_v[:,0]), len(lon_v[0,:])])
                            romsNewVar = romsNewFile.createVariable('v', 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=romsFillVal)
                            romsNewVar.long_name = 'U-wind Component'
                            romsNewVar.units     = 'm s'                     
                            romsNewVar[i,:,:,:] = romsRawVar
                        else: 
                            romsRawVar = romsRawFile.variables['v'][ntimes[0]+i,romsStart:romsStop,:, :]  
                            romsNewVar[i,:,:,:] = romsRawVar  
                bar.next()
            bar.finish()  

        # If ROMS Vertically Integrated U-momentum Component has been chosen.                                               
        if romsUbar == True:
            print('Working on ROMS Vertically Integrated U-momentum Component.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['ubar'][ntimes[0]+i,j0_u:j1_u, i0_u:i1_u]
                        romsNewVar = np.zeros([len(ntimes),len(lat_u), len(lon_u)])
                        romsNewVar = romsNewFile.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertically Integrated U-momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['ubar'][ntimes[0]+i,j0_u:j1_u, i0_u:i1_u]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['ubar'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_u), len(lon_u)])
                        romsNewVar = romsNewFile.createVariable('ubar', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertically Integrated U-momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['ubar'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 

        # If ROMS Vertically Integrated V-momentum Component has been chosen.                                               
        if romsVbar == True:
            print('Working on ROMS Vertically Integrated V-momentum Component.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if selectRomsBox == True:
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['vbar'][ntimes[0]+i,j0_v:j1_v, i0_v:i1_v]
                        romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('vbar', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertically Integrated U-momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['vbar'][ntimes[0]+i,j0_v:j1_v, i0_v:i1_v]
                        romsNewVar[i,:,:] = romsRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        romsRawVar = romsRawFile.variables['vbar'][ntimes[0]+i,:,:]
                        romsNewVar = np.zeros([len(ntimes),len(lat_v), len(lon_v)])
                        romsNewVar = romsNewFile.createVariable('vbar', 'f', ('ocean_time', 'eta_u', 'xi_u'), fill_value=romsFillVal)
                        romsNewVar.long_name = 'Vertically Integrated V-momentum Component'
                        romsNewVar.units     = 'm s-1'
                        romsNewVar[i,:,:] = romsRawVar  
                    else:
                        romsRawVar = romsRawFile.variables['vbar'][ntimes[0]+i,:,:]
                        romsNewVar[i,:,:] = romsRawVar                     
                bar.next()
            bar.finish() 