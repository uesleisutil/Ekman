"""
File name:      runScript.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        07 Apr 2020
Last modified:  28 Dec 2020
Version:        1.2
Python version: 3.3+

TODO:
      - WW3: Comming soon.
      - Select more files, not just one.

"""

# Import libraries.
import os
from configureProject import *
from setOptions import *

def run():
    # Check if project already exists. If not, create a new one.
    checkFolder = os.path.isdir('./Projects/'+projectName) 
    getFolder   = os.getcwd()
    print('Initializing project: '+projectName)

    if checkFolder == False:
        createFolders('./'+projectName)
        print('Project folders created. Now store raw WRF and ROMS output files inside the respective original file folder in: '+getFolder+'/'+projectName+'/netcdf_files')
        raise SystemExit(0)
    else:
        print('Project folder already exists. Continuing.')
        pass

    if selectRomsVars == True:
        from romsSection import romsVars
        print('------------------------')        
        print('ROMS section initialized')
        print('------------------------')  
        romsOriDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+romsOriginalFilename
        romsNewDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+romsNewFilename
        romsVars(romsOriDir,romsNewDir)
    if selectIceVars == True:
        from iceSection import iceVars
        print('---------------------------')  
        print('Sea-Ice section initialized')
        print('---------------------------')  
        iceOriDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+iceOriginalFilename
        iceNewDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+iceNewFilename
        iceVars(iceOriDir,iceNewDir)
    if selectWrfVars == True:
        from wrfSection import wrfVars
        print('-----------------------')  
        print('WRF section initialized')
        print('-----------------------')  
        wrfOriDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+wrfOriginalFilename
        wrfNewDir = getFolder+'/Projects/'+projectName+'/netcdf_files/'+wrfNewFilename
        wrfVars(wrfOriDir,wrfNewDir)        
        
    # Program finished.
    print('Program finished.')

run()