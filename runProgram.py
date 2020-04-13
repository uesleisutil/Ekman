"""
File name:      main.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        07 Apr 2019
Last modified:  13 Apr 2019
Version:        1.2
Python version: 3.3+

This file runs the Ekman Toolbox.
WARNING: Do not change anything in this file.
"""

# Import libraries.
import os
from   configureProject import *
from   setOptions import *
from   romsSection import *
from   iceSection import *
 
# Check if project already exists. If not, create a new one.
checkFolder = os.path.isdir('./'+projectName) 
getFolder   = os.getcwd()
print('Initializing project: '+projectName)

if checkFolder == False:
    createFolders('./'+projectName)
    print('Project folders created. Now store raw WRF and ROMS output files inside the respective original file folder in: '+getFolder+'/'+projectName+'/netcdf_files')
    raise SystemExit(0)
else:
    print('Project folder already exists. Continuing.')
    pass

# ROMS section.
if selectRomsVars or selectRomsBox or selectRomsLevel == True:
    print('ROMS section activated.')
    romsOriDir = getFolder+'/'+projectName+'/netcdf_files/'+romsOriginalFilename
    romsNewDir = getFolder+'/'+projectName+'/netcdf_files/'+romsNewFilename
    romsVars(romsOriDir,romsNewDir)
if selectIceVars == True:
    print('Sea-Ice section activated.')
    iceOriDir = getFolder+'/'+projectName+'/netcdf_files/'+iceOriginalFilename
    iceNewDir = getFolder+'/'+projectName+'/netcdf_files/'+iceNewFilename
    iceVars(iceOriDir,iceNewDir)
# Program finished.
print('Program finished.')