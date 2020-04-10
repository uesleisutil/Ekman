"""
File name:      main.py
Author:         Ueslei Adriano Sutil
Email:          uesleisutil1@gmail.com
Created:        07 Apr 2019
Last modified:  08 Apr 2019
Version:        1.1
Python version: 3.3+

TODO
      - ROMS: Insert zeta, h, sensible and latent variables.
      - SeaIce: Comming soon.
      - WRF: Comming soon.
      - WW3: Comming soon.
"""

# Import libraries.

import os
from   configureProject import *
from   romsSection import *
from   setOptions import *


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
if selectRomsVars == True:
    romsOriDir = getFolder+'/'+projectName+'/netcdf_files/ROMS/'+romsOriginalFilename
    romsNewDir = getFolder+'/'+projectName+'/netcdf_files/ROMS/'+romsNewFilename
    romsVars(romsOriDir,romsNewDir)

# Program finished.
print('Program finished.')