import os

def createFolders(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            os.makedirs(directory+'/netcdf_files')
            os.makedirs(directory+'/netcdf_files/ROMS')         
            os.makedirs(directory+'/netcdf_files/WRF')
            os.makedirs(directory+'/figures/evaluation/WRF')
            os.makedirs(directory+'/figures/evaluation/ROMS')
            os.makedirs(directory+'/analysis/ROMS')
            os.makedirs(directory+'/analysis/WRF')            
    except OSError:
        print ('Error: Creating project folder: ' + directory)
