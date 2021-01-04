import os

def createFolders(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            os.makedirs(directory+'/netcdf_files')       
    except OSError:
        print ('Error: Creating project folder: ' + directory)
