import multiprocessing, logging
logger = multiprocessing.log_to_stderr()
logger.setLevel(multiprocessing.SUBDEBUG)

import math
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.path
import h5py
import os
import json
import fnmatch
from multiprocessing import Pool


from bandproc import bandproc
from iceconc import asi_hdf

import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import cartopy.crs as ccrs

print('hello!')
with open("config_dpr.json", "r") as read_file:
    data = json.load(read_file)
    
if not os.path.exists(data['out_folder']):
    os.mkdir(data['out_folder'])
    
fnsDPR = os.listdir(data['inp_folder_DPR'])
fnsGMI = os.listdir(data['inp_folder_GMI'])    

fns1_DPR = []
for entry in fnsDPR:  
    if fnmatch.fnmatch(entry, data['part_name_DPR']):
        fns1_DPR.append(entry)    

fns1_DPR = []
fns1_GMI = []
for entryDPR in fnsDPR:                  
    for entryGMI in fnsGMI:
        if fnmatch.fnmatch(entryGMI, data['part_name_GMI']):
            if fnmatch.fnmatch(entryDPR, data['part_name_DPR']): 
                print(entryDPR)
                print(entryGMI)
                if entryDPR.split('.')[5] == entryGMI.split('.')[5]:
                    fns1_DPR.append(entryDPR)
                    if data['ice']:
                        fns1_GMI.append(entryGMI)                        
                    else:
                        fns1_GMI.append('')    

    q = np.empty([2,len(fns1_DPR)], dtype=object)
    print(fns1_DPR)
    q[0,:] = fns1_DPR
    q[1,:] = fns1_GMI

    fns = list(np.transpose(q))   
with open("config_dpr.json", "r") as read_file:
    data = json.load(read_file)   
fg = h5py.File(os.path.join(data['inp_folder_GMI'],fnsGMI[0]))
print(fg['S1']['SCstatus'].keys())   