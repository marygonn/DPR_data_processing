import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path

def smooth(ocean,swath,nf,filter):
    nvals = int((nf*nf-1)*0.5)
    dw = int((nf-1)*0.5)

    szs = np.shape(swath)
    swath1ij = swath
    
    for i in range(dw, szs[0]-dw):
        for j in range(dw, szs[1]-dw):            
            arr = [] 
            k = 0         
            if ocean[i,j]:
                for im in range(nf):
                    for jm in range(nf):
                        q = swath[i-dw+im,j-dw+jm]
                        if q < 9999:
                            arr.append(q)
                            k = k +1                  
                if k>=nvals:
                    if filter=='median':
                        swath1ij[i,j]=np.median(np.array(arr))
                    if filter=='mean':                        
                        swath1ij[i,j]=np.mean(np.array(arr))
    return(swath1ij)   