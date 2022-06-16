import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path

def shellin(A0,n):
# Input is ocean mask. The function removes coastal or near ice zone.
    for i in range (n):    
        A0 = shelli(A0)  
    return(A0)

def shelli(A1):
    shp = np.shape(A1)  
    A2 = np.logical_not(A1)
    A = A2.astype(int)
    ww = 0

    for i1 in range(1,shp[0]-1):
        for j1 in range(1,shp[1]-1):  
            M = A[i1-1:i1+2,j1-1:j1+2]           
            if any(q==1 for q in M.flatten()) and not(M[1,1]==1):
                A[i1,j1]=2
    A1 = A.astype(bool)                            
    A = np.logical_not(A1)
    return(A)