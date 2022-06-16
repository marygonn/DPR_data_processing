import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path
import sys
import copy

def kurt_whole(tanth, sigdb, sigpdf, filenm):     
    nal = np.shape(tanth)[0]
    nsc = np.shape(tanth)[1]

    A1 = np.zeros(nsc)
    A2 = np.zeros(nsc)
    kurtosis = np.zeros(nal)

    outp = []
    with open('kout.txt','w') as f:
        #sys.stdout = f
        for j in range (nal):
            if np.logical_not(np.any(sigdb[j,:]<=-9999.)):
                for i in range (nsc):
                    A1[i] = tanth[j,i]
                    if i < nsc/2:
                        A1[i] = -tanth[j,i]       
                    A2[i] = sigpdf[j,i]                  

                S1 = np.sum(np.multiply(A1,A2))
                SS0 = np.sum(A2)
                    
                with np.errstate(divide='raise'):
                    try:
                        S1/SS0   # this gets caught and handled as an exception
                    except FloatingPointError:
                        print('in file '+filenm + ', in str '+str(j))

                if SS0==0: 
                    kurtosis[j] = -9999
                else:            
                    fcm1=S1/SS0
                    
                    S200=np.sum(np.multiply(A2,np.power(A1-fcm1,2)))
                    S4=np.sum(np.multiply(A2,np.power(A1-fcm1,4)))            
                    M4 = S4/SS0

                    if S200/SS0>=0:
                        Bd = np.sqrt(S200/SS0)
                    else:
                        kurtosis[j]=-9999
                    if Bd==0:
                        kurtosis[j]=-9999
                    else:
                        kurtosis[j] = M4/Bd**4 - 3. 
            else:
                kurtosis[j] = -9999 
                print('in kurt -9999'+ filenm)
    
    #sys.stdout = original_stdout 
    #np.savetxt('kurt_out.txt',np.array(outp))               
    return(kurtosis)


def kurt_twosides(tanth, sigdb, sigpdf, filenm):
    sigpdf_rev = np.fliplr(sigpdf)
    A1 = copy.deepcopy(sigpdf)
    nsc1 = len(tanth[0])
    nal = len(tanth)
    half = int(np.floor(nsc1/2))
    
    #A1[:,0:half] = sigpdf[:,0:half]
    A1[:,half:nsc1] = sigpdf_rev[:,half:nsc1]
         
    A2 = copy.deepcopy(sigpdf)
    A2[:,0:half] = sigpdf_rev[:,0:half]
    #A2[:,half:nsc1] = sigpdf[:,half:nsc1] 
    
    excess_left = kurt_whole(tanth, sigdb, A1, filenm)
    excess_right = kurt_whole(tanth, sigdb, A2, filenm)
   
    exc_swath = np.zeros((nal,2))
    
    exc_swath[:,0] = excess_left
    exc_swath[:,1] = excess_right      
    
    return exc_swath

def kurt(tanth, sigdb, sigpdf, filenm):     

    nal = np.shape(tanth)[0]
    nsc = np.shape(tanth)[1]

    A1 = np.zeros(nsc)
    A2 = np.zeros(nsc)
    kurtosis = np.zeros(nal)

    outp = []
    with open('kout.txt','w') as f:
        #sys.stdout = f
        for j in range (nal):
            if np.logical_not(np.any(sigdb[j,:]<=-9999.)):
                for i in range (nsc):
                    A1[i] = tanth[j,i]
                    if i < nsc/2:
                        A1[i] = -tanth[j,i]       
                    A2[i] = sigpdf[j,i]                  

                S1 = np.sum(np.multiply(A1,A2))
                SS0 = np.sum(A2)
                    
                with np.errstate(divide='raise'):
                    try:
                        S1/SS0   # this gets caught and handled as an exception
                    except FloatingPointError:
                        print('in file '+filenm + ', in str '+str(j))

                if SS0==0: 
                    kurtosis[j] = -9999
                else:            
                    fcm1=S1/SS0
                    
                    S200=np.sum(np.multiply(A2,np.power(A1-fcm1,2)))
                    S4=np.sum(np.multiply(A2,np.power(A1-fcm1,4)))            
                    M4 = S4/SS0

                    if S200/SS0>=0:
                        Bd = np.sqrt(S200/SS0)
                    else:
                        kurtosis[j]=-9999
                    if Bd==0:
                        kurtosis[j]=-9999
                    else:
                        kurtosis[j] = M4/Bd**4 - 3. 
            else:
                kurtosis[j] = -9999 
                print('in kurt -9999'+ filenm)
    
    #sys.stdout = original_stdout 
    #np.savetxt('kurt_out.txt',np.array(outp))               
    return(kurtosis)
