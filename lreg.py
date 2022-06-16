import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path

def lreg(tanth,sigpdf,nangles,nsamples,filenm):
    nangles = 4
    nsamples = 4
    # data for regression 
    if np.all(sigpdf>0):
        y = np.log(sigpdf)  
        x = np.power(tanth,2)  
        f = np.zeros(2)    

        larr = len(x)    
        # count different angles, difference up to
        lsetinc = len(set(np.around(x,decimals=3)))
        
        if (lsetinc < nangles) or (larr/lsetinc < nsamples):            
            S0 = 99999
            sxx = 99999
            err_sxx = 99999
            err_S0 = 99999
            nn = 0
            corr = 0
            scoreh = 0
        else:        
            modelh = HuberRegressor().fit(x.reshape((-1, 1)), y)        
            scoreh = modelh.score(x.reshape((-1, 1)), y)
                    
            coef = modelh.coef_
            inter = modelh.intercept_       
            
            sxx = -1 / modelh.coef_ * 0.5 
            S0 = np.exp(modelh.intercept_)  

            # -----errors------------------------------
            not_out = np.logical_not(modelh.outliers_)
            x = x[not_out]
            y = y[not_out]
            y1 = x * modelh.coef_[0] + modelh.intercept_             
            Dlt = (y - y1)
            N = len(x)
            
            if not(N==0):
                s2 = np.sum(np.square(Dlt))/(N-2)
                xav = np.average(x)
                varpx = np.sum(np.square(x-xav))        
                var_coef = s2/varpx
                err_coef = np.sqrt(var_coef)

                var_inter = s2 * (1/N + xav**2/varpx)
                err_inter = np.sqrt(var_inter)

                err_sxx = err_coef / modelh.coef_[0]**2 * 0.5        
                err_S0 = err_inter * S0
                
                not_out = np.logical_not(modelh.outliers_)
                corr = np.corrcoef(x, y)[0,1] 
                scoreh = modelh.score(x.reshape((-1, 1)), y)
            else:
                sxx = 999999
                S0 = 999999
                corr = 999999
                err_sxx = 999999
                err_S0 = 999999
                scoreh = 999999

    else:
        print('in lreg -9999 ' + filenm)
        sxx = 999999
        S0 = 999999
        corr = 999999
        err_sxx = 999999
        err_S0 = 999999
        scoreh = 999999
    return{'sxx':sxx,'S0':S0,'corr':corr,'err_sxx':err_sxx, 'err_S0':err_S0, 'r2': scoreh}
