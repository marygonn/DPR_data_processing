import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path
from pyresample import kd_tree, geometry, image

from shellin import shellin
from swathproc import swathproc
from kurt import kurt_twosides

import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
#from mpl_toolkits.axes_grid1.colorbar import colorbar
import cartopy.crs as ccrs
LAND = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                    edgecolor='face',
                                    facecolor=cfeature.COLORS['land'])

def make_map(bbox, projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(8, 6),
                           subplot_kw=dict(projection=projection))
    ax.set_extent(bbox)
    #ax.add_feature(LAND, facecolor='0.75')
    ax.coastlines(resolution='50m')
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax

def bandproc(d,data,filenm):
    print(filenm)
    # ------ice_conc resampling -----------   
    if data['ice']:
        infl_radius = 15000
        ku_swath_def = geometry.SwathDefinition(lons = d['Lon'], lats = d['Lat'])
        gmi_swath_def = geometry.SwathDefinition(lons = d['LonGMI'], lats = d['LatGMI'])
        ice_conc_dpr = kd_tree.resample_nearest(gmi_swath_def,  d['ice_conc_GMI'],
        ku_swath_def, fill_value = 9999, radius_of_influence = infl_radius, epsilon = 0.5, nprocs=1) 
    
    else:
        ice_conc_dpr = []
    # ---------- data preparing-------------
    lo = d['Lon']
    la = d['Lat']
    th = d['IncAngle']
    sigdb = d['sigma']
    prec = d['precipRate'] 
    surf = d['surftype']

    # --------------- masks ----------------------
    if data['ice']:
        icenan = np.zeros(ice_conc_dpr.shape)
        icenan[ice_conc_dpr==9999] = 1    
        ice_dpr = ice_conc_dpr>0.001
        noice = (np.logical_not(ice_dpr))
    else:
        noice = np.ones(lo.shape)

    norain = prec < 1  
    
    noland = np.logical_and(surf>=0,surf<=99)    
    noland_internal = np.logical_and(surf>=300,surf<=399) 
    noland = np.logical_or(noland,noland_internal)
    noland = np.logical_and(noland,sigdb>-9999.9)
    #landice_pluscoast = shellin(np.logical_and(noland, noice),1)

    landice_pluscoast = np.logical_and(noland, noice)
    ocean = np.logical_and(norain,landice_pluscoast)
    masktot = np.logical_and(ocean, th < data['max_inc'])

    mdfl = int((data['dfl']-1)*0.5)
    masktot[0:mdfl,:] = False
    masktot[-mdfl:,:] = False   
    
    mask_area_lo = np.logical_and(lo>=data['lon_min'],lo<=data['lon_max'])
    mask_area_la = np.logical_and(la>=data['lat_min'],la<=data['lat_max'])
    mask_area = np.logical_and(mask_area_la,mask_area_lo)

    masktot = np.logical_and(masktot,mask_area)   
    '''
    fig,ax = make_map([-180,180, -90,90])
    ax.scatter(lo[masktot],la[masktot],c='r',s = 1)
    plt.savefig('masktot.png')    
    print('f')

    fig,ax = make_map([-180,180, -90,90])
    ax.scatter(lo[mask_area],la[mask_area],c='r',s = 1)
    plt.savefig('mask_area.png')    

    fig,ax = make_map([-180,180, -90,90])
    ax.scatter(lo[ice_dpr],la[ice_dpr],c='r',s = 1)
    plt.savefig('ice.png')    
    '''
    
    indmask = np.argwhere(masktot)
    sort_ind = np.sort(indmask[:,0]) 
    
    if len(sort_ind)>=1:
        imin = sort_ind[0]
        imax = sort_ind[-1]

        masktot[0:imin,:] = False
        masktot[imax:,:] = False
        #indmask = np.argwhere(masktot)
        # ^ TODO: move the mask constructing code to a separate function

        sig = np.power(10,sigdb*0.1)
        if data['mss']:
            # mss        
            nangles = 3
            nsamples = 3
            swath = swathproc(th,sig,masktot,data,nangles,nsamples,filenm)      
        else:
            swath = []

        if data['kurt']:
            print('in kurt')
            tanth = np.tan(th/180*math.pi)
            sigpdf = np.multiply(sig,(np.power(np.cos(th/180*math.pi),4)))        
            kurtosis = kurt_twosides(tanth, sigdb, sigpdf,filenm) 
            print(np.shape(kurtosis))
        else:
            kurtosis = []   
    else:
        print('file error, sort_ind ' + filenm)        
        swath = {'sxx_all':[], 'S0_all': [],'sxx_filtered':[], 'S0_filtered':[], 'err_sxx': [], 'err_S0': [], 'r2':[], 'corr': []}        
        
        kurtosis = []
        ice_conc_dpr = []
        
    rez = {'swath': swath, 'kurt':kurtosis, 'conc': ice_conc_dpr}
    return(rez)    
