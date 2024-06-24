import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import math
import matplotlib.pyplot as plt
import matplotlib.path
import cartopy.crs as ccrs
from shapely.geometry import Point, Polygon, LineString
import h5py
import matplotlib.path

from lreg import lreg
from smooth import smooth

def swathproc(th,sig,mask,data,nangles,nsamples,filenm):
  dfl = data['dfl']
  dsc = data['dsc']

  # size of the array
  nfl = len(th)
  nsc = len(th[0])

  mdfl = int((dfl-1)*0.5)
  mdsc = int((dsc-1)*0.5)

  # windows
  boxa = np.zeros((dfl,dsc))
  boxb = np.zeros((dfl,dsc))

  # the middle part boundaries for arrays A and B, with bad data at low theta
  bleft = int((nsc - 1) / 2 - 1) 
  bright = int(bleft + 3)

  sxx_all = np.ones((nfl,nsc)) * 9999
  S0_all = np.ones((nfl,nsc)) * 9999
  sxx_filtered = np.ones((nfl,nsc)) * 9999
  S0_filtered = np.ones((nfl,nsc)) * 9999

  corr = np.ones((nfl,nsc)) * 9999
  err_sxx = np.ones((nfl,nsc)) * 9999
  err_S0 = np.ones((nfl,nsc)) * 9999
  r2 = np.ones((nfl,nsc)) * 9999

  indmask = np.argwhere(mask) 

  for item in indmask:
    bc_i = item[0]
    bc_j = item[1]

    thc = th[bc_i,bc_j]
    boxa = th[bc_i-mdfl:bc_i+mdfl+1,bc_j-mdsc:bc_j+mdsc+1]
    # ^ angles in the box
    boxb = sig[bc_i-mdfl:bc_i+mdfl+1,bc_j-mdsc:bc_j+mdsc+1]
    # ^ signal values in the box
    boxm = mask[bc_i-mdfl:bc_i+mdfl+1,bc_j-mdsc:bc_j+mdsc+1]
    # ^ mask values in the box
    boxa = np.tan(boxa/180*math.pi)
    boxb = np.multiply(boxb,(np.power(np.cos(boxa/180*math.pi),4)))    		
    vreg = lreg(boxa.flatten()[boxm.flatten()], boxb.flatten()[boxm.flatten()],nangles,nsamples,filenm)
    sxx_all[bc_i,bc_j] = vreg['sxx']
    S0_all[bc_i,bc_j] = vreg['S0']

    if vreg['corr'] < data['corr'] and vreg['r2'] > data['r2']:
      sxx_filtered[bc_i,bc_j] = vreg['sxx']
      S0_filtered[bc_i,bc_j] = vreg['S0']
    corr[bc_i,bc_j] = vreg['corr']
    err_sxx[bc_i,bc_j] = vreg['err_sxx']
    err_S0[bc_i,bc_j] = vreg['err_S0']
    r2[bc_i,bc_j] = vreg['r2']                                                    
      
  hnsc = int(nsc * 0.5)
  # ^ number of scans
  maskS0 = np.ones(mask.shape) * mask[:,:]
  masksxx = np.ones(mask.shape) * mask[:,:]
  
  if data['add_mid']:         
    S0_filtered[:,hnsc-1:hnsc+2][maskS0[:,hnsc-1:hnsc+2]==1] =sig[:,hnsc-1:hnsc+2][maskS0[:,hnsc-1:hnsc+2]==1]
    S0_all[:,hnsc-1:hnsc+2][maskS0[:,hnsc-1:hnsc+2]==1] =sig[:,hnsc-1:hnsc+2][maskS0[:,hnsc-1:hnsc+2]==1]            
    maskS0[:,hnsc-1:hnsc+2] = True
          
  if data['nf'] > 1:                                
    S0_filtered = smooth(maskS0, S0_filtered, data['nf'], data['filter'])
    sxx_filtered = smooth(masksxx, sxx_filtered, data['nf'], data['filter'])    
  rez = {'sxx_all':sxx_all, 'S0_all': S0_all,'sxx_filtered':sxx_filtered, 'S0_filtered': S0_filtered, 'err_sxx': err_sxx, 'err_S0': err_S0, 'r2':r2, 'corr': corr}
  return(rez)
