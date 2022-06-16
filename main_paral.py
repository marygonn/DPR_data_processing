import multiprocessing, logging
logger = multiprocessing.log_to_stderr()
logger.setLevel(multiprocessing.SUBDEBUG)

import math
import numpy as np
import matplotlib.pyplot as plt
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


def f(fns):   
    with open("config_dpr.json", "r") as read_file:
        data = json.load(read_file)
    if data['Ka_mode'] =='Xmode':
        file_length_min = 2
    else:        
        file_length_min = 2*10**6

    for ifn in range(len(fns)):   
        if os.path.getsize(os.path.join(data['inp_folder_DPR'],fns[ifn][0]))>file_length_min:         
            if data['Ka_mode'] =='Xmode':
                #Lat=Lat,Lon=Lon,S0=S0,Inc=Inc,Surf=Surf,Prec=Prec,Secofday=Secofday
                #Lat=Lat,Lon=Lon,S0=S0,Inc=Inc,Surf=Surf,Prec=Prec,Secofday=Secofday
                #np.savez_compressed  
                npz_ka = np.load(os.path.join(data['inp_folder_DPR'],fns[ifn][0]))
                d = dict()                  
                d['Lat'] = npz_ka['Lat']
                d['Lon'] = npz_ka['Lon']
                d['sigma'] = npz_ka['S0']
                d['secOfDay'] = npz_ka['Secofday']
                d['precipRate'] = npz_ka['Prec']
                plt.hist(d['precipRate'][50][d['precipRate'][50]>-9999])
                #plt.imshow(d['precipRate'][0])
                plt.savefig('rain0.jpg')

                d['IncAngle'] = npz_ka['Inc']
                d['surftype'] = npz_ka['Surf']

                fg = h5py.File(os.path.join(data['inp_folder_GMI'],fns[ifn][1]), 'r')                 
                d['LatGMI'] = np.array(fg['S1']['Latitude'])
                d['LonGMI'] = np.array(fg['S1']['Longitude'])               
                Tc = fg['S1']['Tc']  
                d['ice_conc_GMI'] = asi_hdf(Tc)

                rezka = bandproc(d,data,fns[ifn][0])

                Lat = d['Lat']
                Lon = d['Lon']    
                S0_filtered = rezka['swath']['S0_filtered']
                sxx_filtered = rezka['swath']['sxx_filtered']
                S0_all = rezka['swath']['S0_all']
                sxx_all = rezka['swath']['sxx_all']
                err_S0 = rezka['swath']['err_S0']
                err_sxx = rezka['swath']['err_sxx']
                corr = rezka['swath']['corr']
                r2 = rezka['swath']['r2']                        
                Rain = d['precipRate']
                Sec = d['secOfDay']

                if data['kurt']:                    
                    kurt = rezka['kurt'] 
                    ice_conc = rezka['conc']
                else:
                    kurt = [] 
                    ice_conc = []                   

                np.savez_compressed(os.path.join(data['out_folder'],'N'+fns[ifn][0]), 
                Lat=Lat,Lon=Lon,S0_filtered=S0_filtered,sxx_filtered=sxx_filtered,S0_all=S0_all,sxx_all=sxx_all,
                err_S0=err_S0,err_sxx=err_sxx,corr=corr,r2=r2,Rain=Rain,Sec=Sec) 
            else:               
                print(fns[ifn][0])        
                # DPR - 0            
                f = h5py.File(os.path.join(data['inp_folder_DPR'],fns[ifn][0]), 'r')                  
                hf = h5py.File(os.path.join(data['out_folder'],'N'+fns[ifn][0]),'w') 

                d = dict()
                if data['ice']:
                    # GMI - 1              
                    fg = h5py.File(os.path.join(data['inp_folder_GMI'],fns[ifn][1]), 'r')   
                    d['LatGMI'] = np.array(fg['S1']['Latitude'])
                    d['LonGMI'] = np.array(fg['S1']['Longitude'])               
                    Tc = fg['S1']['Tc']  
                    d['ice_conc_GMI'] = asi_hdf(Tc)
                else:
                    d['LatGMI'] = []
                    d['LonGMI'] = []            
                    d['ice_conc_GMI'] = []
                    
                grKu = hf.create_group('Ku')
                grKa = hf.create_group('Ka')
                
                if data['Ku']:
                    l = 'NS'
                    d['Lat'] = np.array(f['/%s/Latitude' % l])
                    d['Lon'] = np.array(f['/%s/Longitude' % l])

                    d['sigma'] = np.array(f['/%s/VER/sigmaZeroNPCorrected' % l])
                    d['secOfDay'] = np.array(f['/%s/ScanTime/SecondOfDay' % l])
                    d['precipRate'] = np.array(f['/%s/SLV/precipRateESurface' % l])
                    d['IncAngle'] = np.array(f['/%s/PRE/localZenithAngle' % l])        
                    d['surftype'] = np.array(f['/%s/PRE/landSurfaceType' % l])     
                    rezku = bandproc(d,data,fns[ifn][0])

                    grKu.create_dataset('Lat', data=d['Lat'])
                    grKu.create_dataset('Lon', data=d['Lon'])     
                    if data['mss']:
                        grKu.create_dataset('S0_filtered', data=rezku['swath']['S0_filtered'])
                        
                        if data['full_output']:
                            grKu.create_dataset('sxx_filtered', data=rezku['swath']['sxx_filtered'])
                            grKu.create_dataset('S0_all', data=rezku['swath']['S0_all'])
                            grKu.create_dataset('sxx_all', data=rezku['swath']['sxx_all'])
                            grKu.create_dataset('err_S0', data=rezku['swath']['err_S0'])
                            grKu.create_dataset('err_sxx', data=rezku['swath']['err_sxx'])
                            grKu.create_dataset('corr', data=rezku['swath']['corr'])
                            grKu.create_dataset('r2', data=rezku['swath']['r2'])

                        grKu.create_dataset('Rain', data = d['precipRate'])
                    if data['kurt']:
                        # replace kurt by 2 columns!                                    
                        grKu.create_dataset('kurt', data=rezku['kurt'])  
                        grKu.create_dataset('ice_conc', data = rezku['conc'])          
                    hf.create_dataset('Sec', data=d['secOfDay'])

                if data['Ka']:                    
                    l = 'MS'
                    d['Lat'] = np.array(f['/%s/Latitude' % l])
                    d['Lon'] = np.array(f['/%s/Longitude' % l])
                    d['sigma'] = np.array(f['/%s/VER/sigmaZeroNPCorrected' % l])
                    d['secOfDay'] = np.array(f['/%s/ScanTime/SecondOfDay' % l])
                    d['precipRate'] = np.array(f['/%s/SLV/precipRateESurface' % l])
                    d['IncAngle'] = np.array(f['/%s/PRE/localZenithAngle' % l])
                    d['surftype'] = np.array(f['/%s/PRE/landSurfaceType' % l])

                rezka = bandproc(d,data,fns[ifn][0])

                grKa.create_dataset('Lat', data=d['Lat'])
                grKa.create_dataset('Lon', data=d['Lon'])     
                if data['mss']:
                    grKa.create_dataset('S0_filtered', data=rezka['swath']['S0_filtered'])

                    if data['full_output']:
                        grKa.create_dataset('sxx_filtered', data=rezka['swath']['sxx_filtered'])
                        grKa.create_dataset('S0_all', data=rezka['swath']['S0_all'])
                        grKa.create_dataset('sxx_all', data=rezka['swath']['sxx_all'])
                        grKa.create_dataset('err_S0', data=rezka['swath']['err_S0'])
                        grKa.create_dataset('err_sxx', data=rezka['swath']['err_sxx'])
                        grKa.create_dataset('corr', data=rezka['swath']['corr'])
                        grKa.create_dataset('r2', data=rezka['swath']['r2'])
                        
                    grKa.create_dataset('Rain', data = d['precipRate'])
                if data['kurt']:                    
                    grKa.create_dataset('kurt', data=rezka['kurt'])  
                    grKa.create_dataset('ice_conc', data = rezka['conc'])     
                if not('Sec' in hf.keys()):     
                    hf.create_dataset('Sec', data=d['secOfDay'])
            
                hf.attrs['inp_folder_DPR'] =  data['inp_folder_DPR']
                hf.attrs['inp_folder_GMI'] =  data['inp_folder_GMI']
                hf.attrs['part_name_DPR'] = data['part_name_DPR']
                hf.attrs['part_name_GMI'] = data['part_name_GMI']
                hf.attrs['out_folder'] = data['out_folder']

                hf.attrs['Ku'] = data['Ku']
                hf.attrs['Ka'] = data['Ka']
                hf.attrs['kurt'] = data['kurt']
                hf.attrs['mss'] = data['mss']
                hf.attrs['add_mid'] = data['add_mid']
                hf.attrs['corr'] = data['corr']
                hf.attrs['r2'] = data['r2']
                
                hf.attrs['dfl'] = data['dfl']
                hf.attrs['dsc'] = data['dsc']
                hf.attrs['nf'] = data['nf']
                hf.attrs['filter'] = data['filter']
                
                hf.attrs['lon_min'] = data['lon_min']
                hf.attrs['lon_max']  = data['lon_max']
                hf.attrs['lat_min'] = data['lat_min']
                hf.attrs['lat_max'] = data['lat_max']   
                hf.attrs['max_inc'] = data['max_inc']

                hf.close()
                f.close()
                if data['ice']:
                    fg.close()
        else:
            print('short file: ' + os.path.join(data['inp_folder_DPR'],fns[ifn][0]))
    
      
def f_mp(fns,nprocs):
    print('fns')
    print(fns)
    chunks = [fns[i::nprocs] for i in range(nprocs)] 
    print('chunks')
    print(chunks)	
    pool = Pool(processes=nprocs)
    result = pool.map(f, chunks)
    return ()

if __name__ == '__main__':    
    print('hello!')
    with open("config_dpr.json", "r") as read_file:
        data = json.load(read_file)
        
    if not os.path.exists(data['out_folder']):
        os.mkdir(data['out_folder'])
        
    fnsDPR = os.listdir(data['inp_folder_DPR'])
    fnsGMI = os.listdir(data['inp_folder_GMI'])    
    
    # print(fnsDPR)
    # print(fnsGMI)

    '''
    fns1_DPR = []
    for entry in fnsDPR:  
        if fnmatch.fnmatch(entry, data['part_name_DPR']):
            fns1_DPR.append(entry)    
    '''

    fns1_DPR = []
    fns1_GMI = []
    for entryDPR in fnsDPR:                  
        for entryGMI in fnsGMI:
            if fnmatch.fnmatch(entryGMI, data['part_name_GMI']):
                if fnmatch.fnmatch(entryDPR, data['part_name_DPR']): 
                    # print(entryDPR)
                    # print(entryGMI)
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
    print(q)
    fns = list(np.transpose(q))
    print('fns')
    print(fns)
    f_mp(fns,data['nprocs'])
    
