import h5py as h5
import pandas as pd
import numpy as np
from os.path import join as pj

class data:
    def __init__(self,metadata,data_directory):
        self.meta = pd.read_csv(metadata)
        self.direc = data_directory

    def load_file(self,file):
        file = h5.File(pj(self.direc,file), 'r')
        gi = tgf(file)
        ts = tsf(file)
        ti = tif(file)
        mt = msf(file)
        am = amf(file)
        sk = skf(file)
        rl = rlf(file)
        return dict(gauge_info=gi,timeseries=ts,tides=ti,mean_sea_level=mt,
                    annual_maxima=am,skew_surge=sk,return_levels=rl)

    def load_group(self,file,group):
        datafile = h5.File(pj(self.direc,file), 'r')
        return group2func(group,datafile)

    def load_dataset(self,file,group,dataset):
        datafile = h5.File(pj(self.direc,file), 'r')
        groupdata = group2func(group,datafile)
        if dataset == 'time' or dataset == 'year':
            out = groupdata.index
        else:
            out = groupdata[dataset]
        return out

def group2func(group,file):
    groupfuncs = dict(gauge_info=tgf,
                  timeseries=tsf,
                  tides=tif,
                  mean_sea_level=msf,
                  annual_maxima=amf,
                  skew_surge=skf,
                  return_levels=rlf)
    return groupfuncs[group](file)

# tide gauge info
def tgf(file):
    f = file['tide_gauge_info']
    l = [f[key][:].flatten() for key in f.keys()]
    for i in [0, 2, 3]:
        l[i] = [st.decode('utf-8') for st in l[i]]
        if i != 0:
            l[i] = str(l[i][0])
    dnames = ['contributor_flag_information', 'contributor_flags', 'datum',
                'gauge_type', 'gesla_flags', 'latitude', 'longitude']
    gi = {}
    for i, n in enumerate(dnames):
        gi[n] = l[i]
    return gi

# times of the timeseries  
def tmf(file):
    return pd.to_datetime([st.decode('utf-8') for st in file[:].flatten()])

# water levels, water levels detrended, tide and surge
def tsf(file):
    t = tmf(file['times']['time'])
    ts = pd.DataFrame({
        'water_level': file['water_levels']['water_level'][:].flatten(), 
        'water_level_detrended': file['water_levels']['water_level_detrended'][:].flatten(),
        'tide': file['tides']['tide'][:].flatten(),
        'non_tidal_residual': file['surge']['non_tidal_residual'][:].flatten()},
        index = t)
    ts.rename_axis('time', inplace=True)
    return ts

# tides: 15mins tide and harmonic analysis  
def t15(file):
    f = file['tides']
    t = tmf(f['tide_15min_res_times'])
    tides15 = pd.DataFrame({'tide': f['tide_15min_res'][:].flatten()}, index=t)
    tides15.rename_axis('time', inplace=True)
    return tides15

def hf(file):
    f = file['tides']
    h = f['harmonic_analysis'][:]
    i = h[0,1:].astype(str)
    c = h[1:,0].astype(str)
    h = h[1:,1:].T.astype(float)
    hdf = pd.DataFrame(h, columns=c, index=i)
    hdf.rename_axis('constituents', inplace=True)
    return hdf

def tif(file):
    t = t15(file)
    h = hf(file)
    return {'tide_15min_res':t, 'harmonic_analysis':h}
    
# mean sea level trend
def msf(file):
    f = file['water_levels']
    mt = pd.DataFrame(f['mean_sea_level_trend'][:].T, columns=('year','msl_trend'))
    mt['year'] = mt['year'].astype(int)
    mt.set_index('year', inplace=True)
    ms = pd.DataFrame({'mean_sea_level':f['mean_sea_level'][:].T.flatten()})
    ms['time'] = tmf(file['times']['time'])
    ms.set_index(ms['time'],inplace=True)
    ms = ms.drop('time',axis=1)
    return {'mean_sea_level':ms,'mean_sea_level_trend':mt}

# annual maxima
def amf(file):
    f = file['annual_maxima']
    l = [f[key][:].flatten() for key in f.keys()]
    li = [ind for ind, key in enumerate(f.keys()) if 'time' in key]
    for lix in li:
        l[lix] = pd.to_datetime([st.decode('utf-8') for st in l[lix]])
    am = pd.DataFrame(list(map(list, zip(*l))), columns=f.keys())
    am['year'] = am['year'].astype(int)
    am.set_index('year', inplace=True)
    am = am.reindex(columns=['time_of_water_level', 'water_level',
                                'time_of_water_level_detrended', 'water_level_detrended',
                                'time_of_non_tidal_residual', 'non_tidal_residual',
                                'time_of_skew_surge', 'skew_surge'])
    return am
        
# skew surge
def skf(file):
    f = file['skew_surge']
    l = [f[key][:].flatten() for key in f.keys()]
    li = [ind for ind, key in enumerate(f.keys()) if 'time_of' in key]
    for lix in li:
        l[lix] = pd.to_datetime([st.decode('utf-8') for st in l[lix]])
    sk = pd.DataFrame(list(map(list, zip(*l))), columns=f.keys())
    sk = sk.reindex(columns=['skew_surge', 'time_to_predicted_from_actual_high_water_(minutes)',
                                'predicted_tidal_high_water', 'time_of_predicted_tidal_high_water',
                                'actual_high_water', 'time_of_actual_high_water'])
    return sk

# return levels (& exceedances)
def rlf(file):
    f = file['return_levels']
    l = [f[key][:].flatten() for key in f.keys()]
    colnames = list(f.keys())
    rls = pd.DataFrame(list(map(list, zip(*l[slice(0,7,3)]))),columns=colnames[slice(0,7,3)])
    for sl in (1,4,7):
        rls[colnames[sl]] = list(l[sl].reshape(-1,2,order='F'))
    yr = pd.to_datetime([st.decode('utf-8') for st in l[10]])
    colnums = (1,5,10,25,50,100)
    colabbr = ['gev','gp','gum']                                           
    z = [[colabbr[i] + str(j) for j in colnums] for i in range(0,len(colabbr))]
    epy = pd.DataFrame(np.hstack(l[slice(2,9,3)]).reshape(-1,18,order='F'),columns=sum(z,[]))  
    epy['year'] = yr
    epy.set_index('year', inplace=True)
    return {'return_levels':rls, 'exceedances_per_year':epy}
