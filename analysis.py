# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:38:19 2020

@author: Reed Cohen
"""

import datetime
import pytz
import numpy as np
import scipy as sp
import pandas as pd
import os
from numpy.polynomial.polynomial import Polynomial
import threading
import zipfile
import fnmatch
import scipy.stats as stats
strptime = datetime.datetime.strptime
timedelta = datetime.timedelta
polyfit = Polynomial.fit

IOLock = threading.Lock()

peak_separation = 5  
same_peak = .005
dt_to_ind = 250  # conversion from microseconds to indices, divide by 1000 to ns, then 4ns per index
binnum = 'fd'

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def rms(array):
    '''
    Obtains the deviation from the mean.
    '''
    avg = 0
    for val in array:
        avg += val**2
    avg = avg/len(array)
    avg = np.sqrt(avg)
    return avg

def blinfit(x,y,err):
    '''
    Mathematically exact best linear fit for data with errors.
    '''
    x = sp.array(x)
    x = x.astype(float)
    y = sp.array(y)
    y = y.astype(float)
    err = sp.array(err)
    w = err**(-2)
    
    wx = w*x
    wx2 = w*(x**2)
    wy = w*y
    wxy = (w*x)*y
    
    D = w.sum()*wx2.sum() - wx.sum()**2
    c = (wx2.sum()*wy.sum() - wx.sum()*wxy.sum())/D
    m = (w.sum()*wxy.sum() - wx.sum()*wy.sum())/D
    
    ac = np.sqrt(wx2.sum()/D)
    am = np.sqrt(w.sum()/D)
    
    return (c,m,ac,am)

def get_peaks(x, y):
    '''
    Analyze waveform and find peaks with filtering then return info
    about peaks.
    '''
    threshold = np.mean(y)
    # (None,None) puts value in props but doesn't add a filter cond.
    peak_inds, props = sp.signal.find_peaks(y, prominence=np.mean(y), height=threshold,
                                            distance=10, width=same_peak * dt_to_ind)
    peak_inds = list(peak_inds)

    # if the waveform has no pulse, take the largest noise peak
    if len(peak_inds) == 0:
        peak_inds = [max(y)]

    # filter out cutoff peaks
    for i in peak_inds:
        if i < peak_separation:
            # print('peak cutoff')
            return None
    
    # irrelevant in current version, we only take clean waveforms (one pulse)
    '''
    #filter out crowded peaks
    for pair in combinations(peak_inds,r=2):
        if abs(pair[0] - pair[1]) < peak_separation * ind_to_dt:
            try:
                peak_inds.remove(pair[0])
                peak_inds.remove(pair[1])
            except ValueError:
                continue
    '''

    # get valid data for fitting baseline
    start = find_start(x, y, peak_inds[0])
    if start is None:
        # print('bad start')
        return None
    bldat = y[:start]

    bl, = polyfit(x[:start], y[:start], deg=0)

    bldat = np.array(bldat)
    bldat -= bl
    bldat = list(bldat)
    noiserms = rms(bldat)

    # identify unwanted pulses occuring right before recording
    if bldat.index(max(bldat)) == 0 and y[1] > y[2]:
        return None

    # for now, no more than one pulse is allowed
    if len(peak_inds) > 1:
        return None

    elif len(peak_inds) == 1:
        peak_dat = dict({'peak_height': props['peak_heights'][0],
                         'baseline': bl,
                         'noise_rms': noiserms,
                         'peak_width': props['right_ips'][0] - props['left_ips'][0],
                         'peak_time': x[peak_inds[0]]})

    return peak_dat

def get_dt(path):
    '''
    A function to get the datetime from the first waveform in a runfile.
    '''
    path = str(path)
    pathlist = path.split('/')
    date = pathlist[3]
    date = date[:4] + '/' + date[4:6] + '/' + date[6:]
    try:
        z = zipfile.ZipFile(path + '/w1.zip', 'r')
        with z.open('w1.txt', 'r') as f:
            lines = f.readlines()
    except IOError:
        with open(path + '/w1.txt','r') as f:
            lines = f.readlines()
    timerow = lines[1]
    time = timerow.split('\t')[1]
    time = time.split('\r')[0]
    dt = date + ' ' + time
    dt = dt.rstrip()
    try:
        dt = strptime(dt, '%Y/%m/%d %H:%M:%S.%f')
    except:
        pass
    return dt

def line_to_dt(line):
    '''
    Turns a m/d/y string into a datetime object.
    '''
    dt = strptime(' '.join(line.split()[:2]), '%m/%d/%Y %H:%M:%S.%f')
    return dt

def process_waveforms(dfname, rundir):
    '''
    Function to analyze all wfs in a given directory. Returns a DataFrame.
    '''
    num = len(fnmatch.filter(os.listdir(rundir), 'w*'))
    wfs = fnmatch.filter(os.listdir(rundir), 'w*')
    
    skipped = 0

    DF = pd.DataFrame(data=None, columns=['time', 'peak_height', 'baseline',
                                          'noise_rms', 'peak_width',
                                          'peak_time'])
    DF.index.name = 'wnum'
    
    
    for wf in wfs:
        wnum = wf.split('.')[0][1:]
        x, y, time = get_wf(wf, rundir)
        
        if x is None:
            # print( 'w{} was skipped'.format(wnum)
            skipped += 1
            continue

        peak_dat = get_peaks(x, y)

        if peak_dat is None:
            # print( 'w{} was skipped'.format(wnum))
            skipped += 1
            continue

        # start time of the waveform, peak_time is when the PEAK happened
        peak_dat.update({'time': time})

        series = pd.Series(peak_dat, name=str(wnum))
        DF = DF.append(series)
        
        
        
    print( '# of waveforms total: ', num)
    print( '# of skipped waveforms: ', skipped)

    return DF

def find_start(x, y, peak):
    '''
    Follows a peak down the right side until the right side reaches
    less than the minima of the left, then returns that distance less
    than the peak index.
    '''
    i = peak
    while True:
        if y[i] < min(y[:peak]):
            break
        elif i == len(y) - 1:
            break
        else:
            i += 1
            continue
    width = abs(peak - i)
    
    start = peak - width
    
    if start < 10:
        return None
    return start

def get_wf(wf, path):
    '''
    Returns x,y,time for a wf in the given path. Returns None triplet
    if there are no wfs or some other error occurs.
    
    This assumes that the time recorded in the waveform and the filename
    is using the EST timezone and affected by DST.
    '''
    
    wpath = path + '/' + str(wf)
    wnum = wf.split('.')[0][1:]
    try:
        z = zipfile.ZipFile(wpath, 'r')
        with z.open('w{}.txt'.format(wnum), 'r') as f:
            lines = f.readlines()
    except (IOError,zipfile.BadZipfile):
        wfpath = path + '/w{}.txt'.format(wnum)
        try:
            with open(wfpath,'r') as f:
                lines = f.readlines()
        except IOError:
            return None,None,None
    r = lines[16:]
    rows = [l.split() for l in r]
    
    timerow = lines[1]
    time = timerow.split('\t')[1]
    timestr = time.split('\r')[0].strip()
    pathlist = path.split('/')
    datestr = pathlist[3].strip()
    dtstr = datestr + ' ' + timestr
    dt = datetime.datetime.strptime(dtstr, '%Y%m%d %H:%M:%S.%f')
    tz = pytz.timezone('EST') #assume text file used EST with DST effects
    dt = tz.localize(dt) 
    utcdt = dt.astimezone(pytz.utc)
    time = utcdt.timestamp()
    
    try:
        x, y = zip(*rows)
        x = [(np.float(dat) + 2e-6)*1e6 for dat in x]
        y = [np.float(dat) for dat in y]
        return x,y,time

    except ValueError:
        print('wnum {} gave an error, no wfs?'.format(wnum))
        return None,None,None

def histogram(peaks, info=False):
    
    peaks = np.array(peaks)
    histogram = np.histogram(peaks, bins=binnum)
    
    bins = histogram[1]
    centers = (bins[1:] + bins[:-1])/2
    counts = histogram[0]
    
    #initial guesses
    p0 = (1,) + stats.norm.fit(peaks)
    
    popt,pcov = sp.optimize.curve_fit(gauss,centers,counts,p0 = p0,
                                      method='lm' ,maxfev = 10000)
    
    popt[2] = abs(popt[2]) #only pos. value for sigma
    
    #Restrict actual fitting domain
    i = -1
    while True:
        i -= 1
        if centers[i] < popt[1] - 2*popt[2]:
            break
    left = i
    
    i = 0
    while True:
        i += 1
        if i == len(centers)-1 and centers[i] < popt[1] + 0.5*popt[2]:
            raise ValueError
        elif centers[i] > popt[1] + 2*popt[2]:
            break
        elif i == len(centers)-1:
            break
    right = i
    
    centers = np.array(centers[left:right])
    counts = np.array(counts[left:right])
    vpeaks = peaks[centers[0] < peaks]
    vpeaks = vpeaks[centers[-1] > vpeaks]
    
    sigma = np.sqrt(counts*(1 - counts/len(peaks))) #if bin heights are binomial
    
    res = counts - gauss(centers,*popt)
    #plt.plot(centers,res)
    res = res**2
    res = res/(sigma**2)
    res = res[np.invert(np.isinf(res))]
    Chi2 = float(res.sum())/(len(counts) - 3)
    
    #Perform second and proper fit
    popt,pcov = sp.optimize.curve_fit(gauss,centers,counts,p0 = popt,
                                  sigma=sigma,
                                  method='lm' ,maxfev = 10000,
                                  absolute_sigma=True)
    popt[2] = abs(popt[2]) #only pos. value for sigma
    valpeaks = counts.sum()
    
    if info == True:
        return popt,pcov,Chi2,valpeaks,centers
    else:
        return popt,pcov,Chi2

def alpha_df(rundir, peaksDF, tempsDF, peakid, makenew=False):
    '''
    Collects info from peak DF and temperature DF to create new
    single row DF which can be appended to a DF in the alphas.h5 hdf.
    '''
#     get label data
    rundir = str(rundir)
    pathlist = rundir.split('/')
    rawdate = pathlist[3]
    sep = pathlist[2].split('_')[2]
    separation = ''.join([char for char in sep if char.isdigit() or char == '.'])
    biasV = pathlist[-1].split('_')[4]
    setpt = pathlist[-1].split('_')[3]

    #get peak data
    peaks = peaksDF['peak_height'] - peaksDF['baseline']
    peaks = np.array(peaks)
    
    try:
        popt,pcov,Chi2 = histogram(peaks)
    except:
        popt = np.zeros((3,3),int)
        pcov = np.zeros((3,3),int)
    
    popt[2] = abs(popt[2]) #only pos. value for sigma
    
    #get temperature data
    temps = tempsDF['temperature']
    tempavg = np.mean(temps)
    temprms = rms([t - tempavg for t in temps])
    
    #collect data in dict
    data = {'date': rawdate, 'separation': separation,
            'biasV': biasV, 'setpoint': setpt,
            'midpoint': popt[1], 'midpt_error': np.sqrt(pcov[1][1]),
            'sigma': popt[2], 'sigma_error': np.sqrt(pcov[2][2]),
            'temperature_avg': tempavg, 'temperature_rms': temprms,
            'wfpath': rundir}
    
    columns = ['date', 'separation', 'biasV',
               'midpoint', 'midpt_error', 'sigma',
               'sigma_error', 'temperature_avg',
               'temperature_rms', 'setpoint', 'wfpath']
    
    DF = pd.DataFrame(data, columns=columns, index=[peakid])
    return DF