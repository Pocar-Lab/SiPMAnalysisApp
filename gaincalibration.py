# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:27:00 2020

@author: lab-341
"""

import datetime
from numpy.polynomial.polynomial import Polynomial
import fnmatch
import pandas as pd
import scipy as sp
import numpy as np
import zipfile
import os
from matplotlib import pyplot as plt
from matplotlib import style
import matplotlib
matplotlib.use("TkAgg", warn=False)
style.use("ggplot")
#import shutil
polyfit = Polynomial.fit
strptime = datetime.datetime.strptime
timedelta = datetime.timedelta
from tempfile import NamedTemporaryFile
import webbrowser

def df_window(df):
    with NamedTemporaryFile(delete = False,suffix = '.html') as f:
        df.to_html(f)
    webbrowser.open(f.name)

def blinfit(x,y,err):
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


def rms(array):
    avg = 0
    for val in array:
        avg += val**2
    avg = avg/len(array)
    avg = np.sqrt(avg)
    return avg

def get_wf(wf, path):
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
    time = time.split('\r')[0]
    
    try:
        x, y = zip(*rows)
        x = [(np.float(dat) + 2e-6)*1e6 for dat in x]
        y = [np.float(dat) for dat in y]
        return x,y,time

    except ValueError:
        print 'wnum {} gave an error, no wfs?'.format(wnum)
        return None,None,None


def process_waveforms(rundir):
    num = len(fnmatch.filter(os.listdir(rundir), 'w*'))
    wfs = fnmatch.filter(os.listdir(rundir), 'w*')
    wfcount = len(wfs)
    yl = []
    
    if wfcount == 0:
        return None
    
    for wf in wfs:
        x, y, time = get_wf(wf, rundir)
        if x == None:
            return None
        #get baseline
        peak = y.index(max(y))
        noisex = np.array(x)
        ind = len(noisex[noisex < 2])
        bl = y[:ind]
        avg = np.mean(bl)
        #-----------#
        output = max(y) - avg
        yl.append(output)
    yl = sp.array(yl)
    avg = yl.mean()
    err = rms(yl - yl.mean())
    return (avg,err,wfcount)



folder = 'D:\Xe\ShAmpCalibrationsNov2019'
dates = os.listdir(folder)

DF = pd.DataFrame(data=None, columns=['date','input','freq','gain','peak_avg','peak_std','wfcount'])

for date in dates[1:-2]:
    print 'On date: ', date
    ext = folder + '\\'+ date + '\Waveforms\\'
    runs = os.listdir(ext)
    
    for run in runs:
        print run
        hold = run.split('_')
        freq = hold[4][:-1]
        gain = hold[3]
        Vin = hold[5]
        
        rundir = ext + '\\' + run
        try:
            avg,err,wfcount = process_waveforms(rundir)
        except TypeError:
            continue
        
        data = {'date':date,'input':Vin,'freq':freq,'gain':gain,'peak_avg':avg,
                'peak_std':err,'wfcount':wfcount}
        series = pd.Series(data)
        DF = DF.append(series,ignore_index=True)
        print 'Run Finished'
        
#solution for the weird ass formatting issue
DF = DF.replace({'14000000000000002':'14mV',
                  '21999999999999996':'22mV',
                  '7000000000000001':'7cV',
                  '41499999999999992':'415mV'})
    
#a solution for the horrible format of '.1V' = '1dV'
def quickreplace(string):
    mag = string[-2]
    try:
        int(mag)
        return float(string[:-1])
    except ValueError:
        pass
    num = float(string[:-2])
    if mag == 'd':
        num = num*0.1
    if mag == 'c':
        num = num*0.01
    if mag == 'm':
        num = num*0.001
    return num

DF.input = DF.input.apply(quickreplace)
##
'''
gains= [0,1,10,11]

for gain in gains:
    plt.figure()
    df = DF[DF.gain == gain]
    df = df.sort_values(by=['input'])
    handles = []
    labels = []
    for date in set(df.date):
        data = df[df.date == date]
        x = data.input
        y = data.peak_avg
        yerr = data.peak_std
        handle = plt.scatter(x,y)
        handles.append(handle)
        labels.append(date)
    
    #now lin fit
    x = df.input
    y = df.peak_avg
    yerr = df.peak_std
    c,m,ac,am = blinfit(x,y,yerr)
    xax = np.linspace(0,max(x))
    handle = plt.plot(xax,m*xax+c)
    handles.append(handle[0])
    labels.append(u'Slope: {:0.2} \u00B1 {:0.1} , \u00B1 Offset: {:0.2} \u00B1 {:0.1}'.format(m,am,c,ac))
    plt.title('Gain Setting {}'.format(gain))
    plt.legend(handles,labels)
    plt.savefig('Gain{}.png'.format(gain))
    plt.xlabel('Input (V)')
    plt.ylabel('Output (V)')
    
    
#method for getting a sample of wf plots
for date in dates[1:-2]:
    #directory stuff
    ext = folder + '\\'+ date + '\Waveforms\\'
    runs = os.listdir(ext)
    run = runs[0]
    #get params
    hold = run.split('_')
    gain = hold[3]
    Vin = hold[5]
    #more dir stuff
    rundir = ext + '\\' + run
    wfs = fnmatch.filter(os.listdir(rundir),'w*')
    wf = wfs[0]
    #get data and plot
    x,y,time = get_wf(wf,rundir)
    plt.figure()
    plt.plot(x,y)
    plt.title('Date: {} , Gain: {}, Input: {}'.format(date,gain,Vin))
    plt.xlabel('Time (us)')
    plt.ylabel('Output (V)')
    plt.savefig('ShaperWF{}.png'.format(date))
#-------------------------------------#
 '''