# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:27:00 2020

@author: lab-341
"""

from scipy import signal  # needs to be here otherwise sp.signal won't be recognized
from tkinter import messagebox
from tkFileDialog import askdirectory, askopenfilename
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import queue
import threading
import datetime
from numpy.polynomial.polynomial import Polynomial
import fnmatch
import pandas as pd
import scipy as sp
import numpy as np
import zipfile
from uuid import uuid4
import re
import os
import ttk
import Tkinter as tk
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
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
        output = max(y)
        yl.append(output)
    
    yl = sp.array(yl)
    avg = yl.mean()
    err = rms(yl - yl.mean())
    return (avg,err,wfcount)

folder = 'D:\Xe\ShAmpCalibrationsNov2019'
dates = os.listdir(folder)

DF = pd.DataFrame(data=None, columns=['date','input','freq','gain','peak_avg','peak_std','wfcount'])

for date in dates:
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