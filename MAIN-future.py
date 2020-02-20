                                                                                     
# =============================================================================
# Currently working on updating this code to work in Python 3.
# Primarily this consists of changing print commands, and updating syntax
# with Tk modules. Hopefully there is minimal syntax change and just a change 
# in the way that the modules are loaded.
# =============================================================================
#GUI imports
import tkinter as tk
from tkinter import messagebox, ttk
from tkinter.filedialog import askdirectory, askopenfilename
import queue
import threading
#Plotting
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib import style
#Analysis
import datetime
import fnmatch
import pandas as pd
import scipy as sp
import numpy as np
from numpy.polynomial.polynomial import Polynomial
#I/O commands and file management
import zipfile
from uuid import uuid4
import re
import os

matplotlib.use("TkAgg", warn=False)
style.use("ggplot")
#import shutil
polyfit = Polynomial.fit
strptime = datetime.datetime.strptime
timedelta = datetime.timedelta
from tempfile import NamedTemporaryFile
import webbrowser

temproot = tk.Tk()
Xepath = askdirectory(
    title='Please select the Xe folder which you would like to work in')
peakspath = Xepath + '/peaks.h5'
alphapath = Xepath + '/alphas.h5'
recordspath = Xepath + '/records.h5'
temppath = Xepath + '/temps.h5'
temproot.destroy()

LARGE_FONT = ('Verdana', 10)


def df_window(df):
    with NamedTemporaryFile(delete = False,suffix = '.html') as f:
        df.to_html(f)
    webbrowser.open(f.name)

# =============================================================================
# Analysis Constants
# =============================================================================
peak_separation = 5  

same_peak = .005
dt_to_ind = 250  # conversion from microseconds to indices, divide by 1000 to ns, then 4ns per index
binnum = 150
# =============================================================================
# Threading Functions
# =============================================================================

killthread = False
IOLock = threading.Lock()
tasks = queue.Queue()
outputs = queue.Queue()



def process_data(task):
    '''
    Thread-safe function to acquire requested data and run the analysis.
    task = dfname, peakid, rundir, TCPath, make_new
    
    dfname is a parent group in the hdf file where the processed data should
    be put
    peakid is a unique uuid for the particular datapoint
    rundir is the full directory source of the raw data
    TCPath is the location of the temperature file for the run
    make_new is a boolean indicating if a new hdf is to be made for the data
    '''
    
    dfname, peakid, rundir, TCPath, make_new = task
    
    print('Processing Waveforms for {}'.format(rundir))
    peaksDF = process_waveforms(dfname, rundir)
    print('Processing TC Data')
    tempsDF = process_TC(TCPath, rundir)
    alphaDF = alpha_df(rundir, peaksDF, tempsDF, peakid,
                       make_new)  # may be series
    print('Analysis Complete, writing to files...')
    with IOLock:
        # write peakdata to peaks.h5
        with pd.HDFStore(peakspath, 'a') as hdf:
            hdf.put(dfname + '/' + peakid, peaksDF)

        # write temperature data to temps.h5
        with pd.HDFStore(temppath, 'a') as hdf:
            hdf.put(dfname + '/' + peakid, tempsDF)

        # append alpha series to alphas.h5 (or write fresh alphaDF)
        with pd.HDFStore(alphapath, mode='a') as hdf:
            # table format allows for future appending
            hdf.append(dfname, alphaDF, format='table',
                       data_columns=['date', 'separation', 'biasV',
                                     'midpoint', 'midpt_error', 'sigma',
                                     'sigma_error', 'temperature_avg',
                                     'temperature_rms', 'setpoint',
                                     'wfpath'], 
                                     min_itemsize={'wfpath': 150, 'setpoint': 16,
                                                                  'separation':8})
    print('Task Complete')



def threader():
    '''
    Passive function running in thread to pass tasks and call the
    data processing function on them.
    '''
    while True:
        if killthread:
            break
        task = tasks.get()  # blocks if the queue is empty, so this doesn't waste CPU
        process_data(task)
        tasks.task_done()
        outputs.put(task) #return just the rundir and if it was new
    print('Thread {} finished'.format(threading.currentThread()))

mythread = threading.Thread(target=threader)
mythread.start()

# =============================================================================
# GUI Functions
# =============================================================================
def treeview_sort_column(tv, col, reverse):
    '''
    A function to be called on a TreeView column to 
    sort the rows by the column value.
    '''
    l = [(tv.set(k, col), k) for k in tv.get_children('')]
    l.sort(reverse=reverse)

    # rearrange items in sorted positions
    for index, (val, k) in enumerate(l):
        tv.move(k, '', index)

    # reverse sort next time
    tv.heading(col, command=lambda:
               treeview_sort_column(tv, col, not reverse))

# WF Import and TC Line Processing

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
    time = time.split('\r')[0]
    
    try:
        x, y = zip(*rows)
        x = [(np.float(dat) + 2e-6)*1e6 for dat in x]
        y = [np.float(dat) for dat in y]
        return x,y,time

    except ValueError:
        print('wnum {} gave an error, no wfs?'.format(wnum))
        return None,None,None

# =============================================================================
# Data Processing
# =============================================================================

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


def process_TC(TCPath, rundir):
    '''
    Function to return DataFrame of Temperature info for a run.
    '''

    times = []
    temps = []
    i = 0
    runtime = get_dt(rundir)
    delta = timedelta(minutes=2)

    while True:
        with open(TCPath, 'r') as f:
            lines = f.readlines()

        for line in lines:
            i += 1
            if line.isspace():
                continue
            linetime = line_to_dt(line)
            if linetime < runtime < linetime + delta:
                temps.append(float(line.split()[-1]))
                times.append(line_to_dt(line))
            else:
                continue

        if temps == []:
            TCPath = askopenfilename(
                title='Select the correct TCTests file for run: {}'.format(rundir))
            continue
        else:
            break

    DF = pd.DataFrame({'time': times, 'temperature': temps})

    return DF


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

#     get peak data
    peaks = peaksDF['peak_height']
    peaks = np.array(peaks)
    

    histogram = np.histogram(peaks, bins=binnum, density=True)

    
    bins = histogram[1]
    centers = (bins[1:] + bins[:-1])/2
    counts = histogram[0]

    p0 = (1,) + sp.stats.norm.fit(peaks)
    
    
    popt, pcov = sp.optimize.curve_fit(gauss, centers, counts, p0=p0,maxfev = 100000)
    
#     get temperature data
    temps = tempsDF['temperature']
    tempavg = np.mean(temps)
    temprms = rms([t - tempavg for t in temps])

#     collect data in dict
    data = {'date': rawdate, 'separation': separation,
            'biasV': biasV, 'setpoint': setpt,
            'midpoint': popt[1], 'midpt_error': pcov[1][1],
            'sigma': popt[2], 'sigma_error': pcov[2][2],
            'temperature_avg': tempavg, 'temperature_rms': temprms,
            'wfpath': rundir}

    columns = ['date', 'separation', 'biasV',
               'midpoint', 'midpt_error', 'sigma',
               'sigma_error', 'temperature_avg',
               'temperature_rms', 'setpoint', 'wfpath']

    DF = pd.DataFrame(data, columns=columns, index=[peakid])
    return DF


# =============================================================================
# GUI Classes
# =============================================================================

class LXeDataManager(tk.Tk):
    '''
    The primary applet. Initializes the start page and creates the
    other pages, sets up the page2page functionality.
    '''
    
    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage,  AlphaViewer, RunFitting): #WaveformViewer,
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky='NSEW')

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class StartPage(tk.Frame):
    '''
    Start Page of the applet
    '''
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Start Page", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        #button1 = ttk.Button(self, text="Visit Waveform Viewer",
        #                    command=lambda: controller.show_frame(WaveformViewer))
        #button1.pack()

        button2 = ttk.Button(self, text="Visit Alpha Data Viewer",
                             command=lambda: controller.show_frame(AlphaViewer))
        button2.pack()

        
        button3 = ttk.Button(self,text = "Visit Run Fit Viewer",
                                 command = lambda: controller.show_frame(RunFitting))
        button3.pack()

'''
class WaveformViewer(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent, bg='LightSteelBlue4')
        
        self.controller = controller
        
        self.browse = ttk.Button(
            self, text='Browse', command=self.browse_button)
        self.rundir = tk.StringVar()
        self.rundir.set('Select a folder of Waveforms...')
        self.folder = ttk.Label(self, textvariable=self.rundir)
        self.browse.grid(row=0, column=0)
        self.folder.grid(row=0, column=1, columnspan=2)
        
        self.wflist = tk.Listbox(self, selectmode=tk.BROWSE)
        self.wflist.grid(row=1, column=0, sticky='NS',padx = 10)
        
        
        self.button1 = ttk.Button(self, text='PLOT WAVEFORM',
                                  command=self.plot_wf)
        self.button1.grid(row=2, column=0)
        
        self.f = Figure(figsize=(5, 5), dpi=100)
        self.ax = self.f.add_subplot(111)

        self.graphframe = tk.Frame(self)
        self.graphframe.grid(row=1, column=1)
        
        self.canvas = FigureCanvasTkAgg(self.f, master=self.graphframe)
        self.canvas.get_tk_widget().pack()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.graphframe)
        self.toolbar.update()
        
        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.grid(row=2, column=1)
        
        button2 = ttk.Button(self, text="Get Missing",
                             command=self.get_missing)
        button2.grid(row=3,column=0)
        
    def browse_button(self):
        filename = askdirectory()
        self.rundir.set(filename)
        wfs = os.listdir(filename)
        wfs = filter(lambda name: re.match('w', name), wfs)
        wfs.sort(key=lambda name: int(name[1:-4]))
        self.wflist.delete(0, 'end')
        for wf in wfs:
            self.wflist.insert('end', wf)
    
    def get_missing(self):
        cont = self.controller
        aview = cont.frames[AlphaViewer]
        dfname = aview.dfname
        peakid, = aview.table.selection()
        
        with IOLock:
            with pd.HDFStore(alphapath,'r') as hdf:
                row = hdf.select(dfname,where='index=={}'.format(peakid))
            
            with pd.HDFStore(peakspath,'r') as hdf:
                DF = hdf[dfname + '/' + peakid]
        
        indices = DF.index
        indices = indices.astype(int)
        
        missing = []
        for i in range(1,int(indices[-1])):
            if i not in indices:
                missing.append(i)
        
        self.rundir.set(row['wfpath'][0])
        wfs = ['w{}.zip'.format(num) for num in missing]
        
        self.wflist.delete(0,'end')
        for wf in wfs:
            self.wflist.insert('end',wf)
        
    def plot_wf(self):
        
        wf = self.wflist.curselection()
        if len(wf) != 0:
            wfname = self.wflist.get(wf[0])
            wnum = int(wfname[1:-4])
        else:
            return
        
        wfname = 'w{}.zip'.format(wnum)
        x, y, time = get_wf(wfname, str(self.rundir.get()))
        
        if x is None:
            print( "WF doesn't exist")
            return

        threshold = np.mean(y)
        peak_inds, props = sp.signal.find_peaks(y, prominence=np.mean(y), height=threshold,
                                                distance=10, width=same_peak * dt_to_ind)

        start = find_start(x, y, peak_inds[0])
        
        self.x = x
        self.y = y
        self.ax.clear()
        self.ax.plot(x, y)
        if start is not None:
            xstart = x[start]
            self.ax.axvline(xstart)
        for peak in peak_inds:
            self.ax.plot(x[peak], y[peak], 'b.')
        self.canvas.draw()
'''

class RunFitting(tk.Frame):
    '''
    Applet page used for plotting runs from the alphas.h5 DF and fitting them.
    '''
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self,parent,bg = 'LightSteelBlue4')
        
        #collect the relevant info in a local df variable
        with pd.HDFStore(alphapath,'r') as hdf:
            i = 0
            for key in hdf.keys():
                if i == 0:
                    df = hdf[key]
                else:
                    df = df.append(hdf[key])
                i += 1
        
        separations = list(set(df.separation))
        biases = list(set(df.biasV))
        store = pd.HDFStore(peakspath,'r')
        i = 0
        
        #Adds a counts column to the DF
        for s in separations:
            keys = df[df.separation == s].index
            s = s[:2]
            top = '/VUV4_{}mm/'.format(s)
            n = []
            for key in keys:
                try:
                    n.append(len(store[top + key]))
                except KeyError:
                    n.append(sp.nan)
                    print( 'Missing {}'.format(top + key))
                    pass
            if i == 0:
                i += 1
                counts = pd.Series(n,keys)
            else:
                counts = counts.append(pd.Series(n,keys))
        store.close()
        counts.name = 'num'
        df = df.join(counts)
        
        self.df = df
        
        self.tef = tk.IntVar()
        self.postb = tk.IntVar()
        
        self.sepmenu = ttk.Combobox(self, values=separations, text='Separation')
        self.biasmenu = ttk.Combobox(self, values=biases, text='BiasVoltage')
        self.teflon = ttk.Checkbutton(self, text='Teflon', variable=self.tef)
        self.baked = ttk.Checkbutton(self, text='PostBaking', variable=self.postb)
        
        self.datelist = ttk.Combobox(self)
        
        self.seplabel = ttk.Label(self, text='Separation')
        self.biaslabel = ttk.Label(self, text='Bias Voltage')
        self.seplabel.grid(row=0,column=0, sticky='NSEW',pady=1)
        self.biaslabel.grid(row=1,column=0, sticky='NSEW',pady=1)
        self.sepmenu.grid(row=0,column=1, sticky='NSEW',pady=1)
        self.biasmenu.grid(row=1,column=1, sticky='NSEW',pady=1)
        self.teflon.grid(row=2,column=0, sticky='NSEW',pady=1)
        self.baked.grid(row=2,column=1, sticky='NSEW',pady=1)
        self.datelabel = ttk.Label(self, text='Date')
        self.datelabel.grid(row=3,column=0,sticky='NSEW',pady=1)
        self.datelist.grid(row=3,column=1, sticky='NSEW',pady=1)
        
        self.refresh_dates()
        
        self.button = ttk.Button(self,text = 'Refresh',
                                 command=self.refresh_dates)
        self.button.grid(row=4,column=0)
        
        self.f = Figure(figsize=(7, 7), dpi=100)
        self.ax = self.f.add_subplot(111)
        
        self.graphframe = tk.Frame(self)
        self.graphframe.grid(row=0, column=2,rowspan = 7, sticky='NSEW',padx=10,pady=10)
        
        self.canvas = FigureCanvasTkAgg(self.f, master=self.graphframe)
        self.canvas.get_tk_widget().pack()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.graphframe)
        self.toolbar.update()
        
        
        self.plotbutton = ttk.Button(self, text='Plot Data',
                                     command=self.plot_data)
        self.plotbutton.grid(row=4,column=1)
        
        homebutton = ttk.Button(self, text='Back to Home',
                                command=lambda: controller.show_frame(StartPage))
        homebutton.grid(row=5, column=0,columnspan=2,pady = 10)
        
        self.rowconfigure(6, weight=1)
        
        
    def refresh_dates(self):
        '''
        Refreshes the dates menu for selection. Only presents dates
        that meets the other given requirements.
        '''
        sep = self.sepmenu.get()
        bv = self.biasmenu.get()
        tef = self.tef.get()
        postb = self.postb.get()
        df = self.df
        
        data = df[df.biasV == bv]
        data = data[data.separation == sep]
        
        if tef:
            data = data[data.date.astype(int) >= 20190424]
        
        if postb:
            data = data[data.date.astype(int) >= 20190605]
        
        dates = list(set(data.date))
        self.datelist['values'] = dates
    
    def plot_data(self):
        '''
        Uses the date given and other parameters to select and plot data
        and a linaer fit.
        '''
        bv = self.biasmenu.get()
        df = self.df
        
        date = int(self.datelist.get())
        dates = [date,date+1,date+2,date-1,date-2]
        dates = [str(item) for item in dates]
        
        data = df[df.biasV == bv]
        data = data[data.date.isin(dates)]
        
        x = data.temperature_avg
        xerr = data.temperature_rms
        y = data.midpoint
        yerr = data.sigma/np.sqrt(data.num)
        
        self.ax.clear()
        self.ax.errorbar(x,y,yerr,xerr,fmt='none')
        
        c,m,ac,am = blinfit(x,y,yerr)
        xax = np.linspace(min(x),max(x),num=5000)
        self.ax.plot(xax,m*xax + c)
        
        legend = [u'Slope: {:0.2} \u00B1 {:0.2} , Intercept: {:0.2} \u00B1 {:0.2}'.format(
                m,am,c,ac)]
        self.ax.legend(legend)
        
        self.ax.set_title('Date: {}, Bias: {}'.format(date,bv))
        self.ax.set_xlabel('Temperature (K)')
        self.ax.set_ylabel('Alpha Peak (V)')
        
        self.canvas.draw()

    
class AlphaViewer(tk.Frame):
    '''
    Applet page used for viewing the contents of a DataFrame and processing
    new data to add to the DF.
    '''
    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent, bg='LightSteelBlue4')

        self.rundir = Xepath #default filebrowser prompt for single run select
        self.tcpath = Xepath #default filebrowser prompt for TC selection
        self.runlist = Xepath #default filebrowser prompt for multiple run select

        self.table = ttk.Treeview(self)
        self.table.heading('#0', text='WFPeaksID')
        self.table.column('#0', width=30, stretch=False)
        self.table.grid(row=0, column=2, rowspan=10, sticky='NSEW',padx = 10,pady = 10)
        self.columnconfigure(2, weight=1)
        #self.rowconfigure(0, weight=1)
        self.rowconfigure(9,weight = 1)
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)
        
        self.namelist = namelist = tk.Listbox(self, selectmode=tk.BROWSE)
        namelist.grid(row=1, column=1, rowspan = 4,sticky='NS')
        
        listlabel = ttk.Label(self,text = 'Alpha Frame List',anchor = 'center')
        listlabel.grid(row=0, column=1, sticky = 'SEW')
        
        homebutton = ttk.Button(self, text='Back to Home',
                                command=lambda: controller.show_frame(StartPage))
        homebutton.grid(row=0, column=0,pady = 10)

        button1 = ttk.Button(self, text='Display Data',
                             command=self.pop_table)
        button1.grid(row=2, column=0,sticky='NSEW')

        button2 = ttk.Button(self, text='New Entry',
                             command=self.new_entry)
        button2.grid(row=6, column=0,sticky = 'SEW',pady = (10,0))

        #label3 = tk.Label(self, text='Create New Alpha DF')
        #label3.grid(row=0, column=0)
        self.make_new = tk.BooleanVar()
        self.button3 = ttk.Checkbutton(self, variable=self.make_new,
                                       onvalue=True, offvalue=False,
                                       text = 'Create New Alpha DF')
        self.button3.grid(row=7, column=0,sticky = 'EW')
        self.entry3 = ttk.Entry(self)
        self.entry3.grid(row=7, column=1,sticky = 'EW',pady = 5)
        
        button4 = ttk.Button(self, text='Refresh',
                             command=self.refresh_list)
        button4.grid(row=1, column=0,sticky = 'NSEW')

        button5 = ttk.Button(self, text='Delete Frame',
                             command=self.del_frame)
        button5.grid(row=4, column=0, sticky='NSEW')
        
        button6 = ttk.Button(self, text='Display Histogram',
                             command=self.display_hist)
        button6.grid(row=8, column=1,sticky = 'EW')
        
        button7 = ttk.Button(self, text='Save Table to Clipboard',
                             command=self.save_to_clipboard)
        button7.grid(row=3, column=0,sticky = 'NSEW')

        button8 = ttk.Button(self, text='Delete Row',
                             command=self.del_row)
        button8.grid(row=6, column=1,sticky = 'SEW',pady = (10,0))
        self.check_queue()
        
        self.multiple_dirs = tk.BooleanVar()
        self.button9 = ttk.Checkbutton(self, variable=self.multiple_dirs,
                                       text='Select Multiple Runs',
                                       onvalue=True, offvalue=False)
        self.button9.grid(row=8, column=0)
        
        
        self.queuelist = ttk.Treeview(self)
        self.queuelist.grid(row=9,column=0,columnspan=2,sticky='NSEW',padx=10,pady=10)
        self.queuelist.tag_configure('completed',background = 'green')
        self.queuelist.config(columns = ['queue'])
        self.queuelist.heading('#0',text = 'Rundir')
        self.queuelist.heading('queue',text = 'Order')
        self.queuelist.column('#0',width = 150)
        self.queuelist.column('queue',width = 35)
        
        self.refresh_list()    
        
    def check_queue(self):
        '''
        Automatically updates the GUI Queue display and refreshes the DF
        whenever an item from the Queue is completed or new items are added
        '''
#       print( 'checking queue')
        try:
            dfname, peakid, rundir, TCPath, make_new = outputs.get(block=False)
            if make_new:
                self.refresh_list() #new table means new item for namelist
            #repopulate the alpha table
            index = self.namelist.get(0,'end').index(dfname)
            self.namelist.selection_set(index)
            print('populating table')
            self.pop_table()
            self.refresh_queue(completed = rundir) #refresh GUI display of queue
            self.master.after(2500,self.check_queue)
        except queue.Empty:
            self.master.after(2500, self.check_queue)
        except ValueError:
            print('DFName {} not found in namelist, couldn\'t update table'.format(dfname))
            dfname, peakid, rundir, TCPath, make_new = outputs.get(block=False)
            self.refresh_queue(completed = rundir)
            self.master.after(2500,self.check_queue)
    
    def pop_table(self):
        '''
        Populates the DataFrame display in the GUI according to selection
        parameters.
        '''
        
        self.table.grid_remove() #remove the table from the GUI
        self.table.delete(*self.table.get_children())
        self.update_dfname()
        dfname = self.dfname
        table = self.table
        with IOLock:
            with pd.HDFStore(alphapath, 'a') as hdf:
                df = hdf[dfname]
        '''
            self.df = df = pd.DataFrame()
            with pd.HDFStore(alphapath,'r') as hdf:
                  for key in keys:
                        df = df.append(hdf[key])
            '''
        columns = [str(header) for header in df.columns]
        table.config(columns=columns)
        
        #default widths for the DF columns
        widths = [60,27,37,56,59,50,59,36,52,50,20]
        table.column('#0',width = 30,stretch = False)
        for i,col in enumerate(columns):
            table.heading(col, text=col, command=lambda _col=col:
                          treeview_sort_column(table, _col, False))
            if i == 10:
                table.column(col, width=widths[i],stretch = True)
            else:
                table.column(col, width=widths[i],stretch = False)

        for peakid in df.index:
            table.insert('', 'end', peakid, text=peakid)
            i = 0
            for col in columns:
                if i < 2:
                    table.set(peakid, col, df.loc[peakid][col])
                    i += 1
                else:
                    try:
                        value = float(df.loc[peakid][col])
                        value = '{:0.3g}'.format(value)
                        table.set(peakid, col, value)
                    except (ValueError, TypeError):
                        table.set(peakid, col, df.loc[peakid][col])
        #put the new updated table back in the GUI
        self.table.grid(row=0, column=2, rowspan=10, sticky='NSEW',padx = 10,pady = 10)
    
    def new_entry(self):
        '''
        Initialize selection windows to process new data and format it
        to be passed along to the Queue.
        '''
        self.update_dfname()
        dfname = self.dfname
        make_new = self.make_new.get()
        multiple_dirs = self.multiple_dirs.get()
        
        if dfname == '':
            print('Do not use an empty DFName')
            return
        
        #Collect a list of the run files if multiple are being used
        if multiple_dirs:
            path = askdirectory(initialdir=self.runlist,
                                title='Select a folder containing multiple runs')
            if not path:
                print('No Directory selected.')
                return
            
            self.runlist = path
            runs = []
            #only selftrig runs with biasV greater than 46
            for run in os.listdir(path):
                if re.match('SelfTrig', run):
                    try:
                        x = int(run.split('_')[4][:-1])
                        if x > 46:
                            runs.append(path + '/' + run)
                    except:
                        print(run)
                        continue
            thisrun = runs[0]
        
        #Otherwise use only the specific file selected.
        elif not multiple_dirs:
            rundir = askdirectory(initialdir=self.rundir,
                                  title='Choose a run to process')
            if not rundir:
                print('No Directory selected.')
                return
            self.rundir = rundir #store last chosen directory
            thisrun = rundir
                
        #get list of paths that have been processed in this table
        with IOLock:
            with pd.HDFStore(alphapath, 'r') as hdf:
                try:
                    paths = list(hdf[dfname]['wfpath'])
                except KeyError:
                    paths = []
        
        #ignore disk name when checking
        queued = [task[2] for task in tasks.queue]
        if multiple_dirs:
            runs = [run[2:] for run in runs if run[2:] not in paths and run[2:] not in queued]
            if len(runs) == 0 or (len(runs) == 1 and runs[0] in self.get_items()):
                print("All of these runs have been processed in this table.")
                return
        
        elif rundir[2:] in paths or rundir[2:] in queued:
            print("This run has already been processed in this table.")
            return
        
        TCPath = askopenfilename(initialdir=self.tcpath,
                                 title='Choose the TC_tests file for run: {}'.format(thisrun))
        
        self.tcpath = '/'.join(TCPath.split('/')[:-1]) #set the most recent location of TCdata
        if self.tcpath == '':
            print("No TC file selected.")
            return
        
        #Get a unique ID# for each new run
        with IOLock:
            with pd.HDFStore(peakspath, mode='a') as hdf:
                if multiple_dirs:
                    idlist = []
                    for _ in runs:
                        while True:
                            peakid = 'h' + uuid4().hex[-8:]
                            keys = [key.split('/')[-1] for key in hdf.keys()]
                            if peakid in keys:
                                continue
                            elif peakid in idlist:
                                continue
                            else:
                                idlist.append(peakid)
                                break
                elif not multiple_dirs:
                    while True:
                        peakid = 'h' + uuid4().hex[-8:]
                        keys = [key.split('/')[-1] for key in hdf.keys()]
                        if peakid in keys:
                            continue
                        else:
                            break
        
        #Give the tasks to the queue
        if multiple_dirs:
            for i in range(len(runs)):
                task = (dfname, idlist[i], runs[i], TCPath, make_new)
                tasks.put(task) #place in threadsafe queue
            self.refresh_queue(new = runs) #now update queuelist (GUI display)
        
        #or just one
        else:
            task = (dfname, peakid, rundir, TCPath, make_new)
            tasks.put(task)
            self.refresh_queue(new = rundir)
    
    def update_dfname(self):
        '''
        Assign current entry to the DFName variable or collect the selected key.
        '''
        if self.make_new.get():
            self.dfname = self.entry3.get()
        else:
            self.dfname = self.namelist.selection_get()
            
    
    def refresh_list(self):
        '''
        Refresh the list of DFNames that are shown.
        '''
        self.namelist.delete(0, tk.END)
        with IOLock:
            with pd.HDFStore(alphapath, 'a') as hdf:
                for dfname in hdf.keys():
                    self.namelist.insert('end', dfname)
        self.refresh_queue()
        
        
    def del_frame(self):
        '''
        Deletes the currently selected DataFrame.
        '''
        MsgBox = messagebox.askquestion ('Delete table','Are you sure you want to delete this *entire* frame?',icon = 'warning')
        if MsgBox == 'no':
            return
        dfname = self.namelist.selection_get()
        with IOLock:
            with pd.HDFStore(alphapath, 'a') as hdf:
                hdf.remove('/' + dfname)
        self.refresh_list()

    def display_hist(self):
        '''
        Plots a histogram corresponding to the datapoint selected in
        the displayed DataFrame.
        '''
        dfname = self.dfname
        peakid, = self.table.selection()
        key = dfname + '/' + peakid
        info = self.table.item(peakid)['values']
        date = info[0]
        sigma = info[5]
        mdpt = info[3]
        biasv = info[2]
        temp = info[7]

        with IOLock:
            with pd.HDFStore(peakspath, 'r') as hdf:
                DF = hdf[key]
        
        #x0 = float(mdpt)
        #b = float(sigma)
        '''
        def visualgauss(x,a):
            return a*np.exp(-(x-x0)**2/(2*b**2))
        '''
        
        peaks = DF['peak_height']
        num = len(peaks)

        legend = '# of Entries: {} \n Midpoint: {}V \n Sigma: {}'.format(
            num, mdpt, sigma)
        
        
        hist = plt.hist(peaks, bins=binnum, density=False, color='b')
        
        bins = hist[1]
        centers = (bins[1:] + bins[:-1])/2
        counts = hist[0]
        
        
        p0 = (100,) + sp.stats.norm.fit(peaks)
        popt,pcov = sp.optimize.curve_fit(gauss,centers,counts,p0 = p0,bounds = (0,np.inf),maxfev = 10000)
        
        i =  -1
        while True:
            i -= 1
            if centers[i] < popt[1] - 2.5*popt[2]:
                break
        left = i
        
        i = 0
        while True:
            i += 1
            if i == len(centers):
                break
            elif centers[i] > popt[1] + 2.5*popt[2]:
                break
        right = i
        
        
        centers = np.array(centers[left:right])
        counts = np.array(counts[left:right])
        res = counts - gauss(centers,*popt)
        plt.plot(centers,res)
        res = res**2
        res = res/gauss(centers,*popt)
        Chi2 = res.sum()/len(centers)
        
        
        legend2 = 'Midpoint of Fit: {:0.3} \n Sigma of Fit: {:0.3} \n Chi Squared:\
        {:0.3}'.format(popt[1],popt[2],Chi2,pcov[1][1],pcov[2][2])
        # \n Midpoint Variance: {:0.3} \n Sigma Variance: {:0.3}\
        
        
        x = np.linspace(0, hist[1][-1], num=1000)
        plt.plot(x, gauss(x,*popt), 'r-')
        plt.title('Alpha Peak Histogram \n Date: {}, BiasV: {}, Temperature: {}K'.format(
            date, biasv, temp))
        plt.xlabel('Peak Amplitude')
        plt.ylabel('Counts')
        plt.legend([legend2])#,legend])
        plt.show()
        
    def save_to_clipboard(self):
        '''
        Saves the current displayed DataFrame to the clipboard.
        '''
        with IOLock:
            with pd.HDFStore(alphapath, 'r') as hdf:
                DF = hdf[self.dfname]
            DF.to_clipboard()
    
    def del_row(self):
        '''
        Deletes the selected row from the displayed DataFrame.
        '''
        MsgBox = messagebox.askquestion ('Entry Deletion','Are you sure you want to delete this row?',icon = 'warning')
        if MsgBox == 'no':
            return
        dfname = self.namelist.selection_get()
        peakid, = self.table.selection()
        dfname = self.dfname
        with IOLock:
            with pd.HDFStore(alphapath, 'a') as hdf:
                hdf.remove(dfname, where='index == {}'.format(peakid))
        self.pop_table()
    
    
    def refresh_queue(self,new = None,completed = None):
        '''
        Constructs and paints the Queuelist displayed in the DF for keeping
        track of what runs have been queued up and processed.
        '''
        queuelist = self.queuelist
        queuelist.column('#0', width=150)
        queuelist.column('queue', width=35)
        
        if isinstance(new,list):
            parents = [item.split('/')[3:]  for item in new]
            parents = zip(*parents)
            children = parents[-1]
            parents = [item[0] for item in parents[:-1]]
            parents = [parents[:i+1] for i in range(len(parents))]
            
            #the iid is it's path in queulist, like '20190101/Waveforms/Alpha'
            #except the rundir's whose iid are FULL path on disk (for easy comparison)
            
            for i in range(len(parents)):
                if queuelist.exists('/'.join(parents[i])): # check if that parent exists already
                    continue
                elif i == 0: # if it's the first entry, put it in the toplevel
                    queuelist.insert('','end',parents[i][0],text = parents[i][-1]) 
                else: # otherwise put it in its parent
                    queuelist.insert('/'.join(parents[i-1]),'end','/'.join(parents[i]),text = parents[i][-1])
            for i in range(len(new)): #the iid here is the FULL rundir, for uniqueness (treeview items have global id regardless of parent)
                queuelist.insert('/'.join(parents[-1]),'end',new[i],text = children[i])
        
        #the same deal but when only one new entry is selected
        elif isinstance(new,str):
            parents = new.split('/')[3:]
            for i in range(len(parents)):
                if i == 0: #top level parent
                    queuelist.insert('','end',parents[i],text = parents[i])
                elif i == len(parents) - 1: #actual item, use new as iid
                    queuelist.insert(parents[i-1],'end',new,text = parents[i])
                else: #midlevel parents
                    queuelist.insert(parents[i-1],'end',parents[i],text = parents[i])
        
        items = self.get_items() #retrieve all items in queuelist
        
        queued = []
        for task in tasks.queue: #list queued tasks
            queued.append(task[2]) #get the directory from the task
        
        #set queue value to match order of processing
        for item in items:
            if item in queued:
                queuelist.item(item,values = [queued.index(item)])
            else:
                queuelist.item(item,values = []) #queue = '' for items not in queue
        
        if completed:
            queuelist.item(completed,tag = 'completed')
        
    def get_items(self,item = ''):
        '''
        Return simplified list of all the items in the queuelist Treeview.
        '''
        queuelist = self.queuelist
        thislist = []
        if len(queuelist.get_children(item)) == 0:
            thislist.append(item)
        else:
            for child in queuelist.get_children(item):
                if child == '':
                    continue
                else:
                    thislist += self.get_items(child)
        return thislist

app = LXeDataManager()
s = ttk.Style()
s.configure('Treeview')
app.mainloop()

killthread = True