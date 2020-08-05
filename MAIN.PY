#%%
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
import matplotlib.cm as cm
#Analysis
import datetime
import pytz
import pandas as pd
import numpy as np
from numpy.polynomial.polynomial import Polynomial
#I/O commands and file management
from uuid import uuid4
import re
import os
import analysis
import copy
#%%
#matplotlib.use("TkAgg", warn=False)
style.use("ggplot")
#import shutil
polyfit = Polynomial.fit
strptime = datetime.datetime.strptime
timedelta = datetime.timedelta
from tempfile import NamedTemporaryFile
import webbrowser

#temproot = tk.Tk()
#Xepath = askdirectory(
#    title='Please select the Xe folder which you would like to work in')
Xepath = 'C:/Xe'
peakspath = Xepath + '/peaks.h5'
alphapath = Xepath + '/alphas.h5'
recordspath = Xepath + '/records.h5'
temppath = Xepath + '/temps.h5'
#temproot.destroy()

LARGE_FONT = ('Verdana', 10)

binnum = 'fd'

def df_window(df):
    with NamedTemporaryFile(delete = False,suffix = '.html') as f:
        df.to_html(f.name)
    webbrowser.open(f.name)

#%%
# =============================================================================
# Threading Functions
# =============================================================================

killthread = False
IOLock = threading.Lock()
tasks = queue.Queue()
outputs = queue.Queue()

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

#%%
# =============================================================================
# I/O Functions
# =============================================================================

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
    peaksDF = analysis.process_waveforms(dfname, rundir)
    print('Processing TC Data')
    tempsDF = process_TC(TCPath, rundir)
    alphaDF = analysis.alpha_df(rundir, peaksDF, tempsDF, peakid,
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

def process_TC(TCPath, rundir):
    '''
    Function to return DataFrame of Temperature info for a run.
    '''

    times = []
    temps = []
    i = 0
    runtime = analysis.get_dt(rundir)
    delta = timedelta(minutes=2)
    
    while True:
        with open(TCPath, 'r') as f:
            lines = f.readlines()

        for line in lines:
            i += 1
            if line.isspace():
                continue
            linetime = analysis.line_to_dt(line)
            if linetime < runtime < linetime + delta:
                temps.append(float(line.split()[-1]))
                times.append(analysis.line_to_dt(line))
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

#%%
# =============================================================================
# GUI Functions
# =============================================================================

norm = matplotlib.colors.Normalize(vmin=46, vmax=54, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)    

def popupmsg(msg):
    popup = tk.Tk()
    
    def leavemini():
        popup.quit()
        popup.destroy()
        
    
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg)
    label.pack(side="top", fill="x", padx=10, pady=10)
    B1 = ttk.Button(popup, text="Okay", command=leavemini)
    B1.pack()
    popup.mainloop()

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

class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        cy = self.widget.winfo_height()
        x = self.widget.winfo_rootx()
        y = cy + self.widget.winfo_rooty()
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                      background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        if widget.winfo_toplevel().tooltips == True:
            toolTip.showtip(text)
        else:
            pass
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

def TempPlotOnClick(event):
    artist = event.artist
    histid = artist.get_xdata().index[event.ind[0]]
    sep = artist.get_xdata().sep
    if sep == '31.07':
        sep = '31'
    key = f'/VUV4_{sep}mm/{histid}'
    
    with IOLock:
            with pd.HDFStore(temppath, 'r') as hdf:
                DF = hdf[key]
            
            with pd.HDFStore(alphapath,'r') as hdf:
                cont = key.split('/')
                dat = hdf[cont[1]].loc[cont[2]]
    
    date = dat.date
    biasv = dat.biasV
    temp = dat.temperature_avg
    
    plt.figure()
    temps = DF['temperature']
    time = DF['time']
    timeaxis = time - time[0]
    plt.plot(timeaxis,temps,color='#e25822')
    time = datetime.datetime.fromtimestamp(int(time[0]))
    time = time.strftime('%m/%d/%Y %H:%M:%S.%f')
    plt.title('Temperature Time Series \n Date: {}, BiasV: {}, Avg Temperature:{:0.3g}\n'
              'True Datetime: {}'.format(date, biasv, temp, time))
    plt.xlabel('Time (S)')
    plt.ylabel('Tempearture (K)')
    plt.show()

# =============================================================================
# GUI Classes
# =============================================================================
#%%
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
        
        self.tooltips = True
        self.sigfigs = {'biasV':3,'midpoint':3,'midpt_error':1,'sigma':3,\
           'sigma_error':1,'temperature_avg':5,'temperature_rms':1}
        self.widths = {'date':60,'separation':27,'biasV':37,'midpoint':56,
                       'midpt_error':59,'sigma':50,'sigma_error':59,
                       'temperature_avg':36,'temperature_rms':52,
                       'setpoint':50,'wfpath':20}
        
        menu = tk.Menu(container)
        self.config(menu=menu)
        filetab = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='File', menu=filetab)
        filetab.add_command(label='Exit', command=self.app_exit)
        
        settingstab = tk.Menu(menu, tearoff=0)
        menu.add_cascade(label='Settings', menu=settingstab)
        settingstab.add_command(label='Toggle Tooltips',
                                command=lambda:self.toggle('tooltips'))
        settingstab.add_command(label='DataFrame Precision',
                                command=self.sigfigprompt)
        settingstab.add_command(label='DataFrame Column Widths',
                                command=self.widthprompt)
        
        self.frames = {}
        for F in (StartPage,  AlphaViewer, RunFitting, SlowControlViewer):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky='NSEW')
        self.show_frame(StartPage)
    
    def toggle(self,attr):
        new = not getattr(self,attr)
        setattr(self,attr,new)
    
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
    
    def app_exit(self):
        self.quit()
        self.destroy()
    
    def sigfigprompt(self):
        prompt = tk.Toplevel()
        prompt.title("Significant Figure Settings")
        
        entries = copy.copy(self.winfo_toplevel().sigfigs)
        
        msg = tk.Message(prompt, text='Enter the desired number of significant'
                         ' figures for each column in the AlphaViewer data'
                         'frames.', aspect=450, anchor='nw')
        msg.grid(row=0,column=0,columnspan=len(entries))
        
        labelrow=1
        entryrow=2
        for i,col in enumerate(entries):
            if i > 3:
                labelrow = 3
                entryrow = 4
            entry = ttk.Entry(prompt)
            entry.insert('end',entries[col])
            entries[col] = entry
            entry.grid(row=entryrow, column=i%4)
            label = ttk.Label(prompt, text=col)
            label.grid(row=labelrow, column=i%4)
        
        sigfigs = self.winfo_toplevel().sigfigs
        def save():
            for col,entry in entries.items():
                sigfigs[col] = entry.get() 
            prompt.destroy()
        
        def cancel():
            prompt.destroy()
        
        savebutton = ttk.Button(prompt, text='Save', command=save)
        cancelbutton = ttk.Button(prompt, text='Cancel', command=cancel)
        savebutton.grid(column=3,row=3, sticky='NSEW')
        cancelbutton.grid(column=3,row=4, sticky='NSEW')
    
    def widthprompt(self):
        prompt = tk.Toplevel()
        prompt.title("Set the Default Width of DataFrame Columns")
        entries = copy.copy(self.winfo_toplevel().widths)
        
        msg = tk.Message(prompt, text='Enter the desired widths of the'
                         ' columns in the AlphaViewer data frames.',
                         aspect=450, anchor='nw')
        msg.grid(row=0,column=0,columnspan=len(entries))
        
        for i,col in enumerate(entries):
            labelrow = 2*(i%3) + 1
            entryrow = labelrow + 1
            labelcol = int(np.floor(i/3))
            entrycol = labelcol
            label = ttk.Label(prompt, text=col)
            label.grid(row=labelrow, column=labelcol)
            entry = ttk.Entry(prompt)
            entry.insert('end',entries[col])
            entries[col] = entry
            entry.grid(row=entryrow, column=entrycol)
        
        widths = self.winfo_toplevel().widths
        def save():
            for col,entry in entries.items():
                widths[col] = entry.get() 
            prompt.destroy()
        def cancel():
            prompt.destroy()
        
        savebutton = ttk.Button(prompt, text='Save', command=save)
        cancelbutton = ttk.Button(prompt, text='Cancel', command=cancel)
        savebutton.grid(column=3,row=5, sticky='NSEW')
        cancelbutton.grid(column=3,row=6, sticky='NSEW')
        
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
        
        button3 = ttk.Button(self, text = "Visit Run Fit Viewer",
                                 command = lambda: controller.show_frame(RunFitting))
        button3.pack()
        
        button4 = ttk.Button(self, text="Visit Slow Control Data Viewer",
                             command = lambda: controller.show_frame(SlowControlViewer))
        button4.pack()


'''WaveformViewer
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
        self.controller = controller
        self.parent = parent
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
        separations.sort()
        self.biases = biases = list(set(df.biasV))
        biases.sort()
        i = 0
        
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
        
        self.refresh = ttk.Button(self,text = 'Refresh',
                                 command=self.refresh_dates)
        
        
        self.refresh.grid(row=4,column=0)
        
        self.f = Figure(figsize=(7, 7), dpi=100)
        self.ax = self.f.add_subplot(111)
        
        self.graphframe = tk.Frame(self)
        self.graphframe.grid(row=0, column=2, rowspan = 10,
                             sticky='NSEW',padx=10,pady=10)
        
        self.canvas = FigureCanvasTkAgg(self.f, master=self.graphframe)
        self.canvas.mpl_connect('pick_event', TempPlotOnClick)
        self.canvas.get_tk_widget().pack()
        
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.graphframe)
        self.toolbar.update()
        
        self.plotbutton = ttk.Button(self, text='Plot Data',
                                     command=self.plot_data)
        self.plotbutton.grid(row=4,column=1,sticky='EW',pady=1)
        
        
        homebutton = ttk.Button(self, text='Back to Home',
                                command=lambda: controller.show_frame(StartPage))
        homebutton.grid(row=5, column=0,pady = 1,sticky='EW')
        
        self.rowconfigure(9, weight=1)
        
        clearbutton = ttk.Button(self, text='Clear Plot',
                                 command = self.clearplot)#self.ax.clear)
        clearbutton.grid(row=5,column=1,pady=1,sticky='EW')        
        self.allbv = tk.IntVar()
        allbiases = ttk.Checkbutton(self, text='Plot All Biases',variable=self.allbv)
        allbiases.grid(row=6, column=0,pady = 1,sticky='EW')
        
        self.handles = dict()
        self.legend = dict()
        
        CreateToolTip(self.teflon, "If selected, only shows dates after"
                      " addition of the Teflon reflector")
        CreateToolTip(self.baked, "If selected, only shows dates after"
                      " the baking incident.")
        CreateToolTip(self.seplabel, "Separation between SiPM and source in mm.")
        CreateToolTip(self.plotbutton, "Plot the data that corresponds to the"
                     " highlighted date and other criteria. Joins data from\n"
                     " consecutive days into one fit. More runs can be put on"
                     " a single plot by changing the selection and plotting\n"
                     " without clearing first.")
        CreateToolTip(self.plotbutton, "Plots data for the selected date"
                     " at all voltages.")
        CreateToolTip(self.refresh, "Updates the datelist to contain the dates"
                     " which meet the other criteria.")
    def clearplot(self):
        self.ax.clear()
        self.legend = dict()
        self.canvas.draw()
    
        
    def refresh_dates(self):
        '''
        Refreshes the dates menu for selection. Only presents dates
        that meets the other given requirements.
        '''
        sep = self.sep = self.sepmenu.get()
        bv = self.biasmenu.get()
        tef = self.tef.get()
        postb = self.postb.get()
        df = self.df
        
        data = df[df.biasV == bv]
        data = data[data.separation == sep]
        
        if tef:
            data = data[data.date.astype(int) >= 20190605]
        
        if postb:
            data = data[data.date.astype(int) >= 20190424]
        
        
        dates = list(set(data.date))
        dates.sort()
        self.datelist['values'] = dates
    
    def plot_data(self):
        '''
        Uses the date given and other parameters to select and plot data
        and a linaer fit.
        '''
        bv = [self.biasmenu.get()]
        df = self.df
        
        if self.allbv.get():
            bv = self.biases
        
        
        date = int(self.datelist.get())
        dates = [date,date+1,date+2,date-1,date-2]
        dates = [str(item) for item in dates]
        
        for v in bv:
            data = df[df.biasV == v]
            data = data[data.date.isin(dates)]
            
            x = data.temperature_avg
            x.sep = self.sep
            xerr = data.temperature_rms
            y = data.midpoint
            yerr = data.midpt_error
            
            if len(x) == 0:
                continue
            
            self.ax.errorbar(x, y, yerr, xerr, marker='.',linestyle='none', color='black',picker=5)
            
            
            nans = y.notnull() & x.notnull()
            x = x[nans]
            y = y[nans]
            yerr = yerr[nans]
            xerr = xerr[nans]
            c,m,ac,am = analysis.blinfit(x,y,yerr)
            xax = np.linspace(min(x),max(x),num=5000)
            handle, = self.ax.plot(xax,m*xax + c,color=mapper.to_rgba(int(v[:-1])))
            self.handles.update({v:handle})
            
            legend = u'{} , Slope: {:0.2} \u00B1 {:0.2}'.format(v,m,am)
            self.legend.update({handle:legend})
            
            self.ax.legend(self.legend.keys(),self.legend.values(),
                           loc='upper center', bbox_to_anchor=(0.5, 1.15),
                           fancybox=True, shadow=True,prop={'size': 7},
                           ncol=2)
            self.ax.set_title('Date: {}, Bias: {}'.format(date,v))
        
        if self.allbv.get():
            self.ax.set_title('Data: {}, All Biases'.format(date))
        
        self.ax.set_xlabel('Temperature (K)')
        self.ax.set_ylabel('Alpha Peak (V)')
        
        
        self.canvas.draw()

class SlowControlViewer(tk.Frame):
    '''
    Applet page used for viewing the slow control data (temperature, flow
    pressure, etc.) in long form time series.
    '''
    
    _handletolabel = {
                        'RTD0 - Coldfinger': 'RTD0',
                        'RTD1 - Cell Bottom': 'RTD1',
                        'RTD2': 'RTD2',
                        'RTD3': 'RTD3',
                        'RTD4 - SiPM/Cell Top': 'RTD4',
                        'Flow Rate': 'Flow',
                        'Total Xenon': 'Flow Integration',
                        'PG1': 'Pressure 1',
                        'PG2': 'Pressure 2',
                        }
    
    _typetogroup = {'Temperature':'Temp_Agg', 'Flow':'FlPr_Agg',
                    'Pressure':'FlPr_Agg'}
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent, bg='LightSteelBlue4')
        
        self.f = Figure(figsize=(7, 7), dpi=100)
        self.temp_ax = self.f.add_subplot(111)
        
        self.graphframe = tk.Frame(self)
        self.graphframe.grid(row=0, column=2, rowspan=10, 
                             sticky='NSEW', padx=10, pady=10)
        
        
        self.canvas = FigureCanvasTkAgg(self.f, master=self.graphframe)
        self.canvas.mpl_connect('pick_event', TempPlotOnClick)
        self.canvas.get_tk_widget().pack()
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.graphframe)
        self.toolbar.update()
        
        homebutton = ttk.Button(self, text='Back to Home',
                                command=lambda: controller.show_frame(StartPage))
        homebutton.grid(row=0, column=0, pady = 5, sticky='NSEW')
        
        
        sensors = self.sensors = {
            'Temperature': ('RTD0 - Coldfinger', 'RTD1 - Cell Bottom',
                            'RTD2', 'RTD3', 'RTD4 - SiPM/Cell Top'),
            'Flow': ('Flow Rate', 'Total Xenon'),
            'Pressure': ('PG1', 'PG2')
            }
        
        
        
        typeio = self.typeio = {}
        
        self.menubutton = ttk.Menubutton(self, text='Data Selection')
        self.menubutton.grid(row=0, column=1, sticky='NSEW', pady=5)
        self.topmenu = tk.Menu(self.menubutton, tearoff=False)
        self.menubutton.configure(menu=self.topmenu)
        for dtype in sensors:
            menu = tk.Menu(self.topmenu, tearoff=False)
            self.topmenu.add_cascade(label=dtype, menu=menu)
            
            for sensor in sensors[dtype]:
                typeio[sensor] = tk.BooleanVar()
                menu.add_checkbutton(label=sensor, variable=typeio[sensor])
        
        self.startdate = ttk.Entry(self, )
        self.enddate = ttk.Entry(self)
        self.startdate.grid(row=1, column=1, sticky='NSEW', pady=5)
        self.enddate.grid(row=2, column=1, sticky='NSEW', pady=5)
        self.startdate.insert(0, '2000-12-30 16:34:56.561683')
        self.enddate.insert(0, '2000-12-31 1:12:21.843155')
        
        
        startlabel = ttk.Label(self, text='Start Date')
        endlabel = ttk.Label(self, text='End Date')
        startlabel.grid(row=1, column=0, pady=5)
        endlabel.grid(row=2, column=0, pady=5)
        
        self.rowconfigure(index=9, weight=1)
        
        self.plotbutton = ttk.Button(self, text='Plot Data', 
                                     command=self.plot_data)
        self.plotbutton.grid(row=3, column=0, columnspan=2, 
                             pady=10, sticky='NS')
    def plot_data(self):
        
        #clear the figure
        self.f.clf()
        
        #Construct list of sensors that are checked to plot
        sensors_to_plot = []
        typeio = self.typeio
        for sensor in typeio:
            if typeio[sensor].get() == True:
                sensors_to_plot.append(sensor)
        
        #Verify which sensor types to create subplots for
        types_to_plot = []
        sensors = self.sensors        
        for sensor_type in sensors:
            for sensor in sensors_to_plot:
                if sensor in sensors[sensor_type]:
                    types_to_plot.append(sensor_type)
                    break
        
        #handle 2 label dictionaries
        self.axes = dict() #{sensor_type1:ax1, ...}
        self.legend = dict() #{ax1:{artist1:label1, ...}, ...}
        
        #build each subplot
        for i,sensor_type in enumerate(types_to_plot):
            if i==0:
                first = ax = self.axes[sensor_type] = self.f.add_subplot(len(types_to_plot), 1, i+1)
                pass
            else:
                ax = self.axes[sensor_type] = self.f.add_subplot(len(types_to_plot), 1, i+1, sharex=first)
            ax.set_title(f'{sensor_type} Data')
            ax.set_xlabel(f'Time')
            ax.set_ylabel(sensor_type)
        
        #Retrieve date entries and convert to Unix time
        startdate = self.startdate.get()
        enddate = self.enddate.get()
        start = datetime.datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S.%f')
        end = datetime.datetime.strptime(enddate, '%Y-%m-%d %H:%M:%S.%f')
        tz = pytz.timezone('US/Eastern')
        start = tz.localize(start).astimezone(pytz.utc).timestamp()
        end = tz.localize(end).astimezone(pytz.utc).timestamp()
        
        #populate each subplot
        for sensor in sensors_to_plot:
            
            #determine sensor_type
            for sensor_type in sensors:
                if sensor in sensors[sensor_type]:
                    sensor_type = sensor_type
                    break
            
            #get data from HDFStore
            with pd.HDFStore('SlowCtrlAgg.h5', 'r') as hdf:
                data = hdf.select(self._typetogroup[sensor_type],
                                  where=f'index > {start} & index < {end}')
                
                x = data.index
                x = pd.to_datetime((x*1e9).astype('int64'), utc=True)
                x = x.tz_convert(pytz.timezone('US/Eastern'))
                x = x.date
                y = data[self._handletolabel[sensor]]
            
            #plot data
            ax = self.axes[sensor_type]
            artist = ax.scatter(x,y)
            
            #track handles to labels
            if sensor_type in self.legend:
                self.legend[sensor_type].update({artist:sensor})
            else: #intialize dict if first line on axes
                self.legend[sensor_type] = {artist:sensor}
        
        for sensor_type in self.axes:
            ax = self.axes[sensor_type]
            ax.legend(self.legend[sensor_type].keys(),
                      self.legend[sensor_type].values(),
                      loc='upper left',
                      fancybox=True, shadow=True, prop={'size': 7})
        
        self.f.tight_layout(pad=2)
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
        #CreateToolTip(self.queuelist, "A display showing all the data that has"
        #              " been queued up for processing. Drop down menus mirror"
        #              " the path location of the data. The DataFrame display"
        #              " will automatically refresh when a run is finished."
        #              " An I7 Intel CPU processed at a speed of about 2 minutes"
        #              " per run, so avoid queueing up too many at a time.")
        self.queuelist.grid(row=9,column=0,columnspan=2,sticky='NSEW',padx=10,pady=10)
        self.queuelist.tag_configure('completed',background = 'green')
        self.queuelist.config(columns = ['queue'])
        self.queuelist.heading('#0',text = 'Rundir')
        self.queuelist.heading('queue',text = 'Order')
        self.queuelist.column('#0',width = 150)
        self.queuelist.column('queue',width = 35)
        
        self.refresh_list()    
        
        
        CreateToolTip(button1, "Displays the DataFrame highlighted in the"
                      " Alpha Frame List.")
        CreateToolTip(button2, "Prompts user to select data for processing."
                      " By default the new data will be appended to the\n"
                      " highlighted DataFrame. The selected folder should"
                      " contain waveforms, like w1.txt or w1.zip.")
        CreateToolTip(self.button3, "If selected, any new entries will be"
                      " placed into a new DataFrame named in the entry to\n"
                      " the right. The forward slash is included automatica"
                      "lly. Refrain from using special characters, spaces, etc.")
        CreateToolTip(button4, "Refresh the list of DataFrames shown in"
                      " the Alpha Frame List.")
        CreateToolTip(button5, "Delete the selected DataFrame from the "
                      "Alpha Frame List. CAUTION: DELETION IS PERMANENT.\n"
                      "CREATE A BACKUP OF THE ALPHAS.H5 FILE BEFORE DELETING.")
        CreateToolTip(button6, "Display a Histogram of the peakdata used to"
                      " create the highlighted row in the DataFrame.")
        CreateToolTip(button7, "Copy the displayed DataFrame to the clipboard."
                      " This can be easily pasted into an excel spreadsheet.")
        CreateToolTip(button8, "Deletes the highlighted row from the DataFrame."
                      " To be used in case of a corrupted or duplicate row. To\n"
                      " be safe, create a backup alphas.h5 before deletion.")
        CreateToolTip(self.button9, "If selected, the New Entry button will"
                      " prompt you to select the parent folder of a set of\n"
                      " runs, subsequently adding all contained runs to the"
                      " Queue list below. This will automatically avoid\n"
                      " processing runs whose paths are already listed in the"
                      " DataFrame, in case you are resuming an interrupted\n"
                      " queue.")
        
    def check_queue(self):
        '''
        Automatically updates the GUI Queue display and refreshes the DF
        whenever an item from the Queue is completed or new items are added
        '''
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
            popupmsg('DFName {} not found in namelist, couldn\'t update table'.format(dfname))
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
            with pd.HDFStore(alphapath, 'r') as hdf:
                df = hdf[dfname]
        '''
            self.df = df = pd.DataFrame()
            with pd.HDFStore(alphapath,'r') as hdf:
                  for key in keys:
                        df = df.append(hdf[key])
            '''
        columns = [str(header) for header in df.columns]
        table.config(columns=columns)
        
        widths = self.winfo_toplevel().widths
        table.column('#0',width = 30,stretch = False)
        for i,col in enumerate(columns):
            table.heading(col, text=col, command=lambda _col=col:
                          treeview_sort_column(table, _col, False))
            if i == 10:
                table.column(col, width=widths[col],stretch = True)
            else:
                table.column(col, width=widths[col],stretch = False)
        
        for peakid in df.index:
            table.insert('', 'end', peakid, text=peakid)
            for i,col in enumerate(columns):
                if i < 2:
                    table.set(peakid, col, df.loc[peakid][col])
                else:
                    try:
                        value = float(df.loc[peakid][col])
                        value = '{0:0.{1}g}'.format(value,
                                        self.winfo_toplevel().sigfigs[col])
                        table.set(peakid, col, value)
                    except (ValueError, TypeError):
                        table.set(peakid, col, df.loc[peakid][col])
        #put the new updated table back in the GUI
        self.table.grid(row=0, column=2, rowspan=10, 
                        sticky='NSEW', padx = 10, pady = 10)
    
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
            popupmsg('Do not use an empty DFName')
            return
        
        #Collect a list of the run files if multiple are being used
        if multiple_dirs:
            path = askdirectory(initialdir=self.runlist,
                                title='Select a folder containing multiple runs')
            if not path:
                popupmsg('No Directory selected.')
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
                popupmsg('No Directory selected.')
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
                popupmsg("All of these runs have been processed in this table.")
                return
        
        elif rundir[2:] in paths or rundir[2:] in queued:
            popupmsg("This run has already been processed in this table.")
            return
        
        TCPath = askopenfilename(initialdir=self.tcpath,
                                 title='Choose the TC_tests file for run: {}'.format(thisrun))
        
        self.tcpath = '/'.join(TCPath.split('/')[:-1]) #set the most recent location of TCdata
        if self.tcpath == '':
            popupmsg("No TC file selected.")
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
        the displayed DataFrame. Also displays a histogram of the temperatures
        for the run.
        '''
        dfname = self.dfname
        peakid, = self.table.selection()
        key = dfname + '/' + peakid
        info = self.table.item(peakid)['values']
        date = info[0]
        biasv = info[2]
        temp = info[7]
        
        with IOLock:
            with pd.HDFStore(peakspath, 'r') as hdf:
                DF = hdf[key]
        
        peaks = DF['peak_height'] - DF['baseline']
        
        error = False
        try:
            popt,pcov,Chi2,valpeaks,centers = analysis.histogram(peaks,info=True)
        except RuntimeError:
            error = True
            popupmsg("Fit could not be found.")
            pass
        except ValueError:
            error = True
            popupmsg("Data is considered cut-off")
        
        plt.figure()
        plt.title('Alpha Peak Histogram \n Date: {}, BiasV: {}, Temperature: {}K'.format(
            date, biasv, temp))
        plt.xlabel('Peak Amplitude')
        plt.ylabel('Counts')
        
        if not error:
            plt.hist(peaks, bins=binnum, density=False, color='b')
            popt[2] = abs(popt[2]) #only pos. value for sigma
            
            legend = u'Midpoint of Fit: {:0.3} \u00B1 {:0.1} \n \
            Sigma of Fit: {:0.3} \u00B1 {:0.1} \n \
            Chi Squared: {:0.3} \n \
            Valid Entries: {}'.format(popt[1],np.sqrt(pcov[1][1]),
                                        popt[2],np.sqrt(pcov[2][2]),
                                        Chi2, valpeaks)
            x = np.linspace(centers[0], centers[-1], num=1000)
            plt.plot(x, analysis.gauss(x,*popt), 'r-')
            plt.legend([legend],loc='upper left')
        
        else:
            plt.hist(peaks, bins=100, density=False, color='b')
        
        
        
        with IOLock:
            with pd.HDFStore(temppath, 'r') as hdf:
                DF = hdf[key]
        plt.figure()
        temps = DF['temperature']
        plt.hist(temps, 20, density=False, color='#e25822')
        plt.title('Temperature Histogram \n Date: {}, BiasV: {}, Temperature:{}\
                  '.format(date, biasv, temp))
        plt.xlabel('Temperature (K)')
        plt.ylabel('Counts')
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
#s = ttk.Style()
#s.configure('Treeview')
app.mainloop()

killthread = True