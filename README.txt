# SiPMAnalysisApp

This is the LXe Applet and accompanying data. Store the parent folder Xe in your root directory C:/
Then run MAIN.PY in a python terminal or IDE.


Instructions:

AlphaViewer - Use this to view the data that has been processed and to process new data

Refresh: refreshes alpha frame list
Display Data: displays selected alpha frame
Save Table to Clip: saves currently displayed table to clipboard (try pasting into excel table)
Delete Frame: delete currently selected frame

New Entry: prompts user to select waveform containing file, puts processed data into currently displayed waveform
Create New Alpha DF: if checked, puts new waveform in new dataframe with name in adjacent box
Select Multiple Runs: if checked, "New Entry" will prompt user to select folder of runs, (folder
containing waveform containing folders- will only select "SelfTrig" runs that have not been processed yet)

Display Histogram: Displays Histogram of currently selected row in table
Delete Row: Deletes currently selected row from the table

----------------------------------------------------------

RunFitViewer - Use this to plot the data from the tables

Separation: Select a separation
Bias V: Select a bias voltage
Teflon: if checked, only shows data after Teflon reflectors were added
Postbaking: if checked, only shows data after the baking incident
Refresh: after choosing previous 4 settings, click refresh to show valid dates
Date: select which data you want to view (if there are multiple consecutive days,
this will plot data from all of them regardless of which you choose)
Plot Data: plots all data matching selections and accompanying fit
