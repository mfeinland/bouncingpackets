import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import sampex
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.signal import find_peaks

'''To be honest I don't think this will work unless you import the other .py functions in the time evolution folder'''
'''but I haven't tested this because I work in Jupyter Notebook. :/ sorry'''

# Read in data
events = pd.read_csv('inner_belt_catalog.csv')

# See if you need to load each HILT object (takes awhile, better if you can avoid it)
if 'nums' in globals():
    old_nums = nums
else:
    old_nums = []
nums = range(0, 45) # events to look at

if (old_nums == nums) & (len(nums)==1):
    load_each = False
else:
    load_each = True

w = []
s = []
t = []

for i in nums:
    print('Microburst', i)
    date = datetime.strptime(events.date[i], "%m/%d/%Y") # change to list iteration for all events
    sstr = events.date[i] + ' ' + events.time[i]
    timestamp = datetime.strptime(sstr, "%m/%d/%Y %H:%M:%S")
    
    if load_each:
        h = sampex.HILT(date)
        h.load()

    obj = ID_Microbursts(h, timestamp)
    
    bad = False # assume everything works at first

    try:
        hd, tr = obj.find_microbursts()
    except:
        bad = True
        print('bad')
        
    if not bad:

        fit = obj.props # data about each curve fit

        # remove poorly-fitted peaks
        good_pks = np.where(fit['r2'] > 0.9)[0]

        # Initialize an empty list to store filtered results
        good_pks_f = []

        # Iterate through good_x and only keep values that are adjacent (spaced by 1 index)
        for i in range(0, len(good_pks)):
            if i == len(good_pks)-1:
                if good_pks[i] - good_pks[i-1] == 1:
                    good_pks_f.append(good_pks[i])
            else:
                if good_pks[i+1] - good_pks[i] == 1:  # Only keep if spaced by 1 index
                    good_pks_f.append(good_pks[i])

        fig, ax = plt.subplots()
        
        bad_2 = False
        try:
            gaus_overlay_plot(fit, hd, ax, tr)
            plt_name = date.strftime('%Y_%m_%d.jpg')
            fig.savefig('fit pics/' + plt_name)
        except:
            bad_2 = True # failed
            print('bad_2')
        
        if not bad_2:
            w_i = fit['fwhm'].values # peak widths
            w_i = np.round(w_i[good_pks_f], 4) # restrict to good pks only, and round values
            w.append(list(w_i)) # add to big array

            # peak spacings
            dt = np.diff(fit.t0)
            s_i = [x.total_seconds() for x in dt]
            s_i.append(None) # add one bc length is len(w)-1 since it's spacing BW peaks
            # need this so that the next line doesn't throw an out-of-bounds index error.
            s_i = [s_i[x] for x in good_pks_f]
            s.append(s_i)
            t.append(timestamp)

w_s = pd.DataFrame({'t': t, 'w': w, 's': s})
# w_s.to_csv('w_s.csv') # uncomment save to .csv if desired
