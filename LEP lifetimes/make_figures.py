# Import packages, functions
import numpy as np
import pandas as pd 
import xarray as xr
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import geopandas as gpd
import sampex
import ast
import warnings
import os
import scipy.optimize
from scipy.signal import find_peaks, peak_widths
import sklearn.metrics
import random

## Functions declared up here

'''Gaussian curve fitting '''
# Most of this code was adapted from Mike Shumko's GitHub repository:
# https://github.com/mshumko/sampex_microburst_widths/tree/main/sampex_microburst_widths/microburst_id

class ID_Microbursts:
    # tell the class the properties you wanna store in here preliminarily.
    def __init__(self, h, packet_timestamp, prominence_rel_height=0.5):
        self.h = h # assign HILT data object to self.h
        self.packet_timestamp = packet_timestamp # assign packet start timestamp
        self.prominence_rel_height=prominence_rel_height # assign prominence height width evaluator
        #(when you find a peak, you measure the width at this % of the full prominence.
        # Default is 0.5, 0 would be the base, 1 would be the top.)
        return

    def obrien(self):
    # Author: Max Feinland, adapted from O'Brien et al. 2003

        def changeCadence(time20, times100):
            # Author: Max Feinland
            # times100 is the list of timestamps for whatever it is: starts, ends, ns
            # times20 is the entire list of time20 for the day in question
            try:
                # try list comprehension first, it's faster
                idx20 = [time20.get_loc(tstmp) for tstmp in times100] 
            except: 
                idx20 = [] # clear and start over
                idx20 = [np.abs((time20 - target).values).argmin() for target in times100] # if a timestamp is missing from time20,
                # you will have to find the minimum timestamp (that is closest in time)
            return idx20
    
        def findActualPeaks(idx20, counts):
            # Author: Max Feinland
            # Date created: 1/9/24
            actual_idx = np.zeros(len(idx20))
            idx20 = np.array(idx20)
            # Find the actual peak in the interval, since it could be anything in the 100-ms window
            start = idx20 - 5
            end = idx20 + 5
    
            for i in range(len(idx20)):
                try:
                    interval = counts[int(start[i]):int(end[i])]
                    actual_idx[i] = int(np.argmax(interval)) + start[i]
                except:
                    actual_idx = []
                    break
            return actual_idx
            
        def qualityCheck(t, x):
            # Author: Max Feinland
            # Date created: 12/24/24
            not_coarse = (len(np.unique(x)) > 5)*1 # are there more than 5 unique values in the interval?
            continuous = (len(np.where(abs(t.to_series().diff()) > timedelta(seconds=1))[0]) == 0)*1 # is the interval more or less continuous?
            return not_coarse, continuous
        
        N20 = self.h['counts'] # count rate sampled every 20 ms
        time20 = self.h['time'] # time every 20 ms
    
        # I'm sure there's a more efficient way to do this, but I don't know it
        df = pd.DataFrame({'time': time20, 'counts': N20})
    
        df.set_index('time', inplace=True) # set time column as the index
    
        # resample the dataframe to 100 ms intervals and sum the counts in each interval
        N100 = df.resample('100ms').sum()
    
        A500 = N100.rolling(5, center=True).mean() # 5-observation centered rolling mean (over 500 ms)
    
        condition = np.divide((N100.counts - A500.counts), np.sqrt(1 + A500.counts)) # O'Brien et al 2003
    
        ns = np.argwhere(condition > 5) # find indices where O'Brien condition is true
        ns = [item[0] for item in ns]
    
        ns20 = changeCadence(time20, N100.index[ns]) # change the cadence from 100ms to 20ms
    
        final_idx = findActualPeaks(ns20, self.h['counts']) # run actual peak finder
        final_idx = np.unique(final_idx).astype(int) # get rid of repeats
    
        qc = {'notcoarse': [], 'continuous': []} # perform quality check
    
        for i in range(len(final_idx)):
            # Create 1-second interval to run quality check
            sidx = int(final_idx[i]-25)
            eidx = int(final_idx[i]+25)
    
            # Pull out variables in this interval
            t = time20[sidx:eidx]
            x = self.h['counts'][sidx:eidx]
    
            # Run quality check
            not_coarse, continuous = qualityCheck(t, x)
    
            # Append results to the dictionary
            qc['notcoarse'].append(not_coarse)
            qc['continuous'].append(continuous)
    
        qc = pd.DataFrame(qc) # convert dictionary to dataframe
        passed_qc = np.where((qc['notcoarse'] == 1) & (qc['continuous'] == 1))[0]
        self.mbidx = final_idx[passed_qc]
        return
        
    
    def find_microbursts(self):
        # Author: Max Feinland
        # Date created: 10/11/2024
        # Last modified:
        # Purpose: finds microbursts, returns indices containing peaks.
        
        self.obrien() # call O'Brien function-- now self contains so, eo, no
        #(start times as identified by O'Brien, end times as identified by O'Brien, 
        # all indices that meet the condition as identified by O'Brien)
        
        t = self.h['times'] # pull out time data
        rate = self.h['counts'] # pull out count rate data
        
        # Returns the index corresponding to the start of the O'Brien ID'd packet.
        # mbidx stands for "microburst index" since I know that's not really clear
        bpidx = np.abs(t[self.mbidx] - self.packet_timestamp).argmin()
        
        interval = rate[self.mbidx[bpidx]-10:self.mbidx[bpidx]+100] # define just one chunk of the rate data, 
        # taken from O'Brien interval
        tint = t[self.mbidx[bpidx]-10:self.mbidx[bpidx]+100]
        
        maxdis = max(interval) - min(interval) # use to generate prominence requirement

        # finding peaks with a high enough prominence
        [pks, props] = find_peaks(interval, prominence=0.25*maxdis, distance=3, wlen=10)
        
        # find instances where a bunch of peaks are evenly spaced (at least 4 peaks)
        indices = np.where(np.convolve(np.abs(np.diff(np.diff(pks))) <= 3, \
                                               np.ones(2), mode='valid') == 2)[0]
        indices = np.append(indices, [indices[-1]+1, indices[-1]+2, indices[-1]+3])

        pks = pks[indices] # pks is a list of indices in the HILT data that has a peak
        # that meets all these conditions. (o'brien, high enough prominence, evenly spaced)

        # Add peak_idx to self-- a list of indices in the HILT dataframe containing peaks
        # (add self.so[mbidx] b/c that's the beginning of the interval)
        self.peak_idx = list(self.mbidx[bpidx]-10 + loc for loc in pks)

        # Convert HILT objects to dataframe
        hilt_df = pd.DataFrame(data={'counts':self.h['counts']}, index=self.h['time'])

        # call an instance of Do_Gaussian class using the day's HILT data, 
        # the list of peak indices, and the prominence width 
        # (how high up the peak, 0 to 1 where 0 is the base and 1 is the top, to evaluate the width)
        gaus = Do_Gaussian(hilt_df, self.peak_idx, self.prominence_rel_height)
        
        # Calculate widths using scipy.signal.peak_widths
        widths_ok = False # assume estimate is bad at first
        counter = 0 # while debugging
        
        while not widths_ok:
            gaus.calc_prominence_widths() # calculate widths with gaussian method
            
            w = gaus.right_peak_base - gaus.left_peak_base # widths for each peak
            widths_ok = all(t <= np.timedelta64(1, 's') for t in w) # are all widths small enough? (if true, exit while loop)

            # if they're not, this line will change the calculation for the next iteration
            gaus.prominence_rel_height = gaus.prominence_rel_height - 0.1
            counter +=1
            
            if counter >= 10:
                widths_ok=True
                print("Errored out")
        
        # Using this initial estimate, calculate widths using Gaussian fit
        fit_df, time_range = gaus.calc_gaus_widths()

        self.props = fit_df # assign attribute 'props' to self
        # (contains information about the fit of the peaks: r^2, amplitude, center, width, etc.)
        # return time range
        return hilt_df, time_range


class Do_Gaussian:
    def __init__(self, hilt_df, peak_idx, prom_h, width_multiplier=3.5, plot_width_s=5):
        """
        Initialize Do_Gaussian class instance with these variables
        """
        self.hilt_df = hilt_df # HILT dataframe
        self.peak_idx = peak_idx # peak indices
        self.prominence_rel_height = prom_h # prominence height to evaluate width
        self.width_multiplier = width_multiplier
        self.plot_width_s = plot_width_s
        self.hilt_times = self.hilt_df.index.to_numpy()
        return

    def calc_prominence_widths(self):
        # Author: Mike Shumko
        # Date created: unknown
        # Last modified: 10/11/2024 by Max Feinland
        # Purpose: Use scipy to find the peak width at self.prominence_rel_height prominence.

        # Check that self.peak_idx correspond to the max values
        peak_check_thresh = 5 # Look 100 ms around the peak count to find the true peak. 
        for i, index_i in enumerate(self.peak_idx):
            self.peak_idx[i] = index_i - peak_check_thresh + \
                np.argmax(self.hilt_df['counts'][index_i-peak_check_thresh:index_i+peak_check_thresh])

        # Use scipy to find the peak width at self.prominence_rel_height prominence
        widths_tuple = peak_widths(self.hilt_df['counts'], self.peak_idx, 
                                            rel_height=self.prominence_rel_height)
        
        self.prom_widths_s = 20E-3*widths_tuple[0]   
        self.width_height = widths_tuple[1]
        self.left_peak_base = self.hilt_times[np.round(widths_tuple[2]).astype(int)]
        self.right_peak_base = self.hilt_times[np.round(widths_tuple[3]).astype(int)]  
            
        return

    def calc_gaus_widths(self, detrend=True):
        # Author: Mike Shumko
        # Date created: unknown
        # Last modified: 10/11/2024 by Max Feinland
        # Purpose: Use prior estimates from scipy estimation to find peak widths, amplitude, 
        # center, opt. linear trend using Gaussian fit.
        
        if not hasattr(self, 'prom_widths_s'):
            raise AttributeError('No prior width estimate exists. Run the '
                                'calc_prominence_widths method first.')

        # Create empty pd.DataFrame for fit variables
        fit_param_names = ['r2', 'adj_r2', 't0_err', 'A', 't0', 'fwhm']
        if detrend:
            fit_param_names.extend(['y-int', 'slope'])
        df = pd.DataFrame(data={key:np.nan*np.ones_like(self.peak_idx) 
                        for key in fit_param_names}, index=self.peak_idx)

        # Loop over every peak
        time_total = []
        for i, (peak_i, width_i, height_i) in enumerate(zip(self.peak_idx, self.prom_widths_s, self.width_height)):
            time_range = [self.left_peak_base[i] - pd.Timedelta(seconds=0.06), 
                          self.right_peak_base[i] + pd.Timedelta(seconds=0.08)]

            width_i = pd.Timedelta(time_range[1] - time_range[0]).total_seconds()/1.5

            t0 = self.hilt_times[peak_i]

            if detrend:
                p0 = [
                    height_i,   # gauss amplitude 
                    t0,         # gauss center time
                    width_i,    # 2x gaus std.
                    self.hilt_df.loc[time_range[0]:time_range[1], 'counts'].median(), # y-intercept
                    0           # Slope
                    ]
            else:
                p0 = [height_i, t0, width_i]
                

            # Catch warnings
            with warnings.catch_warnings(record=True) as w:
                # Catch exceptions
                try:
                    popt, pcov, r2, adj_r2 = self.fit_gaus(time_range, p0)
                    t0_err = pcov[1] # error in location
                            
                except RuntimeError as err:
                    if ('Optimal parameters not found: Number of calls '
                        'to function has reached maxfev') in str(err):
                        continue
                    raise
                if len(w):
                    print(w[0].message, '\n', p0, popt)

            # Save to a pd.DataFrame row.
            if popt[2] < 0:
                popt[2] = popt[2]*-1
                
            area = popt[0]
            df.iloc[i, :2] = r2, adj_r2
            df.iloc[i, 2] = t0_err
            df.iloc[i, 3:] = popt 

            time_total.append(time_range)
        return df, time_total

    def fit_gaus(self, time_range, p0):
        # Author: Mike Shumko
        # Date created: unknown
        # Last modified: 10/11/2024 by Max Feinland
        # Purpose: Fits a gaussian shape with an optional linear detrending term.

        x_data = self.hilt_df.loc[time_range[0]:time_range[1], :].index
        current_date = x_data[0].floor('d')
        x_data_seconds = (x_data-current_date).total_seconds()
        y_data = self.hilt_df.loc[time_range[0]:time_range[1], 'counts']

        if len(x_data) < len(p0):
            raise ValueError('Not enough data points to fit. Increase the '
                            'time_range or self.width_multiplier')

        p0[0] *= 2
        p0[1] = (p0[1] - current_date).total_seconds()
        p0[2] = p0[2]/2 # Convert the microburst width guess to ~std.

        if p0[2] > 1: # width too big
            p0[2] = 0.1 # more reasonable estimate

        popt, pcov = scipy.optimize.curve_fit(Do_Gaussian.gaus_lin_function, 
                                                x_data_seconds, y_data, p0=p0, maxfev=5000, 
                                                bounds=([50, 0, 0.015, -np.inf, -np.inf], 
                                                    [10000, np.inf, 1, np.inf, np.inf]))
        popt_np = -1*np.ones(len(popt), dtype=object)
        popt_np[0] = popt[0]
        popt_np[1] = current_date + pd.Timedelta(seconds=float(popt[1]))
        popt_np[2] = (2*np.sqrt(2*np.log(2)))*popt[2]
        
        if popt_np[2] > 0.2:
            # width TOO BIG. make that boy smaller
            popt_np[2] = 0.15
        
        if len(popt) == 5:
            # If superposed a Gaussian on a linear trend...
            popt_np[3:] = popt[3:]

        y_pred = Do_Gaussian.gaus_lin_function(x_data_seconds, *popt)
        
        try:
            r2, adj_r2 = self.goodness_of_fit(y_data, y_pred, len(popt))
        except ValueError as err:
            if 'Input contains NaN, infinity or a value too large' in str(err):
                print(f'popt={popt}')
                print(f'y-data={y_data}')
                print(f'y_pred={y_pred}')
            raise
        return popt_np, np.sqrt(np.diag(pcov)), r2, adj_r2

    @staticmethod # I think this line just avoids passing "self" to these functions below
    
    def gaus_lin_function(t, *args):
        # Author: Mike Shumko
        # Purpose: Fits a gaussian shape with an optional linear detrending term.

        """
        Args is an array of either 3 or 5 elements. First three elements are
        the Gaussian amplitude, center time, and width. The last two optional
        elements are the y-intercept and slope for a linear trend. 
        """
        
        exp_arg = -(t-args[1])**2/(2*args[2]**2)
        y = args[0]*np.exp(exp_arg)

        if len(args) == 5:
            y += args[3] + t*args[4]
            
        '''Add functionality to calculate area??'''
            
        return y


    def goodness_of_fit(self, y_true, y_pred, n_params):
        # Author: Mike Shumko
        # Purpose: Evaluate goodness of Gaussian fit (returns R^2, adjusted R^2)

        """
        Method to calculate the R^2 coefficient of determination
        and the adjusted R^2 coefficient given the number
        of fit parameters n_params.
        """
        r2 = sklearn.metrics.r2_score(y_true, y_pred)
        n = len(y_true)
        adj_r2 = 1 - (1-r2)*(n-1)/(n-1-n_params)
        return r2, adj_r2

def gaus_overlay_plot(props_df, hilt_df, ax, time_range, good_peaks):
    # Original author: Mike Shumko
    # Date created: unknown
    # Last modified: 10/11/2024 by Max Feinland
    # Purpose: Plots fitted Gaussian curves over raw HILT data.

    cmap = mpl.colormaps['plasma']
    plasma_colors = cmap(np.linspace(0, 1, 21))

    plot_start_time = props_df.t0.iloc[0] - pd.Timedelta(seconds=1)
    u = np.where(hilt_df.index.duplicated()==True)[0]
    df_idx = hilt_df.index[u]
    hilt_df = hilt_df.drop(index=df_idx) # get rid of duplicate indices

    closest_start_idx = hilt_df.index.get_indexer([plot_start_time], method='nearest')[0]
    plot_end_time = props_df.t0.iloc[-1] + pd.Timedelta(seconds=1)
    closest_end_idx = hilt_df.index.get_indexer([plot_end_time], method='nearest')[0]
    plot_time_range = [hilt_df.index[closest_start_idx], hilt_df.index[closest_end_idx]]
    
    time_array = hilt_df.loc[plot_time_range[0]:plot_time_range[-1]].index
    current_date = time_array[0].normalize()
    x_data_seconds = (time_array-current_date).total_seconds()
    y_data = hilt_df.loc[plot_time_range[0]:plot_time_range[1], 'counts']
    ax.plot(time_array, y_data, c=plasma_colors[0], label='HILT data')
    
    solid_counter = 0
    dashed_counter = 0
    
    for i in range(len(props_df)): # Loop through each peak
        popt = props_df.iloc[i, 3:]

        pdate = popt[1] # t0 in datetime terms
        popt[1] = (popt[1] - current_date).total_seconds()
        popt[2] = popt[2]/2.355 # Convert the Gaussian FWHM to std

        gaus_y = Do_Gaussian.gaus_lin_function(x_data_seconds, *popt)

        t_start = np.where(time_array==time_range[i][0])[0][0]
        t_end = np.where(time_array==time_range[i][1])[0][0]


        # If r^2 too low, plot dashed line
        if i in good_peaks:
            line_style = '-'
            solid_counter += 1
            if solid_counter > 1: # Check if one of these has been done already.
                labl=None
            else:
                labl=r'Gaussian fit ($R^2 > 0.9$)'
        else:
            line_style = '--'
            dashed_counter += 1
            if dashed_counter > 1:
                labl = None
            else:
                labl = r'Gaussian fit ($R^2 < 0.9$)'
        ax.plot(time_array[t_start:t_end], gaus_y[t_start:t_end], 
                           c=plasma_colors[10], linestyle=line_style, label=labl)
    
    ax.set_ylim(min(y_data)*0.95, max(y_data)*1.05)
    formatter = mpl.dates.DateFormatter('%H:%M:%S') # 1 decimal place
    ax.xaxis.set_major_formatter(formatter)
    plot_date_string = 'Time on ' + current_date.strftime('%B %-d, %Y') 
    ax.set_xlabel(plot_date_string, fontsize=18)
    ax.set_ylabel('Counts (#/20 ms)', fontsize=18)

    ax.legend(fontsize=12) 
    ax.set_xlim(time_array[0], time_array[-1])
    return 

def get_slotfill_times(param, threshold, bps, random_times):
    # param: dst or lpp
    # threshold: value that the parameter needs to dip below in order for it to be considered a slot-filling event
    # random_times: should you compare the slot-filling times to an equal number of random epochs?

    # Set some constants
    state4_start = datetime(1996, 8, 16)
    state4_end = datetime(2012, 11, 13)
    
    unique_slotfill_times = []
    slot_fill_idx = np.where(param <= threshold)[0] 
    slot_fill = param.index[slot_fill_idx]
    dt = slot_fill.diff()
    changes = np.where(dt > timedelta(hours=1))
    unique_slotfill_times = slot_fill[changes]

    if random_times:
        random_t = [random.random() * (state4_end - state4_start) + state4_start for _ in range(len(unique_slotfill_times))]
        final_slotx = []
        final_randomx = []
        final_L = []
        
        for i in range(len(unique_slotfill_times)):
            slotx = (bps.index - unique_slotfill_times[i])/timedelta(days=1)
            randomx = (bps.index - random_t[i])/timedelta(days=1)
            L = bps['L']
            
            final_slotx.extend(slotx)
            final_randomx.extend(randomx)
            final_L.extend(L.iloc)
    
        return final_slotx, final_randomx, final_L
    else:
        return unique_slotfill_times

# Begin plotting here.

# Text & plot formatting
warnings.filterwarnings("ignore")
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Paired.colors)

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["mathtext.default"] = "regular"

# I used a color theme for this paper, so define that here
cmap = mpl.colormaps['plasma']
plasma_colors = cmap(np.linspace(0, 1, 21))

# Set some constants
state4_start = datetime(1996, 8, 16)
state4_end = datetime(2012, 11, 13)

'''Read in, format, and manipulate some data up here, just for simplicity'''

# Read in microburst data
ibms = pd.read_csv('Data_Files/inner_belt_microbursts.csv', index_col=0) # ibms = inner belt microbursts
ibms.index = pd.to_datetime(ibms.index, format='mixed') # interpret index column as datetime rather than string

# Read in bouncing packet data
all_bps = pd.read_csv('Data_Files/bouncing_packets.csv', index_col=0)
all_bps.index = pd.to_datetime(all_bps.index, format='mixed') # interprets as str otherwise

# limit microbursts to the ones not identified as bouncing packets
bp_idx = []
for t in all_bps.index:
    # find matches between full microburst catalog and bouncing packet catalog
    bp_idx.extend(np.where(np.abs(t - ibms.index) < timedelta(seconds=0.2))[0])
bp_idx = np.unique(bp_idx)
other_microbursts = ibms.drop(ibms.index[bp_idx]) # remove these matches from microburst catalog

# Make dataframe of unique bouncing packets (so there are no duplicate points for the same event)
change_in_id = np.where(all_bps['event_id'].diff() > 0)[0].tolist()
change_in_id.insert(0, 0)
bps = all_bps.iloc[change_in_id]

# OMNI data
omni = pd.read_csv('Data_Files/omni.csv', index_col=0)
omni.index = pd.to_datetime(omni.index, format='%Y-%m-%d %H:%M:%S')

kp_max_24 = omni['Kp'].rolling('24h').max() # maximum Kp value over last 24 hours of recorded data
dst_min_24 = omni['Dst'].rolling('24h').min() # minimum Dst value over last 24 hours of recorded data

# calculate l_ppi from Carpenter & Anderson (1992)
lpp = 5.6 - kp_max_24.mul(0.46)

lpp_dst = pd.DataFrame({'lpp': lpp, 'Dst': dst_min_24}) # make dataframe with minimum plasmapause L and Dst index over last 24 hours


'''Example event time series'''

example_dates = [datetime(1998, 10, 3, 18, 30, 44), datetime(1998, 5, 11, 10, 23, 49), datetime(1999, 10, 27, 7, 50, 13), 
                 datetime(1999, 2, 18, 5, 58, 8), datetime(1999, 12, 25, 7, 48, 39), datetime(2003, 10, 30, 21, 19, 40)]

fig, ax = plt.subplots(3, 2, figsize=(12, 10), layout='tight')
letters = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

for i, ax in enumerate(ax.flat):
    date = example_dates[i]
    
    date_str = date.strftime('%B %d, %Y')
    print(f"Loading HILT data for {date_str}...")
    
    h = sampex.HILT(date) # count rate data
    h.load()

    # Find index of desired microburst
    if i < 5:
        t = ibms.index[abs(ibms.index - date) < timedelta(seconds=5)]
    else:
        t = ibms.index[abs(ibms.index - date) < timedelta(seconds=60)]

    start_time = t[0] - timedelta(seconds=1)
    end_time = t[-1] + timedelta(seconds=2)
        
    tidx = []
    
    try: # Try to find the indices of exact start & end times
        sidx = np.where(h['time'] == start_time)[0][0]
        eidx = np.where(h['time'] == end_time)[0][0]
        for tstmp in t:
            tidx.append(np.where(h['time'] == tstmp)[0][0])
    except: # In case there is no timestamp at that value
        sidx = np.abs((h['time'] - start_time).values).argmin()
        eidx = np.abs((h['time'] - end_time).values).argmin()
        tidx = [] # clear and start over
        for tstmp in t:
            tidx.append(np.abs((h['time'] - tstmp).values).argmin())
    
    # Plot 
    ax.plot(h['time'][sidx:eidx], h['counts'][sidx:eidx], color=plasma_colors[0])
    txt_xpos = [h['time'][sidx+10] if i < 5 else h['time'][sidx+100]]
    txt_ypos = min(h['counts'][sidx:eidx]) + 0.9*np.ptp(h['counts'][sidx:eidx])
    ax.text(txt_xpos, txt_ypos, letters[i], fontsize=18, fontweight='bold')
    if i < 5:
        ax.xaxis.set_major_locator(mpl.dates.SecondLocator(interval=1))
    else:
        ax.xaxis.set_major_locator(mpl.dates.SecondLocator(interval=5))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%-H:%M:%S')) 
    ax.set_xlim(h['time'][sidx], h['time'][eidx])
    ax.set_xlabel(f"Time on {date_str}", fontsize=14)
    ax.set_ylabel("Counts (#/20ms)", fontsize=14)
    ax.tick_params(labelsize=14)

# fig.savefig('Figures/example_time_series.pdf')

'''Gaussian fitting & energy estimation plot'''

timestamp = bps.index[223] # example event I'd like to use

h = sampex.HILT(timestamp) # load HILT data for this event
h.load()

obj = ID_Microbursts(h, timestamp) # call instance of ID_Microbursts class

hd, tr = obj.find_microbursts() # run find_microbursts on this instance

fit = obj.props # data about each Gaussian curve fit

good_pks = fit['r2'] >= 0.9 # find good peaks, with curves w/ R^2 > 0.9

# Find longest possible run of curves w/ R^2 > 0.9
changes = np.diff(good_pks.astype(int))
starts = np.where(changes == 1)[0] + 1
ends = np.where(changes == -1)[0]

if good_pks.iloc[0]:
    starts = np.insert(starts, 0, 0)
if good_pks.iloc[-1]:
    ends = np.append(ends, len(good_pks)-1)

if len(starts) > 0:
    lengths = ends - starts + 1
    longest_idx = np.argmax(lengths)
    good_pks_f = np.arange(starts[longest_idx], ends[longest_idx]+1)

# create figure & subfigures
fig = plt.figure(figsize=(15, 5))
subfigs = fig.subfigures(1, 2, wspace=-0.07)
gaus_ax = subfigs[0].subplots()
energy_ax = subfigs[1].subplots(2, 1, sharex=True, sharey=True)

# call gaussian fitting plot function
gaus_overlay_plot(fit, hd, gaus_ax, tr, good_pks_f)
gaus_ax.tick_params(axis='both', labelsize=16)
gaus_ax.text(datetime(2005, 10, 8, 22, 47, 37), 750, 'a)', fontsize=18) # add subplot label


'''Second/third plots: energy estimates'''

# The distribution in equatorial B-field strength (Beq) divided by local B-field strength (B or Bloc) depends on L, unsurprisingly.
# But you need this quantity to estimate the bounce period in a dipole model. So we fit a curve to the distribution
x = bps['L'] # L-shell
y = np.divide(bps['Beq'], bps['B']) # Ratio between Beq/Bloc

fit = np.polyfit(x, y, 2) # Fit quadratic curve

x0 = np.linspace(1.5, 2.5, 100) # evenly spaced L vector
y0 = fit[0]*x0**2 + fit[1]*x0 + fit[2] # fitted curve of Beq/Bloc

def calc_alpha_eq(b_ratio, alpha_loc): # Calculate equatorial pitch angle from local pitch angle (alpha_loc) and Beq/Bloc
    y = np.sin(np.deg2rad(alpha_loc))*np.sqrt(b_ratio)
    alpha_eq = np.rad2deg(np.arcsin(y))
    return alpha_eq

# Now there is the problem of local pitch angle.
# There are many possible values because SAMPEX has a wide FOV, +-34 degrees in each direction from instrument boresight.
# So to encompass all possible values, I used 90 degrees and the minimum loss cone angle (59 degrees)
possible_alpha_loc = [min(bps['losscone1']), 90]

# Next, the estimates: Schulz & Lanzerotti bounce period equation for dipole model
# (This is Mike Shumko's code)
Re = 6.371E6 # radius of earth
Tsl = lambda L, alpha0, v: 4*Re*np.divide(L, v) * \
       (1.3802 - 0.3198*(np.sin(np.deg2rad(alpha0)) + \
       np.sqrt(np.sin(np.deg2rad(alpha0)))))
beta = lambda Ek: np.sqrt(1-((Ek/511)+1)**(-2)) # fraction of speed of light
c = 3.0E8 # speed of light, m/s
KE_arr = [100, 500, 1000, 2000, 10000] # 100 keV to 10 MeV
L_range = np.linspace(1, 3, 100)

# Step through each possible energy and calculate the relationship to L given two different local pitch angles
for KE in KE_arr:
    tb_min = []
    tb_max = []
    for L in x0:
        b_ratio = fit[0]*L**2 + fit[1]*L + fit[2] # find Beq/Bloc for this L-shell
        alpha_eq_max = calc_alpha_eq(b_ratio, possible_alpha_loc[0]) # find equatorial pitch angle from Beq/Bloc and alpha_loc
        alpha_eq_min = calc_alpha_eq(b_ratio, possible_alpha_loc[1])
    
        tb_min.append(Tsl(L, alpha_eq_max, c*beta(KE)))
        tb_max.append(Tsl(L, alpha_eq_min, c*beta(KE)))

    if KE < 1000:
        label = str(int(KE)) + ' keV'
    else:
        label = str(int(KE/1000)) + ' MeV'

    # Plot the energy lines
    energy_ax[0].plot(x0, tb_min, '--', label=label)
    energy_ax[1].plot(x0, tb_max, '--', label=label)

# Read in spacing data between fitted Gaussians
w_s = pd.read_csv('Data_Files/gaussian_fits.csv', index_col=0)
w_s['s'] = [ast.literal_eval(x) for x in w_s['s']]

# Find minimum spacing for each event
s = []
for i in range(len(bps)):
    s_i = w_s['s'][i]
    if len(s_i) > 0:
        s.append(min(s_i))
    else:
        s.append(None)

letter_list = ['b)', 'c)'] # for subplot labeling

# Making the plots
for i in range(0,2):
    energy_ax[i].scatter(bps['L'], s, s=100, marker='*', color=plasma_colors[20], edgecolor=plasma_colors[5], 
                     zorder=10, label='bouncing packets') # plot the minimum spacing against L-shell
    energy_ax[i].set_xlim(1.5, 2.5) # set x-axis limits
    energy_ax[i].set_ylim(0.1, 0.45) # set y-axis limits
    energy_ax[i].set_yticks(np.arange(0.1, 0.4, 0.1)) # set y-axis ticks
    energy_ax[i].set_ylabel(f'Minimum $t_b$ (s)', fontsize=16) # set y-axis label
    energy_ax[i].tick_params(axis='both', labelsize=14) # make tick labels look nice
    
    alpha_loc_label = r"$\alpha_{loc} = " + str(int(round(possible_alpha_loc[i]))) + r"^{\circ}$"
    energy_ax[i].text(2.3, 0.15, alpha_loc_label, fontsize=16) # add text showing minimum alpha_loc
    energy_ax[i].text(1.53, 0.4, letter_list[i], fontsize=18) # add subplot label

energy_ax[0].legend(fontsize=12, ncol=3, loc='upper right') # add legend
energy_ax[1].set_xlabel('L-shell', fontsize=18) # add x-axis label
plt.subplots_adjust(hspace=0.1)

# fig.savefig('Figures/gaussian_fitting_energy_estimate.pdf')

'''MLT/L polar plot'''
    
# Create the plot
fig, ax = plt.subplots(1, 2, figsize=(12,7), subplot_kw={'projection': 'polar'})

def make_earth(ax): # make the earth
    ax.set_facecolor((0.85, 0.86, 0.9))

    # Nightside
    theta1 = np.linspace(0, 2 * np.pi, 1000) 
    r1 = np.ones_like(theta1)
    ax.fill_between(theta1, 0, r1, color='black')

    # Dayside
    theta2 = np.linspace(np.pi / 2, 3 * np.pi / 2, 500)
    r2 = np.ones_like(theta2)
    ax.fill_between(theta2, 0, r2, color='white') 
    
    # Formatting stuff
    ax.set_rticks(np.arange(0,4), labels=['0', '1', '2', '3'], size=20)  # L shell values
    ax.set_rlabel_position(165)
    ticklocs = np.pi*np.arange(0, 360, 90)/180
    ax.set_xticks(ticklocs.tolist())
    ax.set_xticklabels(['0', '6', '12', '18'], fontsize=20, )#,color="white")
    ax.set_theta_zero_location("E")
    ax.grid(True)
    ax.set_thetagrids(np.linspace(0,360,9))
    ax.text(np.pi, 2, 'L-shell', size=20, rotation=-15)
    _ = plt.text(-0.02, 3.2, 'MLT', color="black", size=20)
    return ax

# call make_earth function
ax0 = make_earth(ax[0])
ax1 = make_earth(ax[1])

# plotting bouncing packets
ax0.scatter(bps['MLT']*np.pi/12, bps['L'], s=100, color=plasma_colors[20], edgecolor=plasma_colors[5],
           marker='*', linewidth=1, label='bouncing packets', zorder=10)
ax0.set_rmax(2.5)

# plotting all microbursts
ax1.scatter(other_microbursts['MLT']*np.pi/12, other_microbursts['L'], s=20, c='aqua', label='other microbursts',
           marker='o', edgecolor=plasma_colors[0], linewidth=1)
ax1.set_rmax(2.5)


# add legend
_ = ax[0].legend(loc='upper right', bbox_to_anchor=(1.2, 1.1), fontsize=14)
_ = ax[1].legend(loc='upper right', bbox_to_anchor=(1.2, 1.1), fontsize=14)

# fig.savefig('Figures/mlt_l_polar_plot.pdf', bbox_inches='tight')

'''Events in space'''

fig, ax = plt.subplots(2, figsize=(15,12))

def make_world(ax):
    ax.set_facecolor('lightgray')
    world = gpd.read_file('Data_Files/ne_110m_admin_0_countries_lakes.shp') # you can get this file from my github or Natural Earth
    world.plot(ax=ax, color='gray', edgecolor='black', zorder=0)
    
    # Formatting shenanigans
    ax.set_xlabel("Longitude", fontsize=16)
    ax.set_ylabel("Latitude", fontsize=16)
    ax.set_ylim(-85, 85)
    ax.set_xlim(-200, 200)

    # Plot L-shell contours
    lgrid = pd.read_csv('Data_Files/Lgrid.dat', delimiter='\t', header=None)
    cols = [1, 3, 15, 17, 29]
    for i in cols:
        min_pos = np.argmin(lgrid.iloc[:,i])
        latl = np.concatenate(([lgrid.iloc[min_pos:,i-1], lgrid.iloc[:min_pos,i-1]]))
        lonl = np.concatenate(([lgrid.iloc[min_pos:,i], lgrid.iloc[:min_pos,i]]))
        ax.plot(lonl, latl, '--', color="black", zorder=0)
    
    # Add L-shell contour labels
    l_labels = {'lat': [48, -40, 59, -55],
               'name': ['L = 2', 'L = 2', 'L = 3', 'L = 3']}
    l_lon = -197
    for i in range(len(l_labels['lat'])):
        ax.text(l_lon, l_labels['lat'][i], l_labels['name'][i], fontsize=12, zorder=50)

    # Set tick labels for x- and y-axes
    ax.tick_params(axis='both', labelsize=12)
    _ = ax.set_xticks(np.arange(-180,181,45))
    _ = ax.set_yticks(np.arange(-90, 90, 30))

    ax.tick_params(labelsize=14)
    return ax

ax0 = make_world(ax[0])
ax1 = make_world(ax[1])

# Make custom plasma diverging colormap
plasma_d_list = np.vstack((mpl.cm.plasma_r(np.linspace(0, 1, 128)), mpl.cm.plasma(np.linspace(0, 1, 128))))
plasma_d = mpl.colors.LinearSegmentedColormap.from_list('plasma_d', plasma_d_list)

# Add bouncing packets
b = ax0.scatter(bps['lon'], bps['lat'], c=bps['pitch'], s=250, label='bouncing packets', marker='*', vmin=0, vmax=180, 
              cmap=plasma_d)

# Add other microbursts
m = ax1.scatter(other_microbursts['lon'], other_microbursts['lat'], c=other_microbursts['pitch'], s=100, label='other microbursts', 
                vmin=0, vmax=180, cmap=plasma_d)

# Add legends
ax0.legend(fontsize=12)
ax1.legend(fontsize=12)

# Add colorbar
cbar_ax = fig.add_axes([0.85, 0.11, 0.02, .77])
cbar = fig.colorbar(m, cax=cbar_ax)
cbar.set_label(r"Center of FOV ($^{\circ}$)", fontsize=16)
cbar.ax.tick_params(labelsize=14)

# fig.savefig('Figures/global_dist.pdf', bbox_inches='tight')

''' Distribution in time '''

# Read in data, get spinning times and missing times
spin = pd.read_csv('Data_Files/spin_times.csv')
spin[['start', 'end']] = spin[['start', 'end']].apply(pd.to_datetime)

missing = pd.read_csv('Data_Files/missing_times.csv')
missing[['start', 'end']] = missing[['start', 'end']].apply(pd.to_datetime)

# Making the plots
fig = plt.figure(layout='constrained', figsize=(12, 15))

# Separate into subfigures
subfigs = fig.subfigures(2, 1, height_ratios=[1, 3])
top_ax = subfigs[0].subplots() # First subfigure: entire time series

# resample monthly so we can pull out monthly-averaged sunspot number
sunspots = omni.resample('ME').mean()

# plot bouncing packets
top_ax.scatter(bps.index, bps['L'], label='bouncing packets', s=100, marker='*', color=plasma_colors[20], edgecolor=plasma_colors[5],
              zorder=20)
# plot other microbursts
top_ax.scatter(other_microbursts.index, other_microbursts['L'], label='other microbursts', color='aqua', edgecolor=plasma_colors[0], 
               s=20, marker='o', zorder=10)

# copy x-axis and plot sunspot number
ssax = top_ax.twinx()
ssax.plot(sunspots.index, sunspots['sunspot'], color=plasma_colors[18], label='sunspot number', linewidth=3)

# add rectangles to show when spacecraft was spinning
spin_color = '#7BDBE6'
for i in range(len(spin)-1):
    ssax.add_patch(mpl.patches.Rectangle((spin['start'][i], 0.5), spin['end'][i]-spin['start'][i], 500, facecolor=spin_color, alpha=0.5))
# add one patch for the legend
ssax.add_patch(mpl.patches.Rectangle((spin['start'].iloc[-1], 0.5), spin['end'].iloc[-1]-spin['start'].iloc[-1], 500, facecolor=spin_color, 
                                     alpha=0.5, label='spin times'))

# add rectangles to show when data was missing or was flagged
missing_color = 'red'
for i in range(len(missing)-1):
    ssax.add_patch(mpl.patches.Rectangle((missing['start'][i], 0.5), missing['end'][i]-missing['start'][i], 500, facecolor=missing_color,
                                         alpha=0.5))
ssax.add_patch(mpl.patches.Rectangle((missing['start'].iloc[-1], 0.5), missing['end'].iloc[-1]-missing['start'].iloc[-1], 500, 
                                     facecolor=missing_color, alpha=0.5, label='flagged/missing data'))

# Format (set axis limits, labels, etc.)
top_ax.set_xlim(datetime(1996,8,15), datetime(2012, 11, 7))
top_ax.set_ylim(1,2.5)
top_ax.tick_params(axis='both', labelsize=16)
top_ax.set_ylabel('L-Shell', fontsize=18)
top_ax.set_xlabel('Year', fontsize=18)
top_ax.set_yticks([1, 1.5, 2, 2.5])

ssax.set_ylim(0, 500)
ssax.tick_params(axis='both', labelsize=16)
ssax.set_ylabel('Sunspot number', fontsize=16)

# Add a legend
lines, labels = top_ax.get_legend_handles_labels()
lines2, labels2 = ssax.get_legend_handles_labels()
_ = ssax.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=12)
ssax.text(datetime(1997, 1, 1), 450, 'a)', fontsize=18) # and a subplot label


''' Distribution in time, zoomed in '''

date_lims = [datetime(1998, 4, 1), datetime(1999, 8, 1), datetime(2000, 11, 1), datetime(2002, 3, 1),
             datetime(2003, 7, 1), datetime(2004, 10, 1), datetime(2005,12,1)]

n_subplots = len(date_lims)-1

alphabet = list(map(chr, range(97, 123)))
alphabet = alphabet[1:n_subplots+1]
letters = [letter + ')' for letter in alphabet]

ax = subfigs[1].subplots(n_subplots)

t = pd.to_datetime(bps.index, format='mixed')

for i in range(n_subplots):
    twinx = ax[i].twinx() # same x, different y (Dst index)
    twiny = ax[i].twiny() # same y, different x (bps)
    
    ax[i].plot(lpp_dst.index, lpp_dst['lpp'], '--', color='midnightblue', label='$L_{ppi}$', linewidth=0.75)
    ax[i].fill_between(lpp_dst.index, lpp_dst['lpp'], 5, color='thistle', alpha=0.5)

    dstl = twinx.plot(lpp_dst.index, lpp_dst['Dst'], color='royalblue', label='Dst index')
    
    bpp = twiny.scatter(bps.index, bps['L'], s =100, c=bps['counts'], label='bouncing packets', cmap='plasma',
                                  marker='*', zorder=50, norm=mpl.colors.LogNorm(vmin=100, vmax=1e+4))

    for j in range(len(spin)-1):
        twinx.add_patch(mpl.patches.Rectangle((spin['start'][j], -500), spin['end'][j]-spin['start'][j], 1000, facecolor=spin_color, 
                                              alpha=0.5))
    # add one patch for the legend
    twinx.add_patch(mpl.patches.Rectangle((spin['start'].iloc[-1], -500), spin['end'].iloc[-1]-spin['start'].iloc[-1], 1000, 
                                          facecolor=spin_color, alpha=0.5, label='spin times'))

    for j in range(len(missing)-1):
        twinx.add_patch(mpl.patches.Rectangle((missing['start'][j], -500), missing['end'][j]-missing['start'][j], 1000, facecolor="red", 
                                              alpha=0.5))
    # add one patch for the legend
    twinx.add_patch(mpl.patches.Rectangle((missing['start'].iloc[-1], -500), missing['end'].iloc[-1]-missing['end'].iloc[-1], 1000, 
                                          facecolor="red", alpha=0.5, label='flagged/missing data'))

    # formatting
    ax[i].set_ylabel('L-shell', fontsize=14)
    ax[i].set_xlim(date_lims[i], date_lims[i+1])
    ax[i].set_ylim(1, 3.5)
    ax[i].set_yticks([1, 2, 3])
    minor_ticks = pd.date_range(date_lims[i], date_lims[i+1], freq='MS')
    major_ticks = pd.date_range(date_lims[i], date_lims[i+1], freq='3MS')
    ax[i].set_xticks(major_ticks)
    ax[i].set_xticks(minor_ticks, minor=True)
    ax[i].xaxis.set_major_formatter(mpl.dates.DateFormatter('%B %Y'))

    twinx.tick_params(axis='x', which='both', bottom=False, top=False)
    twinx.set_ylabel('Dst (nT)', fontsize=14)
    twinx.set_ylim(-400, 100)
    twinx.set_yticks([-400, -200, 0])
    
    twiny.set_xlim(date_lims[i], date_lims[i+1])
    twiny.get_xaxis().set_visible(False)

    # legend formatting
    bp_line, bp_label = twiny.get_legend_handles_labels()
    lpp_line, lpp_label = ax[i].get_legend_handles_labels()
    dst_line, dst_label = twinx.get_legend_handles_labels()
    if i == 0:
        twiny.legend(bp_line + lpp_line + dst_line, bp_label + lpp_label + dst_label, loc='upper right', ncol=2)

    twinx.text(date_lims[i]+timedelta(days=7), 3, letters[i], fontsize=18)

# colorbar
cbar_ax = fig.add_axes([1.01, 0.02, 0.02, 0.72])
cbar = fig.colorbar(bpp, cax=cbar_ax)
cbar.set_label('Maximum count rate', fontsize=14)
cbar.ax.tick_params(labelsize=14)

# fig.savefig('Figures/L_time_series.pdf', bbox_inches='tight')

'''Superposed epoch'''

random.seed(15) # set seed for reproducibility

# Plotting parameters
n_rows = 6
n_cols = 2 # it's actually 4, but don't worry about it
nbins = 20 # number of bins
xmin = -100 # days before slot-filling event
xmax = 100 # days after slot-filling event
ymin = 1.5 # min L-shell
ymax = 2.5 # max L-shell
cmax = 20 # max number of events in colorbar bin
twinxymax = 250 # max on line histogram y-axis

alphabet = list(map(chr, range(97, 123)))
labeled_alphabet = [letter + ')' for letter in alphabet] # subplot labels

# Make the figure
fig = plt.figure(figsize=(14,15))
subfigs = fig.subfigures(1, 2, wspace=0.2)
lpp_subplots = subfigs[0].subplots(n_rows, n_cols, sharex=True, sharey=True)
dst_subplots = subfigs[1].subplots(n_rows, n_cols, sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.1, hspace=0.1)

# plotting dictionary (makes the function very easy)
pdict = {
    'lpp': {
        'thresholds': np.linspace(2, 3, n_rows), # list of thresholds used for this parameter
        'cmap': mpl.cm.plasma, # colormap
        'label_str': r'$L_{ppi}$', # formatted label on plot
        'label_offsets': [0, 6], # for subplot labeling
        'cbar_coords': [.46, 0.11, 0.03, .77] # colorbar coordinates
    },
    'Dst': {
        'thresholds': [int(x) for x in [-200, -150, -130, -100, -80, -60]], # list of thresholds used for this parameter
        'cmap': mpl.cm.viridis, # colormap
        'label_str': r'Dst', # formatted label on plot
        'label_offsets': [12, 18], # for subplot labeling
        'cbar_coords': [1, 0.11, 0.03, .77] # colorbar coordinates
    }
}

def make_hexbin_plots(lpp_dst, param, ax, bps):
    for j in range(n_rows):
        # call slotfill_times function and get out superposed epoch data
        slot_times, random_times, L = get_slotfill_times(lpp_dst[param], pdict[param]['thresholds'][j], bps, random_times=True)

        # making hexagon plots for lpp trigger
        handle = ax[j,0].hexbin(slot_times, L, cmap=pdict[param]['cmap'], gridsize=nbins, 
                                extent=(xmin,xmax,ymin,ymax), vmin=0, vmax=cmax)
        ax[j,1].hexbin(random_times, L, cmap=pdict[param]['cmap'], gridsize=nbins, 
                                extent=(xmin,xmax,ymin,ymax), vmin=0, vmax=cmax)
    
        lines = [slot_times, random_times]
    
        for i in range(n_cols): # step through each column
            twinx = ax[j,i].twinx() # copy the current axes
            twinx.hist(lines[i], bins=np.linspace(xmin,xmax,nbins+1), color='white', histtype='step', linewidth=1.5)
            ax[j,i].vlines(x=0, ymin=1.5, ymax=2.5, color='white', linewidth=1, linestyle='--') # add days=0 line
    
            # formatting
            twinx.set_ylim(0, twinxymax)
            ax[j,i].set_xlim(xmin,xmax)
            ax[j,i].set_ylim(ymin, ymax)
            ax[j,i].set_xticks([-100, -50, 0, 50, 100])
            ax[j,i].set_yticks([1.5, 2, 2.5])
            ax[j,i].tick_params(axis='both', which='major', labelsize=12) 

            ax[j,i].text(-95, 2.4, labeled_alphabet[j+pdict[param]['label_offsets'][i]], fontsize=14, color='white', fontweight='bold')
    
            real_ticks = np.linspace(0, twinxymax, 6)
        
            if i == 0:
                twinx.tick_params(axis='y', which='both', labelright=False)
                ax[j,i].set_xticks([-100, -50, 0, 50, 100])
            else:
                twinx.set_yticks(real_ticks)
                twinx.set_yticklabels([str(int(tick/10)) for tick in real_ticks])
                twinx.set_ylabel('Mean events per day', fontsize=12)
                twinx.tick_params(axis='both', which='major', labelsize=12)
                ax[j,i].set_xticks([-50, 0, 50, 100])
            
            if j == 0:
                ax[j,1].text(10, 2.4, 'random', fontsize=14, color='white')
            elif j == n_rows-1:
                ax[j,i].set_xlabel('Days since slot-filling event', fontsize=12)

        
        param_str = pdict[param]['label_str'] + '$\leq$' + str(pdict[param]['thresholds'][j])

        ax[j,0].text(10, 2.4, param_str, fontsize=14, color='white')
        ax[j,0].set_ylabel('L-shell', fontsize=16)
        
    cbar_ax = fig.add_axes(pdict[param]['cbar_coords'])
    cbar = fig.colorbar(handle, cax=cbar_ax)
    cbar.set_label('Events per day per L', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_ticks(np.arange(0, cmax+1, 5))

# Call the function down here
make_hexbin_plots(lpp_dst, 'lpp', lpp_subplots, bps)
make_hexbin_plots(lpp_dst, 'Dst', dst_subplots, bps)

# fig.savefig('Figures/superposed_epoch.pdf', bbox_inches='tight')

'''Lifetimes'''

fig, ax = plt.subplots(2, 2, figsize=(15, 10))

# thresholds = [-150, -130, -100, -80] # for Dst
thresholds = [2.2, 2.4, 2.6, 2.8] # for L_ppi
letters = ['a)', 'b)', 'c)', 'd)']

for i, ax in enumerate(ax.flat):
    slotfill_times = get_slotfill_times(lpp_dst['lpp'], thresholds[i], bps, random_times=False)
    
    days_since_last_slot = []
    
    for timestamp in bps.index:
        try:
            closest_slot_idx = np.where(np.array(slotfill_times) < timestamp)[0][-1]
            days_since_last_slot.append((timestamp - slotfill_times[closest_slot_idx])/timedelta(days=1))
        except: # sometimes there may not be a SFE before the timestamp, so you can just skip it
            pass

    slot_to_slot = slotfill_times.diff()/timedelta(days=1)
    
    thebins = np.arange(0, 200, 5)
    thebins = np.concatenate((thebins, [2000])) # add overflow bin

    _,_,bphist = ax.hist(days_since_last_slot, bins=thebins, color=plasma_colors[12], edgecolor=plasma_colors[5],
                         label=r'$\Delta t_{slot-bp}$')
    _,_,slothist = ax.hist(slot_to_slot, bins=thebins, color=plasma_colors[0], label=r'$\Delta t_{slot-slot}$')
    if i >= 2:
        ax.set_xlabel(r'$\Delta t$ (days)', fontsize=20)
    if i % 2 == 0:
        ax.set_ylabel('Number of events', fontsize=20)
    ax.tick_params(axis='both', labelsize=16)
    ax.set_ylim(0, 70)


    binscenters = np.array([0.5 * (thebins[i] + thebins[i+1]) for i in range(len(thebins)-1)])
    x0 = np.linspace(0, 200, 1000)
    
    def exponential(x, a, k, b):
        return a*np.exp(-x*k) + b
    
    counts, bins = np.histogram(days_since_last_slot, bins=thebins, density=False)
    fit_params, _ = curve_fit(exponential,  xdata=binscenters, ydata=counts)
    ax.plot(x0, exponential(x0, *fit_params), linestyle='--', color=plasma_colors[3], linewidth=2, 
            label=r'exponential fit ($\Delta t_{slot-bp}$)')
    
    if i == 1:
        ax.legend(fontsize=12)
    
    # adding text with characteristic e-folding time
    bptxt = r'$\tau_{slot-bp} = $' + str(np.round(1/fit_params[1], 2)) + ' days'
    print(thresholds[i], np.mean(days_since_last_slot), np.mean(slot_to_slot)) # for the table

    plotstr = r'$L_{ppi} \leq $' + str(thresholds[i])
    ax.text(20, 62, plotstr, fontsize=14, fontweight='bold')
    ax.text(100, 45, bptxt, fontsize=16)
    ax.text(7, 62, letters[i], fontsize=18)
    ax.set_xticks(np.arange(0, 250, 50))
    ax.set_xticklabels(['0', '50', '100', '150', r'$\geq 200$'])

    ax.set_xlim(0, 200)

# fig.savefig('Figures/dt_dist_mean.pdf', bbox_inches='tight')

'''Lightning data analysis'''
        
# Line plots of mean seasonal lightning for each geographic area
filename = "Data_Files/LISOTD_LRMTS_V2.3.2015.nc" # This is monthly-averaged lightning data from the Lightning Imaging Sensor (LIS)
ds = xr.open_dataset(filename, decode_times=False) # open the dataset
lightning = ds.LRMTS_COM_FR # pull out lightning data
lat = ds.Latitude # pull out lat
lon = ds.Longitude # pull out lon

# I made these up... based on data... but they're still made up
bounding_boxes = {'north america': [-105, -50, 23, 45],
                'south america': [-125, -15, -65, -40],
                'africa': [0, 55, -50, -15],
                'europe': [-13, 42, 36, 60]}

def plot_lightning(region, col, ax):
    lon_min, lon_max, lat_min, lat_max = bounding_boxes[region] 

    # This code from: https://github.com/dervlamk/coding_refs/blob/main/notebooks/working_with_nc_practice/nc_practice_notebook.ipynb
    # (Thank you Kevin & Dervla!)
    l = lightning.sel(Latitude=slice(lat_min,lat_max),Longitude=slice(lon_min,lon_max), Month_since_Jan_95=slice(36, 132))
    
    # determine weights based on latitude value
    weights = np.cos(np.deg2rad(lat))
    weights.name = 'weights'
    
    # calculate area-weighted values
    l_weighted = l.weighted(weights)
    
    # calculate mean of weighted data
    l_weighted_mean = l_weighted.mean(dim=['Latitude','Longitude'], keep_attrs=True)
    
    # Group by month
    month_nums = np.arange(0, 12)
    
    seasons = pd.DataFrame(columns=list_of_months)
    
    for i in range(len(l_weighted_mean)):
        month_num = i % 12
        yr = np.floor(i/12) + 1995
        monthname = list_of_months[month_num]
        seasons.loc[yr, monthname] = l_weighted_mean.values[i]
    
    monthly_mean = seasons.mean()

    labl = 'lightning (' + region + ')'

    r = ax.plot(np.arange(1.5, 13), monthly_mean, color=col, linewidth=2, label=labl)
    return r 

# Plotting
fig, axes = plt.subplots(2, 1, figsize=(12,8))

lpp_thresholds = [2.4, 2.6]
bins = np.arange(1, 14)

americas = bps[bps['lon'] < -10]
africa_europe = bps[bps['lon'] >= -10]
amlines = [americas.index.month]
amlabels = ['bouncing packets (americas)']
afrlines = [africa_europe.index.month]
afrlabels = ['bouncing packets (africa & europe)']

# colorbar
cmap = mpl.colormaps['plasma']
plasma_colors = cmap(np.linspace(0, 1, 21))

paired_colors = plt.cm.Paired.colors
all_colors = [tuple([float(x) for x in plasma_colors[5]])]
all_colors.extend(paired_colors[0:2])

for thres in lpp_thresholds:
    times = get_slotfill_times(lpp_dst['lpp'], thres, bps, random_times=False)
    amlines.append([time.month for time in times])
    afrlines.append([time.month for time in times])
    amlabels.append(r'slot-filling events ($L_{ppi} \leq$' + str(thres) + ')')
    afrlabels.append(r'slot-filling events ($L_{ppi} \leq$' + str(thres) + ')')
    
axes[0].hist(amlines, color=all_colors, bins=bins, label=amlabels, linewidth=2)
axes[1].hist(afrlines, color=all_colors, bins=bins, label=afrlabels, linewidth=2)

list_of_months = list(pd.date_range(start='2024-01-01', periods=12, freq='ME').strftime('%b'))

for i, ax in enumerate(axes.flat):
    major_ticks = np.arange(1.5, 13)
    ax.set_xticks(major_ticks)
    ax.set_xticklabels(list_of_months, fontsize=12)
    ax.set_ylabel('Number of events', fontsize=16)
    ax.set_xlabel('Month', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=12) 

    lightning_ax = ax.twinx()
    if i == 0:
        plot_lightning('north america', 'gold', lightning_ax)
        plot_lightning('south america', 'yellow', lightning_ax)
    else:
        plot_lightning('europe', 'gold', lightning_ax)
        plot_lightning('africa', 'yellow', lightning_ax)
        
    lightning_ax.set_ylabel(f'Mean lightning density (fl/km$^2$/day)', fontsize=12)
    lightning_ax.tick_params(axis='both', which='major', labelsize=12) 
    lightning_ax.set_ylim(0, 0.055)

    slotlines, slotlabels = ax.get_legend_handles_labels()
    llines, llabels = lightning_ax.get_legend_handles_labels()
    _ = ax.legend(slotlines + llines, slotlabels + llabels, loc='upper left', fontsize=12)
    ax.set_ylim(0, 25)
    ax.set_xlim(1, 13)

# fig.savefig('Figures/seasonal_hist.pdf', bbox_inches='tight')
