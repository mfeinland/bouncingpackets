'''Gaussian curve fitting '''
# Most of this code was adapted from Mike Shumko's GitHub repository:
# https://github.com/mshumko/sampex_microburst_widths/tree/main/sampex_microburst_widths/microburst_id

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colormaps
import scipy.optimize
from scipy.signal import find_peaks, peak_widths
import sklearn.metrics
from datetime import datetime, timedelta
import warnings
import os
import random

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

    cmap = colormaps['plasma']
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
    formatter = mdates.DateFormatter('%H:%M:%S') # 1 decimal place
    ax.xaxis.set_major_formatter(formatter)
    plot_date_string = 'Time on ' + current_date.strftime('%B %-d, %Y') 
    ax.set_xlabel(plot_date_string, fontsize=18)
    ax.set_ylabel('Counts (#/20 ms)', fontsize=18)

    ax.legend(fontsize=12) 
    ax.set_xlim(time_array[0], time_array[-1])
    return 