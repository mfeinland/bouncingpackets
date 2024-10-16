'''Find microburst peaks'''
'''MOST OF THIS CODE WAS ADAPTED FROM MIKE SHUMKO'S GITHUB REPOSITORY'''
'''(https://github.com/mshumko/sampex_microburst_widths/tree/main/sampex_microburst_widths/microburst_id)'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import scipy.optimize
import scipy.signal
import sklearn.metrics
import sampex
from scipy.integrate import quad
import warnings

class ID_Microbursts:
    # tell the class the properties you wanna store in here preliminarily.
    def __init__(self, h, packet_timestamp, prominence_rel_height=0.5):
        self.h = h # assign HILT data object to self.h
        self.packet_timestamp = packet_timestamp # assign packet start timestamp
        self.prominence_rel_height=prominence_rel_height # assign prominence height width evaluator
        #(when you find a peak, you measure the width at this % the full prominence.
        # Default is 0.5, 0 would be the base, 1 would be the top.)
        return

    def obrien(self): # O'Brien et al. 2003 criterion, modified to catch lower-count-rate instances
        N20 = self.h['counts'] # count rate sampled every 20 ms
        time20 = self.h['time'] # time every 20 ms

        df = pd.DataFrame({'time': self.h['time'], 'counts': self.h['counts']})

        df.set_index('time', inplace=True) # set time column as the index

        # resample the dataframe to 100 ms intervals and sum the counts in each interval
        N100 = df.resample('100ms').sum()

        A500 = N100.rolling(5, center=True).mean() # 5-observation centered rolling mean (over 500 ms)

        # O'Brien et al 2003
        condition = np.divide((N100.counts - A500.counts), np.sqrt(1 + A500.counts)) 

        ns = np.where(condition > 5)[0] # RHS changed from 10 to 5
        
        epsilon = 10; # if two flagged indices are spaced less than this distance apart, 
        # they are probably part of the same microburst

        starts = [] # initializing
        ends = []

        dn = np.diff(ns) # difference in time between instances of the condition being true

        # finding extended periods of the condition being true
        for i in np.arange(1,len(dn)-10):
            if dn[i] < epsilon and dn[i+1] < epsilon and dn[i-1] >= epsilon: # start condition
                starts.append(ns[i])
                for j in np.arange(i+1, len(dn)-1):
                    if dn[j] < epsilon and dn[j+1] >= epsilon:
                        ends.append(ns[j]) # end condition
                        break
            elif dn[i] <= epsilon and i == 1: # start condition (edge case)
                starts.append(ns[i])
                for j in np.arange(i+1, len(dn-1)):
                    if dn[j] <= epsilon:
                        ends.append(ns[j])
                        break
            elif i == len(dn): # end condition (edge case)
                ends.append(ns[i])

        if len(starts) > len(ends):
            ends.append(starts[len(starts)-1] + 10)

        starts = [x - 2 for x in starts] # pad start times with 0.2 seconds 
        ends = [x + 10 for x in ends] # pad end times with 1 second


        def changeCadence(time20, times100):
            # times100 is the list of timestamps for whatever it is: starts, ends, ns
            # times20 is the entire list of time20 for the day in question
            try:
                # try list comprehension first, it's faster
                idx20 = [time20.get_loc(tstmp) for tstmp in times100] 
            except: 
                # if a timestamp is missing from time20, find the closest timestamp 
                idx20 = [] # clear and start over
                idx20 = [np.abs((time20 - target).values).argmin() for target in times100] 
            return idx20

        # reverting indices to 20 ms cadence
        self.so = changeCadence(time20, N100.index[starts])
        self.eo = changeCadence(time20, N100.index[ends])
        self.no = changeCadence(time20, N100.index[ns])
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
        mbidx =  np.abs((t[self.so] - self.packet_timestamp).values).argmin()
        
        interval = rate[self.so[mbidx]:self.eo[mbidx]] # define just one chunk of the rate data, 
        # taken from O'Brien interval
        
        maxdis = max(interval) - min(interval) # use to generate prominence requirement

        # finding peaks with a high enough prominence
        [pks, _] = find_peaks(interval, prominence=0.25*maxdis, distance=3)

        # find instances where a bunch of peaks are evenly spaced (at least 4 peaks)
        indices = np.where(np.convolve(np.abs(np.diff(np.diff(pks))) <= 3, \
                                               np.ones(2), mode='valid') == 2)[0]
        indices = np.append(indices, [indices[-1]+1, indices[-1]+2, indices[-1]+3])
        
        pks = pks[indices] # pks is a list of indices in the HILT data that has a peak
        # that meets all these conditions. (o'brien, high enough prominence, evenly spaced)

        # Add peak_idx to self-- a list of indices in the HILT dataframe containing peaks
        # (add self.so[mbidx] b/c that's the beginning of the interval)
        self.peak_idx = list(self.so[mbidx] + loc for loc in pks)
        
        # Convert HILT objects to dataframe
        hilt_df = pd.DataFrame(data={'counts':h.counts}, index=self.h['time'])

        # call an instance of Do_Gaussian class using the day's HILT data, 
        # the list of peak indices, and the prominence width 
        # (how high up the peak, 0 to 1 where 0 is the base and 1 is the top, to evaluate the width)
        gaus = Do_Gaussian(hilt_df, self.peak_idx, self.prominence_rel_height)
        
        # Calculate widths using scipy.signal.peak_widths
        gaus.calc_prominence_widths()
        
        # Using this initial estimate, calculate widths using Gaussian fit
        fit_df, tr = gaus.calc_gaus_widths()

        self.props = fit_df # assign attribute 'props' to self
        # (contains information about the fit of the peaks: r^2, amplitude, center, width, etc.)
        # return tr (time range)
        '''Rename tr????'''
        return hilt_df, tr


class Do_Gaussian:
    def __init__(self, hilt_df, peak_idx, prom_h, width_multiplier=2.5, plot_width_s=5):
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
        widths_tuple = scipy.signal.peak_widths(self.hilt_df['counts'], self.peak_idx, 
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

        # Create empty pd.DataFrames for fit variables.
        fit_param_names = ['r2', 'adj_r2', 'A', 't0', 'fwhm']
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
            df.iloc[i, 2:] = popt 

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
       
