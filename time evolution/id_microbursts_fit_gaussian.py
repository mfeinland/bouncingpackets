'''Find microburst peaks'''
'''MOST OF THIS CODE WAS ADAPTED FROM MIKE SHUMKO'S GITHUB REPOSITORY (https://github.com/mshumko/sampex_microburst_widths/tree/main/sampex_microburst_widths/microburst_id)'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import scipy.optimize
import scipy.signal
import sklearn.metrics
import sampex

class ID_Microbursts:
    # tell the class the properties you wanna store in here preliminarily.
    def __init__(self, h, packet_timestamp, prominence_rel_height=0.5):
        self.h = h
        self.packet_timestamp = packet_timestamp
        self.prominence_rel_height=prominence_rel_height
        return

    def obrien(self):
        N20 = self.h['counts'] # count rate sampled every 20 ms
        time20 = self.h['time'] # time every 20 ms

        df = pd.DataFrame({'time': self.h['time'], 'counts': self.h['counts']})

        df.set_index('time', inplace=True) # set time column as the index

        # resample the dataframe to 100 ms intervals and sum the counts in each interval
        N100 = df.resample('100ms').sum()

        A500 = N100.rolling(5, center=True).mean() # 5-observation centered rolling mean (over 500 ms)

        # O'Brien et al 2003
        condition = np.divide((N100.counts - A500.counts), np.sqrt(1 + A500.counts)) 

        ns = np.argwhere(condition > 5)
        ns = [item[0] for item in ns]

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

        starts = [x - 2 for x in starts] # pad with 0.2 seconds 
        ends = [x + 10 for x in ends] # pad with 1 second


        def changeCadence(time20, times100):
            # times100 is the list of timestamps for whatever it is: starts, ends, ns
            # times20 is the entire list of time20 for the day in question
            try:
                # try list comprehension first, it's faster
                idx20 = [time20.get_loc(tstmp) for tstmp in times100] 
            except: 
                idx20 = [] # clear and start over
                idx20 = [np.abs((time20 - target).values).argmin() for target in times100] 
                # if a timestamp is missing from time20,
                # you will have to find the minimum timestamp (that is closest in time)
            return idx20

        # reverting indices to 20 ms cadence
        self.so = changeCadence(time20, N100.index[starts])
        self.eo = changeCadence(time20, N100.index[ends])
        self.no = changeCadence(time20, N100.index[ns])
        return
    
    
    def find_microbursts(self):
        self.obrien()
        
        t = self.h['times']
        rate = self.h['counts']
        
        # Returns the index corresponding to the start of the previously ID'd packet.
        mbidx =  np.abs((t[self.so] - self.packet_timestamp).values).argmin()
        
        interval = rate[self.so[mbidx]:self.eo[mbidx]] # define just one chunk of the rate data, taken from O'Brien
        maxdis = max(interval) - min(interval) # use to generate prominence requirement

        # finding peaks with a high enough prominence
        [pks, _] = find_peaks(interval, prominence=0.25*maxdis, distance=3)
        self.peak_idx = list(self.so[mbidx] + loc for loc in pks)
        # Calculate the microburst widths using the prominence method and
        # the Gaussian fit.
        hilt_df = pd.DataFrame(data={'counts':h.counts}, index=self.h['time'])
        gaus = do_gaussian(hilt_df, self.peak_idx, self.prominence_rel_height)
        gaus.calc_prominence_widths()
        fit_df = gaus.calc_gaus_widths(debug=True)

        # Save to a DataFrame
        df = pd.DataFrame(
            data={
                'dateTime':self.h['time'][self.peak_idx],
                'width_s':gaus.prom_widths_s,
                'width_height':gaus.width_height,
                'peak_counts_20ms':self.h['counts'][self.peak_idx],
                'left_peak_base':gaus.left_peak_base,
                'right_peak_base':gaus.right_peak_base,
                },
            index=self.peak_idx
            )
        self.props = fit_df
        return


class do_gaussian:
    def __init__(self, hilt_df, peak_idx, prom_h, width_multiplier=2.5, plot_width_s=5):
        """
        
        """
        self.hilt_df = hilt_df
        self.hilt_times = self.hilt_df.index.to_numpy()
        self.peak_idx = peak_idx
        self.prominence_rel_height = prom_h
        self.width_multiplier = width_multiplier
        self.plot_width_s = plot_width_s
        return

    def calc_prominence_widths(self):
        """
        Use scipy to find the peak width at self.prominence_rel_height prominence
        """
        # Check that self.stb.peak_idx correspond to the max values
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
        return self.prom_widths_s, self.width_height, self.left_peak_base, self.right_peak_base

    def calc_gaus_widths(self, debug=True, detrend=True):
        
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
        for i, (peak_i, width_i, height_i) in enumerate(zip(self.peak_idx, self.prom_widths_s, self.width_height)):
            peak_ok_check = 0 # see if r^2 is close enough to 1
            incrementer = 0 # number of loops you had to do

            while peak_ok_check == 0:
                time_range = [
                    self.hilt_times[peak_i]-pd.Timedelta(seconds=width_i)*self.width_multiplier,
                    self.hilt_times[peak_i]+pd.Timedelta(seconds=width_i)*self.width_multiplier
                            ]
                # If too little data points, assume a 500 ms fit width.
                if len(self.hilt_df.loc[time_range[0]:time_range[1], :].index) < 5:
                    time_range = [
                                self.hilt_times[peak_i]-pd.Timedelta(seconds=0.25),
                                self.hilt_times[peak_i]+pd.Timedelta(seconds=0.25)
                            ]
                t0 = self.hilt_times[peak_i]

    #             if width_i < 0.1:
    #                 # If the prominence method width is small 
    #                 # change it to a 0.1 s width as a starting guess.
    #                 width_i = 0.1 

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
                        if r2 >= 0.85:
                            peak_ok_check = 1 # fit is acceptable
                            print("This fit is good (took", incrementer, "iterations)")
                        elif incrementer < 5:
                            width_i = width_i*0.9 # make the peak skinnier
                            incrementer = incrementer + 1
                        else:
                            peak_ok_check = 1 # just give up LMAO
                            
                    except RuntimeError as err:
                        if ('Optimal parameters not found: Number of calls '
                            'to function has reached maxfev') in str(err):
                            continue
                        raise
                    if len(w):
                        print(w[0].message, '\n', p0, popt)

            # Save to a pd.DataFrame row.
            df.iloc[i, :2] = r2, adj_r2
            df.iloc[i, 2:] = popt 
            if debug:
                self.fit_test_plot(t0, time_range, popt, r2, adj_r2)
        return df

    def fit_gaus(self, time_range, p0):
        """
        Fits a gausian shape with an optinal linear detrending term.
        """
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

        popt, pcov = scipy.optimize.curve_fit(do_gaussian.gaus_lin_function, 
                                                x_data_seconds, y_data, p0=p0, maxfev=5000)
        popt_np = -1*np.ones(len(popt), dtype=object)
        popt_np[0] = popt[0]
        popt_np[1] = current_date + pd.Timedelta(seconds=float(popt[1]))
        popt_np[2] = (2*np.sqrt(2*np.log(2)))*popt[2]
        if len(popt) == 5:
            # If superposed a Gaussian on a linear trend...
            popt_np[3:] = popt[3:]

        y_pred = SAMPEX_Microburst_Widths.gaus_lin_function(x_data_seconds, *popt)
        try:
            r2, adj_r2 = self.goodness_of_fit(y_data, y_pred, len(popt))
        except ValueError as err:
            if 'Input contains NaN, infinity or a value too large' in str(err):
                print(f'popt={popt}')
                print(f'y-data={y_data}')
                print(f'y_pred={y_pred}')
            raise
        return popt_np, np.sqrt(np.diag(pcov)), r2, adj_r2

    @staticmethod
    def gaus_lin_function(t, *args):
        """
        Args is an array of either 3 or 5 elements. First three elements are
        the Guassian amplitude, center time, and width. The last two optional
        elements are the y-intercept and slope for a linear trend. 
        """
        exp_arg = -(t-args[1])**2/(2*args[2]**2)
        y = args[0]*np.exp(exp_arg)

        if len(args) == 5:
            y += args[3] + t*args[4]
        return y

    def fit_test_plot(self, peak_time, time_range, popt, r2, adj_r2, ax=None):
        """
        Make a test plot of the microburst fit and annotate the fit 
        parameters.
        """
        if ax is None:
            _, ax = plt.subplots()
        plot_time_range = [
            peak_time - pd.Timedelta(seconds=self.plot_width_s/2),
            peak_time + pd.Timedelta(seconds=self.plot_width_s/2)
        ]

        time_array = self.hilt_df.loc[plot_time_range[0]:plot_time_range[-1]].index
        current_date = time_array[0].floor('d')
        x_data_seconds = (time_array-current_date).total_seconds()
        y_data = self.hilt_df.loc[plot_time_range[0]:plot_time_range[1], 'counts']

        popt[1] = (popt[1] - current_date).total_seconds()
        popt[2] = popt[2]/2.355 # Convert the Gaussian FWHM to std

        gaus_y = SAMPEX_Microburst_Widths.gaus_lin_function(x_data_seconds, *popt)
        ax.plot(time_array, y_data, c='k')
        ax.plot(time_array, gaus_y, c='r')

        for t_i in time_range:
            ax.axvline(t_i, c='g')

        s = (f'R^2 = {round(r2, 2)}\n'
            f'adj R^2 = {round(adj_r2, 2)}\n'
            f'A = {round(popt[0])} [counts]\n'
            f't0 = {round(popt[1])} [sec_of_day]\n'
            f'FWHM = {round(popt[2]*2.355, 2)} [s]')
        if len(popt) == 5:
            s += (f'\ny-intercept = {round(popt[3])}\n'
                  f'slope = {round(popt[4])}')
        ax.text(0, 0.95, s, va='top', ha='left', transform=ax.transAxes, fontsize=12)

        ax.set(title=f'SAMPEX microburst fit\n{peak_time}', ylim=(y_data.min(), 1.1*y_data.max()))
        plt.show()
        return

    def goodness_of_fit(self, y_true, y_pred, n_params):
        """
        Method to calculate the R^2 coefficient of determination
        and the adjusted R^2 coefficient given the number
        of fit parameters n_params.
        """
        r2 = sklearn.metrics.r2_score(y_true, y_pred)
        n = len(y_true)
        adj_r2 = 1 - (1-r2)*(n-1)/(n-1-n_params)
        return r2, adj_r2
