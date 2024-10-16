import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def gaus_overlay_plot(props_df, hilt_df, ax, time_range):
    # Original author: Mike Shumko
    # Date created: unknown
    # Last modified: 10/11/2024 by Max Feinland
    # Purpose: Plots fitted Gaussian curves over raw HILT data.

    plot_start_time = props_df.t0.iloc[0] - pd.Timedelta(seconds=1)
    u = np.where(hilt_df.index.duplicated()==True)[0]
    df_idx = hilt_df.index[u]
    hilt_df = hilt_df.drop(index=df_idx) # get rid of duplicate indices

    closest_start_idx = hilt_df.index.get_indexer([plot_start_time], method='nearest')[0]
    plot_end_time = props_df.t0.iloc[-1] + pd.Timedelta(seconds=1)
    closest_end_idx = hilt_df.index.get_indexer([plot_end_time], method='nearest')[0]
    plot_time_range = [hilt_df.index[closest_start_idx], hilt_df.index[closest_end_idx]]
    
    time_array = hd.loc[plot_time_range[0]:plot_time_range[-1]].index
    current_date = time_array[0].normalize()
    x_data_seconds = (time_array-current_date).total_seconds()
    y_data = hd.loc[plot_time_range[0]:plot_time_range[1], 'counts']
    ax.plot(time_array, y_data, c='k', label='HILT data')
    
    solid_counter = 0
    dashed_counter = 0
    
    for i in range(len(props_df)): # Loop through each peak
        popt = [props_df.iloc[i, 2], props_df.iloc[i, 3], props_df.iloc[i, 4], 
                props_df.iloc[i, 5], props_df.iloc[i, 6]]

        pdate = popt[1] # t0 in datetime terms
        popt[1] = (popt[1] - current_date).total_seconds()
        popt[2] = popt[2]/2.355 # Convert the Gaussian FWHM to std

        gaus_y = Do_Gaussian.gaus_lin_function(x_data_seconds, *popt)

        t_start = np.where(time_array==time_range[i][0])[0][0]
        t_end = np.where(time_array==time_range[i][1])[0][0]


        # If r^2 too low, plot dashed line
        if i in good_pks_f:
            line_style = '-'
            alpha_level = 1
            solid_counter += 1
            if solid_counter > 1: # Check if one of these has been done already.
                labl=None
            else:
                labl='Gaussian fit ($r^2 > 0.9$)'
        else:
            line_style = '--'
            alpha_level = 0.5
            dashed_counter += 1
            if dashed_counter > 1:
                labl = None
            else:
                labl = 'Gaussian fit ($r^2 < 0.9$)'
        ax.plot(time_array[t_start:t_end], gaus_y[t_start:t_end], 
                           c='r', linestyle=line_style, alpha=alpha_level, label=labl)

        baseline = hd.counts[time_range[0][0]]

        ax.text(pdate, popt[0]+baseline, str(i+1))
        # end of for loop
    ax.set_ylim(min(y_data)*0.8, max(y_data)*1.2)
    formatter = mdates.DateFormatter('%H:%M:%S') # 1 decimal place
    ax.xaxis.set_major_formatter(formatter)
    ax.grid()
    plot_date_string = 'Time on ' + current_date.strftime('%Y/%m/%d') 
    ax.set_xlabel(plot_date_string)
    ax.set_ylabel('Counts (#/20 ms)')

    ax.legend() 
    ax.set_xlim(time_array[0], time_array[-1])
    plt.show()
    return 
