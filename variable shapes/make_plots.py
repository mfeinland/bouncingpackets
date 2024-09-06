# Make plots to visualize data collected by findbouncingpackets search script
# Date created: 7/2/2024 (or similar)
# Last modified: 9/6/2024
# Author: Max Feinland for Blum Research Group

'''Example event plots'''
def figure_1():

    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import sampex # thanks to Mike Shumko for making this package
    import matplotlib.dates as dates
    from scipy.signal import find_peaks

    plt.rcParams["font.family"] = "Arial" # bc I don't like the default Python font

    def plot_date(ax, j):
        # Function that makes plots of events
        # Inputs: axes object (which subplot to plot on), iterator variable, 
        # dictionary containing specifics for each plot

        b1 = 25 # buffer between start time and start of plotting (in indices)
        b2 = 150 # buffer between stop time and end of plotting (in indices)
        t = plotting_data["t"][j] # timestamp of event start

        cur_t = pd.to_datetime(t) # convert text to DateTime
        h = sampex.HILT(cur_t) # request SAMPEX count rate data
        h.load()

        # find index of where SAMPEX data time matches timestamp
        m = np.where(h.times == cur_t)[0][0]
        d_plot = h.times[m-b1:m+b2]
        r_plot = h.counts[m-b1:m+b2]

        maxdis = max(r_plot) - min(r_plot) # range of interval
        # identify peaks in interval
        [pks, _] = find_peaks(r_plot, prominence=0.25*maxdis, distance=3)
        if j in [0, 2, 4]:
            # The actual process of peaks found is a lot more finicky than this,
            # and can be found in the bouncingpacketfunctions search script (bouncingPackets()). 
            # This version runs a lot faster, albeit less accurately, so these lines are just fixing
            # any discrepancy between this method and the more robust way.
            pks = pks[0:4]
        elif j in [3, 7]:
            pks = pks[1:5]

        # for plotting
        txt_loc_x = h.times[m-b1]
        txt_loc_y = 0.8*(max(r_plot) - min(r_plot)) + min(r_plot)
        xlbel = "Time (UTC) on " + cur_t.strftime('%m/%d/%Y')

        ax.plot(d_plot, r_plot)
        ax.plot(d_plot[pks], r_plot[pks], 'r*', markersize=10)
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S')) 
        ax.set_xlabel(xlbel, fontsize=16)
        ax.grid()
        ax.set_ylabel("Count rate (#/20ms)", fontsize=16)
        ax.text(x=txt_loc_x, y=txt_loc_y, s=plotting_data["letter"][j], fontsize=20)
        ax.tick_params(labelsize=12)

    # read in data
    data = pd.read_csv("all_events_v2.csv", index_col=0)
    good_data = data[data.final_eye<2] # restrict to good events only
    good_data = good_data.reset_index(drop=True)

    idx =  [15, 25, 43, 12, 11, 35, 33, 64] # indices of events I picked as examples
    plotting_data = {"t": [good_data.t[x] for x in idx], 
                     "letter": ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)']}

    fig, ax = plt.subplots(4,2, figsize=(15,12), constrained_layout=True)
    ax_flat = ax.flatten()

    # iterate through axes object and plot each example event
    for j, a in enumerate(ax_flat):
        plot_date(a, j)

    
'''MLT & L-shell histograms'''
def figure_2():

    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from datetime import datetime

    plt.rcParams["font.family"] = "Arial"

    # read in data
    data = pd.read_csv("all_events_v2.csv") # my data
    ref = pd.read_csv("microburst_catalog_00.txt", sep=",") # reference data

    # limit reference data to time surveyed by my search script
    time_needed = pd.to_datetime(ref['dateTime'])
    idx = np.where((time_needed >= datetime(2000, 1, 1)) & (time_needed <= datetime(2003, 12, 31)))[0]
    ref = ref.iloc[idx,:]

    fig, ax = plt.subplots(2, 1, figsize=(14,11))
    def make_hist(ax1, param, rmin, rmax, bin_num):
        # Function that makes histograms of variables of your choosing
        # Inputs: axes object (which subplot to plot on), variable name, 
        # minimum range, maximum range, number of bins

        ax2 = ax1.twinx() # make second axis for reference data

        # plot the events
        ax1.hist(data[param][data.final_eye<2], bins=bin_num, range=(rmin,rmax), alpha=0.5,
                 label='good bouncing packets', facecolor='black', edgecolor='black', zorder=50)
        ax2.hist(ref[param], bins=bin_num, range=(rmin, rmax), alpha=0.3, label='all microbursts', 
                 facecolor='red', edgecolor='red', zorder=20) # plot reference data

        # tick/label stuff (formatting)
        ax1.tick_params(labelsize=12)
        ax1.set_xlabel(param, fontsize=14)
        ax1.set_ylabel('Counts (#), good bouncing packets', fontsize=14)
        ax1.grid()

        ax2.set_ylabel('Counts (#), all microbursts', color="lightcoral", fontsize=14)
        ax2.tick_params(axis='y', labelcolor="lightcoral", labelsize=12)

        # adding subplot label
        if param=="MLT":
            ax1.text(0, 19, "a)", fontsize=20)
            ax1.set_xlabel("Magnetic Local Time", fontsize=14) # just because the variable is named
            # "MLT" in the datasets, but I wanted the x-axis label to be more descriptive
        else:
            ax1.text(1, 40, "b)", fontsize=20)
            ax1.set_xlabel("L-shell", fontsize=14) # just because the variable is named
            # "L" in the datasets, but I wanted the x-axis label to be more descriptive


    # Call function for each variable
    make_hist(ax[0], "MLT", 0, 24, 24)
    make_hist(ax[1], "L", 1, 8, 14)
    

'''Period vs. L-shell'''
def figure_3():
    
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import matplotlib.dates as dates
    from matplotlib.legend_handler import HandlerTuple
    import ast

    plt.rcParams["font.family"] = "Arial"
    
    def determine_match(full_diff, half_diff, full_yn, half_yn, pred):  
        if full_yn and half_yn: # if both half & full-per matches
            if min(full_diff) < min(half_diff):
                return 1, pred
            else:
                return 2, pred/2
        elif full_yn: # if full per match only
            return 1, pred
        elif half_yn: # if half per match only
            return 2, pred/2
        else: # no match, see which is closer
            if min(full_diff) < min(half_diff):
                return 0, pred
            else:
                return 0, pred/2

    def plot_specific_shape(ax, k):
        # Function that plots bounce period vs. L shell for a specific shape
        # Inputs: axes object (which subplot to plot on), iterator variable

        current_shape = txt["shapes"][k] # pull stuff out of dict
        current_letter = txt["letter"][k]

        # find indices containing allowable shapes
        indices = np.where(good_data.shapes==current_shape)[0]

        for j in indices:
            hilt_uncertainty = 0.02 # time resolution of instrument
            err = np.sqrt(pers.spread[j]**2 + hilt_uncertainty**2) # propagate error

            fulldiff = abs(pers.dt[j] - good_data.tb[j])
            halfdiff = abs(pers.dt[j] - good_data.tb[j]/2)

            # are any of the spacings within allowable error bounds? 
            fullyn = any([x <= err for x in fulldiff])
            halfyn = any([x <= err for x in halfdiff])
            
            sort, tb = determine_match(fulldiff, halfdiff, fullyn, halfyn, good_data.tb[j])  
            
            props = [['rosybrown', 'dodgerblue', 'limegreen'], [5, 10, 10]]

            # plotting
            ax.plot(good_data.L[j]*np.ones(len(pers.dt[j])), pers.dt[j], 'o', color=props[0][sort], 
                    alpha = 0.8, markersize=10, zorder=props[1][sort])
            ax.errorbar(good_data.L[j], tb, yerr=err, fmt='s', 
                     markerfacecolor='none', capsize=4, markeredgecolor=props[0][sort],
                     zorder=(props[1][sort]+5), markersize=10, color='black', 
                        markeredgewidth=2, alpha=0.8)

        if current_shape == 'decr':
            shape_label = 'decreasing' # to put on the plot as an ID
        elif current_shape == 'incr':
            shape_label = 'increasing' # to put on the plot as an ID
        else:
            shape_label = current_shape
            
        all_label = current_letter + ' - ' + shape_label
        ax.text(2, 0.9, all_label, fontsize=25) # subplot label
        ax.grid()

        if k in [4, 5, 6, 7]: # put xlabel on bottom 4 plots
            ax.set_xlabel("L-shell", fontsize=18)
        if k in [0, 4]: # put ylabel on leftmost 2 plots
            ax.set_ylabel("Mean period (s)", fontsize=18)
        ax.set_ylim(0, 1)
        ax.set_xlim(1.5, 7)
        ax.xaxis.set_tick_params(labelsize=15)
        ax.yaxis.set_tick_params(labelsize=15)

    # Import data
    data = pd.read_csv("all_events_v2.csv", index_col=0)
    good_data = data[data.final_eye<2] # restrict to good events
    good_data = good_data.reset_index(drop=True) # reset index

    # spread in model predictions
    pers = pd.read_csv("model_preds_and_spacing.csv",index_col=0)
    models  = pers[['T89', 'T05', 'OP', 'SL']]
    pers['spread'] = models.max(axis=1) - models.min(axis=1)

    # spacings between peaks for each event
    pers['dt'] = pers['dt'].apply(ast.literal_eval)

    # Making figures
    fig = plt.figure(layout='constrained', figsize=(20,18))
    subfigs = fig.subfigures(2, 1, wspace=0.07)

    ## First figure: bounce period histograms
    ax1 = subfigs[0].subplots()
    newper = [np.mean(x) for x in pers.dt]

    # Calculate period ratios for each type 
    good_ratio = np.divide(newper, good_data.tb)
    ok_ratio = np.divide(data.per[(data.final_eye>1.5) & (data.final_eye<3)], 
                         data.tb[(data.final_eye>1.5) & (data.final_eye<3)])
    bad_ratio = np.divide(data.per[data.final_eye==3], data.tb[data.final_eye==3])

    ax1.hist(good_ratio, bins=15, range=(-0.05,1.45), alpha=0.5, facecolor='cornflowerblue', edgecolor='navy',
            zorder=10, label='good') # good histogram
    ax1.hist(ok_ratio, bins=15, range=(-0.05,1.45), alpha=0.5, facecolor='goldenrod', edgecolor='orangered', 
             zorder=5, label='okay') # ok histogram
    ax1.hist(bad_ratio, bins=15, range=(-0.05,1.45), alpha=0.5, facecolor='black', edgecolor='black',
             zorder=0, label='bad') # bad histogram

    # axis labels & ticks
    ax1.set_xlabel("Ratio of observed period to predicted period", fontsize=20)
    ax1.set_ylabel("Number of observations", fontsize=20)
    ax1.xaxis.set_tick_params(labelsize=20)
    ax1.yaxis.set_tick_params(labelsize=20)
    ax1.legend(fontsize=20)
    ax1.grid()
    ax1.set_xlim(0, 1.5)
    ax1.set_ylim(0, 60)
    ax1.text(0.05, 54, 'a)', fontsize=30) # subplot label

    # Second figure: bounce period vs. L for each shape
    ax2 = subfigs[1].subplots(2, 4, sharex=True, sharey=True)
    ax_flat = ax2.flatten()

    txt = {"shapes": ["decr", "blake", "crown", "flat", "incr", "half", "smile", "other"],
          "letter": ["b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)"]}

    full_color = 'dodgerblue'
    half_color = 'limegreen'
    else_color = 'rosybrown'

    # iterate through axes object and plot each shape type
    for j, a in enumerate(ax_flat):
        plot_specific_shape(a, j)


    p = ax2[0,3]
    p1, = p.plot(0, 0, 'o', color=full_color, alpha=0.8, markersize=10)
    p2, = p.plot(0, 0, 'o', color=half_color, alpha=0.8, markersize=10)
    p3, = p.plot(0, 0, 'o', color=else_color, alpha=0.8, markersize=10)
    f = p.errorbar(0, 0, yerr=0.2, fmt='s', markeredgecolor=full_color, markerfacecolor='none',
                   color='black', markersize=10, capsize=4, markeredgewidth=2)
    h = p.errorbar(0, 0, yerr=0.2, fmt='s', markeredgecolor=half_color, markerfacecolor='none',
                   color='black', markersize=10, capsize=4, markeredgewidth=2)
    x = p.errorbar(0, 0, yerr=0.2, fmt='s', markeredgecolor=else_color, markerfacecolor='none',
                   color='black', markersize=10, capsize=4, markeredgewidth=2)
    l = p.legend([(p1, p2, p3), f, h, x], # i made the legend look nice
             ['Observed spacing', 'T05 model', 'T05 model (half)', 'T05 model (no match)'], 
             handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize=16) 
    

'''Geographical location'''
def figure_4():

    import numpy as np
    import pandas as pd
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import ast

    plt.rcParams["font.family"] = "Arial" # because I don't like the default Python font
    
    def determine_match(full_diff, half_diff, full_yn, half_yn):  
        if full_yn and half_yn: # if both half & full-per matches
            if min(full_diff) < min(half_diff):
                return 1
            else:
                return 2
        elif full_yn: # if full per match only
            return 1
        elif half_yn: # if half per match only
            return 2
        else: # no match
            return 0


    def make_losscone_map(ax):
        # Function that makes world outline and L-shell contours
        # Inputs: axes object (which subplot to plot on)

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world.plot(ax=ax, color='none', edgecolor='black', zorder=40)

        # read in SAMPEX attitude data. If you don't have this file you can download it from here:
        # https://izw1.caltech.edu/sampex/DataCenter/DATA/PSSet/Text/PSSet_6sec_2000022_2000048.txt
        # Make sure to change the pathing so it works with your machine.
        att = pd.read_csv('C:/Users/maxim/sampex-data/Attitude/PSSet_6sec_2000022_2000048.txt', 
                         sep=' ', header=60, on_bad_lines='skip')

        cols = np.array([7, 8, 35]) # take out lat, lon, losscone 2
        att = att.iloc[:150000, cols] # don't need the whole thing, just about 14 days for full coverage
        att.columns = ['lon', 'lat', 'losscone']

        # change longitude to be -180 to 180
        long_idx = np.where(att.lon > 180)[0]
        att.lon[long_idx] = att.lon[long_idx] - 360

        sc = ax.scatter(att.lon, att.lat, c=att.losscone, s=25, vmin=30, vmax=90, zorder=5, cmap='gray')
        c = plt.colorbar(sc)
        c.set_label("Losscone angle (deg)", fontsize=16)

        lgrid = pd.read_csv('Lgrid.dat', delimiter='\t', header=None)

        for i in np.arange(1, 30, 2):
            min_pos = np.argmin(lgrid.iloc[:,i])
            latl = np.concatenate(([lgrid.iloc[min_pos:,i-1], lgrid.iloc[:min_pos,i-1]]))
            lonl = np.concatenate(([lgrid.iloc[min_pos:,i], lgrid.iloc[:min_pos,i]]))
            ax.plot(lonl, latl, '--', color="white", zorder=45)

        ax.text(-180, 7, "Magnetic Equator", fontsize=16, zorder=50)
        ax.text(-197, 48, "L = 2", fontsize=14, zorder=50)
        ax.text(-197, -40, "L = 2", fontsize=14, zorder=50)
        ax.text(-197, 59, "L = 3", fontsize=14, zorder=50)
        ax.text(-197, -55, "L = 3", fontsize=14, zorder=50)
        ax.text(-197, 64, "L = 5", fontsize=14, zorder=50)
        ax.text(-197, -60, "L = 5", fontsize=14, zorder=50)
        ax.text(-197, 71, "L = 8", fontsize=14, zorder=50)
        ax.text(-197, -66, "L = 8", fontsize=14, zorder=50)

        ax.set_xlabel("Longitude", fontsize=20)
        ax.set_ylabel("Latitude", fontsize=20)
        ax.set_ylim(-85, 85)
        return ax

    # Import data
    data = pd.read_csv("all_events_v2.csv", index_col=0)
    good_data = data[data.final_eye<2] # restrict to good events
    good_data = good_data.reset_index(drop=True) # reset index

    # spread in model predictions
    pers = pd.read_csv("model_preds_and_spacing.csv",index_col=0)
    models  = pers[['T89', 'T05', 'OP', 'SL']]
    pers['spread'] = models.max(axis=1) - models.min(axis=1)

    # spacings between peaks for each event
    pers['dt'] = pers['dt'].apply(ast.literal_eval)
    newper = [np.mean(x) for x in pers.dt] # calculate new mean spacing


    fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(18, 14))

    # call losscone map function
    ax1 = make_losscone_map(ax1)
    ax2 = make_losscone_map(ax2)

    fullper = np.zeros(len(good_data)) # initialize vector containing full/half classification

    # create full/half classification
    # this is kept dynamic so you can change the error condition if you need to
    for j in range(len(good_data.tb)):
            hilt_uncertainty = 0.02 # time resolution of instrument
            err = np.sqrt(pers.spread[j]**2 + hilt_uncertainty**2) # propagate error

            fulldiff = abs(pers.dt[j] - good_data.tb[j])
            halfdiff = abs(pers.dt[j] - good_data.tb[j]/2)

            # are any of the spacings within allowable error bounds? 
            fullyn = any([x <= err for x in fulldiff])
            halfyn = any([x <= err for x in halfdiff])
            
            fullper[j] = determine_match(fulldiff, halfdiff, fullyn, halfyn)

    # subplot 1: half/full/fail
    full = ax1.scatter(good_data.lon[fullper==1], good_data.lat[fullper==1], 
                      s=100,  c='dodgerblue', label='full period', zorder=55, edgecolor='navy')
    half = ax1.scatter(good_data.lon[fullper==2], good_data.lat[fullper==2], 
                      s=100,  c='limegreen', label='half period', zorder=50, edgecolor='darkgreen')
    fail = ax1.scatter(good_data.lon[fullper==0], good_data.lat[fullper==0], 
                      s=100,  c='rosybrown', label='no match', zorder=45, edgecolor='maroon')
    ax1.legend(loc='right', fontsize=14).set_zorder(50)

    # subplot 2: shapes
    decr = ax2.scatter(good_data.lon[good_data.shapes=="decr"], good_data.lat[good_data.shapes=="decr"], 
                      s=100,  c='royalblue', label='decreasing', zorder=50, edgecolor='navy')
    blake = ax2.scatter(good_data.lon[good_data.shapes=="blake"], good_data.lat[good_data.shapes=="blake"],
                       s=100, c='mediumslateblue',  label='blake', zorder=50, edgecolor='rebeccapurple')
    crown = ax2.scatter(good_data.lon[good_data.shapes=="crown"], good_data.lat[good_data.shapes=="crown"],
                       s=100, c='fuchsia', label='crown', zorder=50, edgecolor='purple')
    flat = ax2.scatter(good_data.lon[good_data.shapes=="flat"], good_data.lat[good_data.shapes=="flat"],
                       s=100, c='turquoise', label='flat', zorder=50, edgecolor='darkslategray')
    incr = ax2.scatter(good_data.lon[good_data.shapes=="incr"], good_data.lat[good_data.shapes=="incr"],
                       s=100, c='orangered', label='increasing', zorder=50, edgecolor='maroon')
    smile = ax2.scatter(good_data.lon[good_data.shapes=="smile"], good_data.lat[good_data.shapes=="smile"],
                       s=100, c='purple', label='smile', zorder=50, edgecolor='darkmagenta')
    half = ax2.scatter(good_data.lon[good_data.shapes=="half"], good_data.lat[good_data.shapes=="half"],
                       s=100, c='goldenrod', label='half', zorder=50, edgecolor='darkgoldenrod')
    other = ax2.scatter(good_data.lon[good_data.shapes=="other"], good_data.lat[good_data.shapes=="other"],
                       s=100, c='lawngreen', label='other', zorder=50, edgecolor='olivedrab')

    ax2.legend(loc='right', fontsize=14).set_zorder(50)
