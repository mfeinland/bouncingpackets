import numpy as np
import pandas as pd
import sampex
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib as mpl
import geopandas as gpd
import os.path

# Quick-look plots for events in Feinland et al. 2025 (?)
# Date created: 12/27/24
# Last modified: 1/8/24
# Author: Max Feinland
# Purpose: Runs through each id'ed "good" event, and creates plots showing general trend, location, and portion inside loss cone.
# Option for "lunjin," which selects a few events to show to Lunjin Chen for some potential modeling.

save = False # save the figures to .jpg files? (set to False if you're debugging)
lunjin = True # select only a few events to show Lunjin

# Read in events (catalog from Feinland & Blum, 2025-- you can get this file from the bouncingpackets github in the "variable shapes" folder)
events = pd.read_csv("Data_Files/good_events.csv", index_col=0)
t = pd.to_datetime(events.t)
r = range(len(t))

if lunjin:
    laurens_choice = [2, 6, 10, 20, 55, 70, 89, 93, 97] # indices of events that Lauren selected to show Lunjin
    r = laurens_choice

# How should you plot each shape classification?
plot_dict = {'decr': ['v', 100],
             'half': ['s', 100],
             'crown': ['*', 250],
             'other': ['o', 100]}

# For loop to plot each cataloged event
for i in r:
    ## First plot: show events
    fig, ax = plt.subplots(1, 2, figsize=(10,3), width_ratios=[1, 1.4]) # initialize the figure

    tstmp = t[i]
    date_str = tstmp.strftime('%B %d, %Y') # string with date of event

    # Load the SAMPEX/HILT data (thanks Mike!)
    h = sampex.HILT(tstmp)
    h.load()

    # Find the timestamp of the event for HILT
    sidx = np.where(h['time'] == tstmp)[0][0] - 20 # plotting start index
    eidx = sidx+160 # plotting end index

    hilt_sample_rate = 50 # samples per second
    hilt_geometric_factor = 60 # cm^2 sr

    # Correct for sampling rate and geometric factor
    hilt_corrected = h['counts'][sidx:eidx]*hilt_sample_rate/hilt_geometric_factor

    # Now plot the event in time
    ax[0].plot(h['time'][sidx:eidx], hilt_corrected, linewidth=2, label='HILT')

    if lunjin:
        # Load the SAMPEX/PET data (thanks Mike!)
        p = sampex.PET(tstmp)
        p.load()

        # Find the timestamp of the event for LICA
        sidx_pet = np.argmin(np.abs(p['time'] - tstmp)) - 4 # plotting start index
        eidx_pet = sidx_pet+32 # plotting end index

        pet_sample_rate = 10 # samples per second
        pet_geometric_factor = 10 # cm^2 sr

        # Correct for sampling rate and geometric factor
        pet_corrected = p['counts'][sidx_pet:eidx_pet]*pet_sample_rate/pet_geometric_factor

        # Now plot the event in time
        ax[0].plot(p['time'][sidx_pet:eidx_pet], pet_corrected, linewidth=2, label='PET')

    # Axis label formatting
    ax[0].set_xlabel(f'Time on {date_str}')
    ax[0].set_ylabel('Electrons (#/s cm$^2$ sr)')
    ax[0].xaxis.set_major_locator(dates.SecondLocator(interval=1))
    ax[0].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S')) 
    
    # Add L, MLT, period labels on plot
    L_label = 'L = ' + str(events.L[i]) + '\n' # label of L-shell
    MLT_label = 'MLT = ' + str(events.MLT[i]) + '\n' # label of MLT
    per_label = '$t_{obs}$ = ' + str(events.per[i]) + 's \n' # label of observed period
    T05_label = '$t_{pred}$ = ' + str(np.round(events.T05[i], 2)) + 's' # label of predicted period
    ypos = 0.7*np.ptp(hilt_corrected) + min(hilt_corrected) # y-position of text
    xpos = h['time'][sidx+120] # x-position of text
    ax[0].text(xpos, ypos, L_label + MLT_label + per_label + T05_label) # add text to plot

    if lunjin:
        ax[0].legend() # To differentiate the two count rates

    ## Second plot: event in space
    # Plot world
    world = gpd.read_file('Data_Files/ne_110m_admin_0_countries_lakes.shp') # you can get this file from my github or Natural Earth
    world.plot(ax=ax[1], color='none', edgecolor='black')

    # Initialize the colorbar parameters
    cmap = mpl.cm.plasma # set colormap to plasma 
    norm = mpl.colors.Normalize(vmin=0, vmax=1) # normalize colorbar for each plot to go from 0 to 1 (range of p_BLC)

    # Plot event in space
    e=ax[1].scatter(events['lon'][i], events['lat'][i], c=events['pl1'][i], s=plot_dict[events['shapes'][i]][1], 
                    marker=plot_dict[events['shapes'][i]][0], norm=norm, cmap=cmap, edgecolors='navy')

    # Label formatting
    ax[1].set_xlim(-180, 180)
    ax[1].set_xticks(np.arange(-180, 181, 60))
    ax[1].set_yticks(np.arange(-90, 91, 30))
    ax[1].set_xlabel('Longitude')
    ax[1].set_ylabel('Latitude')

    # Add colorbar
    cbar = plt.colorbar(e, ax=ax[1], orientation='vertical', label='$p_{BLC}$')
    plt.tight_layout()
    plt.show() # Show plot

    # Save plot to .jpg file
    if save:
        if lunjin:
            filename = 'Figures/qlp_pet/' + tstmp.strftime('%Y_%m_%d') + '.jpg'
            fig.savefig(filename)
        else:
            # filename would be something like 'Figures/quicklookplots/2000_01_23.jpg'
            filename = 'Figures/quicklookplots/' + tstmp.strftime('%Y_%m_%d') + '.jpg'
            if os.path.exists(filename): # if date already exists
                filename = 'Figures/quicklookplots/' + tstmp.strftime('%Y_%m_%d') + '_1.jpg' 
            fig.savefig(filename)
