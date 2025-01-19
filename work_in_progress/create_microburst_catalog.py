# Scripts to find microbursts using modified O'Brien et al. (2003) criteria
# Author: Max Feinland for Blum Research Group, LASP

# Housekeeping
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sampex
import os.path as op

def obrien(data):
    # Author: Max Feinland, adapted from O'Brien et al. 2003

    def changeCadence(time20, times100):
        # Author: Max Feinland
        # times100 is the list of timestamps for whatever it is: starts, ends, ns
        # times20 is the entire list of time20 for the day in question
        try:
            # try list comprehension first, it's faster
            idx20 = [time20.get_loc(tstmp) for tstmp in times100] 
            print('i guess this needs to be in here', times20[0], times100[0])
        except: 
            idx20 = [] # clear and start over
            idx20 = [np.abs((time20 - target).values).argmin() for target in times100] # if a timestamp is missing from time20,
            # you will have to find the minimum timestamp (that is closest in time)
        return idx20

    def findActualPeaks(idx20, counts):
        # Author: Max Feinland
        # Date created: 1/9/24
        actual_idx = np.zeros(len(idx20))
        # Find the actual peak in the interval, since it could be anything in the 100-ms window
        start = [x - 5 for x in idx20]
        end = [x + 5 for x in idx20]

        for i in range(len(idx20)):
            interval = counts[start[i]:end[i]]
            actual_idx[i] = int(np.argmax(interval)) + start[i]
        return actual_idx
        
    def qualityCheck(t, x):
        # Author: Max Feinland
        # Date created: 12/24/24
        not_coarse = (len(np.unique(x)) > 5)*1 # are there more than 5 unique values in the interval?
        continuous = (len(np.where(t.to_series().diff() > timedelta(seconds=1))[0]) == 0)*1 # is the interval more or less continuous?
        ok_counts = (max(x) > 100)*1 # is this a higher-count rate interval?
        return not_coarse, continuous, ok_counts
    
    N20 = data['counts'] # count rate sampled every 20 ms
    time20 = data['time'] # time every 20 ms

    # I'm sure there's a more efficient way to do this, but I don't know it
    df = pd.DataFrame({'time': data['time'], 'counts': data['counts']})

    df.set_index('time', inplace=True) # set time column as the index

    # resample the dataframe to 100 ms intervals and sum the counts in each interval
    N100 = df.resample('100ms').sum()

    A500 = N100.rolling(5, center=True).mean() # 5-observation centered rolling mean (over 500 ms)

    condition = np.divide((N100.counts - A500.counts), np.sqrt(1 + A500.counts)) # O'Brien et al 2003

    ns = np.argwhere(condition > 5) # find indices where O'Brien condition is true
    ns = [item[0] for item in ns]

    ns20 = changeCadence(time20, N100.index[ns]) # change the cadence from 100ms to 20ms

    final_idx = findActualPeaks(ns20, h['counts']) # run actual peak finder
    final_idx = np.unique(final_idx) # get rid of repeats

    qc = {'notcoarse': [], 'continuous': [], 'highcounts': []} # perform quality check

    for i in range(len(final_idx)):
        # Create 1-second interval to run quality check
        sidx = int(final_idx[i]-10)
        eidx = int(final_idx[i]+40)

        # Pull out variables in this interval
        t = time20[sidx:eidx]
        x = h['counts'][sidx:eidx]

        # Run quality check
        not_coarse, continuous, high_counts = qualityCheck(t, x)

        # Append results to the dictionary
        qc['notcoarse'].append(not_coarse)
        qc['continuous'].append(continuous)
        qc['highcounts'].append(high_counts)

    qc = pd.DataFrame(qc) # convert dictionary to dataframe
    return final_idx, qc

def isleap(year): # self explanatory?
    return year % 4 == 0

def get_att_rows():
    # Create dataframe with dates of each SAMPEX attitude file grouping (necessary to find the right filename)
    dates = pd.date_range(start='6/8/1996', end='11/13/2012', freq='27d')
    att_rows = pd.DataFrame({'year': dates.year, 'day': dates.day_of_year, 'date': dates}) # put in dataframe
    return att_rows

def load_att(day, current_row):
    # Call get_att_rows to pull out the filename
    att_rows = get_att_rows()

    # Access filename
    y1 = str(att_rows['year'][current_row]) # First year in att filename
    y2 = str(att_rows['year'][current_row+1]) # Second year in att filename
    d1 = str(att_rows['day'][current_row]).rjust(3, '0') # First day-of-year in att filename, zero-padded
    d2 = str(att_rows['day'][current_row+1]-1).rjust(3, '0') # Second day-of-year in att filename, zero-padded
    att_filename = '/Users/maxim/sampex-data/Attitude/PSSet_6sec_' + y1 + d1 + '_' + y2 + d2 + '.txt'
    
    if current_row < 49: # if before the year 2000, file will be zipped
        att_filename = att_filename + '.zip'

    # take out year, day of year, seconds, latitude, longitude, altitude, L-shell, # magnetic field magnitude (B-mag),
    # magnetic local time (MLT), losscone 1, losscone 2, pitch, SAA flag, attitude flag
    cols = np.array([0, 1, 2, 7, 8, 9, 20, 21, 22, 34, 35, 58, 68, 71])

    # Read in file
    att = pd.read_csv(att_filename, sep=' ', skiprows=60, on_bad_lines='skip', header=None, usecols=cols)

    # These lines modified from Mike Shumko's sampex package -- turn year, doy, s columns into Timestamp and set as index
    year_doy = [f"{year}-{doy}" for year, doy in att.iloc[:,[0, 1]].values]
    att_dates = pd.to_datetime(year_doy, format="%Y-%j")
    att.index = att_dates + pd.to_timedelta(att.iloc[:,2], unit="s")

    # Drop year, day-of-year, second-of-day columns
    att.drop([0, 1, 2], axis=1, inplace=True)
    
    # Rename columns
    colnames = ['lon', 'lat', 'alt', 'L', 'B', 'MLT', 'losscone1', 'losscone2', 'SAA', 'pitch',  'att']
    att.columns = colnames
    
    # change longitude to be -180 to 180 (rather than 0 to 360) -- also adapted from sampex package
    att['lon'] = np.mod(att['lon'] + 180, 360) - 180
    return att, colnames

'''The driver script'''
# Runs O'Brien & saves to file.
# Author: Max Feinland for Blum Research Group, LASP
# Last modified: 1/19/25

print("Welcome to the inner_belt_microbursts search script.", flush=True)
print("Coded by Max Feinland for Blum Research Group, 2023-2025.", flush=True)

# Initializing these variables
hilt_already_loaded = None
att_already_loaded = None
save = True # do you want to save to a file or not? (make False if debugging)

filename = "Data_Files/rhs_5_v4.csv" # filename to save to

list_of_days = pd.date_range(start='8/16/1996', end='11/13/2012')
last_t = list_of_days[0] # first date in State 4 (20 ms cadence)

# Skip to most recently run day
try: # Try to open file containing last date run
    last_day_ran = open("Data_Files/last_date.txt").read()
    # Read file and increment day by 1
    start_day = datetime.strptime(last_day_ran, '%Y-%m-%d\n') + timedelta(days=1)
    # Find the index in the list_of_days corresponding to this day
    start_idx = np.where(list_of_days >= start_day)[0][0]
    # Restrict to this date and on
    list_of_days = list_of_days[start_idx:]
except: # there was no last day ran, so no file
    start_day = last_t # start day is first day in dataset; use entire list of days

# Now we are getting to the the function calls, etc
att_rows = get_att_rows() # call list of year-doy pairs
    
for j, day in enumerate(list_of_days):
    formatted_date = day.strftime("%B %d, %Y")
    print(f'\nRunning for {formatted_date}...')

    # find the right row in the ephemeris files for today's date
    current_row = np.where(att_rows['date'] <= day)[0][-1]

    # Check to see if attitude loading needed
    if 'old_row' in globals(): # If already ran something:
        if current_row != old_row: # If your current row is no longer the same as the last one:
            print("Loading attitude data...")
            a, cols = load_att(day, current_row) # load attitude data
            old_row = current_row
    else: # If you didn't run anything
        print("Loading attitude data...")
        a, cols = load_att(day, current_row) # load attitude data
        old_row = current_row

    # Load HILT data if needed
    if hilt_already_loaded is not True:
        try:
            print("Loading HILT data...")
            h = sampex.HILT(day) # count rate data (thanks to Mike Shumko for his SAMPEX package)
            h.load()
            proceed = True
            
        except FileNotFoundError:
            print("Looks like this day isn't in the SAMPEX database. :/")
            print("Skipping this day.")
            proceed = False

    if proceed == True:
        print("Running O'Brien algorithm...")

        idx, qc = obrien(h) # call O'Brien function, return microburst indices and quality control variables
        t = [h.times[int(x)] for x in idx] # go from index to timestamp
        
        a_idx = [] # initializing

        # This section is to provide the attitude data for the identified timestamps.
        # First, go through the microburst timestamps and find the 2 closest attitude timestamps: one before, and one after.
        # Append those to the list a_idx.
        for timestamp in t:
            # Right side: the first index in a after the timestamp
            right_idx = np.where(a.index - timestamp >= timedelta(seconds=0))[0][0]
            # Left side: the index before that
            left_idx = right_idx-1
            a_idx.extend([left_idx, right_idx])

        # Restrict dataframe to these times, so it only contains the necessary data.
        # Interpolate linearly between these values, just need the two closest times for each microburst.
        att_data = a.iloc[a_idx,:]
        
        # Create dataframe of NaNs with the exact microburst timestamps and all the same column names as att_data 
        nans = pd.DataFrame(np.nan, index=t, columns=cols)

        # Now, interpolate!
        att_data = pd.concat([att_data, nans]) # concatenate nan matrix and att data
        att_data.sort_index(inplace=True) # sort values by index
        att_data.interpolate(inplace=True, method='time') # interpolate
        att_data = att_data.loc[t, :] # drop everything but interpolated values

        # add quality control variables
        qc.index = t # set quality control dataframe index to be the microburst timestamps
        att_data = pd.concat([att_data, qc], axis=1) # concatenate 
        att_data = att_data.drop_duplicates() # It is possible that a microburst timestamp coincides with an attitude timestamp,
        # In which case there will be a repeat.

        att_data = att_data.round(3) # Round everything down to 3 decimal places
        att_data['SAA'] = att_data['SAA'].round().astype(int) # Round quality flags to int
        att_data['att'] = att_data['att'].round().astype(int)

        if save: # if saving to a file
            if len(att_data) > 0: # if any microbursts detected:
                now = datetime.now().strftime("%B %d, %Y %H:%M:%S") # current timestamp
                print(att_data) # Print the microbursts and their data
                
                # save to file
                if op.isfile('./' + filename):
                    # if file already exists, append
                    att_data.to_csv(filename, mode='a', sep=',', header=False, encoding='utf-8')
                    print("Added to file at " + now + ".")
                else:
                    # if it does not already exist, create it
                    att_data.to_csv(filename, sep=',', encoding='utf-8')
                    print("Saved to file at " + now + ".")
            else:
                print("No inner belt microbursts found on this day.")

            # Create or overwrite file containing the last date successfully ran
            with open("Data_Files/last_date.txt", "w") as text_file:
                print(day.strftime('%Y-%m-%d'), file=text_file)
                    
    last_t = day

print("Queried dates complete.")
