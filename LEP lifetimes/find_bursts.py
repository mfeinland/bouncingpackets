# Runs O'Brien & saves to file.
# Author: Max Feinland for Blum Research Group, LASP
# Last modified: 9/12/25

import os.path as op
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sampex

# Define all the functions up here

def obrien(data):
    # Author: Max Feinland, adapted from O'Brien et al. 2003
    # Purpose: Find bursts using O'Brien et al. (2003) condition
    # Inputs: HILT data
    # Outputs: timestamp of each burst

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

        for i in range(len(idx20)): # Create an interval of points to analyze for continuity & coarseness
            try:
                interval = counts[int(start[i]):int(end[i])]
                actual_idx[i] = int(np.argmax(interval)) + start[i]
            except: # this may not work, so scrap that interval (not worth analyzing)
                actual_idx = []
                break
        return actual_idx
        
    def qualityCheck(t, x): # Check that interval is not coarse in countrate and is continuous in time
        # Author: Max Feinland
        # Date created: 12/24/24
        not_coarse = (len(np.unique(x)) > 5)*1 # are there more than 5 unique values in the interval?
        continuous = (len(np.where(abs(t.to_series().diff()) > timedelta(seconds=1))[0]) == 0)*1 # is the interval more or less continuous?
        return not_coarse, continuous
    
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
    final_idx = np.unique(final_idx).astype(int) # get rid of repeats
    counts = h['counts'][final_idx]

    qc = {'notcoarse': [], 'continuous': []} # perform quality check

    for i in range(len(final_idx)):
        # Create 1-second interval to run quality check
        sidx = int(final_idx[i]-25)
        eidx = int(final_idx[i]+25)

        # Pull out variables in this interval
        t = time20[sidx:eidx]
        x = h['counts'][sidx:eidx]

        # Run quality check
        not_coarse, continuous = qualityCheck(t, x)

        # Append results to the dictionary
        qc['notcoarse'].append(not_coarse)
        qc['continuous'].append(continuous)

    qc = pd.DataFrame(qc) # convert dictionary to dataframe
    passed_qc = np.where((qc['notcoarse'] == 1) & (qc['continuous'] == 1))[0] # Find spots where quality check was passed
    final_idx = final_idx[passed_qc] # Restrict indices of bursts and counts to times when quality check was passed
    counts = counts[passed_qc]
    return final_idx, counts

def get_att_rows():
    dates = pd.date_range(start='6/8/1996', end='11/13/2012', freq='27d')
    att_rows = pd.DataFrame({'year': dates.year, 'day': dates.day_of_year, 'date': dates}) # put in dataframe
    return att_rows

def load_att(day, current_row):
    # Get attitude rows to pull out the right filename
    att_rows = get_att_rows()

    # Access filename
    y1 = str(att_rows['year'][current_row]) # First year in att filename
    y2 = str(att_rows['year'][current_row+1]) # Second year in att filename
    d1 = str(att_rows['day'][current_row]).rjust(3, '0') # First day-of-year in att filename, zero-padded
    d2 = str(att_rows['day'][current_row+1]-1).rjust(3, '0') # Second day-of-year in att filename, zero-padded
    att_filename = '/Users/maxim/sampex-data/Attitude/PSSet_6sec_' + y1 + d1 + '_' + y2 + d2 + '.txt'
    
    if current_row < 49: # if before the year 2000, file will be zipped
        att_filename = att_filename + '.zip'

    # take out year, day of year, seconds, latitude, longitude, altitude, L-shell, local magnetic field magnitude,
    # magnetic local time (MLT), losscone 1, losscone 2, pitch, equatorial B field magnitude, SAA flag, attitude flag
    cols = np.array([0, 1, 2, 7, 8, 9, 20, 21, 22, 34, 35, 45, 58, 68, 71])

    # Read in file
    try:
        att = pd.read_csv(att_filename, sep=' ', skiprows=60, on_bad_lines='skip', header=None, usecols=cols)
    except FileNotFoundError:
        a = sampex.Attitude(day) # download the file
        att = pd.read_csv(att_filename, sep=' ', skiprows=60, on_bad_lines='skip', header=None, usecols=cols)
        # Uncomment these two lines if you get this far in the data and it throws you an error
        # att = pd.read_csv('/Users/maxim/sampex-data/Attitude/PSSet_6sec_2008340_2008366.txt',
                          # sep=' ', skiprows=60, on_bad_lines='skip', header=None, usecols=cols)

    # These lines modified from Mike Shumko's sampex package -- turn year, doy, s columns into Timestamp and set as index
    year_doy = [f"{year}-{doy}" for year, doy in att.iloc[:,[0, 1]].values]
    att_dates = pd.to_datetime(year_doy, format="%Y-%j")
    att.index = att_dates + pd.to_timedelta(att.iloc[:,2], unit="s")

    # Drop year, day-of-year, second-of-day columns
    att.drop([0, 1, 2], axis=1, inplace=True)
    
    # Rename columns
    colnames = ['lon', 'lat', 'alt', 'L', 'B', 'MLT', 'losscone1', 'losscone2', 'Beq', 'SAA', 'pitch', 'att']
    att.columns = colnames
    
    # change longitude to be -180 to 180 (rather than 0 to 360) -- also adapted from sampex package
    att['lon'] = np.mod(att['lon'] + 180, 360) - 180
    return att

# Now start the search script.

print("Welcome to the find_bursts search script.", flush=True)
print("Coded by Max Feinland for Blum Research Group, 2023-2025.", flush=True)

# Initializing these variables
debug = False
filename = "Data_Files/inner_belt_bursts.csv" # filename to save to

list_of_days = pd.date_range(start='8/16/1996', end='11/13/2012')

# Skip to most recently run day
try: 
    with open(filename, "r", encoding="utf-8") as file:
        last_line = file.readlines()[-1] # Open last line of filename
    # The first chunk of that line has the last date ran
    last_day_ran = datetime.strptime(last_line[:last_line.index(" ")], "%Y-%m-%d")
    # Increment day by 1
    start_day = last_day_ran + timedelta(days=1)
    # Find the index in the list_of_days corresponding to this day
    start_idx = np.where(list_of_days >= start_day)[0][0]
    # Restrict to this date and on
    list_of_days = list_of_days[start_idx:]
except: # there was no last day ran, so no file
    start_day = list_of_days[0] # start day is first day in dataset; use entire list of days

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
            a = load_att(day, current_row) # load attitude data
            old_row = current_row
    else: # If you didn't run anything
        print("Loading attitude data...")
        a = load_att(day, current_row) # load attitude data
        old_row = current_row

    # Load HILT data if needed
    if debug is not True:
        try:
            print("Loading HILT data...")
            h = sampex.HILT(day) # count rate data (thanks to Mike Shumko for his SAMPEX package)
            h.load()
            proceed = True
            
        except FileNotFoundError:
            print("Looks like this day isn't in the SAMPEX database. :/")
            print("Skipping this day.")
            proceed = False

    if proceed:
        print("Running O'Brien algorithm...")
    
        try:
            idx, counts = obrien(h) # call O'Brien function, return microburst indices and quality control variables
            t = [h.times[int(x)] for x in idx] # go from index to timestamp
            proceed2 = True
        except ValueError:
            proceed2 = False
    
        if proceed2 == True:
            a_idx = [] # initializing
    
            # This section is to provide the attitude data for the identified timestamps.
            # First, go through the microburst timestamps and find the 2 closest attitude timestamps: one before, and one after.
            # Append those to the list a_idx.
            for timestamp in t:
                try:
                    # Right side: the first index in a after the timestamp
                    right_idx = np.where(a.index - timestamp >= timedelta(seconds=0))[0][0]
                    # Left side: the index before that
                    left_idx = right_idx-1
                    a_idx.extend([left_idx, right_idx])
                    proceed3 = True
                    
                except: # There is a problem with the attitude timestamps, so give up and don't interpolate here
                    print("These data are missing from the attitude file--skipping this day.")
                    proceed3 = False
                    break
    
            if proceed3:
                interp_cols = ['lon', 'lat', 'alt', 'L', 'B', 'MLT', 'losscone1', 'losscone2', 'Beq', 'pitch']
                # Restrict dataframe to these times, so it only contains the necessary data.
                # Interpolate linearly between these values, just need the two closest times for each burst.
                att_data = a.iloc[a_idx,:]
        
                saa_att_df = att_data[['SAA', 'att']][::2] # Every second index, pull out the SAA flag and attitude quality flag.
                # (We do this here because you can't interpolate integer values that are flags.)
                saa_att_df.index = t
        
                # Turn counts at microburst time into dataframe so you can concatenate later
                counts = pd.DataFrame({'counts': counts})
                counts.index = t
                
                # Create dataframe of NaNs with the exact microburst timestamps and the column names in att_data that can be interpolated
                nans = pd.DataFrame(np.nan, index=t, columns=interp_cols)
        
                # Now, interpolate!
                output = pd.concat([att_data[interp_cols], nans]) # concatenate nan matrix and att data
                output.sort_index(inplace=True) # sort values by index
                output.interpolate(inplace=True, method='time') # interpolate
                output = output.loc[t, :].drop_duplicates() # drop everything but interpolated values and remove duplicates
                output = pd.concat([output, saa_att_df, counts], axis=1) # concatenate the dataframes
    
                output = output.round(4) # Round everything down to 3 decimal places
                output[['SAA', 'att']] = output[['SAA', 'att']].fillna(-1).astype(int) # Replace nans so I can round to integer
                output[['SAA', 'att']] = output[['SAA', 'att']].round().astype('int') # Round quality flags to int
                output[['SAA', 'att']] = output[['SAA', 'att']].replace(-1, np.nan) # Put nans back
        
                if not debug: # if saving to a file
                    if len(output) > 0: # if any bursts detected:
                        now = datetime.now().strftime("%B %d, %Y %H:%M:%S") # current timestamp
                        print(output) # Print the bursts and their data
                        
                        # save to file
                        if op.isfile('./' + filename):
                            # if file already exists, append
                            output.to_csv(filename, mode='a', sep=',', header=False, encoding='utf-8')
                            print("Added to file at " + now + ".")
                        else:
                            # if it does not already exist, create it
                            output.to_csv(filename, sep=',', encoding='utf-8')
                            print("Saved to file at " + now + ".")
                    else:
                        print("No inner belt bursts found on this day.")

    time.sleep(5) # to prevent overheating

print("Queried dates complete.")
