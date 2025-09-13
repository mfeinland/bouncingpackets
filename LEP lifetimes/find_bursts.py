# Runs O'Brien & saves to file.
# Author: Max Feinland for Blum Research Group, LASP
# Last modified: 4/30/25

import os.path as op
import sampex
from datetime import datetime
import time

from find_microburst_functions import *

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
