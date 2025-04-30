## Find bouncing packet microbursts from inner belt microburst catalog
# Author: Max Feinland for Blum Research Group, LASP
# Last modified: 4/30/25

import sampex
import pandas as pd 
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os.path as op

# Read in data
data = pd.read_csv('Data_Files/inner_belt_microbursts.csv', index_col=0) # file with all microburst data
data.index = pd.to_datetime(data.index, format='mixed')

# Restrict to high count rates and inner belt
d = data[(data['counts'] > 100) & (data['counts'] < 10000) & (data['L'] < 2.5)]

t = d.index
dt = t.diff() # Difference in time between two microbursts

filename = 'Data_Files/bp_candidates.csv' # file to save candidate events to

# Find times spaced close together
bp_candidate = np.where(dt <= timedelta(seconds=3))[0]
bp_candidate_copy = bp_candidate
bp_candidate = bp_candidate - 1 # roll back indices by 1
bp = np.unique(np.sort(np.concatenate((bp_candidate, bp_candidate_copy))))

# # Now bp contains indices of potential bouncing packets
d_bp = np.where(np.diff(bp) > 1)[0]

iterator = range(len(t))

if op.isfile(filename):
    # Find location of last date
    with open(filename, "r", encoding="utf-8") as file:
        last_line = file.readlines()[-1] # Open last line of filename
    # The first chunk of that line has the last date ran
    last_time = datetime.strptime(last_line[:last_line.index(",")], "%Y-%m-%d %H:%M:%S.%f")
    print(last_time)

    # Find last microburst ID
    comma_idx = [i for i, ltr in enumerate(last_line) if ltr == ","]
    last_id_idx = comma_idx[-1] + 1
    last_id = int(last_line[last_id_idx:-1])
    
    # Find the index in the list_of_days corresponding to this day
    start_idx = np.where(t > last_time)[0][0]
    iterator = iterator[start_idx:]
    current_event_id = last_id+1 # Start out with event 1

else:
    current_event_id = 1 # start with 1 otherwise
    
event_id = []
timestamp = []
idx = []

debug = False

for i in iterator:
    current_t = t[i]
    date_str = date_str = current_t.strftime('%B %d, %Y')

    # Decide whether to load HILT data
    if not debug:
        if i == iterator[0]: # if first date, load for sure
            print(f"Loading HILT data for {date_str}...")
            h = sampex.HILT(current_t) # count rate data
            h.load()
        else: # otherwise, check with last date
            last_t = t[i-1]
            if current_t.date() != last_t.date():
                print(f"Loading HILT data for {date_str}...")
                h = sampex.HILT(current_t) # count rate data
                h.load()


    # if this is the start of an event:
    if dt[i] >= timedelta(seconds=3):
        current_event_id += 1 # start of new event

    event_id.append(current_event_id)
    timestamp.append(current_t)
    idx.append(i)

    # Check if this is the last timestamp in the event
    if dt[i+1] > timedelta(seconds=3):
        # Locate all the timestamps in the event
        current_event_indices = [i for i, x in enumerate(event_id) if x == current_event_id]
        
        start_time = timestamp[current_event_indices[0]] - timedelta(seconds=1)
        end_time = current_t + timedelta(seconds=2)
        tidx = []
    
        try: # Try to find the indices of exact start & end times
            sidx = np.where(h['time'] == start_time)[0][0]
            eidx = np.where(h['time'] == end_time)[0][0]
            for j in range(len(current_event_indices)):
                tidx.append(np.where(h['time'] == timestamp[current_event_indices[j]])[0][0])
        except: # In case there is no timestamp at that value, find the closest values
            sidx = np.abs((h['time'] - start_time).values).argmin()
            eidx = np.abs((h['time'] - end_time).values).argmin()
            tidx = [] # clear and start over
            for j in range(len(current_event_indices)):
                tidx.append(np.abs((h['time'] - timestamp[current_event_indices[j]]).values).argmin())
        
        # Plot 
        fig, ax = plt.subplots()
        ax.plot(h['time'][sidx:eidx], h['counts'][sidx:eidx])
        ax.plot(h['time'][tidx], h['counts'][tidx], 'r*')
        ax.text(h['time'][eidx-50], min(h['counts'][sidx:eidx]) + 0.8*np.ptp(h['counts'][sidx:eidx]), f"L = {d['L'][i]}")
        ax.xaxis.set_major_locator(mdates.SecondLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S')) 
        ax.set_xlabel(f"Time on {date_str}")
        ax.set_ylabel("Counts (#/20ms)")
        plt.show()

        satisfactory_answer = False # initialize
        # Ask the user if they think this looks like a bouncing packet
        while not satisfactory_answer:
            check_yn = input("Does this look like a bouncing packet candidate? (y/n) ")
            if check_yn.lower() == "y":
                idxx = [idx[x] for x in current_event_indices]
                output = d.iloc[idxx]
                output['event_id'] = [event_id[x] for x in current_event_indices]
                print(output)

                now = datetime.now().strftime("%B %d, %Y %H:%M:%S")
                # save to file
                if op.isfile('./' + filename):
                    # if file already exists, append
                    output.to_csv(filename, mode='a', sep=',', header=False, encoding='utf-8')
                    print("Added to file at " + now + ".\n")
                else:
                    # if it does not already exist, create it
                    output.to_csv(filename, sep=',', encoding='utf-8')
                    print("Saved to file at " + now + ".\n")

                t_str = h['time'][tidx[0]].strftime('%Y%m%d_%H%M%S')
                picture_name = 'Figures/bp_candidates/' + t_str + '.jpg'
                fig.savefig(picture_name)
                satisfactory_answer = True
                
            elif check_yn.lower() == "n":
                satisfactory_answer = True
