# The "main" script, if you will.
# Prompts the user for dates, with a number of options (mostly so my life is easier.)
# Runs O'Brien, then bouncingPackets. Plots the results and saves to a dataframe, which is then
# saved to a file.
# Date created: 01/04/2024
# Last modified: 05/23/2024
# Author: Max Feinland for Blum Research Group, LASP

# Housekeeping
import os.path
import sys
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import sampex
module_directory = 'C:/Users/maxim/IRBEM-main/python/IRBEM'
sys.path.append(module_directory)
from bouncingpacketfunctions import *

print("Welcome to the ChorusWaves search script.", flush=True)
print("Coded by Max Feinland for Blum Research Group, 2023-2024.", flush=True)

## This is all for user input.
print_yn = 1
hilt_already_loaded = None

if ('first_day_d' in globals()) & ('last_day_d' in globals()):
    fdo = first_day_d
    ldo = last_day_d
    f_fd = fdo.strftime("%B %d, %Y")
    f_ld = ldo.strftime("%B %d, %Y")
    if f_fd == f_ld:
        print(f"The most recently queried date was {f_fd}.", flush=True)
    else:
        print(f"The most recently queried dates were {f_fd} to {f_ld}.", flush=True)
        
    print("* Press o to set the first date as 1 day after the last date of the last query.", 
          flush=True)
    print("* Press r to rerun the most recently queried dates.", flush=True)
    first_day = input("Otherwise, enter a start date formatted as \"YYYY, MM, DD\": ")
else:
    first_day = input("Enter a start date formatted as \"YYYY, MM, DD\": " )

if first_day == "r":
    try:
        first_day_d = fdo
        last_day_d = ldo
        if fdo == ldo:
            hilt_already_loaded = True # if you keep running the same day over and over, like I do, 
            # this saves time
    except:
        print("Looks like you cleared your old dates; can't use that command, sorry.")
elif first_day == "o":
    first_day_d = ldo + timedelta(days=1)
    last_day = input("Enter an end date, or enter the number of days as an integer (from 1 to 30): ")
    try:
        last_day = int(last_day)
        last_day_d = first_day_d + timedelta(days = last_day - 1)
    except:
        try:
            last_day = last_day.split(", ")
            last_day_d = datetime(int(last_day[0]), int(last_day[1]), int(last_day[2]))
        except:
            print("Unrecognized command; could not execute.")
else:
    first_day = first_day.split(", ")
    first_day_d = datetime(int(first_day[0]), int(first_day[1]), int(first_day[2]))
    last_day = input("Enter an end date, or enter the number of days as an integer (from 1 to 30): ")
    try:
        last_day = int(last_day)
        if last_day != 0:
            last_day_d = first_day_d + timedelta(days = last_day - 1)
        else:
            print("Zero days is not allowed for obvious reasons.")
    except:
        try:
            last_day = last_day.split(", ")
            last_day_d = datetime(int(last_day[0]), int(last_day[1]), int(last_day[2]))
        except:
            print("Unrecognized command; could not execute.")
date_list = [first_day_d + timedelta(days=i) for i in range((last_day_d - first_day_d).days + 1)]


# Now we are getting to the the function calls, etc

for day in date_list:
    formatted_date = day.strftime("%B %d, %Y")

    if print_yn == 1:
        print(f"\nCollecting data for {formatted_date}...")

    if hilt_already_loaded == True:
        pass
    else:
        h = sampex.HILT(day) # count rate data (thanks Mike!)
        h.load()

    if print_yn == 1:
        print("Running O'Brien algorithm...")

    so, eo, ns = obrien(h) # call O'Brien function

    if print_yn == 1:
        print("Searching for bouncing packets...")

    [stt, ett, pksit,  data] = bouncingPackets(so, eo, h)
    filename = "outerbeltevents.csv"


    if print_yn == 1:
        if len(stt) == 0:
            print("No intervals found.")
        elif len(stt) == 1:
            print("1 interval found. Plotting now.")
        else:
            print(len(stt), "intervals were found. Plotting now.")

    for j in range(len(stt)):
        fig = plt.figure(j)
        ax = fig.add_subplot()
        plt.grid(True)

        ax.plot(h.times[stt[j]-50:ett[j]+75], h.counts[stt[j]-50:ett[j]+75])
        ax.plot(h.times[pksit[j]], h.counts[pksit[j]], marker='d')
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S.%f')) 
        ax.set_xlabel('Time (UTC)')
        ax.set_ylabel('Counts (#/20 ms)')
        ax.set_title('Bouncing Packets on ' + h.times[stt[j]].strftime('%Y-%m-%d'))
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(colors='white', which='both')
        ax.title.set_color('white')
        ax.title.set_color('white')
        plt.show()
        if len(stt) == 1:
            plt_name = h.times[stt[j]].strftime('%Y%m%d') + '.jpg'
        else:
            plt_name = h.times[stt[j]].strftime('%Y%m%d') + "_" + str(j+1) + '.jpg'
        fullplt_name = 'pics of bps/' + plt_name
        fig.savefig(fullplt_name)
    if os.path.isfile('./' + filename):
        data.to_csv(filename, mode='a', sep=',', index=False, header=False, encoding='utf-8')
    else:
        data.to_csv(filename, sep=',', index=False, encoding='utf-8')


print('Script complete.')