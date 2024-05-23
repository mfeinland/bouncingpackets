def obrien(data):
    ## O'Brien

    # Date created: 3/21/23 (in MATLAB)
    # or 01/04/24 in Python
    # Last modified: 05/23/24 to make this runnable as a .py
    # Author: Max Feinland for Blum Research Group, LASP

    # Inputs: sampex package HILT object

    # Outputs: starting & ending time of each microburst, plus time intervals

    # Housekeeping
    
    import sampex
    import pandas as pd
    import numpy as np

    N20 = data['counts'] # count rate sampled every 20 ms
    time20 = data['time'] # time every 20 ms

    # I'm sure there's a more efficient way to do this, but I don't know it
    df = pd.DataFrame({'time': data['time'], 'counts': data['counts']})

    df.set_index('time', inplace=True) # set time column as the index

    # resample the dataframe to 100 ms intervals and sum the counts in each interval
    N100 = df.resample('100ms').sum()

    A500 = N100.rolling(5, center=True).mean() # 5-observation centered rolling mean (over 500 ms)

    condition = np.divide((N100.counts - A500.counts), np.sqrt(1 + A500.counts)) # O'Brien et al 2003
    
    ns = np.argwhere(condition > 10)
    ns = [item[0] for item in ns]

    epsilon = 10; # if two flagged indices are spaced less than this distance apart, 
    # they are probably part of the same microburst

    # initializing
    starts = []
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
    starts20 = changeCadence(time20, N100.index[starts])
    ends20 = changeCadence(time20, N100.index[ends])
    ns20 = changeCadence(time20, N100.index[ns])

    # yay, your output is ready!
    d = dict(); 
    d['st'] = starts20
    d['et'] = ends20
    d['ns'] = ns20
    
    return starts20, ends20, ns20

def tsyganenko05(tstmp):
    # Author: Max Feinland, although Mike Shumko gets a lot of credit for the skeleton
    # Purpose: Calculates bounce period from T05 model.
    
    # If all of the preliminary conditions checked are true, the T05 model
    # will be called to check if the bounce period is within tolerance.
    
    # Also pulls attitude data. This saves some time, since you won't need to download
    # the attitude data otherwise.
    
    # Housekeeping
    import numpy as np
    import spacepy.datamodel
    from datetime import datetime
    import dateutil.parser
    from IRBEM import MagFields
    import sampex

    a = sampex.Attitude(tstmp) # attitude data (thanks Mike!)
    a.load()
    
    # change this for your needs. Of course, this pathing is for my machine.
    omniLoc = 'C:/Users/maxim/.spacepy/data/omnidata.h5'
    omniData = spacepy.datamodel.fromHDF5(omniLoc)
    omniT = np.array([dateutil.parser.parse(i.decode()) for i in omniData['UTC']])
    idx = np.where(tstmp >= omniT)[0][-1] # find relevant omni index
    T05Keys = ['Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    maginput05 = {}
    maginput = {}

    X = {}
    # find relevant attitude index
    a_idx = np.where(tstmp >= a['time'])[0][-1]
    
    # X is the time & location of the particle (lla)
    X['dateTime'] = tstmp.strftime('%Y-%m-%d %H:%M:%S.%f')
    X['x1'] = a['Altitude'][a_idx]
    X['x2'] = a['GEO_Lat'][a_idx]
    X['x3'] = a['GEO_Long'][a_idx]
    # These are returned for data processing reasons (later on). Not necessary for the model
    L = a['L_Shell'][a_idx]
    MLT = a['MLT'][a_idx]
    
    KE = 1000 # keV, true of all particles observed by SAMPEX

    #T05 model
    model05 = MagFields(options = [0,0,0,0,0], kext = 11)
    for i in T05Keys:
        maginput05[i] = float(omniData[i][idx])
    Tb = model05.bounce_period(X, maginput05, KE) # call model
    return Tb, X, L, MLT

def bouncingPackets(so, eo, h):
    # Date created: not sure. Sort of a Ship of Theseus thing. I originally created this in MATLAB
    # 04/25/23, but have no idea when it was created in Python, unfortunately. Likely Jan 2024.
    # Last modified: 4/30/24
    # Purpose: takes intervals specified from the O'Brien algorithm and searches for consistently-spaced,
    # low-width-variance, decreasing, isolated peaks whose bounce periods 
    # match up with the expected values.
    
    # Housekeeping
    import pandas as pd
    import numpy as np
    from scipy.signal import find_peaks, peak_widths
    import time
    import stopit
    
    # Intializing variables
    pksi = []
    st = []
    et = []
    per = []
    L = []
    X = []
    lat = []
    lon = []
    MLT = []
    TbT05 = []
    tb_check = []
    
    # load in HILT data
    t = h.times
    rate = h.counts
    
    # the "with" statement prevents excessive runtimes, at the expense of a couple of days' 
    # worth of data
    with stopit.ThreadingTimeout(60) as context_manager:
        for i in range(len(so)):

            interval = rate[so[i]:eo[i]] # define just one chunk of the rate data, taken from O'Brien
            maxdis = max(interval) - min(interval) # use to generate prominence requirement

            # finding peaks with a high enough prominence
            [pks, _] = find_peaks(interval, prominence=0.25*maxdis, distance=3)
            npks = list(so[i] + loc for loc in pks)
            pksi.extend(npks) # add to "pksi" list, to access later (plotting)

            loc = pd.to_datetime(t[npks])

            dloc = [loc[i + 1] - loc[i] for i in range(len(loc) - 1)] # peak spacing
            dloc = [x.total_seconds() for x in dloc]
            ddloc = np.diff(dloc) # change in peak spacing; must be sufficiently small, i.e.
            # the peaks must be very regularly spaced.

            if len(ddloc) >= 2:
                threshold = 0.05 # lowest allowable change in peak spacing/period

                # find a run of at least 2 indices that meet this criterion. That means 4 peaks
                indices = np.where(np.convolve(np.abs(ddloc) < threshold, \
                                               np.ones(2), mode='valid') == 2)[0]
                if len(indices) == 1: # if there is exactly one microburst
                    tstmp_start = loc[indices[0]] - pd.Timedelta(seconds=0.2)
                    tstmp_final = loc[indices[0]+3] + pd.Timedelta(seconds=0.2)
                    et.extend(np.where(t == tstmp_final)[0])
                    st.extend(np.where(t == tstmp_start)[0])
                elif len(indices) > 1:
                    for j in range(len(indices)-1):
                        # if you are looking at the first index and there is not a jump, that is a start time
                        if j == 0 and indices[j + 1] - indices[j] == 1:
                            tstmp_start = loc[indices[j]] - pd.Timedelta(seconds=0.2)
                            st.extend(np.where(t == tstmp_start)[0])
                        elif j == 0: # otherwise, nothing happens!
                            pass
                        # if previously, you were not in a consecutive streak, 
                        # but now you are, that is a start time
                        elif indices[j+1] - indices[j] == 1 and indices[j] - indices[j-1] > 1:
                            tstmp_start = loc[indices[j]] - pd.Timedelta(seconds=0.2)
                            st.extend(np.where(t == tstmp_start)[0])
                        # if previously, you were in a consecutive streak, but now you are not,
                        # that is an endtime
                        elif indices[j+1] - indices[j] > 1 and indices[j] - indices[j-1] == 1:
                            tstmp_final = loc[indices[j]+3] + pd.Timedelta(seconds=0.2)
                            et.extend(np.where(t == tstmp_final)[0])

                    # end condition if the previous one wasn't met
                    if len(st) > len(et):
                        tstmp_final = loc[indices[j]+3] + pd.Timedelta(seconds=0.2)
                        et.extend(np.where(t == tstmp_final)[0])

    # The variables st and et contain all intervals that satisfy the strict bounce period condition.
    # Now we will check the other conditions.
        
    # Check if timed out
    if context_manager.state == context_manager.EXECUTED:
        st = np.sort(st)
        et = np.sort(et)
    
    elif context_manager.state == context_manager.TIMED_OUT:
        st = []
        et = []
        print("Timed out execution. Skipping this day.")

    # initialize things for DataFrame
    final_st = []
    final_et = []
    finalpksi = []
    dec = []
    MLT = []
    L = []
    TbT05 = []
    lat = []
    lon = []

    for k in range(len(st)):
        # Find all of the indices in the rate data corresponding to this interval.
        ok_range_for_this_index = np.arange(st[k], et[k]+1).tolist()
        # Find all the values in pksi that are contained within these indices.
        these_pk_indices = [index for index, value in enumerate(pksi) if value in ok_range_for_this_index]
        relative_pks = [pksi[x]-st[k] for x in these_pk_indices]
        current_per = np.mean(np.diff(relative_pks))*0.02 # mean period
        
        interval = rate[st[k]:et[k]]
        unq = len(np.unique(interval)) > 5
        if unq == True:
            try:
                widths = peak_widths(interval, relative_pks, rel_height=0.5)
                incrw = np.all(np.diff(widths[0]) >= -3)
            except:
                # sometimes peak_widths errors out :(
                incrw = False

            # Needed to make the condition "isol"
            # Find the index of the time that is 3 seconds before the identified starttime
            desired_starttime = t[st[k]] - pd.Timedelta(seconds=3)
            idxs = np.abs(np.array(t, dtype='datetime64') - np.datetime64(desired_starttime)).argmin()
            interval2_start = idxs

            # Find the index of the time that is 3 seconds after the identified endtime
            desired_endtime = t[et[k]] + pd.Timedelta(seconds=3)
            idxe = np.abs(np.array(t, dtype='datetime64') - np.datetime64(desired_endtime)).argmin()
            interval2_end = idxe
            interval2 = rate[interval2_start:interval2_end]

            interval = rate[st[k]:et[k]]
            maxdis = max(interval) - min(interval)
            
            # Find all peaks within the new, longer interval that satisfy the prominence requirement
            [loc2, _] = find_peaks(interval2, prominence=0.25*maxdis, distance=3)
            intermediate_list = []
            intermediate_list = (interval2_start + loc2).tolist()
            all_pks_for_interval = [pksi[x] for x in these_pk_indices]
            for item in all_pks_for_interval:
                if item in intermediate_list:
                    intermediate_list.remove(item)
            isol = len(intermediate_list) < 12
                
            # this is for the decreasing interpolation
            timestamps = t[st[k]:et[k]]
            timestamps = np.array(timestamps, dtype='datetime64')
            timedelta_index = pd.to_timedelta(timestamps - timestamps.astype('datetime64[D]'))
            intt = timedelta_index.total_seconds()

            peaks, props = find_peaks(-interval, prominence=0.15 * np.ptp(interval))
            peaks = np.concatenate(([0], peaks, [len(interval) - 1]))

            vq = np.interp(intt, [intt[x] for x in peaks], [interval[x] for x in peaks])
            adj = interval - vq + np.min(interval)
            
            # sometimes there's a funny edge condition that outputs the out-of-bound indices as peaks.
            # so this "if" statement is just taking care of that
            if len(adj) in relative_pks:
                relative_pks.remove(len(adj))
                
            max_prom = max(props["prominences"])
            # all peaks, except the first one, have to decrease (or not increase by a whole lot)
            decr = np.all(np.diff(adj[relative_pks[1:]]) < 0.2*max_prom)

            passesMostConditions = isol and unq and incrw and decr
        else:
            passesMostConditions = False


        # This is split up into two cases because the TbT05 calculation takes FOREVER to run.
        # So, the first three conditions must be true before the last one will be checked.
        if passesMostConditions == True:
            print('Calling Tsyganenko 2005 model. This may take a while; please be patient.')
            try:
                [Tb, current_X, current_L, current_MLT] = tsyganenko05(t[st[k]])
            except ValueError:
                Tb = None
                print("There was a problem calling the model, and this interval was not used.")

            if Tb is not None:
                per.append(current_per)
                L.append(current_L)
                MLT.append(current_MLT)
                TbT05.append(Tb)
                final_st.append(st[k])
                final_et.append(et[k])
                finalpksi.append(all_pks_for_interval)
                lat.append(current_X['x2'])
                lon.append(current_X['x3'])
            # dec.append(decr)
    
    data = pd.DataFrame({'t': [h.times[x] for x in final_st], 'per': per, 'tb': TbT05, 'L': L, 
                        'MLT': MLT, 'lat': lat, 'lon': lon})
    # so the dataframe isn't long af
    data = data.round({'per': 4, 'tb': 4, 'L': 3, 'MLT': 3, 'lat': 3, 'lon': 3})
    return final_st, final_et, finalpksi, data