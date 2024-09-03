%% Bouncing Packet Locator Script - Low L-shells Only

% Date Created: 5/15/23 (Original Script)
% Date Created: 6/09/23 (This Script)
% Last Modified: 9/2/24

% Author: Max Feinland for Blum Research Group, LASP

% Inputs: start date, end date, some other user prompts

% Description: complete driver script. The user will be prompted to specify
% the dates they would like to see. This script will then print out
% and plot all of the microbursts found in the specified time period.

% Housekeeping
clc
close all
clear

% desired epoch to look at (user defined)
first_day = input('Enter a start date (formatted as YYYY-MM-DD): ', 's');
last_day = input('Enter an end date (formatted as YYYY-MM-DD): ', 's');

% call plotEpoch with user inputs
[num_bursts, dates, numeric_data] = plotEpoch(first_day, last_day);

if ~isempty(dates)

    fprintf(['Date \t Timestamp \t P (s)\t Lat \t Long\t Alt\t' ...
        'L \t MLT \t Magnetic Lat\n'])

    dates.Format = 'yyyy-MM-dd HH:mm:ss';

    for j = 1:length(dates)
        fprintf('%s %.3f %.4f %.4f %.4f %.4f %.4f %.4f\n', dates(j), numeric_data(j,1), ...
            numeric_data(j,3), numeric_data(j,4), numeric_data(j,5), ...
            numeric_data(j,6), numeric_data(j,7), numeric_data(j,8))
    end

else
    fprintf('No microbursts found.\n')
end

fprintf('Script completed.\n')

%% Functions Needed

%% PlotEpoch
function [num_bursts, out_dts, out_dats] = plotEpoch(first_day, last_day)

% Date Created: 5/15/23
% Last Modified: 6/4/23
% Author: Max Feinland for Blum Research Group, LASP

% Inputs: first day of epoch, last day of epoch, L-shell data (time and 
% value), minimum allowable period, maximum allowable period

% Outputs: plots of each microburst in epoch (by day), total number of
% microbursts, prints the period, std, L-shell, and percent error

% this helps the formats to line up, makes it easier to find datetimes
datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')

dates_ok = 1; % assume by default that these dates are valid

try
    first_day_d = datetime(first_day); % start of epoch, datetime format
    last_day_d = datetime(last_day); % end of epoch, datetime format
    if first_day_d > last_day_d % if nonsequential dates
            fprintf(['Error: non-sequential dates, i.e. ' ...
                'the first date is after the last date. Please try a new input.\n'])
            % do not plot or print anything
        dates_ok = 0; % dates not valid
        num_bursts = 0; % no bursts
        out_dts = [];
        out_dats = [];
    end
catch % nonreal dates
    fprintf(['Error: non-real date. This can happen if a specified date does not exist' ...
        ', e.g. February 31st, \nor if there was a typo during input. Please try a new input.\n'])
    % do not print or plot anything
    dates_ok = 0; % dates not valid
    num_bursts = 0; % no bursts
    out_dts = [];
    out_dats = [];
end

if dates_ok == 1 % if dates exist and are sequential

doy1 = day(first_day_d, 'dayofyear'); % day of start (integer)
doy2 = day(last_day_d, 'dayofyear'); % day of end (integer)
y = year(first_day_d); % year in question
y2 = year(last_day_d); % year of end date (important only if y != y2)

% initializing
num_bursts = 0; % total number of bursts found
multBursts = 0; % multiple bursts (if necessary)

if y == 1996 && doy1 < 220 || y < 1996 || y > 2013 || y == 2013 && doy1 > 312
    % if the data is outside the time range that the data spans
    disp(['Invalid command; please choose dates between August 7th, 1996 ' ...
        'and November 8th, 2013.'])
    num_bursts = 0; out_dts = []; out_dats = [];
    return;
elseif y ~= y2
    disp('Invalid command; please choose two dates in the same year.')
    num_bursts = 0; out_dts = []; out_dats = [];
    return; 
elseif y == 2011 && sum(ismember([doy1 doy2], 297:323)) > 0 || ...
       y == 2005 && sum(ismember([doy1 doy2], 1:9)) > 0
     disp('Request contains data outside ephemeris range. Please try again.')
     num_bursts = 0; out_dts = []; out_dats = [];
     return;
else
    fprintf('Working...\n') % status update for user

    doyvec = doy1:doy2; % array of each day included
    datetimevec = first_day_d + days(doyvec - doyvec(1)); % transform to datetime
    
    % this is for pulling ephemeris data
    dayyearvec = dayarrayfun; % contains day-year pairs

    % search for row that contains correct year and acceptable days
    desired_row = find(dayyearvec(:,2) == year(first_day_d) & dayyearvec(:,1) > doy1);

    if isempty(desired_row)
        % has catch condition in case the day of the year is too late
        desired_row = find(dayyearvec(:,2) == year(first_day_d) & dayyearvec(:,1) <= doy1);
        n = desired_row(end)+1;
    else
        n = desired_row(1);
    end

    % search for row that contains correct year and acceptable days
    desired_row_2 = find(dayyearvec(:,2) == year(last_day_d) & dayyearvec(:,1) > doy2);

    if isempty(desired_row_2)
        % has catch condition in case the day of the year is too late
        desired_row_2 = find(dayyearvec(:,2) == year(last_day_d) & dayyearvec(:,1) <= doy2);

        n2 = desired_row_2(end)+1; % needed row of dayyearvec
    else
        n2 = desired_row_2(1); % needed row of dayyearvec
    end


    % initializing
    l_data = []; 
    ephem_filename = cell(size(n:n2));

    for i = n:n2
        day1var = num2str(dayyearvec(i-1,1));
        day2var = num2str(dayyearvec(i,1)-1);

        % checking if leading zeroes are needed
        if strlength(day1var) == 1
            day1var = ['00' day1var]; %#ok<AGROW> 
        elseif strlength(day1var) == 2
            day1var = ['0' day1var]; %#ok<AGROW> 
        end

        if strlength(day2var) == 1
            day2var = ['00' day2var]; %#ok<AGROW> 
        elseif strlength(day2var) == 2
            day2var = ['0' day2var]; %#ok<AGROW> 
        end

        % status update
        fprintf('Retrieving ephemeris data... (%i/%i)\n', i-n+1, length(n:n2))

        % filename(s) needed
        ephem_filename{i-n+1} = ['PSSet_6sec_' num2str(dayyearvec(i-1,2)) day1var '_' ...
            num2str(dayyearvec(i,2)) day2var '.txt'];

        if i >= 51 % portion for which data is not zipped
            url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/PSSet/Text/', ...
            ephem_filename{i-n+1}]; % url to pull

            file = webread(url);
            start = strfind(file,"BEGIN DATA") + 11; % start of data
            file = file(start:end);

            A = '%*f '; % don't include columns
            B = '%f '; % do include columns
            C = repmat(A, [1 38]);
            D = [B B B A A A A B B B repmat(A, [1 10]) B B B repmat(A, [1 10]) B C];
            % tells textscan to look at columns 1, 2, 3, 8, 9, 10, 21, 22,
            % 23, and 34
            D(end) = [];

            celltab = textscan(file, D, 'Delimiter', '\t');
            data_chunk = [celltab{1}, celltab{2}, celltab{3}, celltab{4} ...
                celltab{5} celltab{6} celltab{7} celltab{8} celltab{9}...
                celltab{10}];

            l_data = [l_data; data_chunk]; %#ok<AGROW> 

        else % data for which it's zipped
            url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/PSSet/Text/', ...
            ephem_filename{i-n+1} '.zip'];
            
            unzip(url) % find and unzip data file
            data_chunk = readmatrix(ephem_filename{i-n+1}, 'HeaderLines', 60); % read in data file
            data_chunk = [data_chunk(:,1:3) data_chunk(:,8:10) data_chunk(:,21:23) ...
                 data_chunk(:,34)];

            l_data = [l_data; data_chunk]; %#ok<AGROW> 
        end
    end

    % status update
    fprintf('Ephemeris data retrieved.\n')

    % extracting relevant data
    secondsvec = l_data(:,3);
    geo_long = l_data(:,4);
    geo_lat = l_data(:,5);
    geo_alt = l_data(:,6);
    lshell_all = l_data(:,7);
    bfield_all = l_data(:,8);
    mlt_all = l_data(:,9);
    maglat_all = l_data(:,10);

    % status update
    fprintf('Searching for bouncing packets...\n')

    for i = 1:length(doyvec)

        existsDate = 1; % assume by default that the date exists in the repository

        % makes sure filename has leading zeros if necessary
        if strlength(num2str(doyvec(i))) == 3
            filename = ['hhrr', num2str(y), num2str(doyvec(i)), '.txt'];
        elseif strlength(num2str(doyvec(i))) == 2
            filename = ['hhrr', num2str(y), '0', num2str(doyvec(i)), '.txt'];
        elseif strlength(num2str(doyvec(i))) == 1
            filename = ['hhrr', num2str(y), '00', num2str(doyvec(i)), '.txt'];
        end

        % pulls file from online SAMPEX data repository
        if doyvec(i) <= 182 && y == 2004 || y < 2004 % data for which it's zipped
            try
                url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
                    filename, '.zip'];
                unzip(url) % find and unzip data file
                data = readmatrix(filename); % read in data file
            catch
                existsDate = 0; % does not exist in the repository
            end
        else % data for which it's not zipped
            try
                url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
                    filename];
                file = webread(url);
                
                % parse the character vector into a table (code adapted from CGPT)
                celltab = textscan(file, '%f %f %f %f %f %f %f', 'Delimiter', '\t', 'HeaderLines', 1);
                
                % converts the parsed data to a table
                data = [celltab{1}, celltab{2}, celltab{3}, celltab{4}, celltab{5}, celltab{6}, celltab{7}];
            catch
                existsDate = 0; % does not exist in the repository
            end

        end

        if existsDate == 1 % if the date does exist
            % call O'Brien 
            [~, d, r, so, eo, ~] = obrien100(data, string(datetimevec(i)));
            
            num = 0;

            % call bouncingPacket
            [st, et, pksi] = bouncingPacketSplits(so, eo, r);

            % Outputting relevant information
            fprintf('%s\n', datetimevec(i)) % header, includes date

            if ~isempty(st)
                % finding closest time
                time_elapsed = seconds(timeofday(d(round(st)))); % seconds since midnight
                allowable_rows = find(l_data(:,2) == doyvec(i)); % days in the list that are this day's date
                dist = pdist2(secondsvec(allowable_rows), time_elapsed, 'euclidean'); % closest second to actual

                % initializing
                yesPlot = zeros(size(st));
                used_lshells = [];
                used_dts = NaT(size(st));
                used_mdl = []; 
                used_mlt = [];
                used_maglats = [];
                used_p_expect = [];
                used_geocoords = [];
                used_bfield = [];

                for j = 1:length(st)

                    [~, n_time] = min(dist, [], 1); % find closest time index
                    long = geo_long(allowable_rows);
                    lat = geo_lat(allowable_rows);
                    alt = geo_alt(allowable_rows);
                    lshell = lshell_all(allowable_rows);
                    bfield = bfield_all(allowable_rows);
                    mlt = mlt_all(allowable_rows);
                    maglat = maglat_all(allowable_rows);

                    ns = pksi >= st(j) & pksi <= et(j); % finds indices with peaks

                    mdl = mean(diff(pksi(ns)/50)); % mean bounce period

                    % calculate expected period
                    p_expect = expectedPeriod(lshell(n_time(j)), maglat(n_time(j)));
             
                    length_int = et(j) - st(j); % number of indices in interval
                    unqvals = unique(r(round(st(j)):round(et(j)))); % unique numbers in interval

                    % (second condition exists since sometimes the data
                    % only shows a few possible values)

                    % to be printed, all these conditions must be true:
                    % short, good resolution, and low l-shell

                    conds = [length_int < 750, length(unqvals) > 10, ...
                         lshell(n_time(j)) <= 3];
        
                    if all(conds) % if they are indeed all met
                        fprintf('[Event %i] %s %.3f (%.3f, %.3f, %.3f) %.4f %.4f %.4f\n', ...
                            num+1, d(st(j)), mdl, lat(n_time(j)), ...
                            long(n_time(j)), alt(n_time(j)), ...
                            lshell(n_time(j)), mlt(n_time(j)), maglat(n_time(j)))

                        num = num + 1; % increment number of bursts found
                        yesPlot(j) = 1; % include this interval in output
                        used_dts(j) = d(st(j));
                        used_mdl(end+1) = mdl; %#ok<AGROW>
                        used_lshells(end+1) = lshell(n_time(j));%#ok<AGROW>
                        used_mlt(end+1) = mlt(n_time(j));%#ok<AGROW>
                        used_maglats(end+1) = maglat(n_time(j));%#ok<AGROW>
                        used_p_expect(end+1) = p_expect;%#ok<AGROW>
                        used_geocoords(end+1,:) = [lat(n_time(j)) long(n_time(j)) alt(n_time(j))];%#ok<AGROW>
                        used_bfield(end+1) = bfield(n_time(j));%#ok<AGROW>

                    end
                    plotting_data = [used_mdl' used_p_expect' used_geocoords ...
                        used_lshells' used_mlt' used_maglats' used_bfield'];
                end
            end

            % plotting
            if num > 0 % plots events of interest and prints number observed
                % remove empty parts of initialized matrices
                used_dts = rmmissing(used_dts);

                % calling PlotInputs
                PlotInputs(num, r, d,  ['Bouncing Packets on ', char(datetimevec(i))], ...
                    st(yesPlot == 1), et(yesPlot == 1), pksi, plotting_data)

                if multBursts == 1 % if multiple bursts found 
                    out_dts = [out_dts; used_dts]; %#ok<AGROW> 
                    out_dats = [out_dats; plotting_data]; %#ok<AGROW> 
                elseif multBursts == 0 && num > 0 % if not multiple, but still some found
                    out_dts = used_dts;
                    out_dats = plotting_data;
                else % if none found
                    out_dts = [];
                    out_dats = [];
                end
                % printing to command window
                fprintf('Number of events observed on %s: %i\n', datetimevec(i), num)
            else
                % if no events, don't plot anything
                fprintf('Number of events observed on %s: 0\n', datetimevec(i))
                if multBursts ~= 1 % if only one burst, or none
                    out_dts = [];
                    out_dats = [];
                end

            end
    
            num_bursts = num_bursts + num; % total microbursts observed
    
            if isfile(filename) % if you had to download the .zip, delete it now
                delete(filename)
            end

        else % if existsDate == 0, i.e. the file does not exist in the repository
            fprintf('This day (%s) does not exist in the SAMPEX data repository.\n', datetimevec(i))
            if ~exist('out_dats', 'var')
                out_dats = []; % no data points to display
                out_dts = [];
            end
        end
    if num_bursts > 0 % so that data isn't deleted if mulitple bursts found
        multBursts = 1;
    end
    end
        fprintf('Total events observed in epoch: %i\n', num_bursts)
end

% deleting downloaded ephemeris files
for k = 1:length(ephem_filename)
    if isfile(ephem_filename{k})
        delete(ephem_filename{k})
    end
end

end
end

%% O'Brien

function [t50, d50, r50, starts, ends, ns] = obrien100(data, day)
% Date created: 3/21/23
% Last modified: 5/15/23
% Author: Max Feinland for Blum Research Group, LASP

% Inputs: name of file, sampling frequency, day, and the date vector from
% formatData if the sampling frequency is 50 Hz

% Outputs: starting & ending time of each microburst, plus time intervals

rates = [data(:,2:5) data(:,7)];
rate = sum(rates, 2); % sum across rows


N100 = rate;
A500 = movmean(N100, 5);   
condition = (N100 - A500)./sqrt(1 + A500); % O'Brien et al 2003
ns = find(condition > 5); % modified from 10 to 5 to find lower-count peaks

epsilon = 10; % if two flagged indices are spaced less than this distance apart,
% they are probably part of the same microburst

% initializing
starts = zeros(size(ns));
ends = zeros(size(ns));
dn = diff(ns);

% finding extended periods of the condition being true
for i = 2:length(dn)-10
    if dn(i) < epsilon && dn(i+1) < epsilon && dn(i-1) >= epsilon % start condition
        starts(i) = ns(i);
        for j = i+1:length(dn)-1
            if dn(j) < epsilon && dn(j+1) >= epsilon
                ends(j) = ns(j); % end condition
                break
            end
        end
    elseif dn(i) <= epsilon && i == 2 % start condition (edge case)
        starts(i) = ns(i);
        for j = i+1:length(dn)
            if dn(j) <= epsilon
                ends(j) = ns(j);
                break
            end
        end
    elseif i == length(dn) % end condition (edge case)
        ends(i) = ns(i);
    end
end

starts = nonzeros(starts);
ends = nonzeros(ends);

if length(starts) > length(ends)
    ends = [ends; starts(end) + 50];
end


starts = starts*5 - 10; % pad with 0.2 seconds 
ends = ends*5 + 50; % pad with 1 second
ns = ns*5;

% Returns the rate and date data
t10 = data(:,1); % time of original matrix (10Hz data)
r1 = data(:,2); % 0 to 20 ms
r2 = data(:,3); % 20 to 40 ms
r3 = data(:,4); % 40 to 60 ms
r4 = data(:,5); % 60 to 80 ms
r5 = data(:,7); % 80 to 100 ms


% 50 Hz data (needs to be created)
t50 = zeros(length(t10)*5, 1); % create/initialize 50Hz time vector 
for i = 1:length(t10)
    t50(5*i - 4) = t10(i); % every 5th index, the time is that of the 10Hz vector
end

% fill in the rest of time vector
for i = 2:length(t50)
    if t50(i) == 0 
        t50(i) = t50(i-1) + 0.02; % add 20 ms after each index
    end
end

d50 = datetime(day) + seconds(t50); % convert from seconds to UTC (nice for plotting)
d50.Format = 'MMM dd, yyyy HH:mm:ss.SSS';

r50 = zeros(length(t10)*5, 1); % create/initialize new rate vector 
for i = 1:length(t10) % fill in each row
    r50(5*i - 4) = r1(i);
    r50(5*i - 3) = r2(i);
    r50(5*i - 2) = r3(i);
    r50(5*i - 1) = r4(i);
    r50(5*i) = r5(i);
end

end

%% BouncingPacketSplits
function [st, et, pksi] = bouncingPacketSplits(st, et, rate)
% Date created: 4/25/23
% Last modified: 5/12/23

% Author: Max Feinland for Blum Research Group, LASP

% Inputs: starting times (indices), ending times(indices), rate data

% Outputs: starting & ending time of each microburst, plus peak locations

% Purpose: This checks for closely-spaced peaks that satisfy a prominence
% requirement and a period consistency requirement.

% initializing
pksi = [];
fs = 50;
st_counter = 1;
et_counter = 1;

keepvar = zeros(size(st)); % this variable tells you which starttimes to delete
st_split = zeros(size(rate));
et_split = zeros(size(rate));

for i = 1:length(st)
    interval = rate(st(i):et(i));
    maxdis = max(interval) - min(interval);
    % finding peaks with a high enough prominence

    [~, loc] = findpeaks(interval, fs, 'MinPeakProminence', 0.25*maxdis, ...
        'MinPeakDistance', 0.06);

    pksi = [pksi; st(i) + fs*loc]; %#ok<AGROW>  % location of each peak
    dloc = diff(loc); % peak spacing
    ddloc = diff(dloc); % spacing between peak spacings
    mdl = mean(dloc);

    % lowest allowable difference between two change-in-periods
    threshold = 0.15*mdl;
    
    % look for runs of three consecutive values less than the threshold
    if length(ddloc) >= 2
        indices = find(conv(double(abs(ddloc) < threshold), ones(1,2), 'valid') == 2);
        if length(indices) == 1 % if there is exactly one microburst
            et(i) = st(i) + fs*loc(indices+3) + fs/5; 
            st(i) = st(i) + fs*loc(indices) - fs/5; % zoom in on it!
            keepvar(i) = 1; % keep this segment
        elseif length(ddloc) >= 4 && length(indices) >= 1 % splits up sequences
            for j = 1:length(indices)-1
                % if you are looking at the first index and there is not a 
                % jump, that is a start time
                if j == 1 && indices(j+1) - indices(j) == 1
                    st_split(st_counter) = st(i) + fs*loc(indices(j)) - fs/5;
                    st_counter = st_counter + 1; % increment counter variable
                elseif j == 1 % otherwise, nothing happens!
                % if previously, you were not in a consecutive streak, but
                % now you are, that is a start time
                elseif indices(j+1) - indices(j) == 1 && indices(j) - indices(j-1) > 1
                    st_split(st_counter) = st(i) + fs*loc(indices(j)) - fs/5;
                    st_counter = st_counter + 1; % increment counter variable
                % if previously, you were in a consecutive streak, but now
                % you are not, that is an endtime
                elseif indices(j+1) - indices(j) > 1 && indices(j) - indices(j-1) == 1
                    et_split(et_counter) = st(i) + fs*loc(indices(j)+3)+fs/5;
                    et_counter = et_counter + 1; % increment counter variable
                end
            end
            
            % end condition if the previous one wasn't met
            if st_counter > et_counter
                et_split(et_counter) = st(i) + fs*loc(indices(j)+3)+fs/5;
                et_counter = et_counter + 1; % increment counter variable
            end
        end
    end

end

st_split = st_split(1:st_counter-1);
et_split = et_split(1:et_counter-1);

st = sort([st(keepvar == 1); st_split]);
et = sort([et(keepvar == 1); et_split]);
end


%% PlotInputs
function PlotInputs(num, r, d, title, st, et, pksi, plotting_data)

% Inputs: number of events to display, rate data, date data, title,
% starting times, ending times, peaks

% Outputs: a plot (or multiple!)

onebynum = 1:2; % if this is how many events, plot nx1
twobynum = 3:8; % if this is how many events, plot nx2
threebynum = 9:15; % if this is how many events, plot nx3
large = 16; % if more than 15, will need to plot a 15x3 and then extras

if num >= large
    plotnum53 = floor(num/15); % this is how many 5x3 plots you will need
    regular = rem(num, 15); % the remainder after division by 15
    for i = 1:plotnum53
        extra_string = [title, ' (', num2str(i), '/', num2str(plotnum53+1), ')'];
        plotFlags(r, d, 5, 3, 1+(15*(i-1)), extra_string, st, et, pksi, plotting_data, 'r*', 0)
    end
    sn = 15*plotnum53 + 1; % starting number for the remainder plot
    num = regular;
    title = [title, ' (', num2str(plotnum53+1), '/', num2str(plotnum53+1), ')'];
else
    sn = 1;
end

if ismember(num, onebynum) % plotting the nx1 plots
    plotFlags(r, d, num, 1, sn, title, st, et, pksi, plotting_data, 'r*', 0)
elseif ismember(num, twobynum) % plotting the nx2s
    plotFlags(r, d, floor(num/2), 2, sn, title, st, et, pksi, plotting_data, 'r*', rem(num, 2))
elseif ismember(num, threebynum) % plotting the nx3s
    plotFlags(r, d, floor(num/3), 3, sn, title, st, et, pksi, plotting_data, 'r*', rem(num, 3))
end

end

%% plotFlags
function plotFlags(rate, d, nr, nc, sn, title, st, et, ns, plotting_data, col, primeyn)
% Inputs: Rate data, date data, number of subplots rows, number of subplot
% columns, starting number of subplots, text, starttime, endtime, color,
% number of subplots at the end that are not necessarily aligned

% Outputs: none. Plots for you :)

mdl = plotting_data(:,1);
p_expect = plotting_data(:,2);
lat = plotting_data(:,3);
long = plotting_data(:,4);
alt = plotting_data(:,5);
L = plotting_data(:,6);
mlt = plotting_data(:,7);

figure('WindowState', 'maximized')
sgtitle(title)

if primeyn == 1 % if the number you are plotting is not divisible by 2 or 3
    for j = sn:((nr*nc)+sn-1)
        subplot(nr+1,nc,j-sn+1)
        ri = rate(st(j):et(j));
        di = d(st(j):et(j));
        nsf = find(ns > st(j) & ns < et(j));
        plot(di, ri, 'k', 'LineWidth', 1)
        hold on
        plot(d(ns(nsf)), rate(ns(nsf)), col)
        xlabel('Time (UTC)'), ylabel('Count')
        xticks(linspace(di(1),di(end),3))
        text(min(xlim), min(ylim), num2str(j), 'BackgroundColor', 'white')
        text(min(xlim) + seconds(0.05), max(ylim) - range(ri)*0.1, ...
            ['P = ', num2str(mdl(j), 3) ' ('...
            num2str(p_expect(j), 3) '), L = ', num2str(L(j)), ', MLT = ' num2str(mlt(j))], 'FontSize', 10)
        text(max(xlim) - seconds(0.6), max(ylim) - range(ri)*0.1, ...
            ['(' num2str(lat(j)) ', ' num2str(long(j)) ', ', num2str(alt(j)) ')'], 'FontSize', 10);
        set(gca,'tag',num2str(j))
    end
    for j = (nr*nc)+sn % add extra at the bottom
        subplot(nr+1, nc, nr*nc + 1 + 0.5*(nc-1))
        ri = rate(st(j):et(j));
        di = d(st(j):et(j));
        nsf = find(ns > st(j) & ns < et(j));
        plot(di, ri, 'k', 'LineWidth', 1)
        hold on
        plot(d(ns(nsf)), rate(ns(nsf)), col)
        xlabel('Time (UTC)'), ylabel('Count')
        xticks(linspace(di(1),di(end),3))
        text(min(xlim), min(ylim), num2str(j), 'BackgroundColor', 'white')
        text(min(xlim) + seconds(0.05), max(ylim) - range(ri)*0.1, ...
            ['P = ', num2str(mdl(j), 3) ' ('...
            num2str(p_expect(j), 3) '), L = ', num2str(L(j)), ', MLT = ' num2str(mlt(j))], 'FontSize', 8)
        text(max(xlim) - seconds(0.6), max(ylim) - range(ri)*0.1, ...
            ['(' num2str(lat(j)) ', ' num2str(long(j)) ', ', num2str(alt(j)) ')'], 'FontSize', 8);
        set(gca,'tag',num2str(j))
    end

elseif primeyn == 2
    for j = sn:((nr*nc)+sn-1)
        subplot(nr+1,nc,j-sn+1)
        ri = rate(st(j):et(j));
        di = d(st(j):et(j));
        nsf = find(ns > st(j) & ns < et(j));
        plot(di, ri, 'k', 'LineWidth', 1)
        hold on
        plot(d(ns(nsf)), rate(ns(nsf)), col)
        xlabel('Time (UTC)'), ylabel('Count')
        xticks(linspace(di(1),di(end),3))
        text(min(xlim), min(ylim), num2str(j), 'BackgroundColor', 'white')
        text(min(xlim) + seconds(0.05), max(ylim) - range(ri)*0.1, ...
            ['P = ', num2str(mdl(j), 3) ' ('...
            num2str(p_expect(j), 3) '), L = ', num2str(L(j)), ', MLT = ' num2str(mlt(j))], 'FontSize', 8)
        text(max(xlim) - seconds(0.6), max(ylim) - range(ri)*0.1, ...
            ['(' num2str(lat(j)) ', ' num2str(long(j)) ', ', num2str(alt(j)) ')'], 'FontSize', 8);
        set(gca,'tag',num2str(j))
    end
    for j = (nr*nc)+sn:(nr*nc)+sn+1 % add extra at the bottom
        subplot(nr+1, nc, j-sn+1.5)
        ri = rate(st(j):et(j));
        di = d(st(j):et(j));
        nsf = find(ns > st(j) & ns < et(j));
        plot(di, ri, 'k', 'LineWidth', 1)
        hold on
        plot(d(ns(nsf)), rate(ns(nsf)), col)
        xlabel('Time (UTC)'), ylabel('Count')
        xticks(linspace(di(1),di(end),3))
        text(min(xlim), min(ylim), num2str(j), 'BackgroundColor', 'white')
        text(min(xlim) + seconds(0.05), max(ylim) - range(ri)*0.1, ...
            ['P = ', num2str(mdl(j), 3) ' ('...
            num2str(p_expect(j), 3) '), L = ', num2str(L(j)), ', MLT = ' num2str(mlt(j))], 'FontSize', 8)
        text(max(xlim) - seconds(0.6), max(ylim) - range(ri)*0.1, ...
            ['(' num2str(lat(j)) ', ' num2str(long(j)) ', ', num2str(alt(j)) ')'], 'FontSize', 8);
        set(gca,'tag',num2str(j))
    end
    
else % if you are just plotting normally
    for j = sn:((nr*nc)+sn-1)
        subplot(nr,nc,j-sn+1)
        ri = rate(st(j):et(j));
        di = d(st(j):et(j));
        nsf = find(ns > st(j) & ns < et(j));
        plot(di, ri, 'k', 'LineWidth', 1)
        hold on
        plot(d(ns(nsf)), rate(ns(nsf)), col)
        xlabel('Time (UTC)'), ylabel('Count')
        xticks(linspace(di(1),di(end),3))
        text(min(xlim), min(ylim), num2str(j), 'BackgroundColor', 'white')
        text(min(xlim) + seconds(0.05), max(ylim) - range(ri)*0.1, ...
            ['P = ', num2str(mdl(j), 3) ' ('...
            num2str(p_expect(j), 3) '), L = ', num2str(L(j)), ', MLT = ' num2str(mlt(j))], 'FontSize', 8)
        text(max(xlim) - seconds(0.6), max(ylim) - range(ri)*0.1, ...
            ['(' num2str(lat(j)) ', ' num2str(long(j)) ', ', num2str(alt(j)) ')'], 'FontSize', 8);
        set(gca,'tag',num2str(j))
    end
end
end

%% dayarrayfun
function daysvec = dayarrayfun
% Author: ChatGPT (prompting by Max Feinland)
startYear = 1996;
startDay = 160;
endYear = 2012;
endDay = 283;
spacing = 27;

daysArray = [];
yearsArray = [];

currentYear = startYear;
currentDay = startDay;

while currentYear < endYear || (currentYear == endYear && currentDay <= endDay)
    if isleap(currentYear)
        numDays = 366;
    else
        numDays = 365;
    end
    
    daysArray = [daysArray, currentDay]; %#ok<AGROW> 
    yearsArray = [yearsArray, currentYear]; %#ok<AGROW> 
    
    currentDay = currentDay + spacing;
    if currentDay > numDays
        currentDay = currentDay - numDays;
        currentYear = currentYear + 1;
    end
end
daysvec = [daysArray' yearsArray'];
end


%% isleap
function answ = isleap(year)
% checks for leap year
if mod(year, 4) == 0
    answ = 1;
else
    answ = 0;
end
end

%% expectedPeriod
function p_expect = expectedPeriod(L, maglat)
% Author: Max Feinland for Blum Research Group

% Inputs: L-shell, magnetic latitude, mirror latitude
% Outputs: expected period

c = 3e+8; % m/s
a = 6371000; % m, radius of Earth
m_e = 9.11e-31; % electron mass, kg
E_rest = m_e*c^2; % MeV
E = 1.6021766e-13; % 1 MeV in J
total_E = E + E_rest; % J
v = c*sqrt(1-(m_e*c^2/total_E)^2); % electron velocity, m/s

m = abs(maglat);

y = (1 + 3*cosd(m).^2).^(-1/4).*sind(m).^3;
T_0 = 1 + (1/(2*sqrt(3)))*log(2+sqrt(3));
T_1 = sqrt(2)*pi/6;
T = T_0 - 0.5*(T_0 - T_1)*(y + sqrt(y));

p_expect =  4*T*L*a/v;

end
