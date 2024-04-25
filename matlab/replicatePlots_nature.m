%% Relativistic Electron Microbursts in the Inner Radiation Belt
% Code to plot
% Author: Max Feinland
% Last modified: 4/23/24

% Housekeeping
clc
close all
clear

% this section allows the user to specify the plots they want to see
functionnums = [1 2 3 4 5];
funs = {@figure_1 @figure_2 @figure_3 @figure_4 @ figure_5};

pickplot = input(['Which plot(s) do you want to generate?\n Please enter '...
    'an array of numbers. ']);

for i = 1:length(pickplot)
    fig = pickplot(i);
    fun = funs{fig};
    fun()
end

%% Figure 1: Examples of Events
function figure_1

% reading in data
data = readtable('lshelldata.csv', 'VariableNamingRule', 'preserve');

list_of_indices = [2 17 22 13]; % index of example events (crown, decreasing, increasing, other)
letters = {'a)' 'b)' 'c)' 'd)'}; % for labelling

figure()

for i = 1:4
    % call plotEpoch to get the start and end times
    [~, ~, out_dats, d, r, st, et] = plotEpoch(data.date(list_of_indices(i)), ...
        data.time(list_of_indices(i)));
    % st = st + 5;

    pk_locs = d(out_dats); % peak locations, datetime
    pk_vals = r(out_dats); % peak locations, rate

    d = d(st-10:et+50); % narrowing down to portion of interest
    r = r(st-10:et+50); % standardizing amount of time makes plotting clearer

    t = posixtime(d); % necessary for the background interpolation

    % this is where the text will go
    ypos = 0.9*range(r) + min(r);
    xpos = 0.05*range(d) + min(d);

    axis tight

    TF = islocalmin(r); % finding "valleys" (local minima, to determine baseline)
    vq = interp1(t(TF == 1), r(TF == 1), t); % interpolate to get background trend

    adj = r - vq + min(r); % background-adjusted rate data

    % plotting
    subplot(2, 2, i)
    plot(d, r, 'k--', 'LineWidth', 1, 'DisplayName', 'Raw')
    hold on
    plot(d, adj,  'b-', 'LineWidth', 1, 'DisplayName', 'Adjusted')
    plot(pk_locs, pk_vals, 'r*', 'DisplayName', 'Identified Peak')
    legend()
    xlabel('Time (UTC)'), ylabel('Count Rate (#/20 ms)')
    % title(['Bouncing Packets on ', char(data.date(list_of_indices(i)))])
    text(xpos, ypos, letters{i}, 'FontSize', 18)
    text(d(end)-seconds(.5), range(r)/3+min(r), ...
        ['L = ' num2str(round(data.L(list_of_indices(i)), 2))])
end
end

%% Figure 2: Properties of Events

function figure_2

% reading in data
data = readtable('lshelldata.csv', 'VariableNamingRule', 'preserve');
load counts.mat num_pks
num_pks(num_pks > 100) = [];
time = data.date + data.time;
figure()

% omni data
dst = readtable("OMNI2_H0_MRG1HR_252725.csv", 'ReadVariableNames', false);
% contains DST information from January 1, 1997 to Dec 31, 2006. Available
% from OMNIweb

% formatting from csv; there are a lot of extraneous headers, cell arrays,
% strings, etc
dst(1:75, :) = [];
dst(end-3:end, :) = [];
dat = cell(height(dst), 1);
tim = table2array(dst(:,1));
for i = 1:height(tim)
    str = tim{i};
    dat{i} = str(1:10);
end
dstval = table2array(dst(:,2));

dsty = zeros(size(time));
% Loop through each value in t2
for i = 1:numel(time)
    [~, idx] = min(abs(time(i) - dat));
    dsty(i) = dstval(idx);
end

tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

% L-shell histogram
nexttile
histogram(data.L, linspace(1.5, 2.5, 9), 'FaceColor', 'y')
xlabel('L-shell'), ylabel('Number of Events')
text(1.5, 18, 'a)', 'FontSize', 18)

% MLT histogram
nexttile
histogram(data.MLT, linspace(0, 24, 9), 'FaceColor', 'cyan')
xlim([0 24])
xticks(linspace(0, 24, 9))
xlabel('Magnetic Local Time'), ylabel('Number of Events')
text(1, 10, 'b)', 'FontSize', 18)

nexttile
histogram(num_pks, 'FaceColor', 'g')
xlabel('Number of Peaks'), ylabel('Number of Events')
text(4, 9, 'c)', 'FontSize', 18)

nexttile
yyaxis left
histogram(dstval, linspace(-250, 50, 20), ...
    "DisplayName", "All Time")
ylabel('Number of Events (all time)')
yyaxis right
histogram(dsty, linspace(-250, 50, 20), ...
    "DisplayName", "Microburst Events")
xlabel('Dst index'), ylabel('Number of Events (microbursts)')
text(-250, 18, 'd)', 'FontSize', 18)
legend("Location", "west")

end

%% Figure 3: Locations of Events
function figure_3

% reading in data
data = readtable('lshelldata.csv', 'VariableNamingRule', 'preserve');
lgrid = readmatrix('Data/Lgrid.dat'); % provided to me by Sergio Vidal-Luengo
lgrid = lgrid(4000:end, :);

% extracting L-shell lat/lon data
glat_l2_nh = lgrid(:,1); % L = 2 lats, northern hemisphere
glon_l2_nh = lgrid(:,2); % L = 2 lons, northern hemisphere
glat_l3_nh = lgrid(:,3); % L = 3 lats, northern hemisphere
glon_l3_nh = lgrid(:,4); % L = 3 lons, northern hemisphere

glat_l2_sh = lgrid(:,15); % L = 2 lats, southern hemisphere
glon_l2_sh = lgrid(:,16); % L = 2 lons, southern hemisphere
glat_l3_sh = lgrid(:,17); % L = 3 lats, southern hemisphere
glon_l3_sh = lgrid(:,18); % L = 3 lons, southern hemisphere

glat_eq = lgrid(:,29); % equator lats
glon_eq = lgrid(:,30); % equator lons

% makes for easier plotting
for i = 1:length(data.long)
    if data.long(i) > 180
        data.long(i) = data.long(i) - 360;
    end
end

% date inputs for magnetic field background
first_day = datetime('2001-01-01');
last_day = datetime('2001-01-10');
[numeric_data] = plotEpoch2(first_day, last_day);

% extracting data
lat = numeric_data(:,1);
long = numeric_data(:,2);
b = numeric_data(:,4);

% also makes for easier plotting
long(long > 180) = long(long > 180) - 360;

% plotting
figure()
geoscatter(lat, long, 15, b, 'HandleVisibility', 'off')
c = colorbar;
c.Label.String = 'B-field intensity (Gauss)';
geolimits([-70, 70],[-180 180])
hold on
load coastlines coastlat coastlon
geoplot(coastlat, coastlon, 'k-', 'HandleVisibility', 'off')
title('Magnetic Field Intensity with Geographic Location')
colormap parula
geoscatter(data.lat, data.long, 30, 'MarkerFaceColor', [1 0.1 0.1], ...
     'MarkerEdgeColor', [0 0 0], 'DisplayName', 'Identified Packets')
hold on

% this is for text locations
midpt = round(length(glat_l2_nh)/2);
l2_col = [1 0.8 1]; % color coding Lshell lines
l3_col = [1 0.7 0.8]; % color coding Lshell lines

% adding L = 2, L = 3, equator lines
geoplot(glat_l2_nh, glon_l2_nh, '--', 'Color', l2_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(glat_l2_nh(midpt), glon_l2_nh(midpt), 'L = 2', 'Color', [1 1 1])

geoplot(glat_l3_nh, glon_l3_nh, '--', 'Color', l3_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(glat_l3_nh(midpt), glon_l3_nh(midpt), 'L = 3', 'Color', [1 1 1])

geoplot(glat_l2_sh, glon_l2_sh, '--', 'Color', l2_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(glat_l2_sh(midpt), glon_l2_sh(midpt), 'L = 2', 'Color', [1 1 1])

geoplot(glat_l3_sh, glon_l3_sh, '--', 'Color', l3_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(glat_l3_sh(midpt)+10, glon_l3_sh(midpt)+10, 'L = 3', 'Color', [1 1 1])

geoplot(glat_eq, glon_eq, '--', 'Color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')
text(glat_eq(1)+5, glon_eq(1)+10, 'Magnetic Equator', 'Color', [1 1 1])

geolimits([-60 60], [-135 90])
geolimits('manual')
legend()
title('Locations of Low L-shell Packets')

end

%% Figure 4: Lightning
function figure_4

% reading in data & formatting
events = readtable("lshelldata.csv", "VariableNamingRule", "preserve"); % events that I found
tstmp = events.date + events.time; % time of event
tstmp.Format = 'yyyy-MM-dd hh:mm:ss';

% conjugate magnetic points, calculated in Python using SpacePy coords
% class
ccor = readtable("interesting/conjlatlon.txt", "VariableNamingRule", "preserve");

% formatting longitude
for i = 1:length(events.long)
    if events.long(i) > 180
        events.long(i) = events.long(i) - 360;
    end
end

% lightning data provided by Dr. Ryan Said at Vaisala (thank you Ryan!!)
lightning = [];
for i = 1:45
    f = readmatrix(['interesting/' num2str(i+1) '.txt']);
    lightning = [lightning; f]; %#ok<AGROW>
end

% reading out data
y = lightning(:,2);
m = lightning(:,3);
d = lightning(:,4);
h = lightning(:,5);
minu = lightning(:,6);
s = lightning(:,7);
ns = lightning(:,8);

% formatting as one datetime
t = datetime(y, m, d, h, minu, s, ns/1e+6);

lat = lightning(:,9); % latitude, degrees
lon = lightning(:,10); % longitude, degrees
pkcur = abs(lightning(:,11)); % peak current, kA

% cell array containing lightning strikes for each event
idx = cell(45,1);
for k = 1:length(idx)
    idx{k} = find(abs(seconds(t - tstmp(k))) < 7);
end

loi = [1 10 35]; % indices of events that I used as examples
numberings1 = {'a)' 'b)' 'c)'};
numberings2 = {'d)' 'e)' 'f)'};
figure()
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact')

for j = 1:3
    if events.lat(loi(j)) < 0
        longdiff = abs(lon(idx{loi(j)}) - ccor.conjlon(loi(j)));
    else
        longdiff = abs(lon(idx{loi(j)}) - events.long(loi(j)));
    end

    nexttile(j)
    if j == 1
        legendyn = 1;
    else
        legendyn = 0;
    end
    plotDate2(events.date(loi(j)), events.time(loi(j)), t(idx{loi(j)}), pkcur(idx{loi(j)}), longdiff, legendyn)
    tx = events.date(loi(j)) + events.time(loi(j)) - seconds(1.5);
    text(tx, 160, numberings1{j}, 'FontSize', 18)

    nexttile(j+3)
    % title(['Lightning on ', char(events.date(loi(j)))])
    geoscatter(events.lat(loi(j)), events.long(loi(j)), 'red', 'filled', ...
        'DisplayName', 'Observed Event')
    hold on
    geobasemap satellite
    geoscatter(lat(idx{loi(j)}), lon(idx{loi(j)}), 'yellow', 'filled', ...
        'DisplayName', 'Lightning Strikes')
    geoscatter(ccor.conjlat(loi(j)), ccor.conjlon(loi(j)), 'green', 'filled', ...
        'DisplayName', 'Conjugate Mirror Point')
    if j == 1
    legend('Location', 'southwest')
    end
    text(35, -150, numberings2{j}, 'Color', 'white', 'FontWeight', 'bold' , ...
        'FontSize', 18)
        % text(50, -90, numberings2{j}, 'Color', 'white', 'FontWeight', 'bold')
end
end

%% Figure 5: DST
function figure_5

figure()
load bsmat.mat b s unqday dstavg timeElapsed % load in data

% reading in data

x = islocalmin(dstavg);
storm_times = find(x == 1);
storm_times = storm_times(dstavg(storm_times) < -50);
dt = diff(unqday(storm_times))*365;

l_vec = 1:0.1:7;

tiledlayout(2, 3)
nexttile([1 3])
[rows, ~] = size(b);
b2 = [];
for i = 1:rows
    b2 = [b2; b(i,:)']; %#ok<AGROW>
end
s2 = repelem(s, numel(l_vec));
y = repmat(l_vec, 1, rows);
scatter(s2, y, 100, log10(b2), 'square', 'filled', ...
    'HandleVisibility', 'off');

x = islocalmin(dstavg);
storm_times = find(x == 1);
storm_times = storm_times(dstavg(storm_times) < -50);

% add low-L events
data = readtable('lshelldata.csv');
L = data.L;
yrs = decyear(data.date);

% Put it all together!
axis tight
yyaxis left
scatter(yrs, L, 50, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], ...
     'DisplayName', 'Observed Events') % plot microburst events
text(1998.1, 4, "a)", "BackgroundColor", "white", "fontsize", 20)

ylabel('L-shell')
xlabel('Decimal year')
ylim([1.5 4])

yyaxis right
plot(unqday, dstavg, '-.', 'Color', [0 0 0], 'LineWidth', 1, ...
    'DisplayName', 'Dst Index'); % plot dst index
hold on
plot(unqday(storm_times), dstavg(storm_times), "pentagram", 'MarkerSize', 10,...
    "MarkerFaceColor", "white",  "MarkerEdgeColor", "magenta", "DisplayName", "Geomagnetic Storms")
xlim([min(s) max(s)])
ylim([-250 50]);
legend()
ylabel('Dst Index')

% making the colorbar look as nice as possible :)
cmap0 = [0 0 0; 0 0 0; 75 0 130; 0 0 255;  0 255 255; 0 255 0; 255 255 0; ...
    255 127 0; 255 0 0];
x0 = linspace(0, 255, 9);
xq = linspace(0, 255, 255);
vq = interp1(x0, cmap0, xq);
vq = vq/255;
c = colorbar('Ticks', [1 1.5 2 2.5 3 3.5 4 4.5], 'Limits', [1 5]);
colormap(vq)

c.Label.String = 'Log_{10} >1 MeV Electron Counts/20ms';
c.Label.FontSize = 11;

nexttile(4, [1 3])
%% Time elapsed figure
mybins = 0:5:70;
histogram(dt, [mybins Inf], "DisplayName", "Delay between storms")
hold on
histogram(timeElapsed, mybins, "DisplayName", ...
    "Delay between microbursts & storms")
xlim([-5 75])
xlabel("Days since storm")
ylabel("Number of events")
legend('FontSize', 10)
text(-3, 25, "b)", "BackgroundColor", "white", "fontsize", 20)
% Initialize the cell array
cell_array = cell(1, 60/10 +1);

% Generate the cell array
for i = 1:numel(cell_array)
    cell_array{i} = num2str((i - 1) * 10);
end
cell_array{end+1} = '>=70';
xticks(0:10:70)
xticklabels(cell_array)
end

%% Other Subfunctions (necessary to get some of the code to run)
function [num_bursts, out_dts, out_dats, d, r, st, et] = plotEpoch(first_day, timestamp)

% Date Created: 5/15/23
% Last Modified: 6/4/23
% Author: Max Feinland for Blum Research Group, LASP

% Inputs: first day of epoch, last day of epoch, L-shell data (time and 
% value), minimum allowable period, maximum allowable period

% Outputs: plots of each microburst in epoch (by day), total number of
% microbursts, prints the period, std, L-shell, and percent error

% this helps the formats to line up, makes it easier to find datetimes
datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')


% initializing
num_bursts = 0; % total number of bursts found


    y = year(first_day);
    doy = day(first_day, 'dayofyear');
    
        % makes sure filename has leading zeros if necessary
        if strlength(num2str(doy)) == 3
            filename = ['hhrr', num2str(y), num2str(doy), '.txt'];
        elseif strlength(num2str(doy)) == 2
            filename = ['hhrr', num2str(y), '0', num2str(doy), '.txt'];
        elseif strlength(num2str(doy)) == 1
            filename = ['hhrr', num2str(y), '00', num2str(doy), '.txt'];
        end

        % pulls file from online SAMPEX data repository
        if doy <= 182 && y == 2004 || y < 2004 % data for which it's zipped
            url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
                filename, '.zip'];
            unzip(url) % find and unzip data file
            data = readmatrix(filename); % read in data file
        else % data for which it's not zipped
            url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
                filename];
            file = webread(url);
            
            % parse the character vector into a table (code adapted from CGPT)
            celltab = textscan(file, '%f %f %f %f %f %f %f', 'Delimiter', '\t', 'HeaderLines', 1);
            
            % converts the parsed data to a table
            data = [celltab{1}, celltab{2}, celltab{3}, celltab{4}, celltab{5}, celltab{6}, celltab{7}];
        end


            % call O'Brien 
            [~, d, r, so, eo, ~] = obrien100(data, string(first_day));
            
            num = 0;

            % call bouncingPacketSplits
            [st, et, pksi] = bouncingPacketSplits(so, eo, r);

            yesPlot = zeros(size(st));
            for j = 1:length(st)

                ns = pksi >= st(j) & pksi <= et(j); % finds indices with peaks
               
                length_int = et(j) - st(j); % number of indices in interval
                unqvals = unique(r(round(st(j)):round(et(j)))); % unique numbers in interval

                % (second condition exists since sometimes the data
                % gets "corrupted", probably from SEUs, and only shows
                % a few possible values)

                % to be printed, all these conditions must be true:
                % short, not corrupted, and low l-shell

                conds = [length_int < 750, length(unqvals) > 10, ...
                    abs(seconds(timeofday(d(st(j))) - timestamp)) < 3];
    
                if all(conds) % if they are indeed all met
                    % fprintf('[Event %i] %s %.3f\n', ...
                    %     num+1, d(st(j)), mdl)

                    num = num + 1; % increment number of bursts found
                    yesPlot(j) = 1; % include this interval in output
                    yesn = ns == 1;
                end
            end

            first_day.Format = 'dd-MMM-yyyy';

            st = st(yesPlot == 1);
            et = et(yesPlot == 1);
            out_dts = first_day;
            out_dats = [pksi(yesn)];
    
            num_bursts = num_bursts + num; % total microbursts observed
    
            if isfile(filename) % if you had to download the .zip, delete it now
                delete(filename)
            end
end

%% O'Brien

function [t50, d50, r50, starts, ends, ns] = obrien100(data, day)
datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')
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

d50 = datetime(day, 'InputFormat', 'MM/dd/uuuu') + seconds(t50); % convert from seconds to UTC (nice for plotting)
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

%% PlotEpoch
function [out_dats] = plotEpoch2(first_day, last_day)

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
        out_dats = [];
    elseif first_day_d >= datetime('2005-01-01') && first_day_d <= datetime('2005-01-09')...
            || last_day_d >= datetime('2005-01-01') && last_day_d <= datetime('2005-01-09')
   
            fprintf(['Error: time range includes data that does not exist in the SAMPEX ' ...
                'ephemeris repository. Please try a new input.\n'])
            % do not plot or print anything
        dates_ok = 0; % dates not valid
        out_dats = [];
    end
catch % nonreal dates
    fprintf(['Error: non-real date. This can happen if a specified date does not exist' ...
        ', e.g. February 31st, \nor if there was a typo during input. Please try a new input.\n'])
    % do not print or plot anything
    dates_ok = 0; % dates not valid
    out_dats = [];
end

if dates_ok == 1 % if dates exist and are sequential

doy1 = day(first_day_d, 'dayofyear'); % day of start (integer)
doy2 = day(last_day_d, 'dayofyear'); % day of end (integer)
y = year(first_day_d); % year in question
y2 = year(last_day_d); % year of end date (important only if y != y2)

if y == 1996 && doy1 < 220 || y < 1996 || y > 2013 || y == 2013 && doy1 > 312
    % if the data is outside the time range that the data spans
    disp(['Invalid command; please choose dates between August 7th, 1996 ' ...
        'and November 8th, 2013.'])
elseif y ~= y2
    disp('Invalid command; please choose two dates in the same year.')
elseif y == 2011 && sum(ismember([doy1 doy2], 297:323)) > 0
     disp('Request contains data outside ephemeris range. Please try again.')
else
    fprintf('Working...\n') % status update for user
    
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
            C = repmat(A, [1 26]);
            D = [B B B A A A A B B B repmat(A, [1 10]) B B B repmat(A, [1 10]) ...
                B B repmat(A, [1 10]) B C];
            % tells textscan to look at columns 1, 2, 3, 8, 9, 10, 21, 22,
            % 23, 34, and 35
            D(end) = [];

            celltab = textscan(file, D, 'Delimiter', '\t');
            data_chunk = [celltab{1}, celltab{2}, celltab{3}, celltab{4} ...
                celltab{5} celltab{6} celltab{7} celltab{8} celltab{9}...
                celltab{10} celltab{11} celltab{12}];

            l_data = [l_data; data_chunk]; %#ok<AGROW> 

        else % data for which it's zipped
            url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/PSSet/Text/', ...
            ephem_filename{i-n+1} '.zip'];
            
            unzip(url) % find and unzip data file
            data_chunk = readmatrix(ephem_filename{i-n+1}, 'HeaderLines', 60); % read in data file
            data_chunk = [data_chunk(:,1:3) data_chunk(:,8:10) data_chunk(:,21:23) ...
                 data_chunk(:,34:35), data_chunk(:,46)];

            l_data = [l_data; data_chunk]; %#ok<AGROW> 
        end
    end

    % status update
    fprintf('Ephemeris data retrieved.\n')

    % extracting relevant data
    geo_long = l_data(:,4);
    geo_lat = l_data(:,5);
    geo_alt = l_data(:,6);
    bfield_all = l_data(:,8);
    loss_cone = l_data(:,11);
    equator_bfield = l_data(:,12);

    out_dats = [geo_lat geo_long geo_alt bfield_all equator_bfield loss_cone];

   
end

% deleting downloaded ephemeris files
for k = 1:length(ephem_filename)
    if isfile(ephem_filename{k})
        delete(ephem_filename{k})
    end
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

function plotDate2(dat, time, lt, lcur, ldiff, legendyn)

doy = day(dat, 'dayofyear');
y = year(dat);
 % makes sure filename has leading zeros if necessary
if strlength(num2str(doy)) == 3
    filename = ['hhrr', num2str(y), num2str(doy), '.txt'];
elseif strlength(num2str(doy)) == 2
    filename = ['hhrr', num2str(y), '0', num2str(doy), '.txt'];
elseif strlength(num2str(doy)) == 1
    filename = ['hhrr', num2str(y), '00', num2str(doy), '.txt'];
end

% pulls file from online SAMPEX data repository
if doy <= 182 && y == 2004 || y < 2004 % data for which it's zipped
    url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
        filename, '.zip'];
    unzip(url) % find and unzip data file
    data = readmatrix(filename); % read in data file
else % data for which it's not zipped
    url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
        filename];
    file = webread(url);
    
    % parse the character vector into a table (code adapted from CGPT)
    celltab = textscan(file, '%f %f %f %f %f %f %f', 'Delimiter', '\t', 'HeaderLines', 1);
    
    % converts the parsed data to a table
    data = [celltab{1}, celltab{2}, celltab{3}, celltab{4}, celltab{5}, celltab{6}, celltab{7}];
end

[t, d, r, so, eo, ~] = obrien100(data, string(dat));
            
[st, et, pksi] = bouncingPacketSplits(so, eo, r);

for j = length(st):-1:1
    if abs(timeofday(d(st(j))) - duration(time)) < seconds(5)
        ns = pksi >= st(j)+10 & pksi <= et(j)-10; % finds indices with peaks

        t1 = t(st(j)-150:st(j)+300); % narrowing down to portion of interest
        d1 = d(st(j)-150:st(j)+300); % narrowing down to portion of interest
        r1 = r(st(j)-150:st(j)+300); % standardizing amount of time makes plotting clearer
        TF = islocalmin(r1); % finding "valleys" (local minima, to determine baseline)
        vq = interp1(t1(TF == 1), r1(TF == 1), t1); % interpolate to get background trend

        adj = r1 - vq + min(r1); % background-adjusted rate data
        plotLightning(adj, d1,  ['Bouncing Packets and Lightning on ', char(dat)], ...
                st(j), pksi(ns), 'r*', lt, lcur, ldiff, legendyn)
    end
end

% cleaning up afterward
if isfile(filename)
    delete(filename)
end

end

%% plotLightning
function plotLightning(rate, d, titl, ~, ~, ~, lt, lcur, ldiff, legendyn)
% Inputs: Rate data, date data, title, starttime, peaks, color,
% lightning timestamps, lightning current, difference in longitude

% Outputs: none. Plots for you :)

% title(titl)

% on left axis, plot electron flux
% ri = rate(st:et);
% di = d(st:et);
% nsf = find(ns > st+b & ns < et-b);
yyaxis left
plot(d, rate, 'k', 'LineWidth', 1, 'DisplayName', 'Adj. electron count rate')
hold on
xlabel('Time (UTC)'), ylabel('Count (#/20 ms)')
tstart = dateshift(d(1),'start','second');
tend = dateshift(d(end),'end','second');
xticks(tstart:seconds(2):tend)

% on right axis, plot lightning current
yyaxis right
plot([lt lt], [zeros(size(lcur)) lcur], 'k--', 'HandleVisibility', 'off')
nc = ldiff/25;
scatter(lt, lcur, 30, nc, 'filled', 'DisplayName', 'Lightning Discharges');
ylim([0 175]) % so that you can see the comparatively low lightning for Africa

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';

% had to get a little creative with the colorbar for some reason
c = colorbar('Ticks', linspace(0, 1, 6));
clim([0 1]);
c.Label.String = '\Delta Longitude (^{\circ})';
c.TickLabels = string(linspace(0, 25, 6));
if year(d) ~= 2004
    colorbar('off')
end

colormap parula % I like this colormap for the longitude, the yellow "feels"
% further away from the blue imo

ylabel('|Peak Current| (kA)')
if year(d) == 1998
    xlim([min(d) + seconds(1) max(d)])
else
    xlim([min(d) max(d)])
end
if legendyn == 1
    legend()
end
end