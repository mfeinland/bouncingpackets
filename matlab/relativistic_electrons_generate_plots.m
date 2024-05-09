%% Relativistic Electron Microbursts from the Inner Radiation Belt
% Code to plot
% Author: Max Feinland for Blum Research Group (THANK YOU LAUREN)
% Last modified: 5/8/24

% Please email me (max.feinland@colorado.edu) if you have any questions!
% I'm graduating so if that email bounces try my personal email:
% maximiliusf@gmail.com

% Housekeeping
clc
close all
clear

% this section allows the user to specify the plots they want to see
functionnums = 1:5;
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
% first the identified event data
events = readtable('SourceData.xlsx', 'Sheet', 'allfigures', 'VariableNamingRule', 'preserve');

% this contains the HILT data for each event, so you don't have to download
% the whole day of data :)
hilt = readtable('SourceData.xlsx', 'Sheet', 'figure1', 'VariableNamingRule', 'preserve');

loi = [2 17 22 13]; % list of indices of example events 
% (crown, decreasing, increasing, other)
letters = {'a)' 'b)' 'c)' 'd)'}; % for subplot labeling

figure('WindowState', 'maximized')

for i = 1:4
    % extracting columns from xlsx
    % they are named like 'datetime1' 'rate1' ... 'datetime2' and so on
    % so this code just extracts the relevant columns for each subplot.
    datecolname = ['datetime', num2str(i)];
    ratecolname = ['rate', num2str(i)];
    pkcolname = ['pk', num2str(i)];

    d = hilt{:,datecolname};
    r = hilt{:,ratecolname};
    pks = hilt{:,pkcolname};


    t = posixtime(d); % necessary for the background interpolation
    TF = islocalmin(r); % finding "valleys" (local minima, to determine baseline)
    vq = interp1(t(TF == 1), r(TF == 1), t); % interpolate to get background trend
    adj = r - vq + min(r); % background-adjusted rate data

    % this is where the subplot label will go
    ypos = 0.9*range(r) + min(r);
    xpos = 0.05*range(d) + min(d);

    axis tight

    % plotting
    subplot(2, 2, i)
    plot(d, r, 'k--', 'LineWidth', 1, 'DisplayName', 'Raw')
    hold on
    plot(d, adj,  'b-', 'LineWidth', 1, 'DisplayName', 'Adjusted')
    plot(d(pks == 1), r(pks == 1), 'r*', 'DisplayName', 'Identified Peak')
    legend()
    xlabel('Time (UTC)'), ylabel('Count Rate (#/20 ms)')
    text(xpos, ypos, letters{i}, 'FontSize', 18) % put subplot label
    % add L-shell
    text(d(end)-seconds(.5), range(r)/3+min(r), ...
        ['L = ' num2str(round(events.L(loi(i)), 2))])
end
end

%% Figure 2: Properties of Events

function figure_2

% identified event data
events = readtable('SourceData.xlsx', 'Sheet', 'allfigures', 'VariableNamingRule', 'preserve');

% contains DST data from January 1, 1997 00:30 to Dec 31, 2006 at 12:30, 
% spaced every half hour. Available from OMNIweb
dst = readtable('SourceData.xlsx', 'Sheet', 'figure2');

figure()
tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

% L-shell histogram
nexttile
histogram(events.L, linspace(1.5, 2.5, 9), 'FaceColor', 'y')
xlabel('L-shell'), ylabel('Number of Events')
text(1.5, 18, 'a)', 'FontSize', 18)

% MLT histogram
nexttile
histogram(events.MLT, linspace(0, 24, 9), 'FaceColor', 'cyan')
xlim([0 24])
xticks(linspace(0, 24, 9))
xlabel('Magnetic Local Time'), ylabel('Number of Events')
text(1, 10, 'b)', 'FontSize', 18)

% Number of peaks histogram
nexttile
histogram(events.num_pks, 'FaceColor', 'g')
xlabel('Number of Peaks'), ylabel('Number of Events')
text(4, 9, 'c)', 'FontSize', 18)

% Dst histogram -- just found nearest times and put into allfigures sheet
nexttile
yyaxis left
% dst of all time (if all time is from 1998 to 2006)
% I was only alive for four years of that time so it's not all time to me!
% Sorry if that made anyone feel old but I sort of doubt anyone is reading
% these comments too closely
histogram(dst.dst, linspace(-250, 50, 20), ...
    "DisplayName", "All Time")
ylabel('Number of Events (all time)')
yyaxis right
% dst index for the events
histogram(events.Dst, linspace(-250, 50, 20), ...
    "DisplayName", "Microburst Events")
xlabel('Dst index'), ylabel('Number of Events (microbursts)')
text(-250, 18, 'd)', 'FontSize', 18)
legend("Location", "west")
end

%% Figure 3: Locations of Events
function figure_3

% identified event data
events = readtable('SourceData.xlsx', 'Sheet', 'allfigures', 'VariableNamingRule', 'preserve');

% magnetic field data, lat, lon from SAMPEX attitude
sampex = readtable('SourceData.xlsx', 'Sheet', 'figure3_bfield', 'VariableNamingRule', 'preserve'); 

% provided to me by Sergio Vidal-Luengo (thanks Sergio!)
% "I calculated them using AACGM coordinates with a python function"
% -Sergio
lgrid = readtable('SourceData.xlsx', 'Sheet', 'figure3_lcontours', 'VariableNamingRule', 'preserve'); 

% makes for nicer plotting
events.long(events.long > 180) = events.long(events.long > 180) - 360;

figure()

% scatterplot of b-field intensity
geoscatter(sampex.lat, sampex.lon, 15, sampex.bfield, 'filled', 'HandleVisibility', 'off')
c = colorbar;
c.Label.String = 'B-field intensity (Gauss)';
hold on

% add world map outline (just to contextualize locations)
load coastlines coastlat coastlon
geoplot(coastlat, coastlon, 'k-', 'HandleVisibility', 'off')
colormap parula

% add identified events
geoscatter(events.lat, events.long, 30, 'MarkerFaceColor', [1 0.1 0.1], ...
     'MarkerEdgeColor', [0 0 0], 'DisplayName', 'Identified Packets')
hold on

% this is for text locations
midpt = round(length(lgrid.lat2N)/2);

 % color coding L-shell lines
l2_col = [1 0.8 1];
l3_col = [1 0.7 0.8];
textcol = [1 1 1];

% adding contour lines (and labels) for L = 2, L = 3, and equator
geoplot(lgrid.lat2N, lgrid.lon2N, '--', 'Color', l2_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(lgrid.lat2N(midpt), lgrid.lon2N(midpt), 'L = 2', 'Color', textcol)

geoplot(lgrid.lat3N, lgrid.lon3N, '--', 'Color', l3_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(lgrid.lat3N(midpt), lgrid.lon3N(midpt), 'L = 3', 'Color', textcol)

geoplot(lgrid.lat2S, lgrid.lon2S, '--', 'Color', l2_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(lgrid.lat2S(midpt), lgrid.lon2S(midpt), 'L = 2', 'Color', textcol)

geoplot(lgrid.lat3S, lgrid.lon3S, '--', 'Color', l3_col, 'LineWidth', 1, 'HandleVisibility', 'off')
text(lgrid.lat3S(midpt)+10, lgrid.lon3S(midpt)+10, 'L = 3', 'Color', textcol)

geoplot(lgrid.lat0, lgrid.lon0, '--', 'Color', [1 1 1], 'LineWidth', 1, 'HandleVisibility', 'off')
text(lgrid.lat0(1)+5, lgrid.lon0(1)+10, 'Magnetic Equator', 'Color', textcol)

geolimits([-60 60], [-135 90])
legend()
end

%% Figure 4: Lightning
function figure_4

% identified event data
events = readtable('SourceData.xlsx', 'Sheet', 'allfigures', 'VariableNamingRule', 'preserve');

% lightning data only for events, provided by Dr. Ryan Said at Vaisala
% (thanks Ryan!)
lightning = readtable('SourceData.xlsx', 'Sheet', 'figure4', 'VariableNamingRule', 'preserve');

% formatting longitude
events.long(events.long > 180) = events.long(events.long > 180) - 360;

loi = [1 10 35]; % list of indices of events that I used as examples
numberings1 = {'a)' 'b)' 'c)'};
numberings2 = {'d)' 'e)' 'f)'};
figure('WindowState', 'maximized')
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact')

for i = 1:3
    % extracting columns from xlsx
    datecolname = ['hilt_t', num2str(i)];
    ratecolname = ['hilt_r', num2str(i)];
    ltcolname = ['lightning_t', num2str(i)];
    pkcurcolname = ['pkcur', num2str(i)];
    longdiffcolname = ['longdiff', num2str(i)];
    latcolname = ['lat', num2str(i)];
    loncolname = ['lon', num2str(i)];

    d = lightning{:,datecolname};
    r = lightning{:,ratecolname};
    lt = lightning{:,ltcolname};
    pkcur = lightning{:,pkcurcolname};
    longdiff = lightning{:,longdiffcolname};
    lat = lightning{:,latcolname};
    lon = lightning{:,loncolname};

    t = posixtime(d); % necessary for the background interpolation
    TF = islocalmin(r); % finding "valleys" (local minima, to determine baseline)
    vq = interp1(t(TF == 1), r(TF == 1), t); % interpolate to get background trend
    adj = r - vq + min(r); % background-adjusted rate data

    nexttile(i)
    % First, plot upper plot: time
    % on left axis, plot adjusted electron flux
    yyaxis left
    plot(d, adj, 'k', 'LineWidth', 1, 'DisplayName', 'Adj. electron count rate')
    hold on
    xlabel('Time (UTC)'), ylabel('Count (#/20 ms)')

    % make sure the ticks are spaced the way I want
    tstart = dateshift(d(1),'start','second');
    tend = dateshift(d(end),'end','second');
    xticks(tstart:seconds(2):tend)
    
    % on right axis, plot lightning current
    yyaxis right
    % to make it look like a stem plot, per Bob Marshall's request
    for j = 1:numel(pkcur)
        line([lt(j), lt(j)], [0, pkcur(j)], 'LineStyle', '--', 'Color', 'k', ...
            'HandleVisibility', 'off');
    end

    % colorbar is 'longdiff'/25 to standardize colorbar across plots
    scatter(lt, pkcur, 30, longdiff/25, 'filled', 'DisplayName', 'Lightning Discharges');
    ylim([0 175]) % to standardize and show all instances have amperage above 100 kA
    ylabel('|Peak Current| (kA)')
    
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    
    % had to get a little creative with the colorbar for some reason
    c = colorbar('Ticks', linspace(0, 1, 6));
    colormap parula
    clim([0 1]);
    c.Label.String = '\Delta Longitude (^{\circ})';
    c.TickLabels = string(linspace(0, 25, 6));
    if i ~= 3 % only put colorbar on last plot
        colorbar('off')
    end

    % Idk why the first plot is weird but if I don't do this, the subplot 
    % labels are all wrong
    if i == 1
        xlim([min(d)+seconds(1) max(d)])
    else
        xlim([min(d) max(d)])
    end
    % aforementioned subplot label
    tx = events.date(loi(i)) + events.time(loi(i)) - seconds(1.5);
    text(tx, 160, numberings1{i}, 'FontSize', 18)

    if i == 1 % only put legend on first plot
        legend()
    end


    % Next, bottom plot: space
    nexttile(i+3)
    % put events
    geoscatter(events.lat(loi(i)), events.long(loi(i)), 'red', 'filled', ...
        'DisplayName', 'Observed Event')
    hold on
    geobasemap satellite
    % lightning data
    geoscatter(lat, lon, 'yellow', 'filled', 'DisplayName', 'Lightning Strikes')
    % magnetic conjugate point
    geoscatter(events.conjlat(loi(i)), events.conjlon(loi(i)), 'green', 'filled', ...
        'DisplayName', 'Conjugate Mirror Point')

    if i == 1 % only put legend on first plot
        legend('Location', 'west')
    end

    % subplot label
    text(35, -150, numberings2{i}, 'Color', 'white', 'FontWeight', 'bold' , ...
        'FontSize', 18)
end
end

%% Figure 5: DST
function figure_5

% identified event data
events = readtable('SourceData.xlsx', 'Sheet', 'allfigures', 'VariableNamingRule', 'preserve'); 
% dst data and some others
dst = readtable('SourceData.xlsx', 'Sheet', 'figure5', 'VariableNamingRule', 'preserve'); 

% find storm times (Dst < -50 nT)
TF = islocalmin(dst.dst);
storm_times = find(TF == 1);
storm_times = storm_times(dst.dst(storm_times) < -50);
dt = diff(dst.omni_t(storm_times))*365; % time between storms

figure('WindowState', 'maximized')
tiledlayout(2, 3)
nexttile([1 3])

% scatterplot of log10 count rate, binned in l-shell
% Modeled after Li et al 2017, but that was using PET data and we got a
% comment from a reviewer that (rightly) pointed out that this would
% generate a completely different scatterplot using HILT
scatter(dst.sampex_t, dst.lshell, 100, log10(dst.counts), 'square', 'filled', ...
    'HandleVisibility', 'off');

yyaxis left
% plot microburst events
scatter(decyear(events.date), events.L, 50, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0 0 0], 'DisplayName', 'Observed Events') 

% subplot label
text(1998.1, 4, "a)", "BackgroundColor", "white", "fontsize", 20)
ylabel('L-shell')
xlabel('Decimal year')
ylim([1.5 4])

% plot timeseries of dst
yyaxis right
plot(dst.omni_t, dst.dst, '-.', 'Color', [0 0 0], 'LineWidth', 1, ...
    'DisplayName', 'Dst Index'); % all dst
hold on
plot(dst.omni_t(storm_times), dst.dst(storm_times), "pentagram", ...
    "MarkerSize", 10, "MarkerFaceColor", "white",  "MarkerEdgeColor", ...
    "magenta", "DisplayName", "Geomagnetic Storms") % storm times
xlim([min(dst.sampex_t) max(dst.sampex_t)])
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

% Time elapsed figure
nexttile(4, [1 3])
bins = 0:2:70; % bins of days elapsed
histogram(dt, [bins Inf], "DisplayName", "Delay between storms")
hold on
histogram(events.deltat, bins, "DisplayName", ...
    "Delay between microbursts & storms")
xlim([-5 72])
xlabel("Days since storm")
ylabel("Number of events")
legend('FontSize', 10)
% subplot label
text(-3, 17, "b)", "BackgroundColor", "white", "fontsize", 20)

% initialize cell array containing tick labels
myticklabels = cell(1, 60/10 +1);

% put tick labels in there
for i = 1:numel(myticklabels)
    myticklabels{i} = num2str((i-1) * 10);
end
myticklabels{end+1} = '\geq70'; % add custom label at the end for overflow
xticks(0:10:70)
xticklabels(myticklabels)
end