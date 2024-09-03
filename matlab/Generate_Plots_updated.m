%% Relativistic Electron Microbursts from the Inner Radiation Belt
% Code to plot
% Author: Max Feinland for Blum Research Group (THANK YOU LAUREN)
% Last modified: 8/8/24

% Please email me (max.feinland@colorado.edu) if you have any questions!
% I just graduated so if that email bounces try my personal email:
% maxfeinland@gmail.com

% Housekeeping
clc
close all
clear

% this section allows the user to specify the plots they want to see
functionnums = 1:6;
funs = {@figure_1 @figure_2 @figure_3 @figure_4 @figure_5 @figure_6};

pickplot = input(['Which plot(s) do you want to generate?\nPlease enter '...
    'an array of numbers, or press 8 for all. ']);
disp(pickplot)

if pickplot==8
    for i = functionnums
    fun = funs{i};
    fun()
    end
else
    for i = 1:length(pickplot)
    fig = pickplot(i);
    fun = funs{fig};
    fun()
    end
end

%% Figure 1: Examples of Events
function figure_1

% reading in data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 1', ...
    'VariableNamingRule', 'preserve');

letters = {'a' 'b' 'c' 'd'}; % for subplot labeling

figure('WindowState', 'maximized')
subplot_id = data.Var2; % tells you which subplot you're on

for i = 1:4
    % find row containing the text 'Figure 1a' or 'Figure 1b', etc.
    this_figname = strcat('Figure 1', letters{i}); 
    this_fig_idx = find(strcmp(subplot_id, this_figname));

    % find row containing next text, or handle edge case if this is fig 1d
    if i < 4 
        next_figname = strcat('Figure 1', letters{i+1});
        next_fig_idx = find(strcmp(subplot_id, next_figname));
        next_fig_idx = next_fig_idx - 5;
    else
        next_fig_idx = length(subplot_id);
    end

    % extracting columns and relevant rows from spreadsheet
    d = data.timestamp(this_fig_idx:next_fig_idx);
    r = data.countrate(this_fig_idx:next_fig_idx);
    pks = data.pks(this_fig_idx:next_fig_idx);
    L = data.L_shell(this_fig_idx);


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
    plot(d(pks == 1), r(pks == 1), 'r*', 'DisplayName', 'Identified Peak', ...
        'LineWidth', 1)
    if i == 2
        legend()
    end
    xlabel('Time (UTC)', 'FontSize', 12) 
    ylabel('Electron count rate (#/20 ms)', 'FontSize', 12)
    xticks(dateshift(min(d), 'start', 'second'):seconds(1):max(d))
    text(xpos, ypos, strcat(letters{i}, ')'), 'FontSize', 20) % put subplot label
    % add L-shell
    text(d(end)-seconds(.6), range(r)/3+min(r), ...
        ['L = ' num2str(round(L, 2))], 'FontSize', 14, ...
        'BackgroundColor', 'white')
end
end


%% Figure 2: Properties of Events
function figure_2

% identified event data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 2', ...
    'VariableNamingRule', 'preserve');

figure()
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

% L-shell histogram
nexttile
histogram(data.L_shell, 1.5:0.1:2.5, 'FaceColor', 'y')
xlabel('L-shell'), ylabel('Number of events')
text(1.5, 18, 'a)', 'FontSize', 18)

% MLT histogram
nexttile
histogram(data.MLT, linspace(0, 24, 9), 'FaceColor', 'c')
xlim([0 24])
xticks(linspace(0, 24, 9))
xlabel('Magnetic local time'), ylabel('Number of events')
text(1, 10, 'b)', 'FontSize', 18)

% Number of peaks histogram
nexttile
histogram(data.num_pks, 4:12, 'FaceColor', 'g')
xlabel('Number of peaks'), ylabel('Number of events')
text(4, 9, 'c)', 'FontSize', 18)

% Mean peak spacing histogram
nexttile
histogram(data.spacing, 0.15:0.025:0.3, 'FaceColor', 'k')
xlabel("Minimum peak spacing per event (s)"), ylabel("Number of events")
ylim([0 20])
text(0.15, 18, 'd)', 'FontSize', 18)

% Dst histogram -- just found nearest times and put into allfigures sheet
nexttile([1 2])
yyaxis left
% dst of all time (if all time is from 1998 to 2006)
histogram(data.all_dst, linspace(-250, 50, 25), ...
    "DisplayName", "All Time", 'FaceColor', 'blue')
ylabel('Number of events (all time)')
yyaxis right
% dst index for the events
histogram(data.event_dst, linspace(-250, 50, 25), ...
    "DisplayName", "Microburst Events", 'FaceColor', 'red')
xlabel('Dst index'), ylabel('Number of events (microbursts)')
text(-250, 13, 'e)', 'FontSize', 18)
legend("Location", "west", 'FontSize', 10)
end

%% Figure 3: Locations of Events
function figure_3

% identified event data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 3', ...
    'NumHeaderLines', 2);

% provided to me by Sergio Vidal-Luengo (thanks Sergio!)
% "I calculated them using AACGM coordinates with a python function"
% -Sergio

figure()

% find magnetic field intensity data
bfield_idx = find(strcmp(data{:,2}, 'Magnetic field intensity'));

% scatterplot of b-field intensity
geoscatter(data{bfield_idx:end, 3}, data{bfield_idx:end, 4}, 25, ...
    data{bfield_idx:end,5}, 'filled', 'HandleVisibility', 'off')
c = colorbar;
c.Label.String = 'Magnetic field intensity (Gauss)';
c.Label.FontSize = 12;
hold on

% add world map outline (just to contextualize locations)
load coastlines coastlat coastlon
geobasemap darkwater
geoplot(coastlat, coastlon, 'k-', 'HandleVisibility', 'off')
colormap parula

% find start index of l-shell contours
contour_s = find(strcmp(data{:,2}, 'L-shell contours'));
% find end index of l-shell contours
contour_e = find(strcmp(data{:,2}, 'Magnetic field intensity'));
contour_e = contour_e - 4;
% this is for text locations
midpt = round((contour_e - contour_s)/2) + contour_s;

 % color coding L-shell lines
l2_col = [1 0.8 1];
l3_col = [1 0.7 0.8];
textcol = [1 1 1];

% adding contour lines (and labels) for L = 2, L = 3, and equator
geoplot(data{contour_s:contour_e,3}, data{contour_s:contour_e,4}, '--', ...
    'Color', l2_col, 'LineWidth', 2, 'HandleVisibility', 'off')
text(data{midpt, 3}, data{midpt, 4}, 'L = 2', 'Color', textcol, ...
    'FontSize', 10, 'FontWeight', 'bold')

geoplot(data{contour_s:contour_e,5}, data{contour_s:contour_e,6}, '--', ...
    'Color', l3_col, 'LineWidth', 2, 'HandleVisibility', 'off')
text(data{midpt, 5}, data{midpt, 6}, 'L = 3', 'Color', textcol, ...
    'FontSize', 10, 'FontWeight', 'bold')

geoplot(data{contour_s:contour_e,7}, data{contour_s:contour_e,8}, '--', ...
    'Color', l2_col, 'LineWidth', 2, 'HandleVisibility', 'off')
text(data{midpt, 7}, data{midpt, 8}, 'L = 2', 'Color', textcol, ...
    'FontSize', 10, 'FontWeight', 'bold')

geoplot(data{contour_s:contour_e,9}, data{contour_s:contour_e,10}, '--', ...
    'Color', l3_col, 'LineWidth', 2, 'HandleVisibility', 'off')
text(data{midpt, 9}+10, data{midpt, 10}+10, 'L = 3', 'Color', textcol, ...
    'FontSize', 10, 'FontWeight', 'bold')

geoplot(data{contour_s:contour_e,11}, data{contour_s:contour_e,12}, '--', ...
    'Color', [1 1 1], 'LineWidth', 2, 'HandleVisibility', 'off')
text(data{midpt, 11}+10, data{midpt, 12}+10, 'Magnetic Equator', ...
    'Color', textcol, 'FontSize', 10, 'FontWeight', 'bold')

% add identified events
events_s = find(strcmp(data{:,2}, 'Event locations'));
events_e = find(strcmp(data{:,2}, 'L-shell contours'));
events_e = events_e - 4;
geoscatter(data{events_s:events_e,3}, data{events_s:events_e,4}, 30, ...
    'MarkerFaceColor', [1 0.1 0.1], 'MarkerEdgeColor', [0 0 0], ...
    'DisplayName', 'Identified Packets')
hold on

geolimits([-50 50], [-120 75])
legend()
end

%% Figure 4: Example Related Events
function figure_4

% reading in data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 4', ...
    'VariableNamingRule', 'preserve');

figure('Position', [100 100 800 400]) % aspect ratio

% extracting columns and relevant rows from spreadsheet
d = data.timestamp;
r = data.countrate;

plot(d, r, 'k-', 'LineWidth', 1.5)
xlabel("Time (UTC)"), ylabel("Electron count rate (#/20 ms)")
grid on
end


%% Figure 5: Lightning
function figure_5

% identified event data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 5', ...
    'VariableNamingRule', 'preserve');

% subplot labels
numberings1 = {'a)' 'b)' 'c)'};
numberings2 = {'d)' 'e)' 'f)'};
figure('WindowState', 'maximized')
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact')

subplot_id = data{:,2};

for i = 1:3
    this_figname = strcat('Col', num2str(i)); 
    this_fig_idx = find(strcmp(subplot_id, this_figname));

    % find row containing next text, or handle edge case if this is fig
    % 4c/f
    if i < 3
        next_figname = strcat('Col', num2str(i+1));
        next_fig_idx = find(strcmp(subplot_id, next_figname));
        next_fig_idx = next_fig_idx - 4;
    else
        next_fig_idx = length(subplot_id);
    end

    % take out data for plot
    d = data.hilt_t(this_fig_idx:next_fig_idx);
    r = data.hilt_rate(this_fig_idx:next_fig_idx);
    lt = data.lightning_t(this_fig_idx:next_fig_idx);
    pkcur = data.peak_current(this_fig_idx:next_fig_idx);
    londiff = data.lon_diff(this_fig_idx:next_fig_idx);
    lat = data.lightning_lat(this_fig_idx:next_fig_idx);
    lon = data.lightning_lon(this_fig_idx:next_fig_idx);

    nexttile(i)
    % First, plot upper plot: time
    % on left axis, plot adjusted electron flux
    yyaxis left
    % plot(d, adj, 'k', 'LineWidth', 1, 'DisplayName', 'Adj. electron count rate')
    plot(d, r, 'k', 'LineWidth', 1, 'DisplayName', 'Electron count rate')
    hold on
    xlabel('Time (UTC)')
    if i == 1
        ylabel('Count (#/20 ms)')
    end

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

    % colorbar is 'longiff'/25 to standardize colorbar across plots
    scatter(lt, pkcur, 30, londiff/25, 'filled', 'DisplayName', 'Lightning Discharges');
    ylim([0 175]) % to standardize and show all instances have amperage above 100 kA
    if i == 3
    ylabel('|Peak Current| (kA)')
    end
    
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

    xlim([min(d) max(d)])
    tx = data.hilt_t(this_fig_idx) + seconds(0.5);
    text(tx, 160, numberings1{i}, 'FontSize', 20) % subplot label

    if i == 3 % only put legend on last plot
        legend()
    end


    % Next, bottom plot: space
    nexttile(i+3)
    % put events
    geoscatter(data.event_lat(this_fig_idx), data.event_lon(this_fig_idx), ...
        'red', 'filled', 'DisplayName', 'Observed Event')
    hold on
    geobasemap colorterrain
    % lightning data
    geoscatter(lat, lon, 'yellow', 'filled', 'DisplayName', 'Lightning Discharges')
    % magnetic conjugate point
    geoscatter(data.conjlat(this_fig_idx), data.conjlon(this_fig_idx), ...
        'blue', 'filled', 'DisplayName', 'Conjugate Mirror Point')

    if i == 3 % only put legend on last plot
        legend('Location', 'east')
    end

    % subplot label
    text(45, -160, numberings2{i}, 'Color', 'white', 'FontWeight', 'bold' , ...
        'FontSize', 20)

    geolimits([-60 60], [-140 -20])
end
end

%% Figure 6: DST
function figure_6

% reading in data
data = readtable('Source Data.xlsx', 'Sheet', 'Figure 6', ...
    'VariableNamingRule', 'preserve'); 

% find storm times (Dst < -50 nT)
TF = islocalmin(data.dst);
storm_times = find(TF == 1);
storm_times = storm_times(data.dst(storm_times) < -50);
dt = diff(data.omni_t(storm_times))*365; % time between storms

figure('WindowState', 'maximized')
tiledlayout(2, 3)
nexttile([1 3])

% scatterplot of log10 count rate, binned in l-shell

% Modeled after Li et al 2017, but that was using PET data and we got a
% comment from a reviewer that (rightly) pointed out that this would
% generate a completely different scatterplot using HILT

% There was a problem with the original data, since I used every 5 days of
% data to get count rate so as not to download and process three thousand 
% files. this meant that some bins in time were not very full or had only
% one day's worth of data. This chunk of code attempts to find obvious
% discontinuities and replace them with data from the previous L-shell bin
lows = find(data.counts < 100);
for i = 1:length(lows)
    current_l = data.lshell(lows(i));
    last_value_in_l_idx = find(data.lshell == current_l & data.sampex_t ...
        < data.sampex_t(lows(i)) & ~isnan(data.counts));
    if length(last_value_in_l_idx) > 1
        last_value_in_l = data.counts(last_value_in_l_idx(end));
    else
        last_value_in_l = data.counts(last_value_in_l_idx);
    end
    if last_value_in_l/data.counts(lows(i)) > 40
        data.counts(lows(i)) = last_value_in_l;
    end
end

scatter(data.sampex_t, data.lshell, 100, ...
    log10(data.counts), 'square', 'filled', ...
    'HandleVisibility', 'off');

% plot microburst events
yyaxis left
scatter(decyear(data.event_dates), data.event_L, 50, 'MarkerFaceColor', [1 1 1], ...
    'MarkerEdgeColor', [0 0 0], 'DisplayName', 'Microburst events') 

% subplot label
text(1998.1, 4, "a)", "BackgroundColor", "white", "fontsize", 20)
ylabel('L-shell')
xlabel('Decimal year')
ylim([1.5 4])

% plot timeseries of dst
yyaxis right
plot(data.omni_t, data.dst, '-', 'Color', [0 0 0], 'LineWidth', 1, ...
    'DisplayName', 'Dst index'); % all dst
hold on
plot(data.omni_t(storm_times), data.dst(storm_times), "pentagram", ...
    "MarkerSize", 10, "MarkerFaceColor", "white",  "MarkerEdgeColor", ...
    [0 0 0.8], "DisplayName", "Geomagnetic storms", "linewidth", 1.5) % storm times
xlim([min(data.sampex_t) max(data.sampex_t)])
ylim([-250 50]);
legend()
ylabel('Dst Index')

c = colorbar('Ticks', [1 1.5 2 2.5 3 3.5 4 4.5], 'Limits', [1 5]);
colormap('parula')
c.Label.String = 'Log_{10} >1 MeV Electron Counts/20ms';
c.Label.FontSize = 11;

% Time elapsed figure
nexttile(4, [1 3])
bins = 0:2:70; % bins of days elapsed
histogram(dt, [bins Inf], "DisplayName", "Delay between storms", ...
    "FaceColor", "blue")
hold on
histogram(data.delta_t, bins, "DisplayName", ...
    "Delay between microbursts & storms", "FaceColor", "red")
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