% clc
% clear
close all

data = readtable('lshelldata.csv', 'VariableNamingRule', 'preserve');
load pitchangledata.mat
alpha_eq = asind(sind(abs(data.maglat)).*sqrt(beq./b));

c = 3e+8; % m/s
a = 6371000; % m, radius of Earth
m_e = 9.11e-31; % electron mass, kg
E_rest = m_e*c^2; % MeV
E = 1.6021766e-13; % 1 MeV in J
total_E = E + E_rest; % J
v_e = c*sqrt(1-(m_e*c^2/total_E)^2); % electron velocity, m/s

t_y = (1.3802 - 0.3198*(sind(15) + sqrt(sind(15))));
myty = 1.3802 - 0.3198*(sind(alpha_eq(2)) + sqrt(sind(alpha_eq(2))));
L = 2;
tb0 = 0.22; % s

v = 4*a*data.L(2)*myty/0.215;
ke = m_e*c^2/sqrt(1-(v/c)^2) - m_e*c^2;
ke_mev = ke/E;

% 
% st = zeros(size(data.time));
% et = zeros(size(data.time));
dt = cell(size(data.time));
loi = [3 5 9 15 16 21 24 28 30 32 37 39 40];
for i = 12
[~, ~, out_dats, d, r, sti, eti] = plotEpoch(data.date(loi(i)), data.time(loi(i)));
st = sti;
et = eti;
pk_locs = d(out_dats); % peak locations, datetime
dt{loi(i)} = seconds(diff(pk_locs));
pk_vals = r(out_dats); % peak locations, rate

% d = d(st-10:et+50); % narrowing down to portion of interest
% r = r(st-10:et+50); % standardizing amount of time makes plotting clearer
% figure()
% plot(d, r, 'k-', 'LineWidth', 1, 'DisplayName', 'Raw')
% hold on
% plot(pk_locs, pk_vals, 'r*', 'DisplayName', 'Identified Peak')
% legend()
% xlabel('Time (UTC)'), ylabel('Count Rate (#/20 ms)')
% hold off
end


% % locations = [.18 .4 .62 .84 1.06 1.29 1.5 1.7 1.92 2.12]; % for 11/12/98
% locations = [.12 .33 .56 .78 .98 1.2 1.4]; % for 7/18/2000
% dl = diff(locations);
% close all
% for i = 1:10
% dl = dt{loi(i)};
% dl = dl/.02;
% dl = round(dl);
% dl = dl*0.02;
% figure()
% min_bound = min(dl)-0.01;
% max_bound = max(dl)+0.01;
% histogram(dl, min_bound:0.001:max_bound)
% myty = 1.3802 - 0.3198*(sind(alpha_eq(loi(i))) + sqrt(sind(alpha_eq(loi(i)))));
% per_pred = 4*a*(data.L(loi(i))/v_e)*myty;
% xline(per_pred, 'r--', "LineWidth", 2)
% xlabel("Spacing between peaks (s)"), ylabel("Count")
% text(per_pred + .001, 3, "Predicted period (dipole)", "Color", "red")
% xlim([min_bound max_bound])
% % ylim([0 4])
% title(['Bounce period analysis for event on ', char(data.date(loi(i)))])
% end



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
