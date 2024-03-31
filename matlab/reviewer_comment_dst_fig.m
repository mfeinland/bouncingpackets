% clear 
% clc
close all


first_day =  '1998-01-01';
last_day = '1998-12-31';

% [t98, data98] = get_sampex_mux(first_day, last_day);

dayyearvec = dayarrayfun;
% dayyearvec = dayyearvec(23:120,:);
% load longasscomptime.mat

load mat2.mat

dates = dayyearvec(:,2) + dayyearvec(:,1)/365;
date_vec = repmat(dates, [1 61]);
l_vec = 1:0.1:7;

rate2 = data2(:,1);
l2 = data2(:,2);
% dy2 = decyear(t98);

% spaced_dy = min(dy):0.005:max(dy);
% dycell = cell(length(spaced_dy), 2);
% for i = 1:length(spaced_dy)
%     idx = find(dy <= spaced_dy(i));
%     dycell{i,1} = l(idx);
%     dycell{i,2} = rate(idx);
% end

spaced_dy2 = min(dy2):0.005:max(dy2);
dycell2 = cell(length(spaced_dy2), 2);
for i = 1:length(spaced_dy2)
    idx = find(dy2 <= spaced_dy2(i));
    dycell2{i,1} = l2(idx);
    dycell2{i,2} = rate2(idx);
end

% binned_l = zeros(length(dayyearvec), 61);
% 
% for i = 1:length(date_vec)-1
%     times_to_consider = dy >= date_vec(i,1) & dy <= date_vec(i+1,1);
%     for j = 1:61
%         lbin = 1 + 0.1*(j-1);
%         vals = find(l(times_to_consider) >= lbin & l(times_to_consider) <= (lbin + 0.1));
%         if length(vals) > 10
%             binned_l(i, j) = mean(rate(vals), 'omitnan');
%         else
%             binned_l(i, j) = NaN;
%         end
%     end
% end



binned_l_cell = zeros(length(dycell2), 61);

for i = 1:length(dycell2)
    for j = 1:61
        lbin = 1 + 0.1*(j-1);
        vals = find(dycell2{i,1} >= lbin & dycell2{i,1} <= (lbin + 0.1));
        cur_rate = dycell2{i,2};
        if length(vals) > 10
            binned_l_cell(i, j) = mean(cur_rate(vals), 'omitnan');
        else
            binned_l_cell(i, j) = NaN;
        end
    end
end


[m, col] = size(binned_l_cell);
for i = 1:m
    scatter(spaced_dy2(i)*ones(length(l_vec),1), l_vec, 70, log10(binned_l_cell(i,:)), 'square', 'filled', ...
    'HandleVisibility', 'off')
    hold on
end
colormap hsv
c = colorbar();
% 
% x = all(isnan(binned_l),2);
% date_vec = date_vec(x == 0, :);
%%% figure out how to bin by day instead of every 27 days
% unique_day_vec = unique()
% binned_l = binned_l(x == 0, :);
% [m, col] = size(binned_l);
% for i = 1:m
%     scatter(date_vec(i,:), l_vec, 70, log10(binned_l(i,:)), 'square', 'filled', ...
%     'HandleVisibility', 'off')
%     hold on
% end
% colormap hsv

% c = colorbar();
% c.Label.String = 'Log_{10} ~2 MeV Electron Flux';
% c.Label.FontSize = 11;

% % reading in data
% dst = readtable("OMNI2_H0_MRG1HR_252725.csv", 'ReadVariableNames', false);
% % contains DST information from January 1, 1997 to Dec 31, 2006. Available
% % from OMNIweb
% 
% % formatting from csv; there are a lot of extraneous headers, cell arrays,
% % strings, etc
% dst(1:75, :) = [];
% dst(end-3:end, :) = [];
% dat = cell(height(dst), 1);
% tim = table2array(dst(:,1));
% for i = 1:height(tim)
%     str = tim{i};
%     dat{i} = str(1:10);
% end
% dstval = table2array(dst(:,2));
% dat = datetime(cell2mat(dat));
% 
% motstmp = dateshift(dat, 'start', 'month', 'nearest');
% unqmos = unique(motstmp);
% dstavg = zeros(size(unqmos));
% 
% % monthly average
% for i = 1:numel(unqmos)
%     month_data = dstval(motstmp == unqmos(i));
%     dstavg(i) = mean(month_data);
% end
% 
% unqmos = decyear(unqmos);
% 
% % figure() % plotting DST and events only
% dstcol = [0.95 0.4 0.8];
% colororder([0.2 0.6 0.8; dstcol])
% 
% % this is plotting events, DST, and 2 MeV electron flux
% axis tight
% yyaxis left
% ylabel('L-shell')
% xlabel('Year')
% % ylim([1.5 3.5])
% title('MeV Electron Population, Dst Index, and Events')
% 
% yyaxis right
% % dstplt = plot(unqmos, dstavg, '-.', 'Color', dstcol, 'LineWidth', 2, ...
% %     'DisplayName', 'Dst Index');
% dstplt = plot(unqmos, dstavg, '-.', 'Color', [1 1 1], 'LineWidth', 2, ...
%     'DisplayName', 'Dst Index');
% % xlim([1999 2002])
% legend()
% ylabel('Dst Index')


%% Functions
function daysvec = dayarrayfun
% Author: ChatGPT (prompting by Max Feinland)
startYear = 1996;
startDay = 160;
endYear = 2006;
endDay = 179;
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

function answ = isleap(year)
% checks for leap year
if mod(year, 4) == 0
    answ = 1;
else
    answ = 0;
end
end

function [truncated_time, datavec] = get_sampex_mux(first_day, last_day)
first_day_d = datetime(first_day);
last_day_d = datetime(last_day);


% % Read in SAMPEX six-second data
doy1 = day(first_day_d, 'dayofyear'); % day of start (integer)
doy2 = day(last_day_d, 'dayofyear'); % day of end (integer)

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
hilt_matrix = [];
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
        D = [B B B repmat(A, [1 17]) B repmat(A, [1 51])];
        % tells textscan to look at columns 1, 2, 3, 21
        D(end) = [];

        celltab = textscan(file, D, 'Delimiter', '\t');
        data_chunk = [celltab{1}, celltab{2}, celltab{3}, celltab{4}];

        l_data = [l_data; data_chunk]; %#ok<AGROW> 

    else % data for which it's zipped
        url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/PSSet/Text/', ...
        ephem_filename{i-n+1} '.zip'];

        ephem_savepath = 'C:/Users/maxim/sampex-data/Attitude/'; % save to data folder

        ephem_pathname = fullfile(ephem_savepath, ephem_filename{i-n+1});
        if ~exist(ephem_pathname, 'file')
            unzip(url, ephem_savepath) % find and unzip data file
        end
        
        data_chunk = readmatrix(ephem_pathname, 'HeaderLines', 60); % read in data file
        col1 = data_chunk(:,1);
        col2 = data_chunk(:,2);
        col3 = data_chunk(:,3);
        col4 = data_chunk(:,21);
        data_chunk = [col1 col2 col3 col4];

        l_data = [l_data; data_chunk]; %#ok<AGROW> 
    end
end

% OK, now that ephemeris data has been successfully downloaded, let's deal
% with the HILT data.
day_list = first_day_d:days(5):last_day_d;

% for each day in the list of days spaced by 5 each time:
for i = 1:length(day_list)
    dat = day_list(i);
    doy = day(dat, 'dayofyear'); %#ok<TNDAY>
    y = year(dat); %#ok<TYEAR>
     % makes sure filename has leading zeros if necessary
    if strlength(num2str(doy)) == 3
        filename = ['hhrr', num2str(y), num2str(doy), '.txt'];
    elseif strlength(num2str(doy)) == 2
        filename = ['hhrr', num2str(y), '0', num2str(doy), '.txt'];
    elseif strlength(num2str(doy)) == 1
        filename = ['hhrr', num2str(y), '00', num2str(doy), '.txt'];
    end


    savepath = 'C:/Users/maxim/sampex-data/HILT/State4/'; % save to data folder

    pathname = fullfile(savepath, filename);
    if ~exist(pathname, 'file')
        url = ['https://izw1.caltech.edu/sampex/DataCenter/DATA/HILThires/State4/',...
            filename, '.zip'];
        unzip(url, savepath) % find and unzip data file
    end
    data = readmatrix(pathname); % read in data file
    fprintf("Downloaded %s (%i of %i)\n", filename, i, length(day_list))

    t = data(:,1); % ok TECHNICALLY this is every 100 ms instead of every 
    % 20 but it doesn't matter for this purpose!
    hilt_time = dat + seconds(t);
    col1 = year(hilt_time);
    col2 = day(hilt_time, 'dayofyear');
    col3 = t;
    rates = [data(:,2:5) data(:,7)];
    rate = sum(rates, 2); % sum across rows
    hilt_matrix = [hilt_matrix; col1 col2 col3 rate]; %#ok<AGROW>
end

    % Now we should have the hilt data. Let's find the L-shells that
    % correspond to it.
    att_time = datetime(l_data(:,1), 1, 1) + days(l_data(:,2)) + seconds(l_data(:,3));

    hilt_time = datetime(hilt_matrix(:,1), 1, 1) + days(hilt_matrix(:,2)) + seconds(hilt_matrix(:,3));

    fprintf("Finding matching times")
    [~, matching_times_hilt, matching_times_att] = intersect(hilt_time, att_time);
    truncated_rate = hilt_matrix(matching_times_hilt, 4);
    truncated_ls = l_data(matching_times_att, 4);
    truncated_time = hilt_time(matching_times_hilt);


    bad_inds = truncated_rate < 0 | truncated_rate > 1e+5;

    datavec = [truncated_rate truncated_ls];
    datavec(bad_inds,:) = [];
    % datavec(315000:end, :) = [];
end
