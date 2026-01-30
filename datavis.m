clear
addpath '/Users/henrysun_1/Desktop/Duke/PhD 2025-2026/Class/ENV876/class/project/data'
suffixes = {'2022_04_to_2022_08', '2022_08_to_2023_04', ...
            '2023_04_to_2023_08', '2023_08_to_2024_04', ...
            '2024_04_to_2024_08'};

for i = 1:length(suffixes)
    load(['sbe56_', suffixes{i}, '.mat']);
end

tem = cell(1, 6);
for s = 1:6
    temp_column = [];
    for i = 1:length(suffixes)
        varName = ['sbe56_', suffixes{i}, '_Rtemp'];
        data = eval(varName); 
        temp_column = [temp_column; data(:, s)];
    end
    tem{s} = temp_column;
end

[tem1, tem2, tem3, tem4, tem5, tem6] = deal(tem{:});
all_data = [tem1, tem2, tem3, tem4, tem5, tem6];
% nanCols = any(isnan(all_data), 1); 
% nanCols;

ptime_full = [];
for i = 1:length(suffixes)
    varNameTime = ['sbe56_', suffixes{i}, '_time'];    
    data_time = eval(varNameTime); 
    ptime_full = [ptime_full; data_time(:, 1)];
end
ptime = ptime_full;

% save matrix with ptime and tem1-6 as a csv
matrix= [ptime, tem1, tem2, tem3, tem4, tem5, tem6];
% writematrix(matrix, "raw.csv")

% save ptime and tem1-6 as .mat file
% save('data/raw.mat', 'ptime', 'tem1', 'tem2', 'tem3', 'tem4', 'tem5', 'tem6');


%%
% calculate min median and max for each site
minValues = [min(tem1), min(tem2), min(tem3), min(tem4), min(tem5), min(tem6)];
medianValues = [median(tem1), median(tem2), median(tem3), median(tem4), median(tem5), median(tem6)];
maxValues = [max(tem1), max(tem2), max(tem3), max(tem4), max(tem5), max(tem6)];
x = 1:6;

% plot min median and max for each site as dots
% make x axis a discrete scale from 1-6
figure;
hold on;
xlim([0.5, 6.5])
xticks([1, 2, 3, 4, 5, 6])
plot(x, minValues, 'g-o', 'DisplayName', 'Min');
plot(x, medianValues, 'b-o', 'DisplayName', 'Median');
plot(x, maxValues, 'm-o', 'DisplayName', 'Max');
ylabel('Temperature (°C)')
xlabel('Restoration Site')
hold off;
legend show;

%%
% make boxplot of temperature at each site
data = [tem1; tem2; tem3; tem4; tem5; tem6];            
group = repelem(1:6, numel(tem1));  

figure
hold on;
b = boxchart(group, data, 'BoxWidth', 0.7);
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
plot(x, minValues, 'g-o', 'DisplayName', 'Min');
plot(x, medianValues, 'b-o', 'DisplayName', 'Median');
plot(x, maxValues, 'm-o', 'DisplayName', 'Max');
xlim([0.5, 6.5])
set(gca, 'XTick', 1:6, 'XTickLabel', {'1','2','3','4','5','6'})
ylabel('Temperature (°C)')
xlabel('Restoration Site')

%% 
% make a matrix of temperature values in tem1 to tem6
temperatureMatrix = [tem1, tem2, tem3, tem4, tem5, tem6];
rho = corr(temperatureMatrix);

% create the heatmap and color limits
h = heatmap(rho);
h.ColorLimits = [0.975 1];
customColors = [
    1, 0, 0;   % Red
    1, 1, 1;   % White
    0, 1, 0    % Blue
];
h.Colormap = interp1(linspace(0,1,3), customColors, linspace(0,1,256));
h.XLabel = 'Restoration Site';
h.YLabel = 'Restoration Site';
h.GridVisible = 'off';

%% 
% alternate method
image(rho,'CDataMapping','scaled');
colorbar

%% 
ptime = matrix(:, 1);
R6 = matrix(:, 7); 

ptime_dt = datetime(ptime, 'ConvertFrom', 'datenum');
startDate = datetime('2023-01-01');
endDate = datetime('2023-12-31');
timeIndex = ptime_dt >= startDate & ptime_dt <= endDate;

R6_filtered = R6(timeIndex);
ptime_filtered = ptime_dt(timeIndex);

dates = dateshift(ptime_filtered, 'start', 'day');
hours = hour(ptime_filtered);

%%
tbl = table(dates, R6_filtered);
dayStats = groupsummary(tbl, 'dates', {'mean', 'max'});

[~, minAvgIdx] = min(dayStats.mean_R6_filtered);
[~, maxAvgIdx] = max(dayStats.mean_R6_filtered);
[~, absMaxIdx] = max(dayStats.max_R6_filtered);

targetDates = [dayStats.dates(minAvgIdx), ...
               dayStats.dates(maxAvgIdx), ...
               dayStats.dates(absMaxIdx)];

figure; 
hold on;
colors = ['b', 'r', 'g'];
labels = {'Lowest Avg Temp - 1/15', 'Highest Avg Temp - 8/11', 'Max Absolute Temp - 8/10'};

for i = 1:3
    dayIdx = (dates == targetDates(i));
    [sortedHours, sIdx] = sort(hours(dayIdx));
    tempValues = R6_filtered(dayIdx);
    plot(sortedHours, tempValues(sIdx), colors(i));
end

xlabel('Hour of Day');
ylabel('Temperature (°C)');
legend(labels, 'Location', 'northwest');
xticks(0:2:24);

uniqueDays = unique(datesOnly);
figure; 
hold on;

for d = 1:length(uniqueDays)
    dayIdx = (datesOnly == uniqueDays(d));
    [sortedHours, sIdx] = sort(hours(dayIdx));
    tempValues = R6_filtered(dayIdx);
    plot(sortedHours, tempValues(sIdx), 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
end

colors = ['b', 'r', 'g'];
% use targetdates from above
labels = {['Lowest Avg: ', datestr(targetDates(1), 'mm/dd')], ...
          ['Highest Avg: ', datestr(targetDates(2), 'mm/dd')], ...
          ['Max Absolute: ', datestr(targetDates(3), 'mm/dd')]};

for i = 1:3
    dayIdx = (datesOnly == targetDates(i));
    [sortedHours, sIdx] = sort(hours(dayIdx));
    tempValues = R6_filtered(dayIdx);
    plot(sortedHours, tempValues(sIdx), colors(i), 'LineWidth', 2);
end

xlabel('Hour of Day');
ylabel('Temperature (°C)');
legend(labels, 'Location', 'northwest');
xticks(0:2:24);
