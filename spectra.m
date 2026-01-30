clear
addpath '/Users/henrysun_1/Desktop/Duke/PhD 2025-2026/Class/ENV876/class/project/data'
load raw.mat

% ptime is a vector in matlab time
matrix = [ptime tem1 tem6];

%% 
% calculate sampling frequency from diff in ptime - ptime in matlab time
target_interval_days = 3 / 24; % 3 hours in days
fs = 1 / target_interval_days;
t_uniform = (datetime(2022,1,1):hours(3):datetime(2024,12,31,21,0,0))';
t_uniform_num = datenum(t_uniform);

tem1_interp = interp1(matrix(:,1), matrix(:,2), t_uniform_num, 'linear');
tem6_interp = interp1(matrix(:,1), matrix(:,3), t_uniform_num, 'linear');
tem1_interp(isnan(tem1_interp)) = mean(tem1_interp, 'omitnan');
tem6_interp(isnan(tem6_interp)) = mean(tem6_interp, 'omitnan');

%%
window_size = floor(121.6 * fs);  % 4 months in days
overlap = floor(window_size / 2); % 50% overlap

[pxx1, f] = pwelch(detrend(tem1_interp), hamming(window_size), overlap, window_size, fs);
[pxx6, ~] = pwelch(detrend(tem6_interp), hamming(window_size), overlap, window_size, fs);

%%
bands = [0.00185, 0.0111;   % Annual
         0.0119, 0.143;    % Seasonal
         0.657, 1.156; % Diurnal
         1.7, 2.353]; % Semidiurnal 
band_names = {'Annual', 'Seasonal', 'Diurnal', 'Semidiurnal'};
results = table();
for i = 1:size(bands, 1)
    % Find the indices that fall within the current band
    idx = (f >= bands(i,1) & f <= bands(i,2));    
    if any(idx)
        results.Band{i,1} = band_names{i};
        results.Tem1_Var(i,1) = trapz(f(idx), pxx1(idx), 1);
        results.Tem6_Var(i,1) = trapz(f(idx), pxx6(idx), 1);
    else
        results.Band{i,1} = band_names{i};
        results.Tem1_Var(i,1) = 0;
        results.Tem6_Var(i,1) = 0;
    end
end

disp(results);
%% 
% plot PSDs
figure;
subplot(211)
loglog(f, pxx6, 'r', 'LineWidth', 1.5);
xlim([0.01 4]);
subplot(212)
loglog(f, pxx1, 'b', 'LineWidth', 1.5); 
xlim([0.01 4]);

xlabel('Frequency (cpd)');
ylabel('Power Spectral Density (°C^{2}/cpd)');

patch([0.657 1.156 1.156 0.657], [ylim 1e-5 1e-5], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([1.7 2.353 1.7 2.353], [ylim 1e-5 1e-5], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

%%
sensors = {tem6_interp, tem1_interp}; 
psds = {pxx6, pxx1};
lineColors = {'#E7B800', '#00AFBB'};
sensorNames = {'Shallow', 'Deep'};
varColNames = {'Tem6_Var', 'Tem1_Var'};
subplotLabels = {'a', 'b'}; 

figure('Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);
t_layout = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:2
    data = sensors{s};
    pxx = psds{s};
    
    alpha = 0.05;
    rho = corr(data(1:end-1), data(2:end)); 
    P_red = mean(pxx) * (1 - rho^2) ./ (1 + rho^2 - 2*rho*cos(2*pi*f/fs));    
    dof = 2 * (length(data) / window_size) * (8/3); 
    chi_val = chi2inv(1-alpha, dof) / dof;
    UCL = P_red * chi_val;
    ax(s) = nexttile; 
    loglog(f, pxx, 'Color',lineColors{s}, 'LineWidth', 1.2); hold on;
    
    % plotting significant peaks
    [pks, locs] = findpeaks(pxx, f, 'MinPeakHeight', 0);
    sig_idx = pks > interp1(f, UCL, locs);
    plot(locs(sig_idx), pks(sig_idx), 'kx', 'MarkerSize', 6, 'LineWidth', 1.2);
    relevant_bands = bands(2:4, :); 
    relevant_names = band_names(2:4);
    
    % plot patches
    for b = 1:3
        bandVar = results.(varColNames{s})(b+1); 
        % shade patch
        p = patch([relevant_bands(b,1) relevant_bands(b,2) relevant_bands(b,2) relevant_bands(b,1)], ...
              [1e-10 1e-10 1e10 1e10], 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
        uistack(p, 'bottom');
        %label text
        txtStr = sprintf('%s\nVar: %.4f', relevant_names{b}, bandVar);
        text(mean(relevant_bands(b,:)), 20, txtStr, ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'Margin', 1);
    end
    text(0.02, 1.1, subplotLabels{s}, 'Units', 'normalized', ...
        'FontSize', 14, 'FontWeight', 'Bold','VerticalAlignment', 'top', 'Clipping','off');
    title(sensorNames{s});
    xlim([0.01 4]);
    ylim([1e-4 100]); 
end

xlabel(t_layout, 'Frequency (cpd)');
ylabel(t_layout, 'Power Spectral Density (°C^{2}/cpd)');
linkaxes(ax, 'xy');

%%
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 12, 6]);
exportgraphics(gcf, 'figs/psd_plot.png', 'Resolution', 300);