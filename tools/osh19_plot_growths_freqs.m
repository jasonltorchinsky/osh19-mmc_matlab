function osh19_plot_growths_freqs(params, mode)

% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Create directory to hold plots
plot_path = fullfile(exp_path, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in growthrates and frequencies
growths_freqs_file = fullfile(component_path, 'growths_freqs.nc');

wavenums     = ncread(growths_freqs_file, 'wavenums');

% Get the mode index for the mode we want to plot
modes         = ncread(growths_freqs_file, 'modes');
[~, mode_idx] = min(abs(modes - mode));

% Read in all of the values, then trim
growth_rates = ncread(growths_freqs_file, 'growth_rates');
growth_rates = growth_rates(:, mode_idx);

freqs = ncread(growths_freqs_file, 'freqs');
freqs = freqs(:, mode_idx);

% Example periods
ext_wavenums = -2:max(wavenums)+1;
per3  = 1/3 * ones(size(ext_wavenums));
per6  = 1/6 * ones(size(ext_wavenums));
per30 = 1/30 * ones(size(ext_wavenums));
per60 = 1/60 * ones(size(ext_wavenums));
per90 = 1/90 * ones(size(ext_wavenums));

pers = [per3
    per6
    per30
    per60
    per90];


% Create tiled layout of plots of evolution
tlo = tiledlayout(1,2); % Outer layout

% Plot growth rates
h(1) = nexttile(tlo);

hold on;

scatter(wavenums, growth_rates, ...
    'k', 'x', ...
    'LineWidth', 1);

for per_idx = 1:size(pers,1)
    per = pers(per_idx, :);
    plot(ext_wavenums, per, ...
        'k--', ...
        'Marker', 'none');
    
    str = strcat([num2str(round(1/per(1), 0)), ' Days']);
    text(min(ext_wavenums)+0.2, per(1), str, ...
        'VerticalAlignment', 'bottom', ...
        'Margin', 1);
end

% Vertical line at wavenumber 0
plot([0,0], [0,10], 'k--');

hold off;

ylabel('Growth Rate (d^{-1})');

% Plot frquencies
h(2) = nexttile(tlo);

hold on;

scatter(wavenums, freqs, ...
    'k', 'x', ...
    'LineWidth', 1);

for per_idx = 1:size(pers,1)
    per = pers(per_idx, :);
    plot(ext_wavenums, per, ...
        'k--', ...
        'Marker', 'none');
    
    str = strcat([num2str(round(1/per(1), 0)), ' Days']);
    text(min(ext_wavenums)+0.2, per(1), str, ...
        'VerticalAlignment', 'bottom', ...
        'Margin', 1);
end

% Vertical line at wavenumber 0
plot([0,0], [0,10], 'k--');

hold off;

ylabel('Frquency (d^{-1})');

% Common elements across all plots
xlabel(tlo, 'Zonal Wavenumber');

for plt = h
    xticks(plt, ext_wavenums);
    
    xlim(plt, [min(ext_wavenums), max(ext_wavenums)]);
    ylim(plt, [0, 0.05]);
end

    
    

%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 6.5; % Figure width in inches.
figHeight = 4; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');

% Save plot.
file_name = strcat([params.component_name, '_growths_freqs.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');

end

