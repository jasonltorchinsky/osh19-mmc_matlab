function eof_plot_exp(eof_params)

% Get path to output data.
out_path       = eof_params.out_path;
exp_path       = fullfile(out_path, eof_params.exp_name);
component_path = fullfile(exp_path, eof_params.component_name);
eof_path       = fullfile(exp_path, 'eofs');

addpath(out_path);
addpath(exp_path);
addpath(component_path);
addpath(eof_path);

% Create directory to hold plots
plot_path = fullfile(exp_path, 'plots');
if ~(isfolder(plot_path))
    mkdir(plot_path);
end

addpath(plot_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

t = 0:out_freq:sim_days;

% Read in EOFs, reconstructed MJOs
eof_file_name = 'eofs.nc';
eof_file = fullfile(eof_path, eof_file_name);

exp1 = ncread(eof_file, 'exp1');
exp2 = ncread(eof_file, 'exp2');

% Get PDFs of expansion coefficients
exp1_hist = histogram(exp1);
exp1_pdf  = histcounts(exp1, 'Normalization', 'probability');
exp1_bin_centers = exp1_hist.BinEdges + (exp1_hist.BinWidth / 2);

exp2_hist = histogram(exp2, 25);
exp2_pdf  = histcounts(exp2, 25, 'Normalization', 'probability');
exp2_bin_centers = exp2_hist.BinEdges + (exp2_hist.BinWidth / 2);

% Create tiled layout 
tlo = tiledlayout(1,3); % Outer layout

tlo_exp_tseries = tiledlayout(tlo, 1, 1);
tlo_exp_tseries.Layout.Tile = 1;
tlo_exp_tseries.Layout.TileSpan = [1, 2];

tlo_exp_pdf = tiledlayout(tlo, 1, 1);
tlo_exp_pdf.Layout.Tile = 3;
tlo_exp_pdf.Layout.TileSpan = [1, 1];

% Create plots of expansion coefficients

% Expansion coefficient time series
h(1) = nexttile(tlo_exp_tseries);

hold on;

plot(t, exp1, 'k-', ...
    'DisplayName', 'EOF 1');
plot(t, exp2, 'k--', ...
    'DisplayName', 'EOF 2');

legend();

xlabel('Time (d)');
ylabel('Expansion Coefficient');

hold off;

% Expansion coefficient PDF
h(2) = nexttile(tlo_exp_pdf);

hold on;


plot(exp1_bin_centers(1:end-1), exp1_pdf, 'k-', ...
    'DisplayName', 'EOF 1');
plot(exp2_bin_centers(1:end-1), exp2_pdf, 'k--', ...
    'DisplayName', 'EOF 2');

legend();

xlabel('Expansion Coefficient');
ylabel('Probability');

hold off;



%~ Figure size
set(gcf, 'Units', 'inches');
figWidth  = 7.5; % Figure width in inches.
figHeight = 4.5; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight],...
    'PaperOrientation', 'portrait');

% Save plot.
file_name = strcat([eof_params.component_name, '_exp.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters', '-fillpage');

end

