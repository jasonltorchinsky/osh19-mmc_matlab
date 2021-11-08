function cmg14_plot_evo(params)

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

% Read the input parameters to find out how long the run is and how
% frequently output was created.
params_file = fullfile(component_path, 'params.nc');

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get arrays for each state variable
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out_idxs = 0:floor(sim_days/out_freq);
u_1  = zeros(size(out_idxs));
u_2  = zeros(size(out_idxs));
v    = zeros(size(out_idxs));
w_u  = zeros(size(out_idxs));

for out_idx = out_idxs
    state_file_name = strcat(['state_', num2str(out_idx,'%04u'),'.nc']);
    state_file = fullfile(component_path, state_file_name);
    u_1(out_idx + 1) = ncread(state_file, 'u_1');
    u_2(out_idx + 1) = ncread(state_file, 'u_2');
    v(out_idx + 1)   = ncread(state_file, 'v');
    w_u(out_idx + 1) = ncread(state_file, 'w_u');
end

% Convert seconds to days
days_to_secs = 3600 * 24;
w_u = w_u * days_to_secs;

% Get PDFs of state variables
nbins = 16;
u1_hist = histogram(u_1, nbins);
u1_pdf  = histcounts(u_1, nbins) / numel(u_1) ./ u1_hist.BinWidth;
u1_bin_centers = u1_hist.BinEdges + (u1_hist.BinWidth / 2);

u2_hist = histogram(u_2, nbins);
u2_pdf  = histcounts(u_2, nbins) / numel(u_2) ./ u2_hist.BinWidth;
u2_bin_centers = u2_hist.BinEdges + (u2_hist.BinWidth / 2);

v_hist = histogram(v, nbins);
v_pdf  = histcounts(v, nbins) / numel(v) ./ v_hist.BinWidth;
v_bin_centers = v_hist.BinEdges + (v_hist.BinWidth / 2);

wu_hist = histogram(w_u, nbins);
wu_pdf  = histcounts(w_u, nbins) / numel(w_u) ./ wu_hist.BinWidth;
wu_bin_centers = wu_hist.BinEdges + (wu_hist.BinWidth / 2);

% Gaussian fit for each
x = linspace(-8, 8, 1600);

u1_mean = mean(u_1);
u1_var  = var(u_1);
u1_gauss = 1/sqrt(2 * pi * u1_var) * exp(-1/2 * (x - u1_mean).^2 / u1_var);

u2_mean = mean(u_2);
u2_var  = var(u_2);
u2_gauss = 1/sqrt(2 * pi * u2_var) * exp(-1/2 * (x - u2_mean).^2 / u2_var);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Make plots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Create master tiled layout
tlo = tiledlayout(3, 2);

tlo_mjo = tiledlayout(tlo, 1, 2);
tlo_mjo.Layout.Tile = 1;
tlo_mjo.Layout.TileSpan = [1, 2];

tlo_v = tiledlayout(tlo, 1, 2);
tlo_v.Layout.Tile = 3;
tlo_v.Layout.TileSpan = [1, 2];

tlo_wu = tiledlayout(tlo, 1, 2);
tlo_wu.Layout.Tile = 5;
tlo_wu.Layout.TileSpan = [1, 2];

% First Row - MJO indices

h(1,1) = nexttile(tlo_mjo, 1);

hold on;
plot(h(1,1), u1_bin_centers(1:end-1), u1_pdf, 'k-');
plot(h(1,1), x, u1_gauss, 'r-');
plot(h(1,1), u2_bin_centers(1:end-1), u2_pdf, 'k--');
plot(h(1,1), x, u2_gauss, 'r--');

xlim(h(1,1), [min([u1_bin_centers, u2_bin_centers], [], 'all'), ...
    max([u1_bin_centers, u2_bin_centers], [], 'all')]);
ylim(h(1,1), [min([u1_pdf, u2_pdf], [], 'all'), ...
    max([u1_pdf, u2_pdf], [], 'all')]);

ylabel(h(1,1), 'Probability')

legend(h(1,1), 'u_{1}', '', 'u_{2}');

hold off;

h(1,2) = nexttile(tlo_mjo, 2);

semilogy(h(1,2), u1_bin_centers(1:end-1), u1_pdf, 'k-');
hold on;
semilogy(h(1,2), x, u1_gauss, 'r-');
semilogy(h(1,2), u2_bin_centers(1:end-1), u2_pdf, 'k--');
semilogy(h(1,2), x, u2_gauss, 'r--');

xlim(h(1,2), [min([u1_bin_centers, u2_bin_centers], [], 'all'), ...
    max([u1_bin_centers, u2_bin_centers], [], 'all')]);
ylim(h(1,2), [min([u1_pdf, u2_pdf], [], 'all'), ...
    max([u1_pdf, u2_pdf], [], 'all')]);

ylabel(h(1,2), 'log(Probability)')

legend(h(1,2), 'u_{1}', '', 'u_{2}');

hold off;

% Second row - stochastic damping

h(2,1) = nexttile(tlo_v, 1);

hold on;
plot(h(2,1), v_bin_centers(1:end-1), v_pdf, 'k-');

ylabel(h(2,1), 'Probability')

legend(h(2,1), 'v');

hold off;

h(2,2) = nexttile(tlo_v, 2);

semilogy(h(2,2), v_bin_centers(1:end-1), v_pdf, 'k-');
hold on;

ylabel(h(2,2), 'log(Probability)')

legend(h(2,2), 'v');

hold off;

% Third plot - stochastic phase

h(3,1) = nexttile(tlo_wu, 1);

plot(h(3,1), wu_bin_centers(1:end-1), wu_pdf, 'k-');
hold on;

ylabel(h(3,1), 'Probability')

legend(h(3,1), '\omega_{u}');

hold off;

h(3,2) = nexttile(tlo_wu, 2);

hold on;
semilogy(h(3,2), wu_bin_centers(1:end-1), wu_pdf, 'k-');

ylabel(h(3,2), 'log(Probability)')

legend(h(3,2), '\omega_{u}');

hold off;


% Common across all tiles


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% General figure styling.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Page and figure size
set(gcf, 'PaperUnits', 'inches');
figWidth  = 6; % Figure width in inches.
figHeight = 6; % Figure height in inches.
set(gcf,...
    'PaperPosition', [0, 0, figWidth, figHeight],...
    'PaperSize', [figWidth, figHeight])



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plot.

file_name = strcat([params.component_name, '_pdfs.pdf']);
plot_file = fullfile(plot_path, file_name);
print(plot_file, '-dpdf', '-painters');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end