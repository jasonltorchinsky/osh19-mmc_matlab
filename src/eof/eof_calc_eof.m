function eofs = eof_calc_eof(params)

eofs = struct();

% Get path to output data.
out_path       = params.out_path;
exp_path       = fullfile(out_path, params.exp_name);
component_path = fullfile(exp_path, params.component_name);

addpath(out_path);
addpath(exp_path);
addpath(component_path);

% Read in run parameters
params_file = fullfile(component_path, 'params.nc');

nx       = ncread(params_file, 'nx');

sim_days = ncread(params_file, 'sim_days');
out_freq = ncread(params_file, 'out_freq');

nt = floor(sim_days/out_freq) + 1;

% Get vertical, meridional projections of u, q 
[u_proj, u_std] = eof_proj_u(params);
[q_proj, q_std] = eof_proj_q(params);

eofs.u_std = u_std;
eofs.q_std = q_std;

% Perform the EOF decompostion via eigenvalue decomposition
comb_proj = [u_proj q_proj]; % Combined projection matrix
cov = transpose(comb_proj) * comb_proj;
[eof_mtx, L] = eig(cov);

% Get non-dimensional EOFs, expansion coefficients of combined system
eofs.eof1 = transpose(eof_mtx(:, end));
eofs.eof2 = transpose(eof_mtx(:, end-1));

exp_coeffs = comb_proj * eof_mtx(:, [end end-1]);
eofs.exp1 = transpose(exp_coeffs(1:nt, end-1));
eofs.exp2 = transpose(exp_coeffs(1:nt, end));

% Get the first two EOFs, and expansion coefficients
eofs.u_eof1 = eofs.eof1(1:nx) * u_std;
eofs.u_eof2 = eofs.eof2(1:nx) * u_std;
eofs.q_eof1 = eofs.eof1(nx+1:end) * q_std;
eofs.q_eof2 = eofs.eof2(nx+1:end) * q_std;

eofs.scf = sort(diag(L) / trace(L), 'descend'); % Scaled covariance fraction

eofs = eof_recon_mjo(params, eofs);

end

