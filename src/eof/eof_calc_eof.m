function eofs = eof_calc_eof(params)

eofs = struct();

% Get vertical, meridional projections of u, q 
u_proj = eof_proj_u(params);
q_proj = eof_proj_q(params);

% Take out time-series means
u_proj = detrend(u_proj, 0);
q_proj = detrend(q_proj, 0);

% Normalize u_proj, q_proj by their variance
u_proj_var = std(u_proj, 0, 'all');
u_proj = u_proj / u_proj_var;

q_proj_var = std(q_proj, 0, 'all');
q_proj = q_proj / q_proj_var;

eofs.u_var = u_proj_var;
eofs.q_var = q_proj_var;

% Perform the EOF decompostion via eigenvalue decomposition
comb_proj = [u_proj q_proj]; % Combined projection matrix
cov = transpose(comb_proj) * comb_proj;
[eof_mtx, L] = eig(cov);

% Get "raw" EOFs, xpansion coefficients of combined system
eofs.raw_eof1 = transpose(eof_mtx(:, end));
eofs.raw_eof2 = transpose(eof_mtx(:, end-1));

[nt, nx] = size(u_proj);

exp_coeffs = comb_proj * eof_mtx(:, [end end-1]);
eofs.raw_exp1 = transpose(exp_coeffs(1:nt, end));
eofs.raw_exp2 = transpose(exp_coeffs(1:nt, end-1));

% Get the first two EOFs, and expansion coefficients
eofs.u_eof1 = transpose(eof_mtx(1:nx, end)) * u_proj_var;
eofs.u_eof2 = transpose(eof_mtx(1:nx, end-1)) * u_proj_var;
eofs.q_eof1 = transpose(eof_mtx(nx+1:end, end)) * q_proj_var;
eofs.q_eof2 = transpose(eof_mtx(nx+1:end, end-1)) * q_proj_var;

eofs.scf = sort(diag(L) / trace(L), 'descend'); % Scaled covariance fraction

% Scale EOFs, expansion coefficients so EOFs have unit norm
eof1_norm = norm(eofs.raw_eof1)^2;
eof2_norm = norm(eofs.raw_eof2)^2;

eofs.u_eof1 = eofs.u_eof1 / eof1_norm;
eofs.u_eof2 = eofs.u_eof2 / eof2_norm;
eofs.q_eof1 = eofs.q_eof1 / eof1_norm;
eofs.q_eof2 = eofs.q_eof2 / eof2_norm;

eofs.exp1 = eofs.raw_exp1 * eof1_norm;
eofs.exp2 = eofs.raw_exp2 * eof2_norm;

eofs = eof_recon_mjo(params, eofs);

end

